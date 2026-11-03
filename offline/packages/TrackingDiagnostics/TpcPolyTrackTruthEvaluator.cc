#include "TpcPolyTrackTruthEvaluator.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <ffaobjects/EventHeader.h>

#include <g4main/PHG4EventHeader.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4MCProcessDefs.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <TFile.h>
#include <TNamed.h>
#include <TObject.h>
#include <TParameter.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
  using HitId = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  struct ContributorMetrics
  {
    std::set<HitId> shared_hits;
    std::set<unsigned int> shared_layers;
    unsigned int exclusive_hits{0};
    double fractional_weight{0.0};
  };

  struct TrackResult
  {
    const Tpc_PolyTrack *track{nullptr};
    std::set<HitId> raw_hits;
    std::set<unsigned int> layers;
    unsigned int associated_hits{0};
    unsigned int ambiguous_hits{0};
    unsigned int missing_cluster_keys{0};
    unsigned int cluster_source_mismatches{0};
    unsigned int missing_raw_hit_references{0};
    std::map<int, ContributorMetrics> contributors;
    int best_truth_track_id{0};
    int overlap_assigned_daughter_role{0};
  };

  struct AssignmentScore
  {
    int matched_daughters{0};
    unsigned int shared_layers{0};
    unsigned int shared_hits{0};
    double fractional_weight{0.0};
    unsigned int plus_track_id{std::numeric_limits<unsigned int>::max()};
    unsigned int minus_track_id{std::numeric_limits<unsigned int>::max()};
  };

  bool better_assignment(const AssignmentScore &lhs, const AssignmentScore &rhs)
  {
    if (lhs.matched_daughters != rhs.matched_daughters)
    {
      return lhs.matched_daughters > rhs.matched_daughters;
    }
    if (lhs.shared_layers != rhs.shared_layers)
    {
      return lhs.shared_layers > rhs.shared_layers;
    }
    if (lhs.shared_hits != rhs.shared_hits)
    {
      return lhs.shared_hits > rhs.shared_hits;
    }
    if (std::fabs(lhs.fractional_weight - rhs.fractional_weight) > 1.e-12)
    {
      return lhs.fractional_weight > rhs.fractional_weight;
    }
    return std::tie(lhs.plus_track_id, lhs.minus_track_id) <
           std::tie(rhs.plus_track_id, rhs.minus_track_id);
  }

  float ratio(const double numerator, const double denominator)
  {
    return denominator > 0.0 ? static_cast<float>(numerator / denominator) : 0.F;
  }

  float transverse_momentum(const double px, const double py)
  {
    return static_cast<float>(std::hypot(px, py));
  }

  float pseudorapidity(const double px, const double py, const double pz)
  {
    const double pt = std::hypot(px, py);
    if (pt <= 0.0)
    {
      return std::numeric_limits<float>::quiet_NaN();
    }
    return static_cast<float>(std::asinh(pz / pt));
  }

  float rapidity(const double e, const double pz)
  {
    if (!(e > std::fabs(pz)))
    {
      return std::numeric_limits<float>::quiet_NaN();
    }
    return static_cast<float>(0.5 * std::log((e + pz) / (e - pz)));
  }

  int signed_charge(const double charge)
  {
    return charge > 0.0 ? 1 : (charge < 0.0 ? -1 : 0);
  }

  int parent_pid(const PHG4TruthInfoContainer *truth_info, const PHG4Particle *particle)
  {
    if (!truth_info || !particle)
    {
      return 0;
    }
    const auto *parent = truth_info->GetParticle(particle->get_parent_id());
    return parent ? parent->get_pid() : 0;
  }
}  // namespace

TpcPolyTrackTruthEvaluator::TpcPolyTrackTruthEvaluator(
    const std::string &name,
    const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

TpcPolyTrackTruthEvaluator::~TpcPolyTrackTruthEvaluator()
{
  if (m_file)
  {
    if (m_file->IsOpen())
    {
      m_file->Close();
    }
    delete m_file;
    m_file = nullptr;
  }
}

int TpcPolyTrackTruthEvaluator::Init(PHCompositeNode * /*topNode*/)
{
  if (!m_primary_pion_mode && (m_parent_pdg == 0 || m_positive_daughter_pdg == 0 ||
      m_negative_daughter_pdg == 0 ||
      m_positive_daughter_pdg == m_negative_daughter_pdg))
  {
    std::cerr << Name() << "::Init - invalid decay-chain PDGs: parent="
              << m_parent_pdg << ", positive daughter="
              << m_positive_daughter_pdg << ", negative daughter="
              << m_negative_daughter_pdg << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_file = TFile::Open(m_filename.c_str(), "CREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cerr << Name() << "::Init - could not create " << m_filename << std::endl;
    if (m_file)
    {
      delete m_file;
      m_file = nullptr;
    }
    return Fun4AllReturnCodes::ABORTRUN;
  }

  create_trees();
  TNamed schema_name("schema_name", "TpcPolyTrackTruthEvaluator");
  TNamed matching_contract(
      "matching_contract",
      "unique_raw_hit_binary_overlap_with_equal_split_ambiguity_diagnostic_v1");
  TNamed assignment_contract(
      "assignment_contract",
      m_primary_pion_mode
          ? "primary_pion_quality_assignment_in_analysis_from_unchanged_matchTree_v1"
          : "distinct_polytracks_any_nonzero_daughter_raw_hit_overlap_no_quality_cut_v1");
  TNamed truth_selection("truth_selection", m_primary_pion_mode
      ? "generated_primary_charged_pions_including_zero_hit_v1"
      : "configured_primary_V0_direct_decay_daughters_v1");
  TParameter<int> schema_version("schema_version", 1);
  TParameter<int> parent_pdg("parent_pdg", m_primary_pion_mode ? 0 : m_parent_pdg);
  TParameter<int> positive_daughter_pdg(
      "positive_daughter_pdg", m_primary_pion_mode ? 0 : m_positive_daughter_pdg);
  TParameter<int> negative_daughter_pdg(
      "negative_daughter_pdg", m_primary_pion_mode ? 0 : m_negative_daughter_pdg);
  if (schema_name.Write() <= 0 || matching_contract.Write() <= 0 ||
      assignment_contract.Write() <= 0 || truth_selection.Write() <= 0 || schema_version.Write() <= 0 ||
      parent_pdg.Write() <= 0 || positive_daughter_pdg.Write() <= 0 ||
      negative_daughter_pdg.Write() <= 0 ||
      m_file->TestBit(TFile::kWriteError))
  {
    std::cerr << Name() << "::Init - failed to write schema metadata to "
              << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void TpcPolyTrackTruthEvaluator::reset_rows()
{
  m_event = EventRow{};
  m_truth = TruthRow{};
  m_track = TrackRow{};
  m_match = MatchRow{};
}

void TpcPolyTrackTruthEvaluator::create_trees()
{
  m_file->cd();
  m_event_tree = new TTree("eventTree", "TPC polytrack truth event summary");
  m_truth_tree = new TTree("truthTree", "Truth particles relevant to TPC polytrack matching");
  m_track_tree = new TTree("trackTree", "TPC polytracks and their best truth match");
  m_match_tree = new TTree("matchTree", "All nonzero TPC polytrack-to-truth overlaps");

#define BRANCH_ULL(tree, row, member) tree->Branch(#member, &(row).member, #member "/l")
#define BRANCH_UINT(tree, row, member) tree->Branch(#member, &(row).member, #member "/i")
#define BRANCH_INT(tree, row, member) tree->Branch(#member, &(row).member, #member "/I")
#define BRANCH_FLOAT(tree, row, member) tree->Branch(#member, &(row).member, #member "/F")

  BRANCH_UINT(m_event_tree, m_event, schema_version);
  BRANCH_ULL(m_event_tree, m_event, sample_id);
  BRANCH_ULL(m_event_tree, m_event, entry);
  BRANCH_INT(m_event_tree, m_event, run);
  BRANCH_INT(m_event_tree, m_event, event);
  BRANCH_INT(m_event_tree, m_event, has_event_header);
  BRANCH_UINT(m_event_tree, m_event, n_truth_particles);
  BRANCH_UINT(m_event_tree, m_event, n_ionizing_g4hits);
  BRANCH_UINT(m_event_tree, m_event, n_raw_hits);
  BRANCH_UINT(m_event_tree, m_event, n_truth_associated_raw_hits);
  BRANCH_UINT(m_event_tree, m_event, n_ambiguous_raw_hits);
  BRANCH_UINT(m_event_tree, m_event, n_polytracks);
  BRANCH_UINT(m_event_tree, m_event, n_valid_polytracks);
  BRANCH_UINT(m_event_tree, m_event, n_generated_primary_particles);
  BRANCH_UINT(m_event_tree, m_event, n_generated_primary_pions);
  BRANCH_INT(m_event_tree, m_event, primary_pion_truth_track_id);
  BRANCH_UINT(m_event_tree, m_event, n_primary_parent);
  BRANCH_UINT(m_event_tree, m_event, n_direct_positive_daughter);
  BRANCH_UINT(m_event_tree, m_event, n_direct_negative_daughter);
  BRANCH_INT(m_event_tree, m_event, parent_truth_track_id);
  BRANCH_INT(m_event_tree, m_event, positive_daughter_truth_track_id);
  BRANCH_INT(m_event_tree, m_event, negative_daughter_truth_track_id);
  BRANCH_INT(m_event_tree, m_event, positive_daughter_overlap_poly_track_id);
  BRANCH_INT(m_event_tree, m_event, negative_daughter_overlap_poly_track_id);
  BRANCH_INT(m_event_tree, m_event, positive_daughter_has_overlap_assignment);
  BRANCH_INT(m_event_tree, m_event, negative_daughter_has_overlap_assignment);
  BRANCH_UINT(m_event_tree, m_event, n_primary_kshort);
  BRANCH_UINT(m_event_tree, m_event, n_direct_piplus);
  BRANCH_UINT(m_event_tree, m_event, n_direct_piminus);
  BRANCH_UINT(m_event_tree, m_event, n_missing_cluster_keys);
  BRANCH_UINT(m_event_tree, m_event, n_cluster_source_mismatches);
  BRANCH_UINT(m_event_tree, m_event, n_missing_raw_hit_references);
  BRANCH_UINT(m_event_tree, m_event, n_missing_g4hits);
  BRANCH_UINT(m_event_tree, m_event, n_missing_truth_particles);
  BRANCH_INT(m_event_tree, m_event, kshort_truth_track_id);
  BRANCH_INT(m_event_tree, m_event, piplus_truth_track_id);
  BRANCH_INT(m_event_tree, m_event, piminus_truth_track_id);
  BRANCH_INT(m_event_tree, m_event, piplus_overlap_poly_track_id);
  BRANCH_INT(m_event_tree, m_event, piminus_overlap_poly_track_id);
  BRANCH_INT(m_event_tree, m_event, piplus_has_overlap_assignment);
  BRANCH_INT(m_event_tree, m_event, piminus_has_overlap_assignment);
  BRANCH_INT(m_event_tree, m_event, both_daughters_present);
  BRANCH_INT(m_event_tree, m_event, both_daughters_have_g4hits);
  BRANCH_INT(m_event_tree, m_event, both_daughters_have_digitized_hits);
  BRANCH_INT(m_event_tree, m_event, both_daughters_overlap_assigned_distinct);

  BRANCH_ULL(m_truth_tree, m_truth, sample_id);
  BRANCH_ULL(m_truth_tree, m_truth, entry);
  BRANCH_INT(m_truth_tree, m_truth, run);
  BRANCH_INT(m_truth_tree, m_truth, event);
  BRANCH_INT(m_truth_tree, m_truth, truth_track_id);
  BRANCH_INT(m_truth_tree, m_truth, pid);
  BRANCH_INT(m_truth_tree, m_truth, parent_id);
  BRANCH_INT(m_truth_tree, m_truth, parent_pid);
  BRANCH_INT(m_truth_tree, m_truth, primary_id);
  BRANCH_INT(m_truth_tree, m_truth, vertex_id);
  BRANCH_INT(m_truth_tree, m_truth, barcode);
  BRANCH_INT(m_truth_tree, m_truth, embed_id);
  BRANCH_INT(m_truth_tree, m_truth, is_primary);
  BRANCH_INT(m_truth_tree, m_truth, is_selected_primary_pion);
  BRANCH_INT(m_truth_tree, m_truth, is_selected_parent);
  BRANCH_INT(m_truth_tree, m_truth, is_selected_kshort);
  BRANCH_INT(m_truth_tree, m_truth, direct_daughter_role);
  BRANCH_FLOAT(m_truth_tree, m_truth, px);
  BRANCH_FLOAT(m_truth_tree, m_truth, py);
  BRANCH_FLOAT(m_truth_tree, m_truth, pz);
  BRANCH_FLOAT(m_truth_tree, m_truth, e);
  BRANCH_FLOAT(m_truth_tree, m_truth, pt);
  BRANCH_FLOAT(m_truth_tree, m_truth, eta);
  BRANCH_FLOAT(m_truth_tree, m_truth, phi);
  BRANCH_FLOAT(m_truth_tree, m_truth, rapidity);
  BRANCH_FLOAT(m_truth_tree, m_truth, vx);
  BRANCH_FLOAT(m_truth_tree, m_truth, vy);
  BRANCH_FLOAT(m_truth_tree, m_truth, vz);
  BRANCH_FLOAT(m_truth_tree, m_truth, vt);
  BRANCH_FLOAT(m_truth_tree, m_truth, vertex_radius);
  BRANCH_UINT(m_truth_tree, m_truth, n_ionizing_g4hits);
  BRANCH_UINT(m_truth_tree, m_truth, n_ionizing_g4hit_layers);
  BRANCH_INT(m_truth_tree, m_truth, first_ionizing_g4hit_layer);
  BRANCH_INT(m_truth_tree, m_truth, last_ionizing_g4hit_layer);
  BRANCH_UINT(m_truth_tree, m_truth, n_available_raw_hits);
  BRANCH_UINT(m_truth_tree, m_truth, n_available_layers);
  BRANCH_INT(m_truth_tree, m_truth, first_layer);
  BRANCH_INT(m_truth_tree, m_truth, last_layer);
  BRANCH_INT(m_truth_tree, m_truth, overlap_assigned_poly_track_id);

  BRANCH_ULL(m_track_tree, m_track, sample_id);
  BRANCH_ULL(m_track_tree, m_track, entry);
  BRANCH_INT(m_track_tree, m_track, run);
  BRANCH_INT(m_track_tree, m_track, event);
  BRANCH_UINT(m_track_tree, m_track, poly_track_id);
  BRANCH_UINT(m_track_tree, m_track, source_assembled_track_id);
  BRANCH_INT(m_track_tree, m_track, fit_status);
  BRANCH_INT(m_track_tree, m_track, charge);
  BRANCH_UINT(m_track_tree, m_track, nclusters);
  BRANCH_FLOAT(m_track_tree, m_track, px);
  BRANCH_FLOAT(m_track_tree, m_track, py);
  BRANCH_FLOAT(m_track_tree, m_track, pz);
  BRANCH_FLOAT(m_track_tree, m_track, pt);
  BRANCH_FLOAT(m_track_tree, m_track, eta);
  BRANCH_FLOAT(m_track_tree, m_track, phi);
  BRANCH_UINT(m_track_tree, m_track, n_raw_hits);
  BRANCH_UINT(m_track_tree, m_track, n_layers);
  BRANCH_UINT(m_track_tree, m_track, n_associated_raw_hits);
  BRANCH_UINT(m_track_tree, m_track, n_unassociated_raw_hits);
  BRANCH_UINT(m_track_tree, m_track, n_ambiguous_raw_hits);
  BRANCH_UINT(m_track_tree, m_track, n_truth_contributors);
  BRANCH_UINT(m_track_tree, m_track, n_missing_cluster_keys);
  BRANCH_UINT(m_track_tree, m_track, n_cluster_source_mismatches);
  BRANCH_UINT(m_track_tree, m_track, n_missing_raw_hit_references);
  BRANCH_INT(m_track_tree, m_track, best_truth_track_id);
  BRANCH_INT(m_track_tree, m_track, best_pid);
  BRANCH_INT(m_track_tree, m_track, best_parent_id);
  BRANCH_INT(m_track_tree, m_track, best_parent_pid);
  BRANCH_UINT(m_track_tree, m_track, best_shared_raw_hits);
  BRANCH_UINT(m_track_tree, m_track, best_shared_layers);
  BRANCH_UINT(m_track_tree, m_track, best_exclusive_raw_hits);
  BRANCH_FLOAT(m_track_tree, m_track, best_fractional_shared_weight);
  BRANCH_FLOAT(m_track_tree, m_track, best_hit_overlap_fraction);
  BRANCH_FLOAT(m_track_tree, m_track, best_associated_hit_overlap_fraction);
  BRANCH_FLOAT(m_track_tree, m_track, best_fractional_purity);
  BRANCH_FLOAT(m_track_tree, m_track, best_completeness);
  BRANCH_FLOAT(m_track_tree, m_track, best_layer_completeness);
  BRANCH_INT(m_track_tree, m_track, best_direct_daughter_role);
  BRANCH_INT(m_track_tree, m_track, overlap_assigned_daughter_role);

  BRANCH_ULL(m_match_tree, m_match, sample_id);
  BRANCH_ULL(m_match_tree, m_match, entry);
  BRANCH_INT(m_match_tree, m_match, run);
  BRANCH_INT(m_match_tree, m_match, event);
  BRANCH_UINT(m_match_tree, m_match, poly_track_id);
  BRANCH_UINT(m_match_tree, m_match, source_assembled_track_id);
  BRANCH_INT(m_match_tree, m_match, truth_track_id);
  BRANCH_INT(m_match_tree, m_match, pid);
  BRANCH_INT(m_match_tree, m_match, parent_id);
  BRANCH_INT(m_match_tree, m_match, parent_pid);
  BRANCH_UINT(m_match_tree, m_match, reco_raw_hits);
  BRANCH_UINT(m_match_tree, m_match, reco_associated_raw_hits);
  BRANCH_UINT(m_match_tree, m_match, truth_available_raw_hits);
  BRANCH_UINT(m_match_tree, m_match, truth_available_layers);
  BRANCH_UINT(m_match_tree, m_match, shared_raw_hits);
  BRANCH_UINT(m_match_tree, m_match, shared_layers);
  BRANCH_UINT(m_match_tree, m_match, exclusive_raw_hits);
  BRANCH_FLOAT(m_match_tree, m_match, fractional_shared_weight);
  BRANCH_FLOAT(m_match_tree, m_match, hit_overlap_fraction);
  BRANCH_FLOAT(m_match_tree, m_match, associated_hit_overlap_fraction);
  BRANCH_FLOAT(m_match_tree, m_match, fractional_purity);
  BRANCH_FLOAT(m_match_tree, m_match, completeness);
  BRANCH_FLOAT(m_match_tree, m_match, layer_completeness);
  BRANCH_INT(m_match_tree, m_match, is_best_match);
  BRANCH_INT(m_match_tree, m_match, direct_daughter_role);
  BRANCH_INT(m_match_tree, m_match, overlap_assigned_daughter_role);

#undef BRANCH_ULL
#undef BRANCH_UINT
#undef BRANCH_INT
#undef BRANCH_FLOAT
}

int TpcPolyTrackTruthEvaluator::process_event(PHCompositeNode *topNode)
{
  auto *truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_truth_info_node);
  auto *g4hits = findNode::getClass<PHG4HitContainer>(topNode, m_g4hit_tpc_node);
  auto *hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, m_hitset_node);
  auto *hit_truth_assoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, m_hit_truth_assoc_node);
  auto *poly_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_poly_cluster_node);
  auto *poly_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_poly_track_node);

  if (!truth_info || !g4hits || !hitsets || !hit_truth_assoc || !poly_clusters || !poly_tracks)
  {
    std::cerr << Name() << "::process_event - required node missing:"
              << " truth_info=" << static_cast<bool>(truth_info)
              << " g4hits=" << static_cast<bool>(g4hits)
              << " hitsets=" << static_cast<bool>(hitsets)
              << " hit_truth_assoc=" << static_cast<bool>(hit_truth_assoc)
              << " poly_clusters=" << static_cast<bool>(poly_clusters)
              << " poly_tracks=" << static_cast<bool>(poly_tracks)
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  reset_rows();
  m_event.sample_id = m_sample_id;
  m_event.entry = m_entry;
  m_event.event = static_cast<int>(m_entry);
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, m_event_header_node))
  {
    m_event.has_event_header = 1;
    m_event.run = event_header->get_RunNumber();
    m_event.event = event_header->get_EvtSequence();
  }
  else if (auto *g4_event_header = findNode::getClass<PHG4EventHeader>(topNode, m_event_header_node))
  {
    m_event.has_event_header = 1;
    m_event.run = 1;
    m_event.event = g4_event_header->get_EvtSequence();
  }

  m_event.n_truth_particles = truth_info->size();
  m_event.n_polytracks = poly_tracks->size();

  // Map each actual digitized TPC hit to distinct contributing truth tracks.
  // The same raw hit may appear in several truth denominators when it is a
  // genuinely merged hit.  Repeated G4 hits from one truth track do not
  // multiply its count.
  std::map<HitId, std::set<int>> raw_hit_truth_ids;
  std::map<int, std::set<HitId>> truth_available_hits;
  std::map<int, std::set<unsigned int>> truth_available_layers;
  std::map<int, unsigned int> truth_g4hit_counts;
  std::map<int, std::set<unsigned int>> truth_g4hit_layers;
  std::set<PHG4HitDefs::keytype> missing_g4hit_keys;
  std::set<int> missing_truth_ids;

  // Pre-digitization coverage counts ionizing G4HIT_TPC entries only.  A hit
  // enters this denominator when eion is finite and strictly positive; zero
  // energy steps and geantino markers are intentionally excluded.
  const auto g4hit_range = g4hits->getHits();
  for (auto g4hit_iter = g4hit_range.first; g4hit_iter != g4hit_range.second; ++g4hit_iter)
  {
    const auto *g4hit = g4hit_iter->second;
    if (!g4hit || !std::isfinite(g4hit->get_eion()) || g4hit->get_eion() <= 0.F)
    {
      continue;
    }
    ++m_event.n_ionizing_g4hits;
    ++truth_g4hit_counts[g4hit->get_trkid()];
    truth_g4hit_layers[g4hit->get_trkid()].insert(g4hit->get_layer());
    if (!truth_info->GetParticle(g4hit->get_trkid()))
    {
      missing_truth_ids.insert(g4hit->get_trkid());
    }
  }

  const auto hitset_range = hitsets->getHitSets(TrkrDefs::tpcId);
  for (auto hitset_iter = hitset_range.first; hitset_iter != hitset_range.second; ++hitset_iter)
  {
    const auto hitsetkey = hitset_iter->first;
    const auto *hitset = hitset_iter->second;
    if (!hitset)
    {
      continue;
    }
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    const auto hit_range = hitset->getHits();
    for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
    {
      const HitId hit_id{hitsetkey, hit_iter->first};
      ++m_event.n_raw_hits;

      TrkrHitTruthAssoc::MMap associations;
      hit_truth_assoc->getG4Hits(hitsetkey, hit_iter->first, associations);
      std::set<PHG4HitDefs::keytype> unique_g4hit_keys;
      for (const auto &association : associations)
      {
        unique_g4hit_keys.insert(association.second.second);
      }

      auto &truth_ids = raw_hit_truth_ids[hit_id];
      for (const auto g4hit_key : unique_g4hit_keys)
      {
        const auto *g4hit = g4hits->findHit(g4hit_key);
        if (!g4hit)
        {
          missing_g4hit_keys.insert(g4hit_key);
          continue;
        }
        const int truth_id = g4hit->get_trkid();
        truth_ids.insert(truth_id);
        if (!truth_info->GetParticle(truth_id))
        {
          missing_truth_ids.insert(truth_id);
        }
      }

      if (!truth_ids.empty())
      {
        ++m_event.n_truth_associated_raw_hits;
      }
      if (truth_ids.size() > 1)
      {
        ++m_event.n_ambiguous_raw_hits;
      }
      for (const int truth_id : truth_ids)
      {
        truth_available_hits[truth_id].insert(hit_id);
        truth_available_layers[truth_id].insert(layer);
      }
    }
  }
  m_event.n_missing_g4hits = missing_g4hit_keys.size();
  m_event.n_missing_truth_particles = missing_truth_ids.size();
  if (!missing_g4hit_keys.empty() || !missing_truth_ids.empty())
  {
    std::cerr << Name() << "::process_event - incomplete truth provenance in run "
              << m_event.run << " event " << m_event.event
              << ": missing G4 hits=" << missing_g4hit_keys.size()
              << ", missing truth particles=" << missing_truth_ids.size()
              << std::endl;
    unsigned int printed = 0;
    for (const auto key : missing_g4hit_keys)
    {
      if (printed++ >= 5)
      {
        break;
      }
      std::cerr << "  missing G4HIT_TPC key " << key << std::endl;
    }
    printed = 0;
    for (const int id : missing_truth_ids)
    {
      if (printed++ >= 5)
      {
        break;
      }
      std::cerr << "  missing G4TruthInfo track ID " << id << std::endl;
    }
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Select the configured generated primary V0 and its direct charged decay
  // daughters. A material-interaction secondary can have the V0 as its direct
  // parent too, so parent ID alone is not a decay classification.
  // Iteration is ordered by truth track ID, making multiple-candidate behavior
  // deterministic while event counters retain the anomaly for QA.
  int kshort_id = 0;
  bool have_kshort = false;
  const auto particle_range = truth_info->GetParticleRange();
  std::set<int> generated_primary_pion_ids;
  for (auto particle_iter = particle_range.first; particle_iter != particle_range.second; ++particle_iter)
  {
    const auto *particle = particle_iter->second;
    if (particle && truth_info->is_primary(particle))
    {
      ++m_event.n_generated_primary_particles;
      if (std::abs(particle->get_pid()) == 211)
      {
        ++m_event.n_generated_primary_pions;
        generated_primary_pion_ids.insert(particle->get_track_id());
        if (m_event.n_generated_primary_pions == 1)
        {
          m_event.primary_pion_truth_track_id = particle->get_track_id();
        }
      }
    }
    if (!m_primary_pion_mode && particle && particle->get_pid() == m_parent_pdg &&
        (particle->get_parent_id() == 0 || truth_info->is_primary(particle)))
    {
      ++m_event.n_primary_parent;
      if (!have_kshort)
      {
        have_kshort = true;
        kshort_id = particle->get_track_id();
      }
    }
  }
  if (have_kshort)
  {
    m_event.parent_truth_track_id = kshort_id;
  }

  int piplus_id = 0;
  int piminus_id = 0;
  bool have_piplus = false;
  bool have_piminus = false;
  if (have_kshort)
  {
    for (auto particle_iter = particle_range.first; particle_iter != particle_range.second; ++particle_iter)
    {
      const auto *particle = particle_iter->second;
      if (!particle || particle->get_parent_id() != kshort_id)
      {
        continue;
      }
      const auto *vertex = truth_info->GetVtx(particle->get_vtx_id());
      if (!vertex || vertex->get_process() < 0 ||
          vertex->get_process() >= kMaxMCProcess ||
          vertex->get_process() == PHG4MCProcess::kPNoProcess)
      {
        std::cerr << Name() << "::process_event - direct configured-parent child track "
                  << particle->get_track_id()
                  << " has missing or invalid production-process metadata"
                  << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      if (vertex->get_process() != PHG4MCProcess::kPDecay)
      {
        continue;
      }
      if (particle->get_pid() == m_positive_daughter_pdg)
      {
        ++m_event.n_direct_positive_daughter;
        if (!have_piplus)
        {
          have_piplus = true;
          piplus_id = particle->get_track_id();
        }
      }
      else if (particle->get_pid() == m_negative_daughter_pdg)
      {
        ++m_event.n_direct_negative_daughter;
        if (!have_piminus)
        {
          have_piminus = true;
          piminus_id = particle->get_track_id();
        }
      }
    }
  }
  if (have_piplus)
  {
    m_event.positive_daughter_truth_track_id = piplus_id;
  }
  if (have_piminus)
  {
    m_event.negative_daughter_truth_track_id = piminus_id;
  }
  // Preserve the original Kshort branch contract for existing consumers.
  if (m_parent_pdg == 310 && m_positive_daughter_pdg == 211 &&
      m_negative_daughter_pdg == -211)
  {
    m_event.n_primary_kshort = m_event.n_primary_parent;
    m_event.n_direct_piplus = m_event.n_direct_positive_daughter;
    m_event.n_direct_piminus = m_event.n_direct_negative_daughter;
    m_event.kshort_truth_track_id = m_event.parent_truth_track_id;
    m_event.piplus_truth_track_id = m_event.positive_daughter_truth_track_id;
    m_event.piminus_truth_track_id = m_event.negative_daughter_truth_track_id;
  }
  m_event.both_daughters_present = have_piplus && have_piminus;
  m_event.both_daughters_have_g4hits =
      have_piplus && have_piminus && truth_g4hit_counts[piplus_id] > 0 &&
      truth_g4hit_counts[piminus_id] > 0;
  m_event.both_daughters_have_digitized_hits =
      have_piplus && have_piminus && !truth_available_hits[piplus_id].empty() &&
      !truth_available_hits[piminus_id].empty();

  std::map<TrkrDefs::cluskey, const Tpc_PolyCluster *> clusters_by_key;
  for (unsigned int cluster_index = 0; cluster_index < poly_clusters->size(); ++cluster_index)
  {
    const auto *cluster = poly_clusters->get_cluster(cluster_index);
    if (cluster)
    {
      clusters_by_key.emplace(cluster->get_trkr_cluster_key(), cluster);
    }
  }

  std::vector<TrackResult> results;
  results.reserve(poly_tracks->size());
  for (unsigned int track_index = 0; track_index < poly_tracks->size(); ++track_index)
  {
    const auto *track = poly_tracks->get_track(track_index);
    if (!track || !track->isValid())
    {
      continue;
    }
    ++m_event.n_valid_polytracks;

    TrackResult result;
    result.track = track;
    std::set<TrkrDefs::cluskey> used_cluster_keys;
    if (track->get_cluster_keys().empty() && track->get_nclusters() > 0)
    {
      // This evaluator requires the current key-bearing polytrack schema.  A
      // source-ID fallback could silently attach clusters not used by a track.
      result.missing_cluster_keys = track->get_nclusters();
    }
    for (const auto cluster_key : track->get_cluster_keys())
    {
      if (!used_cluster_keys.insert(cluster_key).second)
      {
        continue;
      }
      const auto cluster_iter = clusters_by_key.find(cluster_key);
      if (cluster_iter == clusters_by_key.end() || !cluster_iter->second)
      {
        ++result.missing_cluster_keys;
        continue;
      }
      const auto *cluster = cluster_iter->second;
      if (cluster->get_source_assembled_track_id() != track->get_source_assembled_track_id())
      {
        ++result.cluster_source_mismatches;
      }
      for (const auto &hit_id : cluster->get_hit_indices())
      {
        result.raw_hits.insert(hit_id);
        result.layers.insert(TrkrDefs::getLayer(hit_id.first));
      }
    }

    for (const auto &hit_id : result.raw_hits)
    {
      const auto truth_iter = raw_hit_truth_ids.find(hit_id);
      if (truth_iter == raw_hit_truth_ids.end())
      {
        ++result.missing_raw_hit_references;
        continue;
      }
      if (truth_iter->second.empty())
      {
        continue;
      }
      ++result.associated_hits;
      const auto &truth_ids = truth_iter->second;
      if (truth_ids.size() > 1)
      {
        ++result.ambiguous_hits;
      }
      const double fractional_weight = 1.0 / static_cast<double>(truth_ids.size());
      for (const int truth_id : truth_ids)
      {
        auto &metrics = result.contributors[truth_id];
        metrics.shared_hits.insert(hit_id);
        metrics.shared_layers.insert(TrkrDefs::getLayer(hit_id.first));
        metrics.fractional_weight += fractional_weight;
        if (truth_ids.size() == 1)
        {
          ++metrics.exclusive_hits;
        }
      }
    }

    for (const auto &contributor : result.contributors)
    {
      if (result.best_truth_track_id == 0)
      {
        result.best_truth_track_id = contributor.first;
        continue;
      }
      const auto &candidate = contributor.second;
      const auto &best = result.contributors.at(result.best_truth_track_id);
      const auto candidate_score = std::make_tuple(
          candidate.shared_hits.size(), candidate.shared_layers.size(),
          candidate.exclusive_hits, candidate.fractional_weight);
      const auto best_score = std::make_tuple(
          best.shared_hits.size(), best.shared_layers.size(),
          best.exclusive_hits, best.fractional_weight);
      if (candidate_score > best_score ||
          (candidate_score == best_score && contributor.first < result.best_truth_track_id))
      {
        result.best_truth_track_id = contributor.first;
      }
    }

    m_event.n_missing_cluster_keys += result.missing_cluster_keys;
    m_event.n_cluster_source_mismatches += result.cluster_source_mismatches;
    m_event.n_missing_raw_hit_references += result.missing_raw_hit_references;
    results.push_back(std::move(result));
  }

  if (m_event.n_missing_cluster_keys != 0 || m_event.n_cluster_source_mismatches != 0 ||
      m_event.n_missing_raw_hit_references != 0)
  {
    std::cerr << Name() << "::process_event - invalid polytrack-to-cluster provenance in run "
              << m_event.run << " event " << m_event.event
              << ": missing cluster keys=" << m_event.n_missing_cluster_keys
              << ", source-ID mismatches=" << m_event.n_cluster_source_mismatches
              << ", missing raw-hit references=" << m_event.n_missing_raw_hit_references
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Enumerate the small two-daughter assignment problem.  A raw-hit overlap is
  // enough to make a candidate; no physics-efficiency threshold is hidden in
  // this evaluator.  Downstream analysis should select on the stored hit/layer
  // overlap and completeness.  Distinct reconstructed tracks are mandatory.
  int selected_plus_index = -1;
  int selected_minus_index = -1;
  AssignmentScore best_assignment;
  bool have_assignment = false;
  std::vector<int> plus_candidates{ -1 };
  std::vector<int> minus_candidates{ -1 };
  for (std::size_t index = 0; index < results.size(); ++index)
  {
    if (have_piplus && results[index].contributors.count(piplus_id) != 0)
    {
      plus_candidates.push_back(static_cast<int>(index));
    }
    if (have_piminus && results[index].contributors.count(piminus_id) != 0)
    {
      minus_candidates.push_back(static_cast<int>(index));
    }
  }
  for (const int plus_index : plus_candidates)
  {
    for (const int minus_index : minus_candidates)
    {
      if (plus_index >= 0 && plus_index == minus_index)
      {
        continue;
      }
      AssignmentScore score;
      for (const auto &candidate : {std::make_pair(plus_index, piplus_id),
                                    std::make_pair(minus_index, piminus_id)})
      {
        if (candidate.first < 0)
        {
          continue;
        }
        ++score.matched_daughters;
        const auto &metrics = results[candidate.first].contributors.at(candidate.second);
        score.shared_layers += metrics.shared_layers.size();
        score.shared_hits += metrics.shared_hits.size();
        score.fractional_weight += metrics.fractional_weight;
      }
      if (plus_index >= 0)
      {
        score.plus_track_id = results[plus_index].track->get_track_id();
      }
      if (minus_index >= 0)
      {
        score.minus_track_id = results[minus_index].track->get_track_id();
      }
      if (!have_assignment || better_assignment(score, best_assignment))
      {
        have_assignment = true;
        best_assignment = score;
        selected_plus_index = plus_index;
        selected_minus_index = minus_index;
      }
    }
  }

  std::map<int, int> assigned_track_by_truth;
  if (selected_plus_index >= 0)
  {
    auto &result = results[selected_plus_index];
    result.overlap_assigned_daughter_role = 1;
    m_event.positive_daughter_has_overlap_assignment = 1;
    m_event.positive_daughter_overlap_poly_track_id = result.track->get_track_id();
    assigned_track_by_truth[piplus_id] = result.track->get_track_id();
  }
  if (selected_minus_index >= 0)
  {
    auto &result = results[selected_minus_index];
    result.overlap_assigned_daughter_role = -1;
    m_event.negative_daughter_has_overlap_assignment = 1;
    m_event.negative_daughter_overlap_poly_track_id = result.track->get_track_id();
    assigned_track_by_truth[piminus_id] = result.track->get_track_id();
  }
  m_event.both_daughters_overlap_assigned_distinct =
      m_event.positive_daughter_has_overlap_assignment &&
      m_event.negative_daughter_has_overlap_assignment &&
      m_event.positive_daughter_overlap_poly_track_id !=
          m_event.negative_daughter_overlap_poly_track_id;
  if (m_parent_pdg == 310 && m_positive_daughter_pdg == 211 &&
      m_negative_daughter_pdg == -211)
  {
    m_event.piplus_has_overlap_assignment =
        m_event.positive_daughter_has_overlap_assignment;
    m_event.piminus_has_overlap_assignment =
        m_event.negative_daughter_has_overlap_assignment;
    m_event.piplus_overlap_poly_track_id =
        m_event.positive_daughter_overlap_poly_track_id;
    m_event.piminus_overlap_poly_track_id =
        m_event.negative_daughter_overlap_poly_track_id;
  }

  const auto daughter_role = [&](const int truth_id) {
    if (have_piplus && truth_id == piplus_id)
    {
      return 1;
    }
    if (have_piminus && truth_id == piminus_id)
    {
      return -1;
    }
    return 0;
  };

  // Fill one row for each valid reconstructed track and one match row for each
  // distinct truth contributor with nonzero raw-hit overlap.
  for (const auto &result : results)
  {
    m_track = TrackRow{};
    m_track.sample_id = m_sample_id;
    m_track.entry = m_entry;
    m_track.run = m_event.run;
    m_track.event = m_event.event;
    m_track.poly_track_id = result.track->get_track_id();
    m_track.source_assembled_track_id = result.track->get_source_assembled_track_id();
    m_track.fit_status = result.track->get_fit_status();
    m_track.charge = signed_charge(result.track->get_charge());
    m_track.nclusters = result.track->get_nclusters();
    m_track.px = result.track->get_px();
    m_track.py = result.track->get_py();
    m_track.pz = result.track->get_pz();
    m_track.pt = transverse_momentum(m_track.px, m_track.py);
    m_track.eta = pseudorapidity(m_track.px, m_track.py, m_track.pz);
    m_track.phi = std::atan2(m_track.py, m_track.px);
    m_track.n_raw_hits = result.raw_hits.size();
    m_track.n_layers = result.layers.size();
    m_track.n_associated_raw_hits = result.associated_hits;
    m_track.n_unassociated_raw_hits = result.raw_hits.size() - result.associated_hits;
    m_track.n_ambiguous_raw_hits = result.ambiguous_hits;
    m_track.n_truth_contributors = result.contributors.size();
    m_track.n_missing_cluster_keys = result.missing_cluster_keys;
    m_track.n_cluster_source_mismatches = result.cluster_source_mismatches;
    m_track.n_missing_raw_hit_references = result.missing_raw_hit_references;
    m_track.best_truth_track_id = result.best_truth_track_id;
    m_track.overlap_assigned_daughter_role = result.overlap_assigned_daughter_role;

    if (result.best_truth_track_id != 0)
    {
      const auto &metrics = result.contributors.at(result.best_truth_track_id);
      const auto *particle = truth_info->GetParticle(result.best_truth_track_id);
      m_track.best_pid = particle ? particle->get_pid() : 0;
      m_track.best_parent_id = particle ? particle->get_parent_id() : 0;
      m_track.best_parent_pid = parent_pid(truth_info, particle);
      m_track.best_shared_raw_hits = metrics.shared_hits.size();
      m_track.best_shared_layers = metrics.shared_layers.size();
      m_track.best_exclusive_raw_hits = metrics.exclusive_hits;
      m_track.best_fractional_shared_weight = metrics.fractional_weight;
      m_track.best_hit_overlap_fraction = ratio(metrics.shared_hits.size(), result.raw_hits.size());
      m_track.best_associated_hit_overlap_fraction = ratio(metrics.shared_hits.size(), result.associated_hits);
      m_track.best_fractional_purity = ratio(metrics.fractional_weight, result.raw_hits.size());
      m_track.best_completeness = ratio(
          metrics.shared_hits.size(), truth_available_hits[result.best_truth_track_id].size());
      m_track.best_layer_completeness = ratio(
          metrics.shared_layers.size(), truth_available_layers[result.best_truth_track_id].size());
      m_track.best_direct_daughter_role = daughter_role(result.best_truth_track_id);
    }
    m_track_tree->Fill();

    for (const auto &contributor : result.contributors)
    {
      const int truth_id = contributor.first;
      const auto &metrics = contributor.second;
      const auto *particle = truth_info->GetParticle(truth_id);
      m_match = MatchRow{};
      m_match.sample_id = m_sample_id;
      m_match.entry = m_entry;
      m_match.run = m_event.run;
      m_match.event = m_event.event;
      m_match.poly_track_id = result.track->get_track_id();
      m_match.source_assembled_track_id = result.track->get_source_assembled_track_id();
      m_match.truth_track_id = truth_id;
      m_match.pid = particle ? particle->get_pid() : 0;
      m_match.parent_id = particle ? particle->get_parent_id() : 0;
      m_match.parent_pid = parent_pid(truth_info, particle);
      m_match.reco_raw_hits = result.raw_hits.size();
      m_match.reco_associated_raw_hits = result.associated_hits;
      m_match.truth_available_raw_hits = truth_available_hits[truth_id].size();
      m_match.truth_available_layers = truth_available_layers[truth_id].size();
      m_match.shared_raw_hits = metrics.shared_hits.size();
      m_match.shared_layers = metrics.shared_layers.size();
      m_match.exclusive_raw_hits = metrics.exclusive_hits;
      m_match.fractional_shared_weight = metrics.fractional_weight;
      m_match.hit_overlap_fraction = ratio(metrics.shared_hits.size(), result.raw_hits.size());
      m_match.associated_hit_overlap_fraction = ratio(metrics.shared_hits.size(), result.associated_hits);
      m_match.fractional_purity = ratio(metrics.fractional_weight, result.raw_hits.size());
      m_match.completeness = ratio(metrics.shared_hits.size(), truth_available_hits[truth_id].size());
      m_match.layer_completeness = ratio(metrics.shared_layers.size(), truth_available_layers[truth_id].size());
      m_match.is_best_match = truth_id == result.best_truth_track_id;
      m_match.direct_daughter_role = daughter_role(truth_id);
      m_match.overlap_assigned_daughter_role =
          (result.overlap_assigned_daughter_role == daughter_role(truth_id)) ?
              result.overlap_assigned_daughter_role : 0;
      m_match_tree->Fill();
    }
  }

  // Write the selected parent/daughters even when they have no digitized hits, and
  // all other truth particles that contributed to at least one actual raw hit.
  std::set<int> truth_ids_to_write;
  if (m_primary_pion_mode)
  {
    truth_ids_to_write.insert(generated_primary_pion_ids.begin(), generated_primary_pion_ids.end());
  }
  for (const auto &available : truth_available_hits)
  {
    truth_ids_to_write.insert(available.first);
  }
  for (const auto &g4hit_count : truth_g4hit_counts)
  {
    truth_ids_to_write.insert(g4hit_count.first);
  }
  if (have_kshort)
  {
    truth_ids_to_write.insert(kshort_id);
  }
  if (have_piplus)
  {
    truth_ids_to_write.insert(piplus_id);
  }
  if (have_piminus)
  {
    truth_ids_to_write.insert(piminus_id);
  }

  for (const int truth_id : truth_ids_to_write)
  {
    const auto *particle = truth_info->GetParticle(truth_id);
    if (!particle)
    {
      continue;
    }
    m_truth = TruthRow{};
    m_truth.sample_id = m_sample_id;
    m_truth.entry = m_entry;
    m_truth.run = m_event.run;
    m_truth.event = m_event.event;
    m_truth.truth_track_id = truth_id;
    m_truth.pid = particle->get_pid();
    m_truth.parent_id = particle->get_parent_id();
    m_truth.parent_pid = parent_pid(truth_info, particle);
    m_truth.primary_id = particle->get_primary_id();
    m_truth.vertex_id = particle->get_vtx_id();
    m_truth.barcode = particle->get_barcode();
    m_truth.embed_id = truth_info->isEmbeded(truth_id);
    m_truth.is_primary = truth_info->is_primary(particle);
    m_truth.is_selected_primary_pion = m_primary_pion_mode && generated_primary_pion_ids.count(truth_id);
    m_truth.is_selected_parent = have_kshort && truth_id == kshort_id;
    m_truth.is_selected_kshort =
        m_parent_pdg == 310 && m_truth.is_selected_parent;
    m_truth.direct_daughter_role = daughter_role(truth_id);
    m_truth.px = particle->get_px();
    m_truth.py = particle->get_py();
    m_truth.pz = particle->get_pz();
    m_truth.e = particle->get_e();
    m_truth.pt = transverse_momentum(m_truth.px, m_truth.py);
    m_truth.eta = pseudorapidity(m_truth.px, m_truth.py, m_truth.pz);
    m_truth.phi = std::atan2(m_truth.py, m_truth.px);
    m_truth.rapidity = ::rapidity(m_truth.e, m_truth.pz);
    if (const auto *vertex = truth_info->GetVtx(particle->get_vtx_id()))
    {
      m_truth.vx = vertex->get_x();
      m_truth.vy = vertex->get_y();
      m_truth.vz = vertex->get_z();
      m_truth.vt = vertex->get_t();
      m_truth.vertex_radius = std::hypot(m_truth.vx, m_truth.vy);
    }
    m_truth.n_ionizing_g4hits = truth_g4hit_counts[truth_id];
    m_truth.n_ionizing_g4hit_layers = truth_g4hit_layers[truth_id].size();
    if (!truth_g4hit_layers[truth_id].empty())
    {
      m_truth.first_ionizing_g4hit_layer = *truth_g4hit_layers[truth_id].begin();
      m_truth.last_ionizing_g4hit_layer = *truth_g4hit_layers[truth_id].rbegin();
    }
    m_truth.n_available_raw_hits = truth_available_hits[truth_id].size();
    m_truth.n_available_layers = truth_available_layers[truth_id].size();
    if (!truth_available_layers[truth_id].empty())
    {
      m_truth.first_layer = *truth_available_layers[truth_id].begin();
      m_truth.last_layer = *truth_available_layers[truth_id].rbegin();
    }
    const auto assigned_iter = assigned_track_by_truth.find(truth_id);
    if (assigned_iter != assigned_track_by_truth.end())
    {
      m_truth.overlap_assigned_poly_track_id = assigned_iter->second;
    }
    m_truth_tree->Fill();
  }

  m_event_tree->Fill();
  ++m_entry;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyTrackTruthEvaluator::End(PHCompositeNode * /*topNode*/)
{
  if (!m_file)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_file->cd();
  const bool write_failed =
      m_event_tree->Write("", TObject::kOverwrite) <= 0 ||
      m_truth_tree->Write("", TObject::kOverwrite) <= 0 ||
      m_track_tree->Write("", TObject::kOverwrite) <= 0 ||
      m_match_tree->Write("", TObject::kOverwrite) <= 0;
  m_file->Flush();
  const bool file_write_failed = m_file->TestBit(TFile::kWriteError);
  m_file->Close();
  delete m_file;
  m_file = nullptr;
  m_event_tree = nullptr;
  m_truth_tree = nullptr;
  m_track_tree = nullptr;
  m_match_tree = nullptr;
  if (write_failed || file_write_failed)
  {
    std::cerr << Name() << "::End - failed to write " << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
