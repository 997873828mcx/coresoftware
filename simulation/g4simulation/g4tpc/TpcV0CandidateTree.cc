#include "TpcV0CandidateTree.h"

#include <g4main/PHG4EventHeader.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>

namespace
{
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kPionMass = 0.13957039;
  constexpr double kProtonMass = 0.938272088;

  template <class T>
  constexpr T square(const T &value)
  {
    return value * value;
  }

  double unwrap_to_previous(double theta, const double previous)
  {
    while (theta - previous > kPi)
    {
      theta -= 2.0 * kPi;
    }
    while (theta - previous < -kPi)
    {
      theta += 2.0 * kPi;
    }
    return theta;
  }
}  // namespace

TpcV0CandidateTree::TpcV0CandidateTree(const std::string &name,
                                       const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

void TpcV0CandidateTree::set_primary_vertex(const double x, const double y, const double z)
{
  m_fixed_primary_vertex = {x, y, z};
}

int TpcV0CandidateTree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = new TFile(m_filename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": failed to create output file " << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_pair_tree = new TTree("pairTree", "TPC truth-point V0 candidates");
  create_branches();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcV0CandidateTree::process_event(PHCompositeNode *topNode)
{
  auto *truth_points = findNode::getClass<PHG4HitContainer>(topNode, m_truth_point_node);
  if (!truth_points)
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << Name() << ": missing truth point node "
                << m_truth_point_node << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto *truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_truth_info_node);
  if (!truth_info && Verbosity() > 0)
  {
    std::cout << PHWHERE << Name() << ": missing truth info node "
              << m_truth_info_node << "; no charged truth tracklets can be built" << std::endl;
  }

  const int event_number = get_event_number(topNode);
  const Vec3 primary_vertex = get_primary_vertex(truth_info);
  const auto tracklet_map = build_tracklets(truth_points, truth_info);

  std::vector<const Tracklet *> tracklets;
  tracklets.reserve(tracklet_map.size());
  for (const auto &entry : tracklet_map)
  {
    tracklets.push_back(&entry.second);
  }

  for (std::size_t i = 0; i < tracklets.size(); ++i)
  {
    for (std::size_t j = i + 1; j < tracklets.size(); ++j)
    {
      make_pair_row(*tracklets[i], *tracklets[j], primary_vertex, event_number);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcV0CandidateTree::End(PHCompositeNode * /*topNode*/)
{
  if (m_file)
  {
    m_file->cd();
    if (m_pair_tree)
    {
      m_pair_tree->Write();
    }
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": pair counters: raw=" << m_counter_raw_pairs
              << " reject_charge=" << m_counter_reject_charge
              << " reject_preselection=" << m_counter_reject_preselection
              << " reject_pca=" << m_counter_reject_pca
              << " reject_pointing=" << m_counter_reject_pointing
              << " reject_ap=" << m_counter_reject_ap
              << " written=" << m_counter_written << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcV0CandidateTree::get_event_number(PHCompositeNode *topNode) const
{
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    return event_header->get_EvtSequence();
  }
  if (auto *g4_event_header = findNode::getClass<PHG4EventHeader>(topNode, "EventHeader"))
  {
    return g4_event_header->get_EvtSequence();
  }
  return 0;
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::get_primary_vertex(PHG4TruthInfoContainer *truth_info) const
{
  if (!m_use_truth_primary_vertex || !truth_info)
  {
    return m_fixed_primary_vertex;
  }

  const auto vtx_range = truth_info->GetPrimaryVtxRange();
  for (auto iter = vtx_range.first; iter != vtx_range.second; ++iter)
  {
    const auto *vtx = iter->second;
    if (vtx)
    {
      return {vtx->get_x(), vtx->get_y(), vtx->get_z()};
    }
  }

  return m_fixed_primary_vertex;
}

std::map<int, TpcV0CandidateTree::Tracklet> TpcV0CandidateTree::build_tracklets(
    PHG4HitContainer *truth_points,
    PHG4TruthInfoContainer *truth_info) const
{
  std::map<int, Tracklet> tracklets;
  if (!truth_points)
  {
    return tracklets;
  }

  const auto hit_range = truth_points->getHits();
  for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
  {
    const PHG4Hit *hit = hit_iter->second;
    if (!hit)
    {
      continue;
    }

    const int track_id = hit->get_trkid();
    Tracklet &tracklet = tracklets[track_id];
    if (tracklet.points.empty())
    {
      tracklet.track_id = track_id;
      tracklet.shower_id = hit->get_shower_id();
    }

    TruthPoint point;
    point.track_id = track_id;
    point.shower_id = hit->get_shower_id();
    point.layer = static_cast<int>(hit->get_layer());
    point.position = {hit->get_x(0), hit->get_y(0), hit->get_z(0)};
    point.momentum = {hit->get_px(0), hit->get_py(0), hit->get_pz(0)};
    point.t = hit->get_t(0);
    point.path = hit->get_path_length();
    if (finite(point.position) && finite(point.momentum))
    {
      tracklet.points.push_back(point);
    }
  }

  for (auto iter = tracklets.begin(); iter != tracklets.end();)
  {
    Tracklet &tracklet = iter->second;
    std::sort(tracklet.points.begin(), tracklet.points.end(),
              [](const TruthPoint &lhs, const TruthPoint &rhs)
              { return lhs.path < rhs.path; });
    tracklet.npoints = static_cast<int>(tracklet.points.size());

    if (tracklet.npoints < m_min_points)
    {
      iter = tracklets.erase(iter);
      continue;
    }

    if (truth_info)
    {
      const PHG4Particle *particle = truth_info->GetParticle(tracklet.track_id);
      if (particle)
      {
        tracklet.pid = particle->get_pid();
        tracklet.parent_id = particle->get_parent_id();
        tracklet.primary_id = particle->get_primary_id();
        tracklet.vtx_id = particle->get_vtx_id();
        tracklet.barcode = particle->get_barcode();
        tracklet.embed_id = truth_info->isEmbeded(tracklet.track_id);
        tracklet.is_primary = truth_info->is_primary(particle) ? 1 : 0;
        tracklet.charge = pdg_charge(tracklet.pid);
        tracklet.truth_momentum = {particle->get_px(), particle->get_py(), particle->get_pz()};
        tracklet.truth_e = particle->get_e();

        if (const PHG4Particle *parent = truth_info->GetParticle(tracklet.parent_id))
        {
          tracklet.parent_pid = parent->get_pid();
        }

        if (auto *vtx = truth_info->GetVtx(tracklet.vtx_id))
        {
          tracklet.truth_vertex = {vtx->get_x(), vtx->get_y(), vtx->get_z()};
          tracklet.truth_vt = vtx->get_t();
        }
      }
    }

    if (tracklet.charge == 0)
    {
      iter = tracklets.erase(iter);
      continue;
    }

    if (m_fit_helix_tracks)
    {
      tracklet.has_helix = fit_helix(tracklet.points, m_fit_first_points, m_bfield_t, tracklet.helix);
      if (!tracklet.has_helix)
      {
        iter = tracklets.erase(iter);
        continue;
      }
      tracklet.position = helix_point(tracklet.helix, tracklet.helix.theta_first);
      tracklet.momentum = helix_momentum(tracklet.helix, tracklet.helix.theta_first);
    }
    else
    {
      tracklet.position = tracklet.points.front().position;
      tracklet.momentum = tracklet.points.front().momentum;
    }

    if (!finite(tracklet.truth_momentum) || norm(tracklet.truth_momentum) <= 0.0)
    {
      tracklet.truth_momentum = tracklet.momentum;
    }

    ++iter;
  }

  return tracklets;
}

bool TpcV0CandidateTree::make_pair_row(const Tracklet &track1, const Tracklet &track2,
                                       const Vec3 &primary_vertex,
                                       const int event_number)
{
  ++m_counter_raw_pairs;

  if (track1.charge == track2.charge)
  {
    ++m_counter_reject_charge;
    return false;
  }

  if (!passes_preselection(track1, track2, primary_vertex))
  {
    ++m_counter_reject_preselection;
    return false;
  }

  Vec3 pca1;
  Vec3 pca2;
  Vec3 mom1 = track1.momentum;
  Vec3 mom2 = track2.momentum;
  double pair_dca = 0.0;
  double theta1 = quiet_nan();
  double theta2 = quiet_nan();
  std::pair<double, double> dca1;
  std::pair<double, double> dca2;

  if (m_fit_helix_tracks && track1.has_helix && track2.has_helix)
  {
    auto candidates = helix_helix_pca_candidates(
        track1.helix, track2.helix, m_theta_extension, m_coarse_steps,
        m_downstream_margin, m_pca_candidates);
    if (candidates.empty())
    {
      ++m_counter_reject_pca;
      return false;
    }

    auto best = candidates.front();
    if (m_prefer_positive_pointing)
    {
      double best_score = std::numeric_limits<double>::max();
      for (const auto &candidate : candidates)
      {
        const Vec3 cand_mom1 = helix_momentum(track1.helix, candidate.theta1);
        const Vec3 cand_mom2 = helix_momentum(track2.helix, candidate.theta2);
        const Vec3 cand_vertex = scale(add(candidate.pca1, candidate.pca2), 0.5);
        const Vec3 flight = subtract(cand_vertex, primary_vertex);
        const Vec3 total_mom = add(cand_mom1, cand_mom2);
        const double cos_theta = vector_cosine(flight, total_mom);
        const double penalty = (std::isfinite(cos_theta) && cos_theta > 0.0) ? 0.0 : 1000.0;
        const double score = penalty + candidate.dca - 1e-3 * cos_theta;
        if (score < best_score)
        {
          best_score = score;
          best = candidate;
        }
      }
    }

    pca1 = best.pca1;
    pca2 = best.pca2;
    pair_dca = best.dca;
    theta1 = best.theta1;
    theta2 = best.theta2;
    mom1 = helix_momentum(track1.helix, theta1);
    mom2 = helix_momentum(track2.helix, theta2);
    dca1 = helix_dca_to_vertex(track1.helix, primary_vertex);
    dca2 = helix_dca_to_vertex(track2.helix, primary_vertex);
  }
  else
  {
    LinePca pca;
    if (!line_line_pca(track1.position, track1.momentum, track2.position, track2.momentum, pca, true))
    {
      ++m_counter_reject_pca;
      return false;
    }
    pca1 = pca.pca1;
    pca2 = pca.pca2;
    pair_dca = pca.dca;
    dca1 = track_dca_to_vertex(track1.position, track1.momentum, primary_vertex);
    dca2 = track_dca_to_vertex(track2.position, track2.momentum, primary_vertex);
  }

  const Vec3 pair_vertex = scale(add(pca1, pca2), 0.5);
  const Vec3 total_mom = add(mom1, mom2);
  const Vec3 flight = subtract(pair_vertex, primary_vertex);
  const double cos_theta = vector_cosine(flight, total_mom);
  if (!std::isfinite(cos_theta) || norm(flight) <= 0.0 || norm(total_mom) <= 0.0)
  {
    ++m_counter_reject_pointing;
    return false;
  }

  const Vec3 &pplus = (track1.charge > 0) ? mom1 : mom2;
  const Vec3 &pminus = (track1.charge > 0) ? mom2 : mom1;
  double alpha = 0.0;
  double qt = 0.0;
  if (!armenteros(pplus, pminus, alpha, qt))
  {
    ++m_counter_reject_ap;
    return false;
  }

  const Vec3 &truth_pplus = (track1.charge > 0) ? track1.truth_momentum : track2.truth_momentum;
  const Vec3 &truth_pminus = (track1.charge > 0) ? track2.truth_momentum : track1.truth_momentum;
  double truth_alpha = quiet_nan();
  double truth_qt = quiet_nan();
  armenteros(truth_pplus, truth_pminus, truth_alpha, truth_qt);

  Vec3 true_decay{quiet_nan(), quiet_nan(), quiet_nan()};
  double pca_to_true_3d = quiet_nan();
  double pca_to_true_xy = quiet_nan();
  double pca_to_true_z = quiet_nan();
  if (track1.parent_id != 0 && track1.parent_id == track2.parent_id)
  {
    true_decay = scale(add(track1.truth_vertex, track2.truth_vertex), 0.5);
    const Vec3 delta = subtract(pair_vertex, true_decay);
    pca_to_true_3d = norm(delta);
    pca_to_true_xy = std::sqrt(square(delta.x) + square(delta.y));
    pca_to_true_z = std::abs(delta.z);
  }

  const Vec3 positive_mom = (track1.charge > 0) ? mom1 : mom2;
  const Vec3 negative_mom = (track1.charge > 0) ? mom2 : mom1;

  reset_pair_row();
  m_pair.run = 1;
  m_pair.evt = event_number;
  m_pair.cross1 = 0;
  m_pair.cross2 = 0;
  m_pair.px1 = static_cast<float>(mom1.x);
  m_pair.py1 = static_cast<float>(mom1.y);
  m_pair.pz1 = static_cast<float>(mom1.z);
  m_pair.px2 = static_cast<float>(mom2.x);
  m_pair.py2 = static_cast<float>(mom2.y);
  m_pair.pz2 = static_cast<float>(mom2.z);
  m_pair.dca_xy1 = static_cast<float>(dca1.first);
  m_pair.dca_z1 = static_cast<float>(dca1.second);
  m_pair.dca_xy2 = static_cast<float>(dca2.first);
  m_pair.dca_z2 = static_cast<float>(dca2.second);
  m_pair.pairDCA = static_cast<float>(pair_dca);
  m_pair.alpha = static_cast<float>(alpha);
  m_pair.qT = static_cast<float>(qt);
  m_pair.charge1 = static_cast<float>(track1.charge);
  m_pair.charge2 = static_cast<float>(track2.charge);
  m_pair.cosThetaReco = static_cast<float>(cos_theta);
  m_pair.Lproj = static_cast<float>(norm(flight));

  m_pair.pca_x = static_cast<float>(pair_vertex.x);
  m_pair.pca_y = static_cast<float>(pair_vertex.y);
  m_pair.pca_z = static_cast<float>(pair_vertex.z);
  m_pair.pca1_x = static_cast<float>(pca1.x);
  m_pair.pca1_y = static_cast<float>(pca1.y);
  m_pair.pca1_z = static_cast<float>(pca1.z);
  m_pair.pca2_x = static_cast<float>(pca2.x);
  m_pair.pca2_y = static_cast<float>(pca2.y);
  m_pair.pca2_z = static_cast<float>(pca2.z);

  m_pair.v0_px = static_cast<float>(total_mom.x);
  m_pair.v0_py = static_cast<float>(total_mom.y);
  m_pair.v0_pz = static_cast<float>(total_mom.z);
  m_pair.v0_pt = static_cast<float>(pt(total_mom));
  m_pair.mass_Kshort = static_cast<float>(invariant_mass(mom1, kPionMass, mom2, kPionMass));
  m_pair.mass_Lambda = static_cast<float>(invariant_mass(positive_mom, kProtonMass, negative_mom, kPionMass));
  m_pair.mass_AntiLambda = static_cast<float>(invariant_mass(positive_mom, kPionMass, negative_mom, kProtonMass));

  m_pair.true_decay_x = static_cast<float>(true_decay.x);
  m_pair.true_decay_y = static_cast<float>(true_decay.y);
  m_pair.true_decay_z = static_cast<float>(true_decay.z);
  m_pair.pca_to_true_3d = static_cast<float>(pca_to_true_3d);
  m_pair.pca_to_true_xy = static_cast<float>(pca_to_true_xy);
  m_pair.pca_to_true_z = static_cast<float>(pca_to_true_z);
  m_pair.truth_alpha = static_cast<float>(truth_alpha);
  m_pair.truth_qT = static_cast<float>(truth_qt);
  m_pair.delta_alpha = static_cast<float>(alpha - truth_alpha);
  m_pair.delta_qT = static_cast<float>(qt - truth_qt);
  m_pair.truth_px1 = static_cast<float>(track1.truth_momentum.x);
  m_pair.truth_py1 = static_cast<float>(track1.truth_momentum.y);
  m_pair.truth_pz1 = static_cast<float>(track1.truth_momentum.z);
  m_pair.truth_px2 = static_cast<float>(track2.truth_momentum.x);
  m_pair.truth_py2 = static_cast<float>(track2.truth_momentum.y);
  m_pair.truth_pz2 = static_cast<float>(track2.truth_momentum.z);
  m_pair.cos_mom1_truth = static_cast<float>(vector_cosine(mom1, track1.truth_momentum));
  m_pair.cos_mom2_truth = static_cast<float>(vector_cosine(mom2, track2.truth_momentum));
  m_pair.pca_theta1 = static_cast<float>(theta1);
  m_pair.pca_theta2 = static_cast<float>(theta2);
  m_pair.track_id1 = track1.track_id;
  m_pair.track_id2 = track2.track_id;
  m_pair.pid1 = track1.pid;
  m_pair.pid2 = track2.pid;
  m_pair.parent_id1 = track1.parent_id;
  m_pair.parent_id2 = track2.parent_id;
  m_pair.parent_pid = (track1.parent_id != 0 && track1.parent_id == track2.parent_id) ? track1.parent_pid : 0;
  m_pair.npoints1 = static_cast<short>(track1.npoints);
  m_pair.npoints2 = static_cast<short>(track2.npoints);

  m_pair_tree->Fill();
  ++m_counter_written;
  return true;
}

void TpcV0CandidateTree::reset_pair_row()
{
  m_pair = {};
  const float nan = quiet_nan();
  m_pair.dca_xy1 = nan;
  m_pair.dca_z1 = nan;
  m_pair.dca_xy2 = nan;
  m_pair.dca_z2 = nan;
  m_pair.pairDCA = nan;
  m_pair.alpha = nan;
  m_pair.qT = nan;
  m_pair.cosThetaReco = nan;
  m_pair.Lproj = nan;
  m_pair.pca_x = nan;
  m_pair.pca_y = nan;
  m_pair.pca_z = nan;
  m_pair.pca1_x = nan;
  m_pair.pca1_y = nan;
  m_pair.pca1_z = nan;
  m_pair.pca2_x = nan;
  m_pair.pca2_y = nan;
  m_pair.pca2_z = nan;
  m_pair.v0_px = nan;
  m_pair.v0_py = nan;
  m_pair.v0_pz = nan;
  m_pair.v0_pt = nan;
  m_pair.mass_Kshort = nan;
  m_pair.mass_Lambda = nan;
  m_pair.mass_AntiLambda = nan;
  m_pair.true_decay_x = nan;
  m_pair.true_decay_y = nan;
  m_pair.true_decay_z = nan;
  m_pair.pca_to_true_3d = nan;
  m_pair.pca_to_true_xy = nan;
  m_pair.pca_to_true_z = nan;
  m_pair.truth_alpha = nan;
  m_pair.truth_qT = nan;
  m_pair.delta_alpha = nan;
  m_pair.delta_qT = nan;
  m_pair.truth_px1 = nan;
  m_pair.truth_py1 = nan;
  m_pair.truth_pz1 = nan;
  m_pair.truth_px2 = nan;
  m_pair.truth_py2 = nan;
  m_pair.truth_pz2 = nan;
  m_pair.cos_mom1_truth = nan;
  m_pair.cos_mom2_truth = nan;
  m_pair.pca_theta1 = nan;
  m_pair.pca_theta2 = nan;
}

void TpcV0CandidateTree::create_branches()
{
  m_pair_tree->Branch("run", &m_pair.run, "run/I");
  m_pair_tree->Branch("evt", &m_pair.evt, "evt/I");
  m_pair_tree->Branch("cross1", &m_pair.cross1, "cross1/S");
  m_pair_tree->Branch("cross2", &m_pair.cross2, "cross2/S");
  m_pair_tree->Branch("px1", &m_pair.px1, "px1/F");
  m_pair_tree->Branch("py1", &m_pair.py1, "py1/F");
  m_pair_tree->Branch("pz1", &m_pair.pz1, "pz1/F");
  m_pair_tree->Branch("px2", &m_pair.px2, "px2/F");
  m_pair_tree->Branch("py2", &m_pair.py2, "py2/F");
  m_pair_tree->Branch("pz2", &m_pair.pz2, "pz2/F");
  m_pair_tree->Branch("dca_xy1", &m_pair.dca_xy1, "dca_xy1/F");
  m_pair_tree->Branch("dca_z1", &m_pair.dca_z1, "dca_z1/F");
  m_pair_tree->Branch("dca_xy2", &m_pair.dca_xy2, "dca_xy2/F");
  m_pair_tree->Branch("dca_z2", &m_pair.dca_z2, "dca_z2/F");
  m_pair_tree->Branch("pairDCA", &m_pair.pairDCA, "pairDCA/F");
  m_pair_tree->Branch("alpha", &m_pair.alpha, "alpha/F");
  m_pair_tree->Branch("qT", &m_pair.qT, "qT/F");
  m_pair_tree->Branch("charge1", &m_pair.charge1, "charge1/F");
  m_pair_tree->Branch("charge2", &m_pair.charge2, "charge2/F");
  m_pair_tree->Branch("cosThetaReco", &m_pair.cosThetaReco, "cosThetaReco/F");
  m_pair_tree->Branch("Lproj", &m_pair.Lproj, "Lproj/F");
  m_pair_tree->Branch("pca_x", &m_pair.pca_x, "pca_x/F");
  m_pair_tree->Branch("pca_y", &m_pair.pca_y, "pca_y/F");
  m_pair_tree->Branch("pca_z", &m_pair.pca_z, "pca_z/F");
  m_pair_tree->Branch("pca1_x", &m_pair.pca1_x, "pca1_x/F");
  m_pair_tree->Branch("pca1_y", &m_pair.pca1_y, "pca1_y/F");
  m_pair_tree->Branch("pca1_z", &m_pair.pca1_z, "pca1_z/F");
  m_pair_tree->Branch("pca2_x", &m_pair.pca2_x, "pca2_x/F");
  m_pair_tree->Branch("pca2_y", &m_pair.pca2_y, "pca2_y/F");
  m_pair_tree->Branch("pca2_z", &m_pair.pca2_z, "pca2_z/F");
  m_pair_tree->Branch("v0_px", &m_pair.v0_px, "v0_px/F");
  m_pair_tree->Branch("v0_py", &m_pair.v0_py, "v0_py/F");
  m_pair_tree->Branch("v0_pz", &m_pair.v0_pz, "v0_pz/F");
  m_pair_tree->Branch("v0_pt", &m_pair.v0_pt, "v0_pt/F");
  m_pair_tree->Branch("mass_Kshort", &m_pair.mass_Kshort, "mass_Kshort/F");
  m_pair_tree->Branch("mass_Lambda", &m_pair.mass_Lambda, "mass_Lambda/F");
  m_pair_tree->Branch("mass_AntiLambda", &m_pair.mass_AntiLambda, "mass_AntiLambda/F");
  m_pair_tree->Branch("true_decay_x", &m_pair.true_decay_x, "true_decay_x/F");
  m_pair_tree->Branch("true_decay_y", &m_pair.true_decay_y, "true_decay_y/F");
  m_pair_tree->Branch("true_decay_z", &m_pair.true_decay_z, "true_decay_z/F");
  m_pair_tree->Branch("pca_to_true_3d", &m_pair.pca_to_true_3d, "pca_to_true_3d/F");
  m_pair_tree->Branch("pca_to_true_xy", &m_pair.pca_to_true_xy, "pca_to_true_xy/F");
  m_pair_tree->Branch("pca_to_true_z", &m_pair.pca_to_true_z, "pca_to_true_z/F");
  m_pair_tree->Branch("truth_alpha", &m_pair.truth_alpha, "truth_alpha/F");
  m_pair_tree->Branch("truth_qT", &m_pair.truth_qT, "truth_qT/F");
  m_pair_tree->Branch("delta_alpha", &m_pair.delta_alpha, "delta_alpha/F");
  m_pair_tree->Branch("delta_qT", &m_pair.delta_qT, "delta_qT/F");
  m_pair_tree->Branch("truth_px1", &m_pair.truth_px1, "truth_px1/F");
  m_pair_tree->Branch("truth_py1", &m_pair.truth_py1, "truth_py1/F");
  m_pair_tree->Branch("truth_pz1", &m_pair.truth_pz1, "truth_pz1/F");
  m_pair_tree->Branch("truth_px2", &m_pair.truth_px2, "truth_px2/F");
  m_pair_tree->Branch("truth_py2", &m_pair.truth_py2, "truth_py2/F");
  m_pair_tree->Branch("truth_pz2", &m_pair.truth_pz2, "truth_pz2/F");
  m_pair_tree->Branch("cos_mom1_truth", &m_pair.cos_mom1_truth, "cos_mom1_truth/F");
  m_pair_tree->Branch("cos_mom2_truth", &m_pair.cos_mom2_truth, "cos_mom2_truth/F");
  m_pair_tree->Branch("pca_theta1", &m_pair.pca_theta1, "pca_theta1/F");
  m_pair_tree->Branch("pca_theta2", &m_pair.pca_theta2, "pca_theta2/F");
  m_pair_tree->Branch("track_id1", &m_pair.track_id1, "track_id1/I");
  m_pair_tree->Branch("track_id2", &m_pair.track_id2, "track_id2/I");
  m_pair_tree->Branch("pid1", &m_pair.pid1, "pid1/I");
  m_pair_tree->Branch("pid2", &m_pair.pid2, "pid2/I");
  m_pair_tree->Branch("parent_id1", &m_pair.parent_id1, "parent_id1/I");
  m_pair_tree->Branch("parent_id2", &m_pair.parent_id2, "parent_id2/I");
  m_pair_tree->Branch("parent_pid", &m_pair.parent_pid, "parent_pid/I");
  m_pair_tree->Branch("npoints1", &m_pair.npoints1, "npoints1/S");
  m_pair_tree->Branch("npoints2", &m_pair.npoints2, "npoints2/S");
}

int TpcV0CandidateTree::pdg_charge(const int pid)
{
  const int apid = std::abs(pid);
  int charge = 0;
  switch (apid)
  {
  case 11:
  case 13:
    charge = -1;
    break;
  case 211:
  case 321:
  case 2212:
  case 3222:
    charge = 1;
    break;
  case 3112:
  case 3312:
  case 3334:
    charge = -1;
    break;
  default:
    charge = 0;
    break;
  }
  return (pid < 0) ? -charge : charge;
}

float TpcV0CandidateTree::quiet_nan()
{
  return std::numeric_limits<float>::quiet_NaN();
}

bool TpcV0CandidateTree::finite(const Vec3 &value)
{
  return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::add(const Vec3 &lhs, const Vec3 &rhs)
{
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::subtract(const Vec3 &lhs, const Vec3 &rhs)
{
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::scale(const Vec3 &value, const double factor)
{
  return {value.x * factor, value.y * factor, value.z * factor};
}

double TpcV0CandidateTree::dot(const Vec3 &lhs, const Vec3 &rhs)
{
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::cross(const Vec3 &lhs, const Vec3 &rhs)
{
  return {
      lhs.y * rhs.z - lhs.z * rhs.y,
      lhs.z * rhs.x - lhs.x * rhs.z,
      lhs.x * rhs.y - lhs.y * rhs.x};
}

double TpcV0CandidateTree::norm(const Vec3 &value)
{
  return std::sqrt(dot(value, value));
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::unit(const Vec3 &value)
{
  const double length = norm(value);
  if (length <= 0.0 || !std::isfinite(length))
  {
    return {quiet_nan(), quiet_nan(), quiet_nan()};
  }
  return scale(value, 1.0 / length);
}

double TpcV0CandidateTree::pt(const Vec3 &value)
{
  return std::sqrt(square(value.x) + square(value.y));
}

double TpcV0CandidateTree::distance(const Vec3 &lhs, const Vec3 &rhs)
{
  return norm(subtract(lhs, rhs));
}

double TpcV0CandidateTree::vector_cosine(const Vec3 &lhs, const Vec3 &rhs)
{
  const double denom = norm(lhs) * norm(rhs);
  if (denom <= 0.0 || !std::isfinite(denom))
  {
    return quiet_nan();
  }
  return dot(lhs, rhs) / denom;
}

bool TpcV0CandidateTree::fit_helix(const std::vector<TruthPoint> &points,
                                   const int fit_first_points,
                                   const double bfield_t,
                                   HelixFit &helix)
{
  const std::size_t nfit =
      (fit_first_points > 0) ? std::min(points.size(), static_cast<std::size_t>(fit_first_points)) : points.size();
  if (nfit < 3)
  {
    return false;
  }

  Eigen::MatrixXd matrix(nfit, 3);
  Eigen::VectorXd rhs(nfit);
  for (std::size_t i = 0; i < nfit; ++i)
  {
    const auto &pos = points[i].position;
    matrix(static_cast<int>(i), 0) = pos.x;
    matrix(static_cast<int>(i), 1) = pos.y;
    matrix(static_cast<int>(i), 2) = 1.0;
    rhs(static_cast<int>(i)) = -(square(pos.x) + square(pos.y));
  }

  const Eigen::Vector3d solution = matrix.colPivHouseholderQr().solve(rhs);
  const double cx = -0.5 * solution(0);
  const double cy = -0.5 * solution(1);
  const double radius2 = square(cx) + square(cy) - solution(2);
  if (radius2 <= 0.0 || !std::isfinite(radius2))
  {
    return false;
  }
  const double radius = std::sqrt(radius2);
  if (radius <= 0.0 || !std::isfinite(radius))
  {
    return false;
  }

  Eigen::MatrixXd z_matrix(nfit, 2);
  Eigen::VectorXd z_rhs(nfit);
  std::vector<double> theta_values;
  theta_values.reserve(nfit);
  for (std::size_t i = 0; i < nfit; ++i)
  {
    double theta = std::atan2(points[i].position.y - cy, points[i].position.x - cx);
    if (!theta_values.empty())
    {
      theta = unwrap_to_previous(theta, theta_values.back());
    }
    theta_values.push_back(theta);
    z_matrix(static_cast<int>(i), 0) = theta;
    z_matrix(static_cast<int>(i), 1) = 1.0;
    z_rhs(static_cast<int>(i)) = points[i].position.z;
  }

  const auto minmax_theta = std::minmax_element(theta_values.begin(), theta_values.end());
  if (minmax_theta.first == theta_values.end() ||
      *minmax_theta.second - *minmax_theta.first < 1e-4)
  {
    return false;
  }

  const Eigen::Vector2d z_solution = z_matrix.colPivHouseholderQr().solve(z_rhs);
  double direction = (theta_values.back() > theta_values.front()) ? 1.0 : -1.0;
  if (std::abs(theta_values.back() - theta_values.front()) <= 0.0)
  {
    direction = 1.0;
  }

  helix.cx = cx;
  helix.cy = cy;
  helix.radius = radius;
  helix.z0 = z_solution(1);
  helix.pitch = z_solution(0);
  helix.theta_first = theta_values.front();
  helix.theta_last = theta_values.back();
  helix.theta_min = *minmax_theta.first;
  helix.theta_max = *minmax_theta.second;
  helix.direction = direction;
  helix.bfield_t = bfield_t;
  return true;
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::helix_point(const HelixFit &helix, const double theta)
{
  return {
      helix.cx + helix.radius * std::cos(theta),
      helix.cy + helix.radius * std::sin(theta),
      helix.z0 + helix.pitch * theta};
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::helix_tangent(const HelixFit &helix, const double theta)
{
  return {
      -helix.radius * std::sin(theta),
      helix.radius * std::cos(theta),
      helix.pitch};
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::helix_momentum(const HelixFit &helix, const double theta)
{
  double p_t = 0.3 * std::abs(helix.bfield_t) * (helix.radius / 100.0);
  if (p_t <= 0.0 || !std::isfinite(p_t))
  {
    p_t = 1.0;
  }

  return {
      helix.direction * p_t * (-std::sin(theta)),
      helix.direction * p_t * std::cos(theta),
      helix.direction * p_t * helix.pitch / helix.radius};
}

std::pair<double, double> TpcV0CandidateTree::theta_search_range(const HelixFit &helix,
                                                                 const double theta_extension,
                                                                 const double downstream_margin)
{
  const double upstream = helix.theta_first - helix.direction * theta_extension;
  const double downstream = helix.theta_first + helix.direction * downstream_margin;
  return {std::min(upstream, downstream), std::max(upstream, downstream)};
}

bool TpcV0CandidateTree::line_line_pca(const Vec3 &pos1, const Vec3 &dir1,
                                       const Vec3 &pos2, const Vec3 &dir2,
                                       LinePca &pca, const bool normalize_dirs)
{
  Vec3 u1 = normalize_dirs ? unit(dir1) : dir1;
  Vec3 u2 = normalize_dirs ? unit(dir2) : dir2;
  if (!finite(u1) || !finite(u2))
  {
    return false;
  }

  const Vec3 w0 = subtract(pos1, pos2);
  const double a = dot(u1, u1);
  const double b = dot(u1, u2);
  const double c = dot(u2, u2);
  const double d = dot(u1, w0);
  const double e = dot(u2, w0);
  const double denom = a * c - b * b;
  if (std::abs(denom) < 1e-12)
  {
    return false;
  }

  const double s = (b * e - c * d) / denom;
  const double t = (a * e - b * d) / denom;
  pca.pca1 = add(pos1, scale(u1, s));
  pca.pca2 = add(pos2, scale(u2, t));
  pca.dca = distance(pca.pca1, pca.pca2);
  pca.step1 = s;
  pca.step2 = t;
  return true;
}

TpcV0CandidateTree::HelixPca TpcV0CandidateTree::refine_helix_pair(
    const HelixFit &helix1, const HelixFit &helix2,
    double theta1, double theta2,
    const double min1, const double max1,
    const double min2, const double max2,
    double max_step)
{
  double best_dca2 = square(distance(helix_point(helix1, theta1), helix_point(helix2, theta2)));

  for (int iter = 0; iter < 30; ++iter)
  {
    LinePca line_pca;
    if (!line_line_pca(helix_point(helix1, theta1), helix_tangent(helix1, theta1),
                       helix_point(helix2, theta2), helix_tangent(helix2, theta2),
                       line_pca, false))
    {
      break;
    }

    const double step1 = std::clamp(line_pca.step1, -max_step, max_step);
    const double step2 = std::clamp(line_pca.step2, -max_step, max_step);
    if (std::abs(step1) < 1e-5 && std::abs(step2) < 1e-5)
    {
      break;
    }

    const double candidate_theta1 = std::clamp(theta1 + step1, min1, max1);
    const double candidate_theta2 = std::clamp(theta2 + step2, min2, max2);
    const double candidate_dca2 = square(distance(helix_point(helix1, candidate_theta1),
                                                  helix_point(helix2, candidate_theta2)));
    if (candidate_dca2 < best_dca2)
    {
      theta1 = candidate_theta1;
      theta2 = candidate_theta2;
      best_dca2 = candidate_dca2;
    }
    else
    {
      max_step *= 0.5;
      if (max_step < 1e-4)
      {
        break;
      }
    }
  }

  HelixPca output;
  output.theta1 = theta1;
  output.theta2 = theta2;
  output.pca1 = helix_point(helix1, theta1);
  output.pca2 = helix_point(helix2, theta2);
  output.dca = distance(output.pca1, output.pca2);
  return output;
}

std::vector<TpcV0CandidateTree::HelixPca> TpcV0CandidateTree::helix_helix_pca_candidates(
    const HelixFit &helix1, const HelixFit &helix2,
    const double theta_extension,
    const int coarse_steps,
    const double downstream_margin,
    const int max_candidates)
{
  const auto range1 = theta_search_range(helix1, theta_extension, downstream_margin);
  const auto range2 = theta_search_range(helix2, theta_extension, downstream_margin);
  const int n_steps = std::max(8, coarse_steps);

  std::vector<double> theta1_values;
  std::vector<double> theta2_values;
  theta1_values.reserve(n_steps);
  theta2_values.reserve(n_steps);
  for (int i = 0; i < n_steps; ++i)
  {
    const double fraction = (n_steps == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n_steps - 1);
    theta1_values.push_back(range1.first + fraction * (range1.second - range1.first));
    theta2_values.push_back(range2.first + fraction * (range2.second - range2.first));
  }

  std::vector<std::tuple<double, int, int>> coarse;
  coarse.reserve(static_cast<std::size_t>(n_steps) * static_cast<std::size_t>(n_steps));
  for (int i = 0; i < n_steps; ++i)
  {
    const Vec3 point1 = helix_point(helix1, theta1_values[i]);
    for (int j = 0; j < n_steps; ++j)
    {
      const Vec3 point2 = helix_point(helix2, theta2_values[j]);
      coarse.emplace_back(square(distance(point1, point2)), i, j);
    }
  }

  std::sort(coarse.begin(), coarse.end(),
            [](const auto &lhs, const auto &rhs)
            { return std::get<0>(lhs) < std::get<0>(rhs); });

  const double max_step = std::max(range1.second - range1.first, range2.second - range2.first) /
                          static_cast<double>(n_steps);
  const int n_candidates = std::min({std::max(1, max_candidates), static_cast<int>(coarse.size())});

  std::vector<HelixPca> candidates;
  candidates.reserve(n_candidates);
  for (int index = 0; index < n_candidates; ++index)
  {
    const int i = std::get<1>(coarse[index]);
    const int j = std::get<2>(coarse[index]);
    candidates.push_back(refine_helix_pair(
        helix1, helix2, theta1_values[i], theta2_values[j],
        range1.first, range1.second, range2.first, range2.second, max_step));
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const HelixPca &lhs, const HelixPca &rhs)
            { return lhs.dca < rhs.dca; });
  return candidates;
}

std::pair<double, double> TpcV0CandidateTree::track_dca_to_vertex(const Vec3 &pos,
                                                                  const Vec3 &mom,
                                                                  const Vec3 &vertex)
{
  const Vec3 rel = subtract(pos, vertex);
  const double pt2 = square(mom.x) + square(mom.y);
  if (pt2 <= 0.0)
  {
    return {quiet_nan(), quiet_nan()};
  }
  const double dca_xy = std::abs(rel.x * mom.y - rel.y * mom.x) / std::sqrt(pt2);
  const double sxy = -(rel.x * mom.x + rel.y * mom.y) / pt2;
  const Vec3 closest = add(pos, scale(mom, sxy));
  return {dca_xy, std::abs(closest.z - vertex.z)};
}

std::pair<double, double> TpcV0CandidateTree::helix_dca_to_vertex(const HelixFit &helix,
                                                                  const Vec3 &vertex)
{
  const double vx = vertex.x - helix.cx;
  const double vy = vertex.y - helix.cy;
  const double distance_to_center = std::sqrt(square(vx) + square(vy));

  double theta_raw = helix.theta_first;
  double dca_xy = helix.radius;
  if (distance_to_center > 0.0)
  {
    theta_raw = std::atan2(vy, vx);
    dca_xy = std::abs(distance_to_center - helix.radius);
  }

  const double wraps = std::round((helix.theta_first - theta_raw) / (2.0 * kPi));
  const double theta = theta_raw + 2.0 * kPi * wraps;
  const Vec3 closest = helix_point(helix, theta);
  return {dca_xy, std::abs(closest.z - vertex.z)};
}

bool TpcV0CandidateTree::armenteros(const Vec3 &pplus, const Vec3 &pminus,
                                    double &alpha, double &qt)
{
  const Vec3 v0p = add(pplus, pminus);
  const Vec3 direction = unit(v0p);
  if (!finite(direction))
  {
    return false;
  }

  const double pl_plus = dot(pplus, direction);
  const double pl_minus = dot(pminus, direction);
  const double denom = pl_plus + pl_minus;
  if (std::abs(denom) < 1e-10)
  {
    return false;
  }

  alpha = (pl_plus - pl_minus) / denom;
  qt = norm(subtract(pplus, scale(direction, pl_plus)));
  return true;
}

double TpcV0CandidateTree::invariant_mass(const Vec3 &mom1, const double mass1,
                                          const Vec3 &mom2, const double mass2)
{
  const double e1 = std::sqrt(dot(mom1, mom1) + square(mass1));
  const double e2 = std::sqrt(dot(mom2, mom2) + square(mass2));
  const Vec3 total_mom = add(mom1, mom2);
  const double mass2_total = square(e1 + e2) - dot(total_mom, total_mom);
  return (mass2_total > 0.0) ? std::sqrt(mass2_total) : 0.0;
}

bool TpcV0CandidateTree::passes_preselection(const Tracklet &track1, const Tracklet &track2,
                                             const Vec3 &primary_vertex) const
{
  if (m_pre_track_pt_min > 0.0 &&
      (pt(track1.momentum) < m_pre_track_pt_min || pt(track2.momentum) < m_pre_track_pt_min))
  {
    return false;
  }

  const auto dca1 = (track1.has_helix) ? helix_dca_to_vertex(track1.helix, primary_vertex)
                                       : track_dca_to_vertex(track1.position, track1.momentum, primary_vertex);
  const auto dca2 = (track2.has_helix) ? helix_dca_to_vertex(track2.helix, primary_vertex)
                                       : track_dca_to_vertex(track2.position, track2.momentum, primary_vertex);
  if (!std::isfinite(dca1.first) || !std::isfinite(dca1.second) ||
      !std::isfinite(dca2.first) || !std::isfinite(dca2.second))
  {
    return false;
  }

  if (m_pre_track_dca_xy_min >= 0.0 &&
      (dca1.first < m_pre_track_dca_xy_min || dca2.first < m_pre_track_dca_xy_min))
  {
    return false;
  }
  if (m_pre_track_dca_z_min >= 0.0 &&
      (dca1.second < m_pre_track_dca_z_min || dca2.second < m_pre_track_dca_z_min))
  {
    return false;
  }
  if (m_pre_track_dca_xy_max >= 0.0 &&
      (dca1.first > m_pre_track_dca_xy_max || dca2.first > m_pre_track_dca_xy_max))
  {
    return false;
  }
  if (m_pre_track_dca_z_max >= 0.0 &&
      (dca1.second > m_pre_track_dca_z_max || dca2.second > m_pre_track_dca_z_max))
  {
    return false;
  }

  if (m_pre_pair_dca_max < 0.0 && m_pre_lproj_min < 0.0 && m_pre_cos_theta_min < -1.0)
  {
    return true;
  }

  LinePca rough_pca;
  if (!line_line_pca(track1.position, track1.momentum, track2.position, track2.momentum, rough_pca, true))
  {
    return false;
  }

  if (m_pre_pair_dca_max >= 0.0 && rough_pca.dca > m_pre_pair_dca_max)
  {
    return false;
  }

  const Vec3 rough_vertex = scale(add(rough_pca.pca1, rough_pca.pca2), 0.5);
  const Vec3 flight = subtract(rough_vertex, primary_vertex);
  const Vec3 total_mom = add(track1.momentum, track2.momentum);
  const double lproj = norm(flight);
  const double cos_theta = vector_cosine(flight, total_mom);

  if (m_pre_lproj_min >= 0.0 && (!std::isfinite(lproj) || lproj < m_pre_lproj_min))
  {
    return false;
  }
  if (m_pre_cos_theta_min >= -1.0 &&
      (!std::isfinite(cos_theta) || cos_theta < m_pre_cos_theta_min))
  {
    return false;
  }

  return true;
}
