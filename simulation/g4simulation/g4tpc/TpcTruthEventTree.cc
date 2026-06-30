#include "TpcTruthEventTree.h"

#include <g4main/PHG4EventHeader.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>

TpcTruthEventTree::TpcTruthEventTree(const std::string &name,
                                     const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

void TpcTruthEventTree::set_primary_vertex(const double x, const double y, const double z)
{
  m_fixed_primary_vertex = {x, y, z};
}

int TpcTruthEventTree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = new TFile(m_filename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": failed to create output file " << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("truthEventTree", "Generator-level event QA summary");
  create_branches();
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTruthEventTree::process_event(PHCompositeNode *topNode)
{
  auto *truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_truth_info_node);
  if (!truth_info)
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << Name() << ": missing truth info node " << m_truth_info_node << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  reset_row();
  m_row.run = get_run_number(topNode);
  m_row.evt = get_event_number(topNode);
  m_row.charged_eta_max = static_cast<float>(m_charged_eta_max);
  m_row.charged_pt_min = static_cast<float>(m_charged_pt_min);
  m_row.v0_abs_y_max = static_cast<float>(m_v0_abs_y_max);
  m_row.v0_pt_min = static_cast<float>(m_v0_pt_min);

  std::map<int, std::vector<const PHG4Particle *>> daughters_by_parent;
  const auto all_range = truth_info->GetParticleRange();
  for (auto iter = all_range.first; iter != all_range.second; ++iter)
  {
    const PHG4Particle *particle = iter->second;
    if (!particle)
    {
      continue;
    }

    ++m_row.n_particles_all;
    if (particle->get_track_id() > 0)
    {
      ++m_row.n_particles_primary;
    }
    else
    {
      ++m_row.n_particles_secondary;
    }

    daughters_by_parent[particle->get_parent_id()].push_back(particle);
  }

  const auto vtx_range = truth_info->GetPrimaryVtxRange();
  for (auto iter = vtx_range.first; iter != vtx_range.second; ++iter)
  {
    if (iter->second)
    {
      ++m_row.n_primary_vertices;
    }
  }

  if (const PHG4VtxPoint *primary_vertex = get_primary_vertex(truth_info))
  {
    m_row.primary_vtx_id = primary_vertex->get_id();
    m_row.primary_embed_id = truth_info->isEmbededVtx(primary_vertex->get_id());
    m_row.primary_x = static_cast<float>(primary_vertex->get_x());
    m_row.primary_y = static_cast<float>(primary_vertex->get_y());
    m_row.primary_z = static_cast<float>(primary_vertex->get_z());
    m_row.primary_t = static_cast<float>(primary_vertex->get_t());
  }
  else
  {
    m_row.primary_x = static_cast<float>(m_fixed_primary_vertex.x);
    m_row.primary_y = static_cast<float>(m_fixed_primary_vertex.y);
    m_row.primary_z = static_cast<float>(m_fixed_primary_vertex.z);
  }

  const auto primary_range = truth_info->GetPrimaryParticleRange();
  for (auto iter = primary_range.first; iter != primary_range.second; ++iter)
  {
    const PHG4Particle *particle = iter->second;
    if (!particle)
    {
      continue;
    }
    const auto daughter_iter = daughters_by_parent.find(particle->get_track_id());
    const bool has_daughters = daughter_iter != daughters_by_parent.end() && !daughter_iter->second.empty();
    count_primary_particle(particle, has_daughters);
  }

  for (auto iter = all_range.first; iter != all_range.second; ++iter)
  {
    const PHG4Particle *particle = iter->second;
    if (!particle || !is_v0_parent(particle->get_pid()))
    {
      continue;
    }
    const auto daughter_iter = daughters_by_parent.find(particle->get_track_id());
    static const std::vector<const PHG4Particle *> empty;
    const auto &daughters = (daughter_iter != daughters_by_parent.end()) ? daughter_iter->second : empty;
    count_v0_parent(particle, daughters);
  }

  m_tree->Fill();
  ++m_counter_events;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTruthEventTree::End(PHCompositeNode * /*topNode*/)
{
  if (m_file)
  {
    m_file->cd();
    if (m_tree)
    {
      m_tree->Write();
    }
    m_file->Close();
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": wrote event QA rows=" << m_counter_events << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTruthEventTree::get_run_number(PHCompositeNode *topNode) const
{
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    return event_header->get_RunNumber();
  }
  return 1;
}

int TpcTruthEventTree::get_event_number(PHCompositeNode *topNode) const
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

const PHG4VtxPoint *TpcTruthEventTree::get_primary_vertex(PHG4TruthInfoContainer *truth_info) const
{
  if (!m_use_truth_primary_vertex || !truth_info)
  {
    return nullptr;
  }

  const auto vtx_range = truth_info->GetPrimaryVtxRange();
  for (auto iter = vtx_range.first; iter != vtx_range.second; ++iter)
  {
    if (iter->second)
    {
      return iter->second;
    }
  }
  return nullptr;
}

bool TpcTruthEventTree::has_expected_charged_decay(
    const PHG4Particle *parent,
    const std::vector<const PHG4Particle *> &daughters) const
{
  if (!parent)
  {
    return false;
  }

  for (std::size_t i = 0; i < daughters.size(); ++i)
  {
    for (std::size_t j = i + 1; j < daughters.size(); ++j)
    {
      if (daughters[i] && daughters[j] &&
          is_expected_daughter_pair(parent->get_pid(), daughters[i]->get_pid(), daughters[j]->get_pid()))
      {
        return true;
      }
    }
  }
  return false;
}

bool TpcTruthEventTree::passes_v0_fiducial(const PHG4Particle *parent) const
{
  if (!parent)
  {
    return false;
  }

  const Vec3 mom = momentum(parent);
  const double parent_pt = pt(mom);
  const double parent_y = rapidity(parent);
  return std::isfinite(parent_pt) && std::isfinite(parent_y) &&
         parent_pt >= m_v0_pt_min && std::abs(parent_y) < m_v0_abs_y_max;
}

void TpcTruthEventTree::count_primary_particle(const PHG4Particle *particle, const bool has_daughters)
{
  if (!particle)
  {
    return;
  }

  const int pid = particle->get_pid();
  if (pid == 211)
  {
    ++m_row.n_pi_plus_primary;
  }
  else if (pid == -211)
  {
    ++m_row.n_pi_minus_primary;
  }
  else if (pid == 321)
  {
    ++m_row.n_k_plus_primary;
  }
  else if (pid == -321)
  {
    ++m_row.n_k_minus_primary;
  }
  else if (pid == 2212)
  {
    ++m_row.n_p_primary;
  }
  else if (pid == -2212)
  {
    ++m_row.n_pbar_primary;
  }

  if (!is_charged(pid))
  {
    return;
  }

  const Vec3 mom = momentum(particle);
  const double particle_pt = pt(mom);
  const double particle_eta = eta(mom);
  const bool in_eta = std::isfinite(particle_eta) && std::abs(particle_eta) < m_charged_eta_max;
  const bool in_pt = std::isfinite(particle_pt) && particle_pt >= m_charged_pt_min;

  ++m_row.n_charged_primary;
  if (!has_daughters)
  {
    ++m_row.n_charged_primary_no_daughters;
  }
  if (in_eta)
  {
    ++m_row.n_charged_primary_eta;
    if (!has_daughters)
    {
      ++m_row.n_charged_primary_no_daughters_eta;
    }
  }
  if (in_eta && in_pt)
  {
    ++m_row.n_charged_primary_eta_pt;
    if (!has_daughters)
    {
      ++m_row.n_charged_primary_no_daughters_eta_pt;
    }
  }
}

void TpcTruthEventTree::count_v0_parent(
    const PHG4Particle *parent,
    const std::vector<const PHG4Particle *> &daughters)
{
  if (!parent)
  {
    return;
  }

  const int pid = parent->get_pid();
  const bool charged_decay = has_expected_charged_decay(parent, daughters);
  const bool fiducial = charged_decay && passes_v0_fiducial(parent);
  const bool primary_track = parent->get_track_id() > 0;

  if (pid == 310)
  {
    ++m_row.n_kshort;
    if (primary_track)
    {
      ++m_row.n_kshort_primary;
    }
    if (charged_decay)
    {
      ++m_row.n_kshort_charged_decay;
      if (primary_track)
      {
        ++m_row.n_kshort_primary_charged_decay;
      }
    }
    if (fiducial)
    {
      ++m_row.n_kshort_fiducial;
    }
  }
  else if (pid == 3122)
  {
    ++m_row.n_lambda;
    if (primary_track)
    {
      ++m_row.n_lambda_primary;
    }
    if (charged_decay)
    {
      ++m_row.n_lambda_charged_decay;
      if (primary_track)
      {
        ++m_row.n_lambda_primary_charged_decay;
      }
    }
    if (fiducial)
    {
      ++m_row.n_lambda_fiducial;
    }
  }
  else if (pid == -3122)
  {
    ++m_row.n_antilambda;
    if (primary_track)
    {
      ++m_row.n_antilambda_primary;
    }
    if (charged_decay)
    {
      ++m_row.n_antilambda_charged_decay;
      if (primary_track)
      {
        ++m_row.n_antilambda_primary_charged_decay;
      }
    }
    if (fiducial)
    {
      ++m_row.n_antilambda_fiducial;
    }
  }
}

void TpcTruthEventTree::create_branches()
{
  m_tree->Branch("run", &m_row.run, "run/I");
  m_tree->Branch("evt", &m_row.evt, "evt/I");

  m_tree->Branch("primary_vtx_id", &m_row.primary_vtx_id, "primary_vtx_id/I");
  m_tree->Branch("primary_embed_id", &m_row.primary_embed_id, "primary_embed_id/I");
  m_tree->Branch("n_primary_vertices", &m_row.n_primary_vertices, "n_primary_vertices/I");
  m_tree->Branch("primary_x", &m_row.primary_x, "primary_x/F");
  m_tree->Branch("primary_y", &m_row.primary_y, "primary_y/F");
  m_tree->Branch("primary_z", &m_row.primary_z, "primary_z/F");
  m_tree->Branch("primary_t", &m_row.primary_t, "primary_t/F");

  m_tree->Branch("n_particles_all", &m_row.n_particles_all, "n_particles_all/I");
  m_tree->Branch("n_particles_primary", &m_row.n_particles_primary, "n_particles_primary/I");
  m_tree->Branch("n_particles_secondary", &m_row.n_particles_secondary, "n_particles_secondary/I");

  m_tree->Branch("n_charged_primary", &m_row.n_charged_primary, "n_charged_primary/I");
  m_tree->Branch("n_charged_primary_no_daughters", &m_row.n_charged_primary_no_daughters, "n_charged_primary_no_daughters/I");
  m_tree->Branch("n_charged_primary_eta", &m_row.n_charged_primary_eta, "n_charged_primary_eta/I");
  m_tree->Branch("n_charged_primary_eta_pt", &m_row.n_charged_primary_eta_pt, "n_charged_primary_eta_pt/I");
  m_tree->Branch("n_charged_primary_no_daughters_eta", &m_row.n_charged_primary_no_daughters_eta, "n_charged_primary_no_daughters_eta/I");
  m_tree->Branch("n_charged_primary_no_daughters_eta_pt", &m_row.n_charged_primary_no_daughters_eta_pt, "n_charged_primary_no_daughters_eta_pt/I");

  m_tree->Branch("n_pi_plus_primary", &m_row.n_pi_plus_primary, "n_pi_plus_primary/I");
  m_tree->Branch("n_pi_minus_primary", &m_row.n_pi_minus_primary, "n_pi_minus_primary/I");
  m_tree->Branch("n_k_plus_primary", &m_row.n_k_plus_primary, "n_k_plus_primary/I");
  m_tree->Branch("n_k_minus_primary", &m_row.n_k_minus_primary, "n_k_minus_primary/I");
  m_tree->Branch("n_p_primary", &m_row.n_p_primary, "n_p_primary/I");
  m_tree->Branch("n_pbar_primary", &m_row.n_pbar_primary, "n_pbar_primary/I");

  m_tree->Branch("n_kshort", &m_row.n_kshort, "n_kshort/I");
  m_tree->Branch("n_lambda", &m_row.n_lambda, "n_lambda/I");
  m_tree->Branch("n_antilambda", &m_row.n_antilambda, "n_antilambda/I");
  m_tree->Branch("n_kshort_primary", &m_row.n_kshort_primary, "n_kshort_primary/I");
  m_tree->Branch("n_lambda_primary", &m_row.n_lambda_primary, "n_lambda_primary/I");
  m_tree->Branch("n_antilambda_primary", &m_row.n_antilambda_primary, "n_antilambda_primary/I");
  m_tree->Branch("n_kshort_charged_decay", &m_row.n_kshort_charged_decay, "n_kshort_charged_decay/I");
  m_tree->Branch("n_lambda_charged_decay", &m_row.n_lambda_charged_decay, "n_lambda_charged_decay/I");
  m_tree->Branch("n_antilambda_charged_decay", &m_row.n_antilambda_charged_decay, "n_antilambda_charged_decay/I");
  m_tree->Branch("n_kshort_primary_charged_decay", &m_row.n_kshort_primary_charged_decay, "n_kshort_primary_charged_decay/I");
  m_tree->Branch("n_lambda_primary_charged_decay", &m_row.n_lambda_primary_charged_decay, "n_lambda_primary_charged_decay/I");
  m_tree->Branch("n_antilambda_primary_charged_decay", &m_row.n_antilambda_primary_charged_decay, "n_antilambda_primary_charged_decay/I");
  m_tree->Branch("n_kshort_fiducial", &m_row.n_kshort_fiducial, "n_kshort_fiducial/I");
  m_tree->Branch("n_lambda_fiducial", &m_row.n_lambda_fiducial, "n_lambda_fiducial/I");
  m_tree->Branch("n_antilambda_fiducial", &m_row.n_antilambda_fiducial, "n_antilambda_fiducial/I");

  m_tree->Branch("charged_eta_max", &m_row.charged_eta_max, "charged_eta_max/F");
  m_tree->Branch("charged_pt_min", &m_row.charged_pt_min, "charged_pt_min/F");
  m_tree->Branch("v0_abs_y_max", &m_row.v0_abs_y_max, "v0_abs_y_max/F");
  m_tree->Branch("v0_pt_min", &m_row.v0_pt_min, "v0_pt_min/F");
}

void TpcTruthEventTree::reset_row()
{
  m_row = Row();
}

bool TpcTruthEventTree::is_charged(const int pid)
{
  const auto *particle = TDatabasePDG::Instance()->GetParticle(pid);
  return particle && std::abs(particle->Charge()) > 0.5;
}

bool TpcTruthEventTree::is_v0_parent(const int pid)
{
  return pid == 310 || std::abs(pid) == 3122;
}

bool TpcTruthEventTree::is_expected_daughter_pair(const int parent_pid, const int pid1, const int pid2)
{
  const int lo = std::min(pid1, pid2);
  const int hi = std::max(pid1, pid2);
  if (parent_pid == 310)
  {
    return lo == -211 && hi == 211;
  }
  if (parent_pid == 3122)
  {
    return (pid1 == 2212 && pid2 == -211) || (pid2 == 2212 && pid1 == -211);
  }
  if (parent_pid == -3122)
  {
    return (pid1 == -2212 && pid2 == 211) || (pid2 == -2212 && pid1 == 211);
  }
  return false;
}

double TpcTruthEventTree::pt(const Vec3 &mom)
{
  return std::sqrt(square(mom.x) + square(mom.y));
}

double TpcTruthEventTree::eta(const Vec3 &mom)
{
  const double pt_value = pt(mom);
  if (pt_value <= 0.0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::asinh(mom.z / pt_value);
}

double TpcTruthEventTree::rapidity(const PHG4Particle *particle)
{
  if (!particle)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double e = particle->get_e();
  const double pz = particle->get_pz();
  if (!std::isfinite(e) || !std::isfinite(pz) || e <= std::abs(pz))
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return 0.5 * std::log((e + pz) / (e - pz));
}

TpcTruthEventTree::Vec3 TpcTruthEventTree::momentum(const PHG4Particle *particle)
{
  if (!particle)
  {
    return {};
  }
  return {particle->get_px(), particle->get_py(), particle->get_pz()};
}
