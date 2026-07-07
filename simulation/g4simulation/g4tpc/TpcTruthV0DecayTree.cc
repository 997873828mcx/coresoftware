#include "TpcTruthV0DecayTree.h"

#include <g4main/PHG4EventHeader.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
  constexpr double kPionMass = 0.13957039;
  constexpr double kProtonMass = 0.938272088;

  TpcTruthV0DecayTree::Vec3 momentum(const PHG4Particle *particle)
  {
    return {particle->get_px(), particle->get_py(), particle->get_pz()};
  }

  TpcTruthV0DecayTree::Vec3 subtract(const TpcTruthV0DecayTree::Vec3 &lhs,
                                     const TpcTruthV0DecayTree::Vec3 &rhs)
  {
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
  }

  double dot(const TpcTruthV0DecayTree::Vec3 &lhs,
             const TpcTruthV0DecayTree::Vec3 &rhs)
  {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
  }

  double norm(const TpcTruthV0DecayTree::Vec3 &value)
  {
    return std::sqrt(dot(value, value));
  }
}  // namespace

TpcTruthV0DecayTree::TpcTruthV0DecayTree(const std::string &name,
                                         const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

void TpcTruthV0DecayTree::set_primary_vertex(const double x, const double y, const double z)
{
  m_fixed_primary_vertex = {x, y, z};
}

int TpcTruthV0DecayTree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = new TFile(m_filename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": failed to create output file " << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("truthV0Tree", "Generator-level K0S/Lambda V0 decays");
  create_branches();
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTruthV0DecayTree::process_event(PHCompositeNode *topNode)
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

  const int event_number = get_event_number(topNode);
  const Vec3 primary_vertex = get_primary_vertex(truth_info);

  const auto range = truth_info->GetParticleRange();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    const PHG4Particle *particle = iter->second;
    if (!particle || !is_v0_parent(particle->get_pid()))
    {
      continue;
    }

    ++m_counter_parents;
    if (fill_decay_row(truth_info, particle, event_number, primary_vertex))
    {
      m_tree->Fill();
      ++m_counter_written;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTruthV0DecayTree::End(PHCompositeNode * /*topNode*/)
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
    std::cout << Name() << ": truth V0 parents=" << m_counter_parents
              << " written charged two-body decays=" << m_counter_written << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcTruthV0DecayTree::get_event_number(PHCompositeNode *topNode) const
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

TpcTruthV0DecayTree::Vec3 TpcTruthV0DecayTree::get_primary_vertex(PHG4TruthInfoContainer *truth_info) const
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

const PHG4VtxPoint *TpcTruthV0DecayTree::get_decay_vertex(PHG4TruthInfoContainer *truth_info,
                                                          const PHG4Particle *daughter1,
                                                          const PHG4Particle *daughter2) const
{
  if (!truth_info || !daughter1 || !daughter2)
  {
    return nullptr;
  }

  const PHG4VtxPoint *vtx1 = truth_info->GetVtx(daughter1->get_vtx_id());
  const PHG4VtxPoint *vtx2 = truth_info->GetVtx(daughter2->get_vtx_id());
  if (vtx1)
  {
    return vtx1;
  }
  return vtx2;
}

bool TpcTruthV0DecayTree::fill_decay_row(PHG4TruthInfoContainer *truth_info,
                                         const PHG4Particle *parent,
                                         const int event_number,
                                         const Vec3 &primary_vertex)
{
  reset_row();
  m_row.run = 1;
  m_row.evt = event_number;
  m_row.parent_id = parent->get_track_id();
  m_row.parent_pid = parent->get_pid();
  m_row.parent_vtx_id = parent->get_vtx_id();
  m_row.parent_embed_id = truth_info ? truth_info->isEmbeded(parent->get_track_id()) : 0;
  m_row.parent_is_primary = truth_info && truth_info->is_primary(parent) ? 1 : 0;
  m_row.parent_is_sphenix_primary = truth_info && truth_info->is_sPHENIX_primary(parent) ? 1 : 0;

  std::vector<const PHG4Particle *> daughters;
  const auto range = truth_info->GetParticleRange();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    const PHG4Particle *daughter = iter->second;
    if (daughter && daughter->get_parent_id() == parent->get_track_id())
    {
      daughters.push_back(daughter);
    }
  }

  m_row.n_daughters = static_cast<int>(daughters.size());

  const PHG4Particle *daughter1 = nullptr;
  const PHG4Particle *daughter2 = nullptr;
  for (std::size_t i = 0; i < daughters.size(); ++i)
  {
    for (std::size_t j = i + 1; j < daughters.size(); ++j)
    {
      if (is_expected_daughter_pair(parent->get_pid(), daughters[i]->get_pid(), daughters[j]->get_pid()))
      {
        daughter1 = daughters[i];
        daughter2 = daughters[j];
        break;
      }
    }
    if (daughter1 && daughter2)
    {
      break;
    }
  }

  if (!daughter1 || !daughter2)
  {
    return false;
  }

  const Vec3 mom1 = momentum(daughter1);
  const Vec3 mom2 = momentum(daughter2);

  const PHG4Particle *positive = daughter1;
  const PHG4Particle *negative = daughter2;
  if (daughter1->get_pid() < 0 && daughter2->get_pid() > 0)
  {
    positive = daughter2;
    negative = daughter1;
  }
  const Vec3 pplus = momentum(positive);
  const Vec3 pminus = momentum(negative);

  double alpha = std::numeric_limits<double>::quiet_NaN();
  double qt = std::numeric_limits<double>::quiet_NaN();
  if (!armenteros(pplus, pminus, alpha, qt))
  {
    return false;
  }

  m_row.daughter1_id = daughter1->get_track_id();
  m_row.daughter1_pid = daughter1->get_pid();
  m_row.daughter2_id = daughter2->get_track_id();
  m_row.daughter2_pid = daughter2->get_pid();

  m_row.primary_x = static_cast<float>(primary_vertex.x);
  m_row.primary_y = static_cast<float>(primary_vertex.y);
  m_row.primary_z = static_cast<float>(primary_vertex.z);

  if (const PHG4VtxPoint *decay_vtx = get_decay_vertex(truth_info, daughter1, daughter2))
  {
    m_row.decay_vtx_id = decay_vtx->get_id();
    m_row.decay_x = static_cast<float>(decay_vtx->get_x());
    m_row.decay_y = static_cast<float>(decay_vtx->get_y());
    m_row.decay_z = static_cast<float>(decay_vtx->get_z());
    m_row.decay_t = static_cast<float>(decay_vtx->get_t());
    const Vec3 decay_pos{decay_vtx->get_x(), decay_vtx->get_y(), decay_vtx->get_z()};
    const Vec3 flight = subtract(decay_pos, primary_vertex);
    m_row.Lxy = static_cast<float>(std::sqrt(square(flight.x) + square(flight.y)));
    m_row.Lxyz = static_cast<float>(norm(flight));
  }

  const Vec3 parent_mom = momentum(parent);
  m_row.parent_px = static_cast<float>(parent_mom.x);
  m_row.parent_py = static_cast<float>(parent_mom.y);
  m_row.parent_pz = static_cast<float>(parent_mom.z);
  m_row.parent_e = static_cast<float>(parent->get_e());
  m_row.parent_pt = static_cast<float>(pt(parent_mom));
  m_row.parent_eta = static_cast<float>(eta(parent_mom));
  m_row.parent_phi = static_cast<float>(std::atan2(parent_mom.y, parent_mom.x));

  m_row.px1 = static_cast<float>(mom1.x);
  m_row.py1 = static_cast<float>(mom1.y);
  m_row.pz1 = static_cast<float>(mom1.z);
  m_row.e1 = static_cast<float>(daughter1->get_e());
  m_row.pt1 = static_cast<float>(pt(mom1));
  m_row.eta1 = static_cast<float>(eta(mom1));
  m_row.phi1 = static_cast<float>(std::atan2(mom1.y, mom1.x));

  m_row.px2 = static_cast<float>(mom2.x);
  m_row.py2 = static_cast<float>(mom2.y);
  m_row.pz2 = static_cast<float>(mom2.z);
  m_row.e2 = static_cast<float>(daughter2->get_e());
  m_row.pt2 = static_cast<float>(pt(mom2));
  m_row.eta2 = static_cast<float>(eta(mom2));
  m_row.phi2 = static_cast<float>(std::atan2(mom2.y, mom2.x));

  m_row.alpha = static_cast<float>(alpha);
  m_row.qT = static_cast<float>(qt);
  m_row.mass_Kshort = static_cast<float>(invariant_mass(mom1, kPionMass, mom2, kPionMass));

  const Vec3 positive_mom = momentum(positive);
  const Vec3 negative_mom = momentum(negative);
  m_row.mass_Lambda = static_cast<float>(invariant_mass(positive_mom, kProtonMass, negative_mom, kPionMass));
  m_row.mass_AntiLambda = static_cast<float>(invariant_mass(positive_mom, kPionMass, negative_mom, kProtonMass));

  return true;
}

void TpcTruthV0DecayTree::create_branches()
{
  m_tree->Branch("run", &m_row.run, "run/I");
  m_tree->Branch("evt", &m_row.evt, "evt/I");
  m_tree->Branch("parent_id", &m_row.parent_id, "parent_id/I");
  m_tree->Branch("parent_pid", &m_row.parent_pid, "parent_pid/I");
  m_tree->Branch("parent_vtx_id", &m_row.parent_vtx_id, "parent_vtx_id/I");
  m_tree->Branch("decay_vtx_id", &m_row.decay_vtx_id, "decay_vtx_id/I");
  m_tree->Branch("parent_embed_id", &m_row.parent_embed_id, "parent_embed_id/I");
  m_tree->Branch("parent_is_primary", &m_row.parent_is_primary, "parent_is_primary/I");
  m_tree->Branch("parent_is_sphenix_primary", &m_row.parent_is_sphenix_primary, "parent_is_sphenix_primary/I");
  m_tree->Branch("daughter1_id", &m_row.daughter1_id, "daughter1_id/I");
  m_tree->Branch("daughter1_pid", &m_row.daughter1_pid, "daughter1_pid/I");
  m_tree->Branch("daughter2_id", &m_row.daughter2_id, "daughter2_id/I");
  m_tree->Branch("daughter2_pid", &m_row.daughter2_pid, "daughter2_pid/I");
  m_tree->Branch("n_daughters", &m_row.n_daughters, "n_daughters/I");

  m_tree->Branch("primary_x", &m_row.primary_x, "primary_x/F");
  m_tree->Branch("primary_y", &m_row.primary_y, "primary_y/F");
  m_tree->Branch("primary_z", &m_row.primary_z, "primary_z/F");
  m_tree->Branch("decay_x", &m_row.decay_x, "decay_x/F");
  m_tree->Branch("decay_y", &m_row.decay_y, "decay_y/F");
  m_tree->Branch("decay_z", &m_row.decay_z, "decay_z/F");
  m_tree->Branch("decay_t", &m_row.decay_t, "decay_t/F");
  m_tree->Branch("Lxy", &m_row.Lxy, "Lxy/F");
  m_tree->Branch("Lxyz", &m_row.Lxyz, "Lxyz/F");

  m_tree->Branch("parent_px", &m_row.parent_px, "parent_px/F");
  m_tree->Branch("parent_py", &m_row.parent_py, "parent_py/F");
  m_tree->Branch("parent_pz", &m_row.parent_pz, "parent_pz/F");
  m_tree->Branch("parent_e", &m_row.parent_e, "parent_e/F");
  m_tree->Branch("parent_pt", &m_row.parent_pt, "parent_pt/F");
  m_tree->Branch("parent_eta", &m_row.parent_eta, "parent_eta/F");
  m_tree->Branch("parent_phi", &m_row.parent_phi, "parent_phi/F");

  m_tree->Branch("px1", &m_row.px1, "px1/F");
  m_tree->Branch("py1", &m_row.py1, "py1/F");
  m_tree->Branch("pz1", &m_row.pz1, "pz1/F");
  m_tree->Branch("e1", &m_row.e1, "e1/F");
  m_tree->Branch("pt1", &m_row.pt1, "pt1/F");
  m_tree->Branch("eta1", &m_row.eta1, "eta1/F");
  m_tree->Branch("phi1", &m_row.phi1, "phi1/F");

  m_tree->Branch("px2", &m_row.px2, "px2/F");
  m_tree->Branch("py2", &m_row.py2, "py2/F");
  m_tree->Branch("pz2", &m_row.pz2, "pz2/F");
  m_tree->Branch("e2", &m_row.e2, "e2/F");
  m_tree->Branch("pt2", &m_row.pt2, "pt2/F");
  m_tree->Branch("eta2", &m_row.eta2, "eta2/F");
  m_tree->Branch("phi2", &m_row.phi2, "phi2/F");

  m_tree->Branch("alpha", &m_row.alpha, "alpha/F");
  m_tree->Branch("qT", &m_row.qT, "qT/F");
  m_tree->Branch("mass_Kshort", &m_row.mass_Kshort, "mass_Kshort/F");
  m_tree->Branch("mass_Lambda", &m_row.mass_Lambda, "mass_Lambda/F");
  m_tree->Branch("mass_AntiLambda", &m_row.mass_AntiLambda, "mass_AntiLambda/F");
}

void TpcTruthV0DecayTree::reset_row()
{
  m_row = Row();
}

bool TpcTruthV0DecayTree::is_v0_parent(const int pid)
{
  return pid == 310 || std::abs(pid) == 3122;
}

bool TpcTruthV0DecayTree::is_expected_daughter_pair(const int parent_pid, const int pid1, const int pid2)
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

double TpcTruthV0DecayTree::pt(const Vec3 &mom)
{
  return std::sqrt(square(mom.x) + square(mom.y));
}

double TpcTruthV0DecayTree::eta(const Vec3 &mom)
{
  const double pt_value = pt(mom);
  if (pt_value <= 0.0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::asinh(mom.z / pt_value);
}

double TpcTruthV0DecayTree::invariant_mass(const Vec3 &mom1, const double mass1,
                                           const Vec3 &mom2, const double mass2)
{
  const double e1 = std::sqrt(square(mass1) + square(mom1.x) + square(mom1.y) + square(mom1.z));
  const double e2 = std::sqrt(square(mass2) + square(mom2.x) + square(mom2.y) + square(mom2.z));
  const double px = mom1.x + mom2.x;
  const double py = mom1.y + mom2.y;
  const double pz = mom1.z + mom2.z;
  const double m2 = square(e1 + e2) - square(px) - square(py) - square(pz);
  return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
}

bool TpcTruthV0DecayTree::armenteros(const Vec3 &pplus, const Vec3 &pminus,
                                     double &alpha, double &qt)
{
  const Vec3 v0p{pplus.x + pminus.x, pplus.y + pminus.y, pplus.z + pminus.z};
  const double v0_norm = norm(v0p);
  if (v0_norm <= 0.0 || !std::isfinite(v0_norm))
  {
    return false;
  }

  const Vec3 direction{v0p.x / v0_norm, v0p.y / v0_norm, v0p.z / v0_norm};
  const double pl_plus = dot(pplus, direction);
  const double pl_minus = dot(pminus, direction);
  const double denom = pl_plus + pl_minus;
  if (std::abs(denom) <= 0.0)
  {
    return false;
  }

  const Vec3 transverse{
      pplus.x - pl_plus * direction.x,
      pplus.y - pl_plus * direction.y,
      pplus.z - pl_plus * direction.z};

  alpha = (pl_plus - pl_minus) / denom;
  qt = norm(transverse);
  return std::isfinite(alpha) && std::isfinite(qt);
}
