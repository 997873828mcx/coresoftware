#include "PHG4TpcTruthPointTree.h"

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

#include <cmath>
#include <iostream>
#include <limits>

namespace
{
  template <class T>
  constexpr T square(const T &x)
  {
    return x * x;
  }

  float quiet_nan()
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float safe_eta(const double px, const double py, const double pz)
  {
    const double p = std::sqrt(square(px) + square(py) + square(pz));
    if (p <= std::abs(pz))
    {
      return quiet_nan();
    }
    return static_cast<float>(0.5 * std::log((p + pz) / (p - pz)));
  }
}  // namespace

PHG4TpcTruthPointTree::PHG4TpcTruthPointTree(const std::string &name,
                                             const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

int PHG4TpcTruthPointTree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = new TFile(m_filename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": failed to create output file " << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event_tree = new TTree("events", "event summary for TPC truth points");
  m_event_tree->Branch("event", &m_evt_event, "event/I");
  m_event_tree->Branch("event_index", &m_evt_event_index, "event_index/I");
  m_event_tree->Branch("n_truth_points", &m_evt_n_truth_points, "n_truth_points/I");
  m_event_tree->Branch("n_truth_particles", &m_evt_n_truth_particles, "n_truth_particles/I");
  m_event_tree->Branch("n_primary_particles", &m_evt_n_primary_particles, "n_primary_particles/I");

  m_point_tree = new TTree("tpc_truth_points", "ideal TPC layer truth points");
  m_point_tree->Branch("event", &m_point_event, "event/I");
  m_point_tree->Branch("event_index", &m_point_event_index, "event_index/I");
  m_point_tree->Branch("hit_key", &m_point_hit_key, "hit_key/L");
  m_point_tree->Branch("track_id", &m_point_track_id, "track_id/I");
  m_point_tree->Branch("shower_id", &m_point_shower_id, "shower_id/I");
  m_point_tree->Branch("layer", &m_point_layer, "layer/I");
  m_point_tree->Branch("side", &m_point_side, "side/I");
  m_point_tree->Branch("pid", &m_point_pid, "pid/I");
  m_point_tree->Branch("parent_id", &m_point_parent_id, "parent_id/I");
  m_point_tree->Branch("primary_id", &m_point_primary_id, "primary_id/I");
  m_point_tree->Branch("vtx_id", &m_point_vtx_id, "vtx_id/I");
  m_point_tree->Branch("barcode", &m_point_barcode, "barcode/I");
  m_point_tree->Branch("embed_id", &m_point_embed_id, "embed_id/I");
  m_point_tree->Branch("is_primary", &m_point_is_primary, "is_primary/I");
  m_point_tree->Branch("x", &m_point_x, "x/F");
  m_point_tree->Branch("y", &m_point_y, "y/F");
  m_point_tree->Branch("z", &m_point_z, "z/F");
  m_point_tree->Branch("t", &m_point_t, "t/F");
  m_point_tree->Branch("r", &m_point_r, "r/F");
  m_point_tree->Branch("phi", &m_point_phi, "phi/F");
  m_point_tree->Branch("path", &m_point_path, "path/F");
  m_point_tree->Branch("px", &m_point_px, "px/F");
  m_point_tree->Branch("py", &m_point_py, "py/F");
  m_point_tree->Branch("pz", &m_point_pz, "pz/F");
  m_point_tree->Branch("truth_px", &m_point_truth_px, "truth_px/F");
  m_point_tree->Branch("truth_py", &m_point_truth_py, "truth_py/F");
  m_point_tree->Branch("truth_pz", &m_point_truth_pz, "truth_pz/F");
  m_point_tree->Branch("truth_e", &m_point_truth_e, "truth_e/F");
  m_point_tree->Branch("vx", &m_point_vx, "vx/F");
  m_point_tree->Branch("vy", &m_point_vy, "vy/F");
  m_point_tree->Branch("vz", &m_point_vz, "vz/F");
  m_point_tree->Branch("vt", &m_point_vt, "vt/F");

  m_cluster_tree = new TTree("combined_clusters", "event-display compatible TPC truth points");
  add_display_branches(m_cluster_tree);

  m_hit_tree = new TTree("combined_hits", "event-display compatible TPC truth points");
  add_display_branches(m_hit_tree);

  if (m_write_particle_tree)
  {
    m_particle_tree = new TTree("truth_particles", "G4 truth particle summary");
    m_particle_tree->Branch("event", &m_particle_event, "event/I");
    m_particle_tree->Branch("event_index", &m_particle_event_index, "event_index/I");
    m_particle_tree->Branch("track_id", &m_particle_track_id, "track_id/I");
    m_particle_tree->Branch("pid", &m_particle_pid, "pid/I");
    m_particle_tree->Branch("parent_id", &m_particle_parent_id, "parent_id/I");
    m_particle_tree->Branch("primary_id", &m_particle_primary_id, "primary_id/I");
    m_particle_tree->Branch("vtx_id", &m_particle_vtx_id, "vtx_id/I");
    m_particle_tree->Branch("barcode", &m_particle_barcode, "barcode/I");
    m_particle_tree->Branch("embed_id", &m_particle_embed_id, "embed_id/I");
    m_particle_tree->Branch("is_primary", &m_particle_is_primary, "is_primary/I");
    m_particle_tree->Branch("px", &m_particle_px, "px/F");
    m_particle_tree->Branch("py", &m_particle_py, "py/F");
    m_particle_tree->Branch("pz", &m_particle_pz, "pz/F");
    m_particle_tree->Branch("e", &m_particle_e, "e/F");
    m_particle_tree->Branch("pt", &m_particle_pt, "pt/F");
    m_particle_tree->Branch("eta", &m_particle_eta, "eta/F");
    m_particle_tree->Branch("phi", &m_particle_phi, "phi/F");
    m_particle_tree->Branch("vx", &m_particle_vx, "vx/F");
    m_particle_tree->Branch("vy", &m_particle_vy, "vy/F");
    m_particle_tree->Branch("vz", &m_particle_vz, "vz/F");
    m_particle_tree->Branch("vt", &m_particle_vt, "vt/F");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcTruthPointTree::process_event(PHCompositeNode *topNode)
{
  int event_number = m_event_index;
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    event_number = event_header->get_EvtSequence();
  }
  else if (auto *g4_event_header = findNode::getClass<PHG4EventHeader>(topNode, "EventHeader"))
  {
    event_number = g4_event_header->get_EvtSequence();
  }

  auto *truth_points = findNode::getClass<PHG4HitContainer>(topNode, m_truth_point_node);
  auto *truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_truth_info_node);

  int n_truth_points = 0;
  if (truth_points)
  {
    PHG4HitContainer::ConstRange range = truth_points->getHits();
    for (auto hit_iter = range.first; hit_iter != range.second; ++hit_iter)
    {
      const PHG4Hit *hit = hit_iter->second;
      if (!hit)
      {
        continue;
      }

      reset_point_row();
      m_point_event = event_number;
      m_point_event_index = m_event_index;
      m_point_hit_key = static_cast<Long64_t>(hit_iter->first);
      m_point_track_id = hit->get_trkid();
      m_point_shower_id = hit->get_shower_id();
      m_point_layer = static_cast<int>(hit->get_layer());
      m_point_x = hit->get_x(0);
      m_point_y = hit->get_y(0);
      m_point_z = hit->get_z(0);
      m_point_t = hit->get_t(0);
      m_point_r = std::sqrt(square(m_point_x) + square(m_point_y));
      m_point_phi = std::atan2(m_point_y, m_point_x);
      m_point_side = (m_point_z >= 0.0F) ? 1 : 0;
      m_point_path = hit->get_path_length();
      m_point_px = hit->get_px(0);
      m_point_py = hit->get_py(0);
      m_point_pz = hit->get_pz(0);

      if (truth_info)
      {
        const PHG4Particle *particle = truth_info->GetParticle(m_point_track_id);
        if (particle)
        {
          m_point_pid = particle->get_pid();
          m_point_parent_id = particle->get_parent_id();
          m_point_primary_id = particle->get_primary_id();
          m_point_vtx_id = particle->get_vtx_id();
          m_point_barcode = particle->get_barcode();
          m_point_embed_id = truth_info->isEmbeded(m_point_track_id);
          m_point_is_primary = truth_info->is_primary(particle) ? 1 : 0;
          m_point_truth_px = particle->get_px();
          m_point_truth_py = particle->get_py();
          m_point_truth_pz = particle->get_pz();
          m_point_truth_e = particle->get_e();

          if (auto *vtx = truth_info->GetVtx(m_point_vtx_id))
          {
            m_point_vx = vtx->get_x();
            m_point_vy = vtx->get_y();
            m_point_vz = vtx->get_z();
            m_point_vt = vtx->get_t();
          }
        }
      }

      m_point_tree->Fill();
      fill_display_trees();
      ++n_truth_points;
    }
  }

  int n_truth_particles = 0;
  int n_primary_particles = 0;
  if (truth_info)
  {
    PHG4TruthInfoContainer::ConstRange particle_range = truth_info->GetParticleRange();
    for (auto particle_iter = particle_range.first; particle_iter != particle_range.second; ++particle_iter)
    {
      const PHG4Particle *particle = particle_iter->second;
      if (!particle)
      {
        continue;
      }
      ++n_truth_particles;
      if (truth_info->is_primary(particle))
      {
        ++n_primary_particles;
      }

      if (!m_particle_tree)
      {
        continue;
      }

      reset_particle_row();
      m_particle_event = event_number;
      m_particle_event_index = m_event_index;
      m_particle_track_id = particle->get_track_id();
      m_particle_pid = particle->get_pid();
      m_particle_parent_id = particle->get_parent_id();
      m_particle_primary_id = particle->get_primary_id();
      m_particle_vtx_id = particle->get_vtx_id();
      m_particle_barcode = particle->get_barcode();
      m_particle_embed_id = truth_info->isEmbeded(m_particle_track_id);
      m_particle_is_primary = truth_info->is_primary(particle) ? 1 : 0;
      m_particle_px = particle->get_px();
      m_particle_py = particle->get_py();
      m_particle_pz = particle->get_pz();
      m_particle_e = particle->get_e();
      m_particle_pt = std::sqrt(square(m_particle_px) + square(m_particle_py));
      m_particle_eta = safe_eta(m_particle_px, m_particle_py, m_particle_pz);
      m_particle_phi = std::atan2(m_particle_py, m_particle_px);

      if (auto *vtx = truth_info->GetVtx(m_particle_vtx_id))
      {
        m_particle_vx = vtx->get_x();
        m_particle_vy = vtx->get_y();
        m_particle_vz = vtx->get_z();
        m_particle_vt = vtx->get_t();
      }

      m_particle_tree->Fill();
    }
  }

  reset_event_row();
  m_evt_event = event_number;
  m_evt_event_index = m_event_index;
  m_evt_n_truth_points = n_truth_points;
  m_evt_n_truth_particles = n_truth_particles;
  m_evt_n_primary_particles = n_primary_particles;
  m_event_tree->Fill();

  ++m_event_index;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcTruthPointTree::End(PHCompositeNode * /*topNode*/)
{
  if (!m_file)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_file->cd();
  if (m_event_tree)
  {
    m_event_tree->Write();
  }
  if (m_point_tree)
  {
    m_point_tree->Write();
  }
  if (m_particle_tree)
  {
    m_particle_tree->Write();
  }
  if (m_cluster_tree)
  {
    m_cluster_tree->Write();
  }
  if (m_hit_tree)
  {
    m_hit_tree->Write();
  }
  m_file->Close();
  delete m_file;
  m_file = nullptr;

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcTruthPointTree::reset_event_row()
{
  m_evt_event = 0;
  m_evt_event_index = 0;
  m_evt_n_truth_points = 0;
  m_evt_n_truth_particles = 0;
  m_evt_n_primary_particles = 0;
}

void PHG4TpcTruthPointTree::reset_point_row()
{
  m_point_event = 0;
  m_point_event_index = 0;
  m_point_hit_key = 0;
  m_point_track_id = 0;
  m_point_shower_id = 0;
  m_point_layer = 0;
  m_point_side = 0;
  m_point_pid = 0;
  m_point_parent_id = 0;
  m_point_primary_id = 0;
  m_point_vtx_id = 0;
  m_point_barcode = 0;
  m_point_embed_id = 0;
  m_point_is_primary = 0;
  m_point_x = quiet_nan();
  m_point_y = quiet_nan();
  m_point_z = quiet_nan();
  m_point_t = quiet_nan();
  m_point_r = quiet_nan();
  m_point_phi = quiet_nan();
  m_point_path = quiet_nan();
  m_point_px = quiet_nan();
  m_point_py = quiet_nan();
  m_point_pz = quiet_nan();
  m_point_truth_px = quiet_nan();
  m_point_truth_py = quiet_nan();
  m_point_truth_pz = quiet_nan();
  m_point_truth_e = quiet_nan();
  m_point_vx = quiet_nan();
  m_point_vy = quiet_nan();
  m_point_vz = quiet_nan();
  m_point_vt = quiet_nan();
  m_display_adc = 1.0F;
  m_display_tdriftmax = 0.0F;
  m_display_drift_velocity = 0.0F;
  m_display_zdriftlength = 0.0F;
  m_display_pad = 0;
  m_display_tbin = 0;
  m_display_cluskey = 0;
  m_display_hitkeykey = 0;
  m_display_hitsetkey = 0;
}

void PHG4TpcTruthPointTree::add_display_branches(TTree *tree)
{
  if (!tree)
  {
    return;
  }

  tree->Branch("event", &m_point_event, "event/I");
  tree->Branch("event_index", &m_point_event_index, "event_index/I");
  tree->Branch("gx", &m_point_x, "gx/F");
  tree->Branch("gy", &m_point_y, "gy/F");
  tree->Branch("gz", &m_point_z, "gz/F");
  tree->Branch("x", &m_point_x, "x/F");
  tree->Branch("y", &m_point_y, "y/F");
  tree->Branch("z", &m_point_z, "z/F");
  tree->Branch("t", &m_point_t, "t/F");
  tree->Branch("r", &m_point_r, "r/F");
  tree->Branch("phi", &m_point_phi, "phi/F");
  tree->Branch("side", &m_point_side, "side/I");
  tree->Branch("layer", &m_point_layer, "layer/I");
  tree->Branch("track_id", &m_point_track_id, "track_id/I");
  tree->Branch("shower_id", &m_point_shower_id, "shower_id/I");
  tree->Branch("pid", &m_point_pid, "pid/I");
  tree->Branch("parent_id", &m_point_parent_id, "parent_id/I");
  tree->Branch("primary_id", &m_point_primary_id, "primary_id/I");
  tree->Branch("vtx_id", &m_point_vtx_id, "vtx_id/I");
  tree->Branch("barcode", &m_point_barcode, "barcode/I");
  tree->Branch("embed_id", &m_point_embed_id, "embed_id/I");
  tree->Branch("is_primary", &m_point_is_primary, "is_primary/I");
  tree->Branch("hit_key", &m_point_hit_key, "hit_key/L");
  tree->Branch("path", &m_point_path, "path/F");
  tree->Branch("px", &m_point_px, "px/F");
  tree->Branch("py", &m_point_py, "py/F");
  tree->Branch("pz", &m_point_pz, "pz/F");
  tree->Branch("truth_px", &m_point_truth_px, "truth_px/F");
  tree->Branch("truth_py", &m_point_truth_py, "truth_py/F");
  tree->Branch("truth_pz", &m_point_truth_pz, "truth_pz/F");
  tree->Branch("truth_e", &m_point_truth_e, "truth_e/F");
  tree->Branch("vx", &m_point_vx, "vx/F");
  tree->Branch("vy", &m_point_vy, "vy/F");
  tree->Branch("vz", &m_point_vz, "vz/F");
  tree->Branch("vt", &m_point_vt, "vt/F");

  tree->Branch("adc", &m_display_adc, "adc/F");
  tree->Branch("tdriftmax", &m_display_tdriftmax, "tdriftmax/F");
  tree->Branch("driftVelocity", &m_display_drift_velocity, "driftVelocity/F");
  tree->Branch("zdriftlength", &m_display_zdriftlength, "zdriftlength/F");
  tree->Branch("pad", &m_display_pad, "pad/I");
  tree->Branch("tbin", &m_display_tbin, "tbin/I");
  tree->Branch("cluskey", &m_display_cluskey, "cluskey/L");
  tree->Branch("hitkeykey", &m_display_hitkeykey, "hitkeykey/L");
  tree->Branch("hitsetkey", &m_display_hitsetkey, "hitsetkey/L");
}

void PHG4TpcTruthPointTree::fill_display_trees()
{
  m_display_adc = 1.0F;
  m_display_tdriftmax = 0.0F;
  m_display_drift_velocity = 0.0F;
  m_display_zdriftlength = 0.0F;
  m_display_pad = m_point_layer;
  m_display_tbin = std::isfinite(m_point_t) ? static_cast<int>(std::lround(m_point_t)) : 0;
  m_display_cluskey = m_point_hit_key;
  m_display_hitkeykey = m_point_hit_key;
  m_display_hitsetkey = m_point_layer;

  if (m_cluster_tree)
  {
    m_cluster_tree->Fill();
  }
  if (m_hit_tree)
  {
    m_hit_tree->Fill();
  }
}

void PHG4TpcTruthPointTree::reset_particle_row()
{
  m_particle_event = 0;
  m_particle_event_index = 0;
  m_particle_track_id = 0;
  m_particle_pid = 0;
  m_particle_parent_id = 0;
  m_particle_primary_id = 0;
  m_particle_vtx_id = 0;
  m_particle_barcode = 0;
  m_particle_embed_id = 0;
  m_particle_is_primary = 0;
  m_particle_px = quiet_nan();
  m_particle_py = quiet_nan();
  m_particle_pz = quiet_nan();
  m_particle_e = quiet_nan();
  m_particle_pt = quiet_nan();
  m_particle_eta = quiet_nan();
  m_particle_phi = quiet_nan();
  m_particle_vx = quiet_nan();
  m_particle_vy = quiet_nan();
  m_particle_vz = quiet_nan();
  m_particle_vt = quiet_nan();
}
