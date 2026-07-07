#include "TpcPatternTrackPointTree.h"

#include <g4main/PHG4EventHeader.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <map>

#include "/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/FinalTrack.h"
#include "/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/FinalTrackContainer.h"
#include "/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/TpcPolyClusterTrack.h"
#include "/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/TpcPolyClusterTrackContainer.h"
#define G4TPC_HAS_INMODULETRACKS 1

namespace
{
  [[maybe_unused]] double momentum_energy(const double px, const double py, const double pz, const double mass = 0.13957039)
  {
    return std::sqrt(px * px + py * py + pz * pz + mass * mass);
  }
}  // namespace

TpcPatternTrackPointTree::TpcPatternTrackPointTree(const std::string &name,
                                                   const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

void TpcPatternTrackPointTree::set_primary_vertex(const double x, const double y, const double z)
{
  m_primary_vertex = {x, y, z};
}

int TpcPatternTrackPointTree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = new TFile(m_filename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": failed to create output file " << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  create_branches();
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPatternTrackPointTree::process_event(PHCompositeNode *topNode)
{
  const int run_number = get_run_number(topNode);
  const int event_number = get_event_number(topNode);

#if !G4TPC_HAS_INMODULETRACKS
  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << Name()
              << ": built without inmoduletracks headers; pattern track-point dump disabled" << std::endl;
  }
  (void) run_number;
  (void) event_number;
  ++m_event_index;
  return Fun4AllReturnCodes::EVENT_OK;
#else
  auto *cluster_track_object = findNode::getClass<PHObject>(topNode, m_pattern_cluster_track_node);
  auto *cluster_tracks = static_cast<TpcPolyClusterTrackContainer *>(cluster_track_object);
  auto *final_track_object = findNode::getClass<PHObject>(topNode, m_pattern_final_track_node);
  auto *final_tracks = static_cast<FinalTrackContainer *>(final_track_object);

  if (!cluster_tracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << Name() << ": missing pattern cluster track node "
                << m_pattern_cluster_track_node << std::endl;
    }
    ++m_event_index;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  if (!final_tracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << Name() << ": missing pattern final track node "
                << m_pattern_final_track_node << "; cannot assign charge" << std::endl;
    }
    ++m_event_index;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  reset_event_row();
  m_evt_run = run_number;
  m_evt_event = event_number;
  m_evt_event_index = m_event_index;

  reset_particle_row();
  m_particle_event = event_number;
  m_particle_event_index = m_event_index;
  m_particle_track_id = 0;
  m_particle_pid = 0;
  m_particle_is_primary = 1;
  m_particle_vx = static_cast<float>(m_primary_vertex.x);
  m_particle_vy = static_cast<float>(m_primary_vertex.y);
  m_particle_vz = static_cast<float>(m_primary_vertex.z);
  m_particle_tree->Fill();

  std::map<unsigned int, const FinalTrack *> final_by_source_full_track_id;
  std::map<unsigned int, const FinalTrack *> final_by_track_id;
  for (unsigned int ifinal = 0; ifinal < final_tracks->size(); ++ifinal)
  {
    const FinalTrack *final_track = final_tracks->get_track(ifinal);
    if (!final_track || final_track->get_fit_status() == 0)
    {
      continue;
    }
    final_by_source_full_track_id[final_track->get_source_full_track_id()] = final_track;
    final_by_track_id[final_track->get_track_id()] = final_track;
  }

  for (unsigned int itrack = 0; itrack < cluster_tracks->size(); ++itrack)
  {
    const TpcPolyClusterTrack *cluster_track = cluster_tracks->get_track(itrack);
    if (!cluster_track || !cluster_track->isValid())
    {
      continue;
    }
    if (cluster_track->size_clusters() < static_cast<unsigned int>(m_min_points))
    {
      continue;
    }

    const unsigned int source_full_track_id = cluster_track->get_source_full_track_id();
    const unsigned int pattern_track_id = cluster_track->get_track_id();
    const int track_id = static_cast<int>(pattern_track_id != 0 ? pattern_track_id : itrack + 1);

    const FinalTrack *final_track = nullptr;
    auto final_iter = final_by_source_full_track_id.find(source_full_track_id);
    if (final_iter != final_by_source_full_track_id.end())
    {
      final_track = final_iter->second;
    }
    else
    {
      final_iter = final_by_track_id.find(pattern_track_id);
      if (final_iter != final_by_track_id.end())
      {
        final_track = final_iter->second;
      }
    }

    if (!final_track)
    {
      ++m_counter_missing_final;
      continue;
    }

    const int charge = sign_to_charge(final_track->get_charge());
    if (charge == 0)
    {
      ++m_counter_bad_charge;
      continue;
    }
    const int pid = charge_to_pion_pid(charge);
    const double px = final_track->get_px();
    const double py = final_track->get_py();
    const double pz = final_track->get_pz();

    reset_particle_row();
    m_particle_event = event_number;
    m_particle_event_index = m_event_index;
    m_particle_track_id = track_id;
    m_particle_pid = pid;
    m_particle_parent_id = 0;
    m_particle_is_primary = 0;
    m_particle_px = static_cast<float>(px);
    m_particle_py = static_cast<float>(py);
    m_particle_pz = static_cast<float>(pz);
    m_particle_e = static_cast<float>(momentum_energy(px, py, pz));
    m_particle_pt = static_cast<float>(std::sqrt(px * px + py * py));
    m_particle_phi = static_cast<float>(std::atan2(py, px));
    m_particle_vx = static_cast<float>(m_primary_vertex.x);
    m_particle_vy = static_cast<float>(m_primary_vertex.y);
    m_particle_vz = static_cast<float>(m_primary_vertex.z);
    m_particle_tree->Fill();

    ++m_evt_n_tracks;
    ++m_counter_tracks;

    for (unsigned int icluster = 0; icluster < cluster_track->size_clusters(); ++icluster)
    {
      const double x = cluster_track->get_cluster_x(icluster);
      const double y = cluster_track->get_cluster_y(icluster);
      const double z = cluster_track->get_cluster_z(icluster);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
      {
        continue;
      }

      reset_point_row();
      m_point_event = event_number;
      m_point_event_index = m_event_index;
      m_point_hit_key = static_cast<long long>(icluster);
      m_point_track_id = track_id;
      m_point_shower_id = static_cast<int>(source_full_track_id);
      m_point_layer = static_cast<int>(cluster_track->get_cluster_layer(icluster));
      m_point_side = cluster_track->get_side();
      m_point_pid = pid;
      m_point_parent_id = 0;
      m_point_x = static_cast<float>(x);
      m_point_y = static_cast<float>(y);
      m_point_z = static_cast<float>(z);
      m_point_r = static_cast<float>(std::sqrt(x * x + y * y));
      m_point_phi = static_cast<float>(std::atan2(y, x));
      m_point_path = static_cast<float>(m_point_layer + 1.0e-3 * static_cast<double>(icluster));
      m_point_px = static_cast<float>(px);
      m_point_py = static_cast<float>(py);
      m_point_pz = static_cast<float>(pz);
      m_point_truth_px = m_point_px;
      m_point_truth_py = m_point_py;
      m_point_truth_pz = m_point_pz;
      m_point_truth_e = static_cast<float>(momentum_energy(px, py, pz));
      m_point_vx = static_cast<float>(m_primary_vertex.x);
      m_point_vy = static_cast<float>(m_primary_vertex.y);
      m_point_vz = static_cast<float>(m_primary_vertex.z);
      m_point_tree->Fill();

      ++m_evt_n_points;
      ++m_counter_points;
    }
  }
  m_event_tree->Fill();
  ++m_event_index;
  return Fun4AllReturnCodes::EVENT_OK;
#endif
}

int TpcPatternTrackPointTree::End(PHCompositeNode * /*topNode*/)
{
  if (m_file)
  {
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
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": wrote tracks=" << m_counter_tracks
              << " points=" << m_counter_points
              << " missing_final=" << m_counter_missing_final
              << " bad_charge=" << m_counter_bad_charge << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPatternTrackPointTree::get_event_number(PHCompositeNode *topNode) const
{
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    return event_header->get_EvtSequence();
  }
  if (auto *g4_event_header = findNode::getClass<PHG4EventHeader>(topNode, "EventHeader"))
  {
    return g4_event_header->get_EvtSequence();
  }
  return m_event_index;
}

int TpcPatternTrackPointTree::get_run_number(PHCompositeNode *topNode) const
{
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    return event_header->get_RunNumber();
  }
  return 1;
}

void TpcPatternTrackPointTree::create_branches()
{
  m_event_tree = new TTree("events", "event summary for TPC pattern track points");
  m_event_tree->Branch("run", &m_evt_run, "run/I");
  m_event_tree->Branch("event", &m_evt_event, "event/I");
  m_event_tree->Branch("event_index", &m_evt_event_index, "event_index/I");
  m_event_tree->Branch("n_tracks", &m_evt_n_tracks, "n_tracks/I");
  m_event_tree->Branch("n_points", &m_evt_n_points, "n_points/I");

  m_point_tree = new TTree("tpc_truth_points", "TPC pattern cluster-track points");
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

  m_particle_tree = new TTree("truth_particles", "pattern track particle summary");
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

void TpcPatternTrackPointTree::reset_event_row()
{
  m_evt_run = 0;
  m_evt_event = 0;
  m_evt_event_index = 0;
  m_evt_n_tracks = 0;
  m_evt_n_points = 0;
}

void TpcPatternTrackPointTree::reset_point_row()
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
  m_point_x = 0.0F;
  m_point_y = 0.0F;
  m_point_z = 0.0F;
  m_point_t = 0.0F;
  m_point_r = 0.0F;
  m_point_phi = 0.0F;
  m_point_path = 0.0F;
  m_point_px = 0.0F;
  m_point_py = 0.0F;
  m_point_pz = 0.0F;
  m_point_truth_px = 0.0F;
  m_point_truth_py = 0.0F;
  m_point_truth_pz = 0.0F;
  m_point_truth_e = 0.0F;
  m_point_vx = 0.0F;
  m_point_vy = 0.0F;
  m_point_vz = 0.0F;
  m_point_vt = 0.0F;
}

void TpcPatternTrackPointTree::reset_particle_row()
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
  m_particle_px = 0.0F;
  m_particle_py = 0.0F;
  m_particle_pz = 0.0F;
  m_particle_e = 0.0F;
  m_particle_pt = 0.0F;
  m_particle_eta = 0.0F;
  m_particle_phi = 0.0F;
  m_particle_vx = 0.0F;
  m_particle_vy = 0.0F;
  m_particle_vz = 0.0F;
  m_particle_vt = 0.0F;
}

int TpcPatternTrackPointTree::sign_to_charge(const double value)
{
  if (!std::isfinite(value) || value == 0.0)
  {
    return 0;
  }
  return value > 0.0 ? 1 : -1;
}

int TpcPatternTrackPointTree::charge_to_pion_pid(const int charge)
{
  if (charge > 0)
  {
    return 211;
  }
  if (charge < 0)
  {
    return -211;
  }
  return 0;
}
