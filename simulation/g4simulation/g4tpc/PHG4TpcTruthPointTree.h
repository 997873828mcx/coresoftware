// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_PHG4TPCTRUTHPOINTTREE_H
#define G4TPC_PHG4TPCTRUTHPOINTTREE_H

#include <fun4all/SubsysReco.h>

#include <RtypesCore.h>

#include <string>

class PHCompositeNode;
class TFile;
class TTree;

class PHG4TpcTruthPointTree : public SubsysReco
{
 public:
  PHG4TpcTruthPointTree(const std::string &name = "PHG4TpcTruthPointTree",
                        const std::string &filename = "TpcTruthPoints.root");
  ~PHG4TpcTruthPointTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_output_file(const std::string &filename) { m_filename = filename; }
  void set_truth_point_node(const std::string &name) { m_truth_point_node = name; }
  void set_truth_info_node(const std::string &name) { m_truth_info_node = name; }
  void set_write_particle_tree(const bool value) { m_write_particle_tree = value; }

 private:
  void reset_event_row();
  void reset_point_row();
  void reset_particle_row();
  void add_display_branches(TTree *tree);
  void fill_display_trees();

  std::string m_filename;
  std::string m_truth_point_node{"G4HIT_TPC_TRUECLUSTER"};
  std::string m_truth_info_node{"G4TruthInfo"};

  TFile *m_file{nullptr};
  TTree *m_event_tree{nullptr};
  TTree *m_point_tree{nullptr};
  TTree *m_particle_tree{nullptr};
  TTree *m_cluster_tree{nullptr};
  TTree *m_hit_tree{nullptr};

  bool m_write_particle_tree{true};
  int m_event_index{0};

  int m_evt_event{0};
  int m_evt_event_index{0};
  int m_evt_n_truth_points{0};
  int m_evt_n_truth_particles{0};
  int m_evt_n_primary_particles{0};

  int m_point_event{0};
  int m_point_event_index{0};
  Long64_t m_point_hit_key{0};
  int m_point_track_id{0};
  int m_point_shower_id{0};
  int m_point_layer{0};
  int m_point_side{0};
  int m_point_pid{0};
  int m_point_parent_id{0};
  int m_point_primary_id{0};
  int m_point_vtx_id{0};
  int m_point_barcode{0};
  int m_point_embed_id{0};
  int m_point_is_primary{0};
  float m_point_x{0};
  float m_point_y{0};
  float m_point_z{0};
  float m_point_t{0};
  float m_point_r{0};
  float m_point_phi{0};
  float m_point_path{0};
  float m_point_px{0};
  float m_point_py{0};
  float m_point_pz{0};
  float m_point_truth_px{0};
  float m_point_truth_py{0};
  float m_point_truth_pz{0};
  float m_point_truth_e{0};
  float m_point_vx{0};
  float m_point_vy{0};
  float m_point_vz{0};
  float m_point_vt{0};
  float m_display_adc{1};
  float m_display_tdriftmax{0};
  float m_display_drift_velocity{0};
  float m_display_zdriftlength{0};
  int m_display_pad{0};
  int m_display_tbin{0};
  Long64_t m_display_cluskey{0};
  Long64_t m_display_hitkeykey{0};
  Long64_t m_display_hitsetkey{0};

  int m_particle_event{0};
  int m_particle_event_index{0};
  int m_particle_track_id{0};
  int m_particle_pid{0};
  int m_particle_parent_id{0};
  int m_particle_primary_id{0};
  int m_particle_vtx_id{0};
  int m_particle_barcode{0};
  int m_particle_embed_id{0};
  int m_particle_is_primary{0};
  float m_particle_px{0};
  float m_particle_py{0};
  float m_particle_pz{0};
  float m_particle_e{0};
  float m_particle_pt{0};
  float m_particle_eta{0};
  float m_particle_phi{0};
  float m_particle_vx{0};
  float m_particle_vy{0};
  float m_particle_vz{0};
  float m_particle_vt{0};
};

#endif
