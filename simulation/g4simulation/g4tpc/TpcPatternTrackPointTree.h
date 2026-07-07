// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_TPCPATTERNTRACKPOINTTREE_H
#define G4TPC_TPCPATTERNTRACKPOINTTREE_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TTree;

class TpcPatternTrackPointTree : public SubsysReco
{
 public:
  TpcPatternTrackPointTree(const std::string &name = "TpcPatternTrackPointTree",
                           const std::string &filename = "TpcPatternTrackPoints.root");
  ~TpcPatternTrackPointTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_output_file(const std::string &filename) { m_filename = filename; }
  void set_pattern_cluster_track_node(const std::string &name) { m_pattern_cluster_track_node = name; }
  void set_pattern_final_track_node(const std::string &name) { m_pattern_final_track_node = name; }
  void set_primary_vertex(const double x, const double y, const double z);
  void set_min_points(const int value) { m_min_points = value; }

 private:
  struct Vec3
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  int get_event_number(PHCompositeNode *topNode) const;
  int get_run_number(PHCompositeNode *topNode) const;
  void create_branches();
  void reset_event_row();
  void reset_point_row();
  void reset_particle_row();

  static int sign_to_charge(double value);
  static int charge_to_pion_pid(int charge);

  std::string m_filename;
  std::string m_pattern_cluster_track_node{"TPCPOLYCLUSTERTRACKS"};
  std::string m_pattern_final_track_node{"FINALTRACKS"};
  Vec3 m_primary_vertex{0.0, 0.0, 0.0};

  int m_min_points{5};
  int m_event_index{0};

  TFile *m_file{nullptr};
  TTree *m_event_tree{nullptr};
  TTree *m_point_tree{nullptr};
  TTree *m_particle_tree{nullptr};

  int m_evt_run{0};
  int m_evt_event{0};
  int m_evt_event_index{0};
  int m_evt_n_tracks{0};
  int m_evt_n_points{0};

  int m_point_event{0};
  int m_point_event_index{0};
  long long m_point_hit_key{0};
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
  float m_point_x{0.0F};
  float m_point_y{0.0F};
  float m_point_z{0.0F};
  float m_point_t{0.0F};
  float m_point_r{0.0F};
  float m_point_phi{0.0F};
  float m_point_path{0.0F};
  float m_point_px{0.0F};
  float m_point_py{0.0F};
  float m_point_pz{0.0F};
  float m_point_truth_px{0.0F};
  float m_point_truth_py{0.0F};
  float m_point_truth_pz{0.0F};
  float m_point_truth_e{0.0F};
  float m_point_vx{0.0F};
  float m_point_vy{0.0F};
  float m_point_vz{0.0F};
  float m_point_vt{0.0F};

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
  float m_particle_px{0.0F};
  float m_particle_py{0.0F};
  float m_particle_pz{0.0F};
  float m_particle_e{0.0F};
  float m_particle_pt{0.0F};
  float m_particle_eta{0.0F};
  float m_particle_phi{0.0F};
  float m_particle_vx{0.0F};
  float m_particle_vy{0.0F};
  float m_particle_vz{0.0F};
  float m_particle_vt{0.0F};

  unsigned long long m_counter_tracks{0};
  unsigned long long m_counter_points{0};
  unsigned long long m_counter_missing_final{0};
  unsigned long long m_counter_bad_charge{0};
};

#endif
