// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_TPCTRUTHV0DECAYTREE_H
#define G4TPC_TPCTRUTHV0DECAYTREE_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;
class TFile;
class TTree;

class TpcTruthV0DecayTree : public SubsysReco
{
 public:
  TpcTruthV0DecayTree(const std::string &name = "TpcTruthV0DecayTree",
                      const std::string &filename = "TruthV0Decays.root");
  ~TpcTruthV0DecayTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_output_file(const std::string &filename) { m_filename = filename; }
  void set_truth_info_node(const std::string &name) { m_truth_info_node = name; }
  void set_primary_vertex(const double x, const double y, const double z);
  void set_use_truth_primary_vertex(const bool value) { m_use_truth_primary_vertex = value; }

  struct Vec3
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

 private:
  struct Row
  {
    int run{0};
    int evt{0};

    int parent_id{0};
    int parent_pid{0};
    int parent_vtx_id{0};
    int decay_vtx_id{0};
    int parent_embed_id{0};
    int parent_is_primary{0};
    int parent_is_sphenix_primary{0};

    int daughter1_id{0};
    int daughter1_pid{0};
    int daughter2_id{0};
    int daughter2_pid{0};
    int n_daughters{0};

    float primary_x{0.0F};
    float primary_y{0.0F};
    float primary_z{0.0F};
    float decay_x{0.0F};
    float decay_y{0.0F};
    float decay_z{0.0F};
    float decay_t{0.0F};
    float Lxy{0.0F};
    float Lxyz{0.0F};

    float parent_px{0.0F};
    float parent_py{0.0F};
    float parent_pz{0.0F};
    float parent_e{0.0F};
    float parent_pt{0.0F};
    float parent_eta{0.0F};
    float parent_phi{0.0F};

    float px1{0.0F};
    float py1{0.0F};
    float pz1{0.0F};
    float e1{0.0F};
    float pt1{0.0F};
    float eta1{0.0F};
    float phi1{0.0F};

    float px2{0.0F};
    float py2{0.0F};
    float pz2{0.0F};
    float e2{0.0F};
    float pt2{0.0F};
    float eta2{0.0F};
    float phi2{0.0F};

    float alpha{0.0F};
    float qT{0.0F};
    float mass_Kshort{0.0F};
    float mass_Lambda{0.0F};
    float mass_AntiLambda{0.0F};
  };

  int get_event_number(PHCompositeNode *topNode) const;
  Vec3 get_primary_vertex(PHG4TruthInfoContainer *truth_info) const;
  const PHG4VtxPoint *get_decay_vertex(PHG4TruthInfoContainer *truth_info,
                                       const PHG4Particle *daughter1,
                                       const PHG4Particle *daughter2) const;
  bool fill_decay_row(PHG4TruthInfoContainer *truth_info,
                      const PHG4Particle *parent,
                      const int event_number,
                      const Vec3 &primary_vertex);
  void create_branches();
  void reset_row();

  static bool is_v0_parent(int pid);
  static bool is_expected_daughter_pair(int parent_pid, int pid1, int pid2);
  static double square(double value) { return value * value; }
  static double pt(const Vec3 &mom);
  static double eta(const Vec3 &mom);
  static double invariant_mass(const Vec3 &mom1, double mass1, const Vec3 &mom2, double mass2);
  static bool armenteros(const Vec3 &pplus, const Vec3 &pminus, double &alpha, double &qt);

  std::string m_filename{"TruthV0Decays.root"};
  std::string m_truth_info_node{"G4TruthInfo"};
  bool m_use_truth_primary_vertex{true};
  Vec3 m_fixed_primary_vertex{0.0, 0.0, 0.0};

  TFile *m_file{nullptr};
  TTree *m_tree{nullptr};
  Row m_row;

  long long m_counter_parents{0};
  long long m_counter_written{0};
};

#endif
