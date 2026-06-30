// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_TPCTRUTHEVENTTREE_H
#define G4TPC_TPCTRUTHEVENTTREE_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;
class TFile;
class TTree;

class TpcTruthEventTree : public SubsysReco
{
 public:
  TpcTruthEventTree(const std::string &name = "TpcTruthEventTree",
                    const std::string &filename = "TruthEventQA.root");
  ~TpcTruthEventTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_output_file(const std::string &filename) { m_filename = filename; }
  void set_truth_info_node(const std::string &name) { m_truth_info_node = name; }
  void set_primary_vertex(const double x, const double y, const double z);
  void set_use_truth_primary_vertex(const bool value) { m_use_truth_primary_vertex = value; }
  void set_charged_eta_max(const double value) { m_charged_eta_max = value; }
  void set_charged_pt_min(const double value) { m_charged_pt_min = value; }
  void set_v0_abs_y_max(const double value) { m_v0_abs_y_max = value; }
  void set_v0_pt_min(const double value) { m_v0_pt_min = value; }

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

    int primary_vtx_id{0};
    int primary_embed_id{0};
    int n_primary_vertices{0};
    float primary_x{0.0F};
    float primary_y{0.0F};
    float primary_z{0.0F};
    float primary_t{0.0F};

    int n_particles_all{0};
    int n_particles_primary{0};
    int n_particles_secondary{0};

    int n_charged_primary{0};
    int n_charged_primary_no_daughters{0};
    int n_charged_primary_eta{0};
    int n_charged_primary_eta_pt{0};
    int n_charged_primary_no_daughters_eta{0};
    int n_charged_primary_no_daughters_eta_pt{0};

    int n_pi_plus_primary{0};
    int n_pi_minus_primary{0};
    int n_k_plus_primary{0};
    int n_k_minus_primary{0};
    int n_p_primary{0};
    int n_pbar_primary{0};

    int n_kshort{0};
    int n_lambda{0};
    int n_antilambda{0};
    int n_kshort_primary{0};
    int n_lambda_primary{0};
    int n_antilambda_primary{0};
    int n_kshort_charged_decay{0};
    int n_lambda_charged_decay{0};
    int n_antilambda_charged_decay{0};
    int n_kshort_primary_charged_decay{0};
    int n_lambda_primary_charged_decay{0};
    int n_antilambda_primary_charged_decay{0};
    int n_kshort_fiducial{0};
    int n_lambda_fiducial{0};
    int n_antilambda_fiducial{0};

    float charged_eta_max{0.0F};
    float charged_pt_min{0.0F};
    float v0_abs_y_max{0.0F};
    float v0_pt_min{0.0F};
  };

  int get_run_number(PHCompositeNode *topNode) const;
  int get_event_number(PHCompositeNode *topNode) const;
  const PHG4VtxPoint *get_primary_vertex(PHG4TruthInfoContainer *truth_info) const;
  bool has_expected_charged_decay(const PHG4Particle *parent,
                                  const std::vector<const PHG4Particle *> &daughters) const;
  bool passes_v0_fiducial(const PHG4Particle *parent) const;
  void count_primary_particle(const PHG4Particle *particle, bool has_daughters);
  void count_v0_parent(const PHG4Particle *parent,
                       const std::vector<const PHG4Particle *> &daughters);
  void create_branches();
  void reset_row();

  static bool is_charged(int pid);
  static bool is_v0_parent(int pid);
  static bool is_expected_daughter_pair(int parent_pid, int pid1, int pid2);
  static double square(double value) { return value * value; }
  static double pt(const Vec3 &mom);
  static double eta(const Vec3 &mom);
  static double rapidity(const PHG4Particle *particle);
  static Vec3 momentum(const PHG4Particle *particle);

  std::string m_filename{"TruthEventQA.root"};
  std::string m_truth_info_node{"G4TruthInfo"};
  bool m_use_truth_primary_vertex{true};
  Vec3 m_fixed_primary_vertex{0.0, 0.0, 0.0};

  double m_charged_eta_max{1.1};
  double m_charged_pt_min{0.2};
  double m_v0_abs_y_max{1.1};
  double m_v0_pt_min{0.0};

  TFile *m_file{nullptr};
  TTree *m_tree{nullptr};
  Row m_row;

  long long m_counter_events{0};
};

#endif
