// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_TPCV0CANDIDATETREE_H
#define G4TPC_TPCV0CANDIDATETREE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class TFile;
class TTree;

class TpcV0CandidateTree : public SubsysReco
{
 public:
  TpcV0CandidateTree(const std::string &name = "TpcV0CandidateTree",
                     const std::string &filename = "TpcV0Candidates.root");
  ~TpcV0CandidateTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_output_file(const std::string &filename) { m_filename = filename; }
  void set_truth_point_node(const std::string &name) { m_truth_point_node = name; }
  void set_truth_info_node(const std::string &name) { m_truth_info_node = name; }
  void set_use_truth_primary_vertex(const bool value) { m_use_truth_primary_vertex = value; }
  void set_primary_vertex(const double x, const double y, const double z);

  void set_min_points(const int value) { m_min_points = value; }
  void set_fit_helix(const bool value) { m_fit_helix_tracks = value; }
  void set_fit_first_points(const int value) { m_fit_first_points = value; }
  void set_bfield(const double value) { m_bfield_t = value; }
  void set_theta_extension(const double value) { m_theta_extension = value; }
  void set_coarse_steps(const int value) { m_coarse_steps = value; }
  void set_pca_candidates(const int value) { m_pca_candidates = value; }
  void set_downstream_margin(const double value) { m_downstream_margin = value; }
  void set_prefer_positive_pointing(const bool value) { m_prefer_positive_pointing = value; }

  void set_pre_track_pt_min(const double value) { m_pre_track_pt_min = value; }
  void set_pre_track_dca_xy_min(const double value) { m_pre_track_dca_xy_min = value; }
  void set_pre_track_dca_z_min(const double value) { m_pre_track_dca_z_min = value; }
  void set_pre_track_dca_xy_max(const double value) { m_pre_track_dca_xy_max = value; }
  void set_pre_track_dca_z_max(const double value) { m_pre_track_dca_z_max = value; }
  void set_pre_pair_dca_max(const double value) { m_pre_pair_dca_max = value; }
  void set_pre_lproj_min(const double value) { m_pre_lproj_min = value; }
  void set_pre_cos_theta_min(const double value) { m_pre_cos_theta_min = value; }

 private:
  struct Vec3
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  struct TruthPoint
  {
    int track_id{0};
    int shower_id{0};
    int layer{0};
    Vec3 position;
    Vec3 momentum;
    double t{0.0};
    double path{0.0};
  };

  struct HelixFit
  {
    double cx{0.0};
    double cy{0.0};
    double radius{0.0};
    double z0{0.0};
    double pitch{0.0};
    double theta_first{0.0};
    double theta_last{0.0};
    double theta_min{0.0};
    double theta_max{0.0};
    double direction{1.0};
    double bfield_t{1.4};
  };

  struct Tracklet
  {
    int track_id{0};
    int shower_id{0};
    int pid{0};
    int parent_id{0};
    int parent_pid{0};
    int primary_id{0};
    int vtx_id{0};
    int barcode{0};
    int embed_id{0};
    int is_primary{0};
    int charge{0};
    int npoints{0};
    Vec3 position;
    Vec3 momentum;
    Vec3 truth_momentum;
    double truth_e{0.0};
    Vec3 truth_vertex;
    double truth_vt{0.0};
    std::vector<TruthPoint> points;
    bool has_helix{false};
    HelixFit helix;
  };

  struct HelixPca
  {
    Vec3 pca1;
    Vec3 pca2;
    double dca{0.0};
    double theta1{0.0};
    double theta2{0.0};
  };

  struct PairRow
  {
    int run{0};
    int evt{0};
    short cross1{0};
    short cross2{0};

    float px1{0.0F};
    float py1{0.0F};
    float pz1{0.0F};
    float px2{0.0F};
    float py2{0.0F};
    float pz2{0.0F};

    float dca_xy1{0.0F};
    float dca_z1{0.0F};
    float dca_xy2{0.0F};
    float dca_z2{0.0F};
    float pairDCA{0.0F};

    float alpha{0.0F};
    float qT{0.0F};
    float charge1{0.0F};
    float charge2{0.0F};
    float cosThetaReco{0.0F};
    float Lproj{0.0F};

    float pca_x{0.0F};
    float pca_y{0.0F};
    float pca_z{0.0F};
    float pca1_x{0.0F};
    float pca1_y{0.0F};
    float pca1_z{0.0F};
    float pca2_x{0.0F};
    float pca2_y{0.0F};
    float pca2_z{0.0F};

    float v0_px{0.0F};
    float v0_py{0.0F};
    float v0_pz{0.0F};
    float v0_pt{0.0F};
    float mass_Kshort{0.0F};
    float mass_Lambda{0.0F};
    float mass_AntiLambda{0.0F};

    float true_decay_x{0.0F};
    float true_decay_y{0.0F};
    float true_decay_z{0.0F};
    float pca_to_true_3d{0.0F};
    float pca_to_true_xy{0.0F};
    float pca_to_true_z{0.0F};
    float truth_alpha{0.0F};
    float truth_qT{0.0F};
    float delta_alpha{0.0F};
    float delta_qT{0.0F};
    float truth_px1{0.0F};
    float truth_py1{0.0F};
    float truth_pz1{0.0F};
    float truth_px2{0.0F};
    float truth_py2{0.0F};
    float truth_pz2{0.0F};
    float cos_mom1_truth{0.0F};
    float cos_mom2_truth{0.0F};
    float pca_theta1{0.0F};
    float pca_theta2{0.0F};

    int track_id1{0};
    int track_id2{0};
    int pid1{0};
    int pid2{0};
    int parent_id1{0};
    int parent_id2{0};
    int parent_pid{0};
    short npoints1{0};
    short npoints2{0};
  };

  struct LinePca
  {
    Vec3 pca1;
    Vec3 pca2;
    double dca{0.0};
    double step1{0.0};
    double step2{0.0};
  };

  int get_event_number(PHCompositeNode *topNode) const;
  Vec3 get_primary_vertex(PHG4TruthInfoContainer *truth_info) const;
  std::map<int, Tracklet> build_tracklets(PHG4HitContainer *truth_points,
                                          PHG4TruthInfoContainer *truth_info) const;
  bool make_pair_row(const Tracklet &track1, const Tracklet &track2,
                     const Vec3 &primary_vertex, const int event_number);
  void reset_pair_row();
  void create_branches();

  static int pdg_charge(int pid);
  static float quiet_nan();
  static bool finite(const Vec3 &value);
  static Vec3 add(const Vec3 &lhs, const Vec3 &rhs);
  static Vec3 subtract(const Vec3 &lhs, const Vec3 &rhs);
  static Vec3 scale(const Vec3 &value, double factor);
  static double dot(const Vec3 &lhs, const Vec3 &rhs);
  static Vec3 cross(const Vec3 &lhs, const Vec3 &rhs);
  static double norm(const Vec3 &value);
  static Vec3 unit(const Vec3 &value);
  static double pt(const Vec3 &value);
  static double distance(const Vec3 &lhs, const Vec3 &rhs);
  static double vector_cosine(const Vec3 &lhs, const Vec3 &rhs);

  static bool fit_helix(const std::vector<TruthPoint> &points, int fit_first_points,
                        double bfield_t, HelixFit &helix);
  static Vec3 helix_point(const HelixFit &helix, double theta);
  static Vec3 helix_tangent(const HelixFit &helix, double theta);
  static Vec3 helix_momentum(const HelixFit &helix, double theta);
  static std::pair<double, double> theta_search_range(const HelixFit &helix,
                                                      double theta_extension,
                                                      double downstream_margin);
  static bool line_line_pca(const Vec3 &pos1, const Vec3 &dir1,
                            const Vec3 &pos2, const Vec3 &dir2,
                            LinePca &pca, bool normalize_dirs);
  static HelixPca refine_helix_pair(const HelixFit &helix1, const HelixFit &helix2,
                                    double theta1, double theta2,
                                    double min1, double max1,
                                    double min2, double max2,
                                    double max_step);
  static std::vector<HelixPca> helix_helix_pca_candidates(const HelixFit &helix1,
                                                         const HelixFit &helix2,
                                                         double theta_extension,
                                                         int coarse_steps,
                                                         double downstream_margin,
                                                         int max_candidates);
  static std::pair<double, double> track_dca_to_vertex(const Vec3 &pos,
                                                       const Vec3 &mom,
                                                       const Vec3 &vertex);
  static std::pair<double, double> helix_dca_to_vertex(const HelixFit &helix,
                                                       const Vec3 &vertex);
  static bool armenteros(const Vec3 &pplus, const Vec3 &pminus,
                         double &alpha, double &qt);
  static double invariant_mass(const Vec3 &mom1, double mass1,
                               const Vec3 &mom2, double mass2);

  bool passes_preselection(const Tracklet &track1, const Tracklet &track2,
                           const Vec3 &primary_vertex) const;

  std::string m_filename;
  std::string m_truth_point_node{"G4HIT_TPC_TRUECLUSTER"};
  std::string m_truth_info_node{"G4TruthInfo"};

  TFile *m_file{nullptr};
  TTree *m_pair_tree{nullptr};
  PairRow m_pair;

  Vec3 m_fixed_primary_vertex{0.0, 0.0, 0.0};
  bool m_use_truth_primary_vertex{true};
  int m_min_points{5};
  bool m_fit_helix_tracks{true};
  int m_fit_first_points{8};
  double m_bfield_t{1.4};
  double m_theta_extension{2.0};
  int m_coarse_steps{64};
  int m_pca_candidates{32};
  double m_downstream_margin{0.2};
  bool m_prefer_positive_pointing{false};

  double m_pre_track_pt_min{0.2};
  double m_pre_track_dca_xy_min{0.03};
  double m_pre_track_dca_z_min{-1.0};
  double m_pre_track_dca_xy_max{-1.0};
  double m_pre_track_dca_z_max{-1.0};
  double m_pre_pair_dca_max{5.0};
  double m_pre_lproj_min{0.2};
  double m_pre_cos_theta_min{-2.0};

  std::uint64_t m_counter_raw_pairs{0};
  std::uint64_t m_counter_reject_charge{0};
  std::uint64_t m_counter_reject_preselection{0};
  std::uint64_t m_counter_reject_pca{0};
  std::uint64_t m_counter_reject_pointing{0};
  std::uint64_t m_counter_reject_ap{0};
  std::uint64_t m_counter_written{0};
};

#endif
