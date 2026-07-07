// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_TPCV0CANDIDATETREE_H
#define G4TPC_TPCV0CANDIDATETREE_H

#include "TpcTrackFit.h"

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
class FinalTrackContainer;
class TpcPolyClusterTrackContainer;

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
  void set_pattern_cluster_track_node(const std::string &name) { m_pattern_cluster_track_node = name; }
  void set_pattern_final_track_node(const std::string &name) { m_pattern_final_track_node = name; }
  void use_pattern_cluster_tracks(const bool value = true) { m_use_pattern_cluster_tracks = value; }
  void set_use_truth_primary_vertex(const bool value) { m_use_truth_primary_vertex = value; }
  void set_primary_vertex(const double x, const double y, const double z);

  void set_min_points(const int value) { m_min_points = value; }
  void set_fit_helix(const bool value)
  {
    m_fit_helix_tracks = value;
    if (value)
    {
      m_fit_kalman_tracks = false;
    }
  }
  void set_fit_kalman(const bool value)
  {
    m_fit_kalman_tracks = value;
    if (value)
    {
      m_fit_helix_tracks = false;
      m_use_final_track_helix = false;
    }
  }
  bool set_track_fit_method(const std::string &mode);
  void set_use_final_track_helix(const bool value) { m_use_final_track_helix = value; }
  bool set_point_order(const std::string &mode);
  void set_fit_first_points(const int value) { m_fit_first_points = value; }
  void set_bfield(const double value)
  {
    m_bfield_t = value;
    m_kalman_config.bfield_t = value;
  }
  void set_theta_extension(const double value) { m_theta_extension = value; }
  void set_coarse_steps(const int value) { m_coarse_steps = value; }
  void set_pca_candidates(const int value) { m_pca_candidates = value; }
  void set_downstream_margin(const double value) { m_downstream_margin = value; }
  void set_kalman_search(const double max_upstream_cm, const double downstream_margin_cm)
  {
    m_kalman_max_upstream_cm = max_upstream_cm;
    m_kalman_downstream_margin_cm = downstream_margin_cm;
  }
  void set_kalman_measurement_sigmas(const double xy_cm, const double z_cm)
  {
    set_kalman_measurement_sigmas(xy_cm, xy_cm, z_cm);
  }
  void set_kalman_measurement_sigmas(const double rphi_cm,
                                     const double r_cm,
                                     const double z_cm)
  {
    m_kalman_config.meas_sigma_rphi_cm = rphi_cm;
    m_kalman_config.meas_sigma_r_cm = r_cm;
    m_kalman_config.meas_sigma_z_cm = z_cm;
  }
  void set_kalman_process_sigmas(const double pos_cm,
                                 const double phi,
                                 const double qop_t,
                                 const double tanl)
  {
    m_kalman_config.process_sigma_pos_cm = pos_cm;
    m_kalman_config.process_sigma_phi = phi;
    m_kalman_config.process_sigma_qop_t = qop_t;
    m_kalman_config.process_sigma_tanl = tanl;
  }
  void set_kalman_material(const double x0_per_cm,
                           const double multiple_scattering_scale,
                           const double energy_loss_gev_per_cm,
                           const double energy_loss_sigma_fraction)
  {
    m_kalman_config.material_x0_per_cm = x0_per_cm;
    m_kalman_config.multiple_scattering_scale = multiple_scattering_scale;
    m_kalman_config.energy_loss_gev_per_cm = energy_loss_gev_per_cm;
    m_kalman_config.energy_loss_sigma_fraction = energy_loss_sigma_fraction;
  }
  void set_prefer_positive_pointing(const bool value) { m_prefer_positive_pointing = value; }

  void set_pre_track_pt_min(const double value) { m_pre_track_pt_min = value; }
  void set_pre_track_dca_xy_min(const double value) { m_pre_track_dca_xy_min = value; }
  void set_pre_track_dca_z_min(const double value) { m_pre_track_dca_z_min = value; }
  void set_pre_track_dca_xy_max(const double value) { m_pre_track_dca_xy_max = value; }
  void set_pre_track_dca_z_max(const double value) { m_pre_track_dca_z_max = value; }
  void set_pre_pair_dca_max(const double value) { m_pre_pair_dca_max = value; }
  void set_pre_lproj_min(const double value) { m_pre_lproj_min = value; }
  void set_pre_cos_theta_min(const double value) { m_pre_cos_theta_min = value; }
  void set_write_same_sign_pairs(const bool value) { m_write_same_sign_pairs = value; }
  void set_write_cluster_residual_tree(const bool value) { m_write_cluster_residual_tree = value; }

 private:
  using Vec3 = TpcTrackVec3;
  using TruthPoint = TpcTrackPoint;
  using HelixFit = TpcTrackHelix;
  using HelixPca = TpcTrackHelixPca;
  using LinePca = TpcTrackLinePca;
  using PointOrder = TpcTrackPointOrder;

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
    bool has_kalman{false};
    TpcKalmanResult kalman;
    double fit_chi2{0.0};
    int fit_ndf{0};
    double fit_chi2_ndf{0.0};
  };

  struct KalmanPca
  {
    Vec3 pca1;
    Vec3 pca2;
    double dca{0.0};
    double s1{0.0};
    double s2{0.0};
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
    float kalman_chi2_1{0.0F};
    float kalman_chi2_2{0.0F};
    float kalman_chi2_ndf1{0.0F};
    float kalman_chi2_ndf2{0.0F};
    float quality1{0.0F};
    float quality2{0.0F};

    int track_id1{0};
    int track_id2{0};
    int pid1{0};
    int pid2{0};
    int parent_id1{0};
    int parent_id2{0};
    int parent_pid{0};
    int kalman_ndof1{0};
    int kalman_ndof2{0};
    short npoints1{0};
    short npoints2{0};
  };

  struct TrackRow
  {
    int run{0};
    int evt{0};
    int track_id{0};
    int shower_id{0};
    int pid{0};
    int parent_id{0};
    int parent_pid{0};
    int charge{0};
    int npoints{0};
    int has_helix{0};
    int has_kalman{0};
    int is_primary{0};

    float px{0.0F};
    float py{0.0F};
    float pz{0.0F};
    float pt{0.0F};
    float p{0.0F};
    float x{0.0F};
    float y{0.0F};
    float z{0.0F};

    float first_x{0.0F};
    float first_y{0.0F};
    float first_z{0.0F};
    float first_r{0.0F};
    float last_x{0.0F};
    float last_y{0.0F};
    float last_z{0.0F};
    float last_r{0.0F};

    float dca_xy{0.0F};
    float dca_z{0.0F};
    float vertex_x{0.0F};
    float vertex_y{0.0F};
    float vertex_z{0.0F};

    float helix_cx{0.0F};
    float helix_cy{0.0F};
    float helix_radius{0.0F};
    float helix_z0{0.0F};
    float helix_pitch{0.0F};
    float helix_theta_first{0.0F};
    float helix_theta_last{0.0F};
    float helix_direction{0.0F};

    float kalman_chi2{0.0F};
    int kalman_ndof{0};
    float kalman_qop_t{0.0F};
    float kalman_omega{0.0F};
    float kalman_cx{0.0F};
    float kalman_cy{0.0F};
    float kalman_radius{0.0F};
    float fit_chi2{0.0F};
    int fit_ndf{0};
    float quality{0.0F};

    float truth_px{0.0F};
    float truth_py{0.0F};
    float truth_pz{0.0F};
    float cos_mom_truth{0.0F};

    std::vector<int> cluster_index_vec;
    std::vector<int> cluster_side_vec;
    std::vector<int> cluster_layer_vec;
    std::vector<float> cluster_z_vec;
    std::vector<float> cluster_r_vec;
    std::vector<float> cluster_phi_vec;
    std::vector<float> residual_z_vec;
    std::vector<float> residual_r_vec;
    std::vector<float> residual_rphi_vec;
  };

  struct ClusterResidualRow
  {
    int run{0};
    int evt{0};
    int track_id{0};
    int charge{0};
    int side{0};
    int layer{0};
    int cluster_index{0};
    int ntp_cluster{0};
    int npoints{0};
    int has_helix{0};
    int has_kalman{0};

    float cluster_x{0.0F};
    float cluster_y{0.0F};
    float cluster_z{0.0F};
    float cluster_r{0.0F};
    float cluster_phi{0.0F};

    float fit_x{0.0F};
    float fit_y{0.0F};
    float fit_z{0.0F};
    float fit_r{0.0F};
    float fit_phi{0.0F};

    float residual_x{0.0F};
    float residual_y{0.0F};
    float residual_z{0.0F};
    float residual_r{0.0F};
    float residual_rphi{0.0F};

    float fit_chi2{0.0F};
    int fit_ndf{0};
    float fit_chi2_ndf{0.0F};
  };

  int get_event_number(PHCompositeNode *topNode) const;
  int get_run_number(PHCompositeNode *topNode) const;
  Vec3 get_primary_vertex(PHG4TruthInfoContainer *truth_info) const;
  std::map<int, Tracklet> build_tracklets(PHG4HitContainer *truth_points,
                                          PHG4TruthInfoContainer *truth_info) const;
  std::map<int, Tracklet> build_pattern_tracklets(TpcPolyClusterTrackContainer *cluster_tracks,
                                                  FinalTrackContainer *final_tracks) const;
  bool make_pair_row(const Tracklet &track1, const Tracklet &track2,
                     const Vec3 &primary_vertex, const int run_number,
                     const int event_number);
  void fill_track_row(const Tracklet &tracklet, const Vec3 &primary_vertex,
                      int run_number, int event_number);
  void fill_cluster_residual_rows(const Tracklet &tracklet, const Vec3 &primary_vertex,
                                  int run_number, int event_number);
  void assign_fit_quality(Tracklet &tracklet) const;
  void reset_pair_row();
  void reset_track_row();
  void reset_cluster_residual_row();
  void create_branches();

  static int pdg_charge(int pid);
  static bool parse_point_order(const std::string &mode, PointOrder &order);
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
  static bool fit_circle_least_squares(const std::vector<TruthPoint> &points,
                                       std::size_t nfit,
                                       double &cx,
                                       double &cy,
                                       double &radius);
  static void order_track_points(std::vector<TruthPoint> &points, PointOrder order);

  static bool fit_helix(const std::vector<TruthPoint> &points, int fit_first_points,
                        int charge, double bfield_t, HelixFit &helix);
  bool fit_kalman(const std::vector<TruthPoint> &points,
                  int charge,
                  TpcKalmanResult &kalman) const;
  static bool helix_from_state(const Vec3 &position, const Vec3 &momentum,
                               int charge, double bfield_t, HelixFit &helix);
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
  static Vec3 kalman_point(const TpcKalmanResult &kalman,
                           double s_cm,
                           const TpcKalmanConfig &config,
                           const Vec3 &reference_vertex);
  static Vec3 kalman_tangent(const TpcKalmanResult &kalman,
                             double s_cm,
                             const TpcKalmanConfig &config,
                             const Vec3 &reference_vertex);
  static Vec3 kalman_momentum(const TpcKalmanResult &kalman,
                              double s_cm,
                              const TpcKalmanConfig &config,
                              const Vec3 &reference_vertex);
  static KalmanPca refine_kalman_pair(const TpcKalmanResult &kalman1,
                                      const TpcKalmanResult &kalman2,
                                      const TpcKalmanConfig &config,
                                      const Vec3 &reference_vertex,
                                      double s1, double s2,
                                      double min1, double max1,
                                      double min2, double max2,
                                      double max_step);
  static std::vector<KalmanPca> kalman_pca_candidates(const TpcKalmanResult &kalman1,
                                                     const TpcKalmanResult &kalman2,
                                                     const TpcKalmanConfig &config,
                                                     const Vec3 &reference_vertex,
                                                     double max_upstream_cm,
                                                     double downstream_margin_cm,
                                                     int coarse_steps,
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
  std::string m_pattern_cluster_track_node{"TPCPOLYCLUSTERTRACKS"};
  std::string m_pattern_final_track_node{"FINALTRACKS"};
  bool m_use_pattern_cluster_tracks{false};

  TFile *m_file{nullptr};
  TTree *m_pair_tree{nullptr};
  TTree *m_track_tree{nullptr};
  TTree *m_cluster_residual_tree{nullptr};
  PairRow m_pair;
  TrackRow m_track;
  ClusterResidualRow m_cluster_residual;

  Vec3 m_fixed_primary_vertex{0.0, 0.0, 0.0};
  bool m_use_truth_primary_vertex{true};
  int m_min_points{5};
  bool m_fit_helix_tracks{true};
  bool m_fit_kalman_tracks{false};
  bool m_use_final_track_helix{false};
  PointOrder m_point_order{PointOrder::Path};
  int m_fit_first_points{8};
  double m_bfield_t{1.4};
  TpcKalmanConfig m_kalman_config;
  double m_kalman_max_upstream_cm{80.0};
  double m_kalman_downstream_margin_cm{5.0};
  double m_theta_extension{2.0};
  int m_coarse_steps{64};
  int m_pca_candidates{32};
  double m_downstream_margin{0.2};
  bool m_prefer_positive_pointing{false};
  bool m_write_cluster_residual_tree{false};

  double m_pre_track_pt_min{0.2};
  double m_pre_track_dca_xy_min{0.03};
  double m_pre_track_dca_z_min{-1.0};
  double m_pre_track_dca_xy_max{-1.0};
  double m_pre_track_dca_z_max{-1.0};
  double m_pre_pair_dca_max{5.0};
  double m_pre_lproj_min{0.2};
  double m_pre_cos_theta_min{-2.0};
  bool m_write_same_sign_pairs{false};

  std::uint64_t m_counter_raw_pairs{0};
  std::uint64_t m_counter_reject_charge{0};
  std::uint64_t m_counter_reject_preselection{0};
  std::uint64_t m_counter_reject_pca{0};
  std::uint64_t m_counter_reject_pointing{0};
  std::uint64_t m_counter_reject_ap{0};
  std::uint64_t m_counter_written{0};
  std::uint64_t m_counter_tracks_written{0};
  std::uint64_t m_counter_cluster_residuals_written{0};
};

#endif
