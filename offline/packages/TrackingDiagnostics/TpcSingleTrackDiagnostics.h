#ifndef TRACKINGDIAGNOSTICS_TPCSINGLETRACKDIAGNOSTICS_H
#define TRACKINGDIAGNOSTICS_TPCSINGLETRACKDIAGNOSTICS_H

#include <fun4all/SubsysReco.h>
#include <tpctrackreco/TpcTrackDiagnosticFitter.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class PHField;
class TFile;
class TTree;
class Tpc_PolyCluster;
class Tpc_PolyTrack;

// Read-only diagnostic consumer of fixed pattern associations. Neither truth
// nor a vertex measurement enters any nominal or held-out Kalman fit.
class TpcSingleTrackDiagnostics : public SubsysReco
{
 public:
  explicit TpcSingleTrackDiagnostics(const std::string &name = "TpcSingleTrackDiagnostics",
                                     const std::string &filename = "single_track_diagnostics.root");
  ~TpcSingleTrackDiagnostics() override;
  int Init(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode *) override;

  void set_sample_id(std::uint64_t value) { m_sample_id = value; }
  void set_magnetic_field(const PHField *field) { m_config = TpcTrackDiagnosticFitter::yuri_config(field); }
  void set_write_leave_one_out(bool value) { m_write_leave_one_out = value; }
  void set_write_truth_diagnostics(bool value) { m_write_truth = value; }
  void set_provenance(const std::string &value) { m_provenance = value; }
  void set_inner_tpc_radius_cm(double value) { m_inner_radius_cm = value; }
  void set_readout_z_cm(double value) { m_readout_z_cm = value; }
  const TpcKalmanConfig &get_kalman_config() const { return m_config; }
  bool get_write_leave_one_out() const { return m_write_leave_one_out; }
  bool get_write_truth_diagnostics() const { return m_write_truth; }
  double get_readout_z_cm() const { return m_readout_z_cm; }
  // Synthetic tests may use uniform field; real campaigns require a field map.
  void set_require_field_map(bool value) { m_require_field_map = value; }
  void set_poly_cluster_node(const std::string &value) { m_cluster_node = value; }
  void set_poly_track_node(const std::string &value) { m_track_node = value; }

 private:
  struct Identity
  {
    std::uint64_t sample_id{0};
    std::uint64_t entry{0};
    int run{0};
    int event{0};
  };
  struct EventRow : Identity
  {
    unsigned int n_input_tracks{0}, n_input_clusters{0}, n_valid_polytracks{0};
    unsigned int n_kalman_success{0}, n_kalman_failure{0}, n_pca_failure{0};
    unsigned int n_missing_cluster_keys{0}, n_duplicate_track_ids{0}, n_duplicate_cluster_ids{0};
    unsigned int n_unassociated_clusters{0}, n_residual_rows{0}, n_loo_success{0}, n_loo_failure{0};
    int has_geometry{0}, has_crossing_decisions{0}, has_truth{0};
    double inner_reference_radius_cm{0.0};
  };
  struct TrackRow : Identity
  {
    unsigned int track_id{0}, source_assembled_track_id{0}, nclusters{0};
    int charge{0}, poly_valid{0}, poly_fit_status{0}, kalman_success{0};
    int crossing_valid{0}, crossing{0}, crossing_status{-1}, crossing_tier{-1};
    double crossing_score{0.0}, poly_chi2{0.0}, poly_ndof{0.0};
    std::array<double, 3> poly_position{}, poly_momentum{};
    std::array<double, 36> poly_covariance{};
    std::array<double, 6> kalman_endpoint_state{}, kalman_pca_state{}, kalman_inner_state{};
    std::array<double, 36> kalman_endpoint_covariance{}, kalman_pca_covariance{}, kalman_inner_covariance{};
    int pca_success{0}, inner_success{0};
    double pca_path_cm{0.0}, inner_path_cm{0.0}, inner_reference_radius_cm{0.0};
    unsigned int inner_intersections{0}, naccepted{0}, nrejected{0};
    double chi2{0.0};
    int ndof{0};
    std::string fit_message, pca_message, inner_message;
    std::vector<std::uint64_t> cluster_ids;
    std::vector<unsigned char> measurement_used, measurement_rejection_reason;
    std::vector<double> measurement_chi2_raw, measurement_weight_scale;
  };
  struct ClusterRow : Identity
  {
    std::uint64_t cluster_id{0}, trkr_cluster_key{0};
    unsigned int source_assembled_track_id{0}, n_associated_tracks{0};
    int valid{0}, layer{-1}, sector{-1}, side{-1}, geometry_valid{0};
    double x{0.0}, y{0.0}, z{0.0}, layer_radius_cm{0.0}, drift_length_cm{0.0};
    std::array<double, 3> world_position{};
    // Rows are local envelope (x,y,z) axes expressed in world coordinates.
    std::array<double, 9> envelope_basis_world{};
    std::array<double, 3> envelope_origin_world{};
    double adc{0.0}, rms_x{0.0}, rms_y{0.0}, rms_z{0.0};
    unsigned int nhits{0}, phi_width{0}, time_width{0};
    std::vector<unsigned int> hitset_keys, hit_keys;
  };
  struct ResidualRow : Identity
  {
    unsigned int track_id{0};
    std::uint64_t cluster_id{0}, trkr_cluster_key{0};
    int residual_type{0};  // 1 included smoother; 2 independent held-out fit.
    int layer{-1}, sector{-1}, side{-1}, charge{0};
    int success{0}, nominal_used{0}, nominal_rejection_reason{-1}, pull_valid{0};
    int independent_fit_charge{0}, independent_charge_ambiguous{0};
    unsigned int excluded_measurements{0}, excluded_shared_raw_hit_measurements{0}, intersections{0};
    double layer_radius_cm{0.0}, path_cm{0.0}, tangent_angle{0.0}, longitudinal_angle{0.0};
    std::array<double, 3> measurement{}, residual{}, pull{};
    std::array<double, 6> prediction_state{};
    std::array<double, 36> prediction_covariance{};
    std::array<double, 9> measurement_covariance_local{}, prediction_covariance_local{}, residual_covariance_local{};
    std::vector<std::uint64_t> retained_cluster_ids;
    std::string message;
  };
  struct TruthReferenceRow : Identity
  {
    std::uint64_t cluster_id{0};  // zero for generated-primary reference row.
    int reference_type{0};  // 1 primary at inner cylinder; 2 cluster-associated truth.
    int truth_track_id{0}, pid{0}, is_primary{0}, layer{-1}, success{0};
    unsigned int shared_raw_hits{0}, ambiguous_raw_hits{0}, intersections{0};
    int interpolation_method{-1};  // 0 exact endpoint; 1 G4-segment interpolation.
    std::uint64_t g4hit_key{0};
    double reference_radius_cm{0.0}, interpolation_span_cm{0.0}, time_ns{0.0};
    std::array<double, 3> position{}, momentum{}, residual{};
    std::string message;
  };

  void create_trees();
  void fill_residual(const TpcTrackDiagnosticMeasurement &, const ClusterRow &,
                     unsigned int, int, int, bool, int,
                     const TpcKalmanSurfaceResult &, const TpcTrackLeaveOneOutResult *);
  int fill_truth(PHCompositeNode *, const std::vector<ClusterRow> &);
  Identity identity() const;

  std::string m_filename, m_provenance;
  std::string m_cluster_node{"TPC_POLYCLUSTERS"}, m_track_node{"TPC_POLYTRACKS"};
  std::uint64_t m_sample_id{0}, m_entry{0};
  bool m_write_leave_one_out{true}, m_write_truth{false}, m_require_field_map{true};
  double m_inner_radius_cm{0.0}, m_readout_z_cm{105.5};
  TpcKalmanConfig m_config{TpcTrackDiagnosticFitter::yuri_config(nullptr)};
  TFile *m_file{nullptr};
  TTree *m_event_tree{nullptr}, *m_track_tree{nullptr}, *m_cluster_tree{nullptr};
  TTree *m_residual_tree{nullptr}, *m_truth_reference_tree{nullptr};
  EventRow m_event;
  TrackRow m_track;
  ClusterRow m_cluster;
  ResidualRow m_residual;
  TruthReferenceRow m_truth_reference;
};

#endif
