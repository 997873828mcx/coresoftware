// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H

#include "TpcClusterBuilder.h"

#include <trackbase/ActsGeometry.h>

#include <g4main/PHG4HitContainer.h>

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <cstddef>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class PHG4TpcPadPlane;
class PHG4TpcDistortion;
class PHCompositeNode;
class TH1;
class TH2;
class TNtuple;
class TTree;
class TFile;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class TrkrTruthTrackContainer;
class TrkrClusterContainer;
class TrkrTruthTrack;
class DistortedTrackContainer;
class TpcClusterBuilder;
class PHG4TpcGeomContainer;
class SvtxTrackMap;
class ClusHitsVerbose;
class ActsGeometry;

class PHG4TpcElectronDrift : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcElectronDrift(const std::string &name = "PHG4TpcElectronDrift");
  ~PHG4TpcElectronDrift() override = default;
  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode *) override;

  void SetDefaultParameters() override;

  //! detector name
  void Detector(const std::string &d)
  {
    detector = d;
  }

  //! detector name
  const std::string &Detector() const
  {
    return detector;
  }

  //! random seed
  void set_seed(const unsigned int iseed);

  //! setup TPC distortion
  void setTpcDistortion(PHG4TpcDistortion *);

  //
  void set_flag_threshold_distortion(bool setflag, float setthreshold);

  //! setup readout plane
  void registerPadPlane(PHG4TpcPadPlane *padplane);

  // cluster the PHG4Tracks individually
  TpcClusterBuilder truth_clusterer{};
  void set_pixel_thresholdrat(double val) { truth_clusterer.set_pixel_thresholdrat(val); };
  void set_max_g4hitstep(float f) { max_g4hitstep = f; };
  void set_ClusHitsVerbose(bool set = true) { record_ClusHitsVerbose = set; };
  void set_zero_bfield_flag(bool flag) { zero_bfield = flag; };
  void set_zero_bfield_diffusion_factor(double f) { zero_bfield_diffusion_factor = f; };
  void use_PDG_gas_params() { m_use_PDG_gas_params = true; }
  void set_force_min_trans_drift_length(double length)
  {
    force_min_trans_drift_length = (length > 0.) ? length : 0.;
    set_double_param("force_min_trans_drift_length", force_min_trans_drift_length);
  }
  void set_enable_cluster_size_fluctuations(bool b) { m_enable_cluster_size_fluctuations = b; }
  // Backward-compatible alias; use set_enable_cluster_size_fluctuations instead.
  void set_enable_laser_clustering(bool b) { set_enable_cluster_size_fluctuations(b); }
  void set_use_pai_cluster_seeds(bool b)
  {
    m_use_pai_cluster_seeds = b;
    set_int_param("use_pai_cluster_seeds", b ? 1 : 0);
  }
  void set_pai_cluster_node_name(const std::string &name) { m_pai_cluster_node_name = name; }
  void set_qa_output_file(const std::string &f) { m_qa_output_file = f; }
  void set_step_qa_output_file(const std::string &f) { m_step_qa_output_file = f; }
  void set_ionization_qa_output_file(const std::string &f) { m_ionization_qa_output_file = f; }
  void set_write_legacy_qa_histograms(bool b) { m_write_legacy_qa_histograms = b; }
  void set_qa_write_only_with_secondaries(bool b) { m_qa_write_only_with_secondaries = b; }
  void set_avg_output_file(const std::string &f) { m_avg_output_file = f; }
  void set_do_ElectronDriftQAHistos(bool b) { do_ElectronDriftQAHistos = b; }
  ClusHitsVerbosev1 *mClusHitsVerbose{nullptr};

 private:
  TrkrHitSetContainer *hitsetcontainer{nullptr};
  TrkrHitTruthAssoc *hittruthassoc{nullptr};
  TrkrTruthTrackContainer *truthtracks{nullptr};
  TrkrTruthTrack *truth_track{nullptr};
  TrkrClusterContainer *truthclustercontainer{nullptr};  // the TrkrClusterContainer for truth clusters
  ActsGeometry *m_tGeometry{nullptr};
  PHG4TpcGeomContainer *seggeo{nullptr};
  PHG4HitContainer *m_truth_intersection_hits{nullptr};
  SvtxTrackMap *m_track_map{nullptr};

  struct TpcLayerRadius
  {
    unsigned int layer{0};
    double radius{0.0};
  };
  std::vector<TpcLayerRadius> m_tpc_layer_radii;

  TNtuple *nt{nullptr};
  TNtuple *nthit{nullptr};
  TNtuple *ntfinalhit{nullptr};
  TNtuple *ntpad{nullptr};
  TTree *m_avgXResidualTree{nullptr};
  float m_avgTree_trackid{0};
  float m_avgTree_residual{0};
  float m_avgTree_avgx{0};
  float m_avgTree_xint{0};
  float m_avgTree_residual_primary{0};
  float m_avgTree_avgx_primary{0};
  float m_avgTree_count_primary{0};
  float m_avgTree_residual_end{0};
  float m_avgTree_avgx_end{0};
  float m_avgTree_count_end{0};
  float m_avgTree_residual_end_primary{0};
  float m_avgTree_avgx_end_primary{0};
  float m_avgTree_count_end_primary{0};
  float m_avgTree_residual_x{0};
  float m_avgTree_avgx_cart{0};
  float m_avgTree_residual_x_primary{0};
  float m_avgTree_avgx_cart_primary{0};
  float m_avgTree_residual_x_end{0};
  float m_avgTree_avgx_cart_end{0};
  float m_avgTree_residual_x_end_primary{0};
  float m_avgTree_avgx_cart_end_primary{0};

  ///@name evaluation histograms
  //@{
  TH1 *dlong{nullptr};
  TH1 *dtrans{nullptr};
  TH1 *ratioElectronsRR{nullptr};
  TH1 *diffDistance{nullptr};
  TH1 *diffDX{nullptr};
  TH1 *diffDY{nullptr};
  TH1 *nPrimaryClusters{nullptr};
  TH1 *primaryClusterMean{nullptr};
  TH1 *diffPerSqrtL{nullptr};
  TH1 *diffDXPerSqrtL{nullptr};
  TH1 *diffDYPerSqrtL{nullptr};
  TH1 *nPrimaryClustersPerCm{nullptr};
  TH1 *nIonizationElectrons{nullptr};
  TH1 *nIonizationElectronsPerCm{nullptr};
  TH1 *g4StepDedx{nullptr};
  TH1 *ionizedElectronsDedx{nullptr};
  TH1 *electronDensityProfile{nullptr};
  TH2 *hitmapstart{nullptr};
  TH2 *hitmapend{nullptr};
  TH2 *hitmapstart_z{nullptr};
  TH2 *hitmapend_z{nullptr};
  TH2 *z_startmap{nullptr};
  TH2 *deltaphi{nullptr};
  TH2 *deltar{nullptr};
  TH2 *deltaphinodiff{nullptr};
  TH2 *deltaRphinodiff{nullptr};
  TH2 *deltaphivsRnodiff{nullptr};
  TH2 *deltaphinodist{nullptr};
  TH2 *deltarnodiff{nullptr};
  TH2 *deltarnodist{nullptr};
  TH2 *deltaz{nullptr};
  TH2 *nPrimaryClustersVsMean{nullptr};
  TH2 *nIonizationElectronsVsPrimaryClusters{nullptr};
  TH2 *ionizedElectronsDedxVsG4StepDedx{nullptr};
  TH2 *diffVsDrift{nullptr};
  TH2 *diffPerSqrtLVsDrift{nullptr};
  TH2 *diffDXVsDrift{nullptr};
  TH2 *diffDYVsDrift{nullptr};
  TH2 *diffDXPerSqrtLVsDrift{nullptr};
  TH2 *diffDYPerSqrtLVsDrift{nullptr};
  TTree *driftXY{nullptr};
  TTree *driftStepQA{nullptr};
  TTree *ionizationEventQA{nullptr};
  TTree *ionizationClusterQA{nullptr};
  TTree *ionizationSegmentQA{nullptr};
  float m_drift_start_x{0.F};
  float m_drift_start_y{0.F};
  float m_drift_end_x{0.F};
  float m_drift_end_y{0.F};
  int m_stepqa_event{0};
  int m_stepqa_track_id{0};
  float m_stepqa_step_length_cm{0.F};
  float m_stepqa_edep_kev{0.F};
  float m_stepqa_eion_kev{0.F};
  float m_stepqa_dedx_kev_per_cm{0.F};
  float m_stepqa_deiondx_kev_per_cm{0.F};
  float m_stepqa_mean_primary_clusters{0.F};
  int m_stepqa_n_primary_clusters{0};
  float m_stepqa_avg_ionization_energy_kev_per_electron{0.F};
  float m_stepqa_mean_ionization_electrons{0.F};
  int m_stepqa_n_ionization_electrons{0};
  float m_stepqa_ionized_electrons_edep_kev{0.F};
  float m_stepqa_ionized_electrons_dedx_kev_per_cm{0.F};
  int m_eventqa_event{0};
  int m_eventqa_use_pai_cluster_seeds{0};
  int m_eventqa_source_hits{0};
  float m_eventqa_path_length_cm{0.F};
  int m_eventqa_n_primary_clusters{0};
  int m_eventqa_n_ionization_electrons{0};
  float m_eventqa_primary_clusters_per_cm{0.F};
  float m_eventqa_ionization_electrons_per_cm{0.F};
  int m_clusterqa_event{0};
  int m_clusterqa_use_pai_cluster_seeds{0};
  int m_clusterqa_track_id{0};
  int m_clusterqa_g4hit_id{0};
  int m_clusterqa_cluster_index{0};
  int m_clusterqa_segment_index{0};
  int m_clusterqa_cluster_size{0};
  float m_clusterqa_track_path_cm{0.F};
  float m_clusterqa_x_cm{0.F};
  float m_clusterqa_y_cm{0.F};
  float m_clusterqa_z_cm{0.F};
  float m_clusterqa_t_ns{0.F};
  int m_segmentqa_event{0};
  int m_segmentqa_use_pai_cluster_seeds{0};
  int m_segmentqa_track_id{0};
  int m_segmentqa_segment_index{0};
  int m_segmentqa_n_primary_clusters{0};
  int m_segmentqa_n_ionization_electrons{0};
  float m_segmentqa_segment_start_cm{0.F};
  float m_segmentqa_segment_length_cm{0.F};
  float m_segmentqa_primary_clusters_per_cm{0.F};
  float m_segmentqa_ionization_electrons_per_cm{0.F};
  //@}

  int event_num{0};

  float max_g4hitstep{7.};
  float thresholdforreachesreadout{0.5};

  double diffusion_trans = std::numeric_limits<double>::quiet_NaN();
  double added_smear_sigma_trans = std::numeric_limits<double>::quiet_NaN();
  double diffusion_long = std::numeric_limits<double>::quiet_NaN();
  double added_smear_sigma_long = std::numeric_limits<double>::quiet_NaN();
  double drift_velocity = std::numeric_limits<double>::quiet_NaN();
  double tpc_length = std::numeric_limits<double>::quiet_NaN();
  double electrons_per_gev = std::numeric_limits<double>::quiet_NaN();
  double primary_clusters_per_gev = std::numeric_limits<double>::quiet_NaN();
  double primary_clusters_per_cm = std::numeric_limits<double>::quiet_NaN();
  double ionization_electrons_per_cm = std::numeric_limits<double>::quiet_NaN();
  double ionization_electrons_per_gev = std::numeric_limits<double>::quiet_NaN();
  double avg_ionization_energy_kev_per_electron = std::numeric_limits<double>::quiet_NaN();
  double fixed_primary_clusters_per_cm = -1.0;
  double min_active_radius = std::numeric_limits<double>::quiet_NaN();
  double max_active_radius = std::numeric_limits<double>::quiet_NaN();
  double min_time = std::numeric_limits<double>::quiet_NaN();
  double max_time = std::numeric_limits<double>::quiet_NaN();
  double zero_bfield_diffusion_factor{3.5};  // at drift voltage of 400 V
  double force_min_trans_drift_length{0.0};
  bool m_density_enabled{false};
  int m_density_layer{-1};
  double m_density_bin_width_cm{0.05};
  double m_density_window_cm{0.5};
  double m_density_layer_radius{0.0};
  double m_density_layer_rlow{0.0};
  double m_density_layer_rhigh{0.0};
  double m_density_radial_margin_cm{0.0};
  double m_density_capture_rlow{0.0};
  double m_density_capture_rhigh{0.0};
  double m_density_hist_half_range{0.0};
  bool m_avg_x_enabled{false};
  int m_avg_layer{-1};
  double m_avg_layer_radius{0.0};
  double m_avg_layer_rlow{0.0};
  double m_avg_layer_rhigh{0.0};
  std::unordered_map<int, double> m_track_path_offset;
  std::unordered_map<int, bool> m_track_anchor_set;
  std::unordered_map<int, double> m_track_anchor_length;
  struct TrackLayerData
  {
    double sum_rphi{0.0};
    double sum_rphi_primary{0.0};
    double sum_rphi_end{0.0};
    double sum_rphi_end_primary{0.0};
    double sum_x{0.0};
    double sum_x_primary{0.0};
    double sum_x_end{0.0};
    double sum_x_end_primary{0.0};
    std::size_t count{0};
    std::size_t count_primary{0};
    std::size_t count_end{0};
    std::size_t count_end_primary{0};
    double base_x{0.0};
    double base_y{0.0};
    double base_z{0.0};
    double dir_x{0.0};
    double dir_y{0.0};
    double dir_z{0.0};
    bool have_line{false};
  };
  std::unordered_map<int, TrackLayerData> m_track_layer_data;
  bool m_uniform_density_test{false};
  std::vector<double> cluster_size_cdf;
  bool m_enable_cluster_size_fluctuations{false};
  bool m_use_pai_cluster_seeds{false};
  std::string m_pai_cluster_node_name{"G4HIT_TPC_PAI_CLUSTER"};
  std::string m_qa_output_file{"ElectronDriftQA.root"};
  std::string m_step_qa_output_file;
  std::string m_ionization_qa_output_file;
  std::string m_avg_output_file{"avgXResidual.root"};
  bool m_qa_write_only_with_secondaries{false};
  bool m_write_legacy_qa_histograms{true};
  bool m_seen_event_with_secondaries{false};

  bool record_ClusHitsVerbose{false};
  bool do_ElectronDriftQAHistos{true};
  bool do_getReachReadout{false};
  bool zero_bfield{false};
  bool m_use_PDG_gas_params{false};

  std::unique_ptr<TrkrHitSetContainer> temp_hitsetcontainer;
  std::unique_ptr<TrkrHitSetContainer> single_hitsetcontainer;
  std::unique_ptr<PHG4TpcPadPlane> padplane;
  std::unique_ptr<PHG4TpcDistortion> m_distortionMap;
  std::unique_ptr<TFile> m_outf;
  std::unique_ptr<TFile> EDrift_outf;
  std::unique_ptr<TFile> m_stepQAOutf;
  std::unique_ptr<TFile> m_ionizationQAOutf;
  std::unique_ptr<TFile> m_avgOutf;

  std::string detector;
  std::string hitnodename;
  std::string seggeonodename;

  //! rng de-allocator
  class Deleter
  {
   public:
    //! deletion operator
    void operator()(gsl_rng *rng) const { gsl_rng_free(rng); }
  };
  std::unique_ptr<gsl_rng, Deleter> RandomGenerator;

  bool cylinder_intersection(double x0, double y0, double z0,
                             double vx, double vy, double vz,
                             double R, double &t_out,
                             double &xi, double &yi, double &zi) const;
};

#endif  // G4TPC_PHG4TPCELECTRONDRIFT_H
