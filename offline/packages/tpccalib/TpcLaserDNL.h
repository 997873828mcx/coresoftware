#ifndef TPCCALIB_TPCLASERDNL_H
#define TPCCALIB_TPCLASERDNL_H
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <string>
#include <memory>
#include <vector>
#include <limits>
#include <map>
#include <Rtypes.h>

class SvtxTrackMap;
class TrkrHitSetContainer;
class PHG4TpcGeomContainer;
class ActsGeometry;
class TrkrClusterContainer;
class TrkrCluster;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class PHG4HitContainer;
class TFile;
class TTree;

class TpcLaserDNL : public SubsysReco
{
public:
  explicit TpcLaserDNL(const std::string& name = "TpcLaserDNL");
  ~TpcLaserDNL() override = default;

int Init(PHCompositeNode* topNode) override;
int InitRun(PHCompositeNode* topNode) override;
int process_event(PHCompositeNode* topNode) override;
int End(PHCompositeNode* topNode) override;

// config
void set_outputfile(const std::string& fn) { m_outfile = fn; }
void set_max_dca(double v) { m_max_dca = v; }
void set_max_dz(double v) { m_max_dz = v; }
  void set_pedestal(double v) { m_pedestal = v; }
  // threshold on the chosen weight (ADC or charge)
  void set_min_adc(double v) { m_min_adc = v; }
  void set_use_pedestal(bool v) { m_use_pedestal = v; }
  // choose weighting mode: true = ADC (default), false = charge (hit energy)
  void set_weight_by_adc(bool v) { m_weight_by_adc = v; }
  // choose source: false = use hits (default), true = use clusters
  void set_use_clusters(bool v) { m_use_clusters = v; }
  // when true, enforce single-hitset (single-sector) contribution per layer result
  void set_restrict_to_single_hitset(bool v) { m_restrict_to_single_hitset = v; }
  // allow disabling use of SvtxTrackMap seeds
  void set_use_reco_seeds(bool v) { m_use_reco_seeds = v; }
  // enable dumping intersection/hit info in event-display format
  void set_write_display_ntuple(bool v) { m_write_display_ntuple = v; }
  // control how many samples are taken along each G4 hit segment for the display
  void set_display_hit_subsamples(unsigned int v) { m_display_hit_subsamples = v; }
  // optionally include extrapolated intersections (no G4 support) in the display output
  void set_include_fallback_intersections(bool v) { m_include_fallback_intersections = v; }
  // enable/disable fitted-layer-44 residual calculation
  void set_enable_fit_layer44_residuals(bool v) { m_enable_fit_layer44_residuals = v; }
  // when true, only use reco hits/clusters that are truth-associated to
  // positive truth-track IDs (trkid > 0)
  void set_primary_hits_only(bool v) { m_primary_hits_only = v; }

private:
// nodes
  SvtxTrackMap* m_track_map{nullptr};
  TrkrHitSetContainer* m_hitsets{nullptr};
  PHG4TpcGeomContainer* m_geom{nullptr};
  ActsGeometry* m_acts{nullptr};
  TrkrClusterContainer* m_clusters{nullptr};

// io
std::string m_outfile{"laser_dnl.root"};
std::unique_ptr<TFile> m_tf;
TTree* m_tt{nullptr};
TTree* m_tt_display_intersections{nullptr};
TTree* m_tt_display_g4hits{nullptr};

// cuts
double m_max_dca{0.3}; // cm
double m_max_dz{1.0}; // cm
  double m_pedestal{74.4};
  double m_min_adc{1.0};
  bool m_use_pedestal{true};
  bool m_weight_by_adc{true};
  bool m_use_clusters{false};
  bool m_restrict_to_single_hitset{false};
  bool m_use_reco_seeds{true};
  bool m_write_display_ntuple{false};
  unsigned int m_display_hit_subsamples{0};
  bool m_include_fallback_intersections{false};
  bool m_enable_fit_layer44_residuals{true};
  bool m_primary_hits_only{true};

  struct LayerPoint
  {
    unsigned int layer{0};
    double radius{0.};
    double x{0.};
    double y{0.};
    double z{0.};
    double dirx{0.};
    double diry{0.};
    double dirz{0.};
    int side{0};
    double path{0.};
    bool from_g4hit{false};
  };

  struct TrackSeed
  {
    int id{0};
    double pt{std::numeric_limits<double>::quiet_NaN()};
    double origin[3]{0., 0., 0.};
    double dir[3]{0., 0., 0.};
    bool dir_valid{false};
    std::vector<LayerPoint> layers;
  };

// tree vars
int m_event{0};
int m_trkid{0};
double m_pt{std::numeric_limits<double>::quiet_NaN()};
unsigned int m_layer{0};
int m_side{0};
int m_sector{-1};
double m_r{0};
double m_phi_true{0};
double m_phi_reco{0};
double m_dphi{0};
double m_dRphi{0};
int m_nused{0};
int m_nhit_scanned{0};
double m_adcsum{0};
double m_xtrue{0}, m_ytrue{0}, m_ztrue{0};
double m_xreco{0}, m_yreco{0}, m_zreco{0};
int m_npad_used{0};
  int m_ntbin_used{0};
  int m_nbins_used{0};
double m_phi_pad_max{std::numeric_limits<double>::quiet_NaN()};
double m_phase{std::numeric_limits<double>::quiet_NaN()};
double m_phase_reco{std::numeric_limits<double>::quiet_NaN()};
std::vector<double> m_pad_phi_centers;
  // fitted straight-line intersection with layer 44 (NaN for other layers or failed fits)
  double m_xfit_layer44{std::numeric_limits<double>::quiet_NaN()};
  double m_yfit_layer44{std::numeric_limits<double>::quiet_NaN()};
  double m_zfit_layer44{std::numeric_limits<double>::quiet_NaN()};
  double m_phi_fit_layer44{std::numeric_limits<double>::quiet_NaN()};
  double m_dphi_fit_layer44{std::numeric_limits<double>::quiet_NaN()};
  double m_dRphi_fit_layer44{std::numeric_limits<double>::quiet_NaN()};

// debug vectors: store per-hit info for used hits
std::vector<ULong64_t> m_hitkeys;
std::vector<ULong64_t> m_hitsetkeys;
  std::vector<unsigned int> m_iphi;
  std::vector<unsigned int> m_tbin;
  std::vector<double> m_hit_charge;
  double m_total_charge_layer{0.0};
  double m_max_charge_layer{0.0};

// helpers
static bool cylinder_intersection(double x0,double y0,double z0,
double vx,double vy,double vz,
double R, double& t_out,
double& xi,double& yi,double& zi);
static double wrap_dphi(double d);

  void build_reco_seeds(std::vector<TrackSeed>& seeds) const;
  void build_truth_seeds(std::vector<TrackSeed>& seeds) const;

  TrkrHitTruthAssoc* m_hittruthassoc{nullptr};
  TrkrClusterHitAssoc* m_cluster_hit_assoc{nullptr};
  PHG4HitContainer* m_g4hits_tpc{nullptr};
  PHG4HitContainer* m_g4hits{nullptr};

  struct DisplayIntersection
  {
    int event{0};
    int trackid{0};
    int layer{0};
    int side{0};
    double gx{0.};
    double gy{0.};
    double gz{0.};
    double r{0.};
    double phi{0.};
    double path{0.};
    int used_in_seed{0};
  };

  struct DisplayG4Hit
  {
    int event{0};
    int trackid{0};
    int layer{0};
    int side{0};
    double gx{0.};
    double gy{0.};
    double gz{0.};
    double r{0.};
    double phi{0.};
    int sample{0};
    double sample_frac{0.};
    double step_path{0.};
    double edep{0.};
    double eion{0.};
    double t0{0.};
    double t1{0.};
    int used_in_track{0};
    ULong64_t hitid{0};
  };

  DisplayIntersection m_display_intersection;
  DisplayG4Hit m_display_g4hit;
};

#endif
