#ifndef G4TPC_PHG4TPCPADPLANEREADOUT_H
#define G4TPC_PHG4TPCPADPLANEREADOUT_H

#include "PHG4TpcPadPlane.h"
#include "TpcClusterBuilder.h"

#include <g4main/PHG4HitContainer.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <climits>
#include <cmath>
#include <string>  // for string
#include <vector>
#include <map>
typedef std::map<TrkrDefs::hitsetkey, std::vector<TrkrDefs::hitkey>> hitMaskTpc;

class PHCompositeNode;
class PHG4TpcGeomContainer;
class PHG4TpcGeom;
class TH2;
class TF1;
class TNtuple;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcPadPlaneReadout : public PHG4TpcPadPlane
{
 public:
  PHG4TpcPadPlaneReadout(const std::string &name = "PHG4TpcPadPlaneReadout");

  ~PHG4TpcPadPlaneReadout() override;

  int InitRun(PHCompositeNode *topNode) override;

  void UseGain(const int flagToUseGain);
  void SetUseModuleGainWeights(const int flag) {m_use_module_gain_weights = flag;}
  void SetModuleGainWeightsFileName(const std::string &name) {m_tpc_module_gain_weights_file = name;}
  void ReadGain();
  void SetUsePolyaGEMGain(const int flagPolya) {m_usePolya = flagPolya;}
  void SetUseLangauGEMGain(const int flagLangau) {m_useLangau = flagLangau;}
  void SetLangauParsFileName(const std::string &name) {m_tpc_langau_pars_file = name;}
  // Pad-sharing method selection: true → SERF polygon overlap; false → analytic triangle
  void UseSerfPadSharing(bool use_serf) { m_use_serf_padsharing = use_serf; }
  // If true and SERF polygons are unavailable, abort InitRun with an error
  void RequireSerfPadSharing(bool require) { m_require_serf = require; }
  // If true, also require polygons for all readout layers (partial coverage aborts)
  void RequireSerfFullCoverage(bool require) { m_require_serf_full = require; }
  // When using analytic sharing (no SERF polygons), choose rectangular response instead of triangular
  void UseRectangularPadResponse(bool use_rectangular) { m_use_rectangular_pad_response = use_rectangular; }

  void SetDriftVelocity(double vd) override { drift_velocity = vd; }
  void SetReadoutTime(float t) override { extended_readout_time = t; }
  // otherwise warning of inconsistent overload since only one MapToPadPlane methow is overridden
  using PHG4TpcPadPlane::MapToPadPlane;

  void MapToPadPlane(TpcClusterBuilder &tpc_clustbuilder, TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc * /*hittruthassoc*/, const double x_gem, const double y_gem, const double t_gem, const unsigned int side, PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/ ) override;
  //void MapToPadPlane(TpcClusterBuilder &tpc_clustbuilder, TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc * /*hittruthassoc*/, const double x_gem, const double y_gem, const double t_gem, const unsigned int side, PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/, TH2* h_adc_ref, TH2* h_adc_serf  ) override;

  void SetDefaultParameters() override;
  void UpdateInternalParameters() override;
void SetMaskChannelsFromFile() 
  {
    m_maskFromFile = true;
  } 
  void SetDeadChannelMapName(const std::string& dcmap) 
  {
    m_maskDeadChannels = true;
    m_deadChannelMapName = dcmap;
  }
  void SetHotChannelMapName(const std::string& hmap) 
  {
    m_maskHotChannels = true;
    m_hotChannelMapName = hmap;
  }
  // Debug printing helpers
  // If set >= 0, limit PadHit prints to a single layer number; otherwise prints for all layers.
  void SetDebugPadHitLayer(int layer) { m_dbg_pad_hit_layer = layer; }
  // Enable a one-shot visualization of a single avalanche cloud overlap with zigzag pads.
  // Passing target_side/target_layer < 0 matches the first cloud encountered.
  // grid_step <= 0 defaults to sigma/30 sampling.
  void EnableSingleCloudVisualization(bool enable,
                                      const std::string &output_file = "AvalancheCloudOverlap.png",
                                      int target_side = -1,
                                      int target_layer = -1,
                                      double grid_step = -1.0);
  void SetVisualizationDumpFile(const std::string &file);
  void SetVisualizeAllClouds(bool enable);

  protected: 
  double Ts = 80.0; // SAMPA peaking time

  double sampaShapingResponseFunction(double tzero, double t) const;

  void sampaTimeDistribution(double tzero, 
                             std::vector<int> &adc_tbin, 
                             std::vector<double> &adc_tbin_share);

 private:

  //  void populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void SERF_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double rad_gem, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  //void populate_tbins(const double t, const std::array<double, 2> &cloud_sig_tt, std::vector<int> &adc_tbin, std::vector<double> &adc_tbin_share);

  double check_phi(const unsigned int side, const double phi, const double radius);
void makeChannelMask(hitMaskTpc& aMask, const std::string& dbName, const std::string& totalChannelsToMask);
  // utility: pick layers whose annulus intersects a radial window around rad
  std::vector<unsigned int> layersInRadialWindow(double rad, double sigma, double nsig) const;

  // utility: find geometry for a given layer
  PHG4TpcGeom* getGeomForLayer(unsigned int layer) const;

  // utility: determine sector for (x,y) and rotate to a canonical frame
  void rotatePointToSector(double x, double y, unsigned int side, int& sectorFound, double& xNew, double& yNew);

  PHG4TpcGeomContainer *GeomContainer = nullptr;
  PHG4TpcGeom *LayerGeom = nullptr;

  // Minimum effective electrons per (pad,tbin) to create a hit.
  // Default to 0.0 so all contributions are kept unless configured otherwise.
  double neffelectrons_threshold = 0.0;

  std::array<double, 3> MinRadius{};
  std::array<double, 3> MaxRadius{};

  static constexpr int NSides = 2;
  static constexpr int NSectors = 12;
  static const int NRSectors = 3;

  double sigmaT = std::numeric_limits<double>::signaling_NaN();
  std::array<double, 2> sigmaL{};
  double phi_bin_width{};
  double drift_velocity = 8.0e-03;  // default value, override from macro
  float extended_readout_time = 0;  // ns
  int NTBins = std::numeric_limits<int>::max();
  int m_NHits = 0;
  // Using Gain maps is turned off by default
  int m_flagToUseGain = 0;

  // Optionally apply a module-by-module weight to the GEM gain
  // Weights are input from a file for all 72 TPC modules
  bool m_use_module_gain_weights = false;
  std::string m_tpc_module_gain_weights_file = "";

  // gaussian sampling
  static constexpr double _nsigmas = 2.5;

  double averageGEMGain = std::numeric_limits<double>::signaling_NaN();
  double polyaTheta = std::numeric_limits<double>::signaling_NaN();

  std::array<std::vector<double>, NSides> sector_min_Phi;
  std::array<std::vector<double>, NSides> sector_max_Phi;

  // return random distribution of number of electrons after amplification of GEM for each initial ionizing electron
  double getSingleEGEMAmplification();
  double getSingleEGEMAmplification(double weight);
  double getSingleEGEMAmplification(TF1 *f);
  bool m_usePolya = false;

  bool m_useLangau = false;
  std::string m_tpc_langau_pars_file = "";

  gsl_rng *RandomGenerator = nullptr;

  std::array<TH2 *, 2> h_gain{nullptr};

  double m_module_gain_weight[2][3][12] = { 
    { {1,1,1,1,1,1,1,1,1,1,1,1},
      {1,1,1,1,1,1,1,1,1,1,1,1},
      {1,1,1,1,1,1,1,1,1,1,1,1} },
    { {1,1,1,1,1,1,1,1,1,1,1,1},
      {1,1,1,1,1,1,1,1,1,1,1,1},
      {1,1,1,1,1,1,1,1,1,1,1,1} } 
  };

  TF1 *flangau[2][3][12] = {};
    hitMaskTpc m_deadChannelMap;
  hitMaskTpc m_hotChannelMap; 
bool m_maskDeadChannels {false};
  bool m_maskHotChannels {false};
  bool m_maskFromFile {false};
  std::string m_deadChannelMapName; 
  std::string m_hotChannelMapName; 
  struct Point { double x, y; };

  struct PadInfo {
    std::string          name;       // pad name
    int                  pad_number; // pad number (number in module)
    int                  pad_bin;    // pad phi bin (number according to get_phi_bin)
    double               cx, cy;     // centroid coords
    double               rad, phi;   // pad radius and phi
    std::vector<Point>   vertices;   // pad polygon
    bool                isedge = false; // whether to keep this pad signal
    void clear() {
      name.clear();
      pad_number = -1;
      pad_bin    = -1;
      cx = cy = rad = phi = 0.0;
      vertices.clear();
      isedge = false;
    }
  };
  
  struct DebugSample
  {
    double x = 0.0;
    double y = 0.0;
    double density = 0.0;
  };

std::array<std::vector<PadInfo>,3*16+7> Pads;
bool pointInPolygon( double x, double y,const std::vector<Point>& poly); 
  double integratedDensityOfCircleAndPad(double hitX,double hitY, double sigma, const std::vector<Point>& pad,double gridStep = 0.0, std::vector<DebugSample>* debug_samples = nullptr);
  // hard‑coded list of input .brd files
  static const std::vector<std::string> brdMaps_;
void loadPadPlanes();
int ntpc_phibins_sector[3] = { 94, 128, 192 };

int findPadForPoint( double x, double y, int tpc_module);
  const std::array<double, 5> Thickness =
      {{
          0.56598621677629212,
          1.0206889851687158,
          1.0970475085472556,
          0.5630547309825637,
          0.56891770257002054,
      }};
double min_radii_module[3]={314.9836110818037, 416.59202613529567, 589.1096495597712};
double max_radii_module[3]={399.85222874031024, 569.695373910603, 753.6667758418596};

  // choose between SERF polygon overlap (default) and triangle response
  bool m_use_serf_padsharing = true;
  bool m_require_serf = false;
  bool m_require_serf_full = false;
  bool m_serf_polygons_present = false; // summarized availability after InitRun
  bool m_warned_serf_fallback = false;  // printed once if per-hit fallback occurs
  // analytic sharing option: rectangular (flat) pad response instead of triangular
  bool m_use_rectangular_pad_response = false;

  // Debug: restrict PadHit prints to a single layer (or all if < 0)
  int m_dbg_pad_hit_layer = 49;

  struct DebugPadContribution
  {
    int pad_bin = -1;
    double charge = 0.0;
    std::vector<Point> polygon;
    double pad_phi = 0.0;
    double pad_r = 0.0;
  };
  struct VisualizationCircle
  {
    double x = 0.0;
    double y = 0.0;
    double radius = 0.0;
  };

  void maybeVisualizeAvalanche(unsigned int side,
                               unsigned int layernum,
                               double phi,
                               double rad_gem,
                               double cloud_sig_rp,
                               double x_center,
                               double y_center,
                               const std::vector<DebugPadContribution> &contribs,
                               const std::vector<DebugSample> &samples,
                               const std::vector<VisualizationCircle> &circles);

  bool m_visualize_single_cloud = false;
  bool m_visualization_done = false;
  int  m_visualization_target_layer = -1;
  int  m_visualization_target_side = -1;
  std::string m_visualization_output = "AvalancheCloudOverlap.png";
  double m_visualization_grid_step = -1.0;
  bool m_visualize_all_matches = false;
  std::string m_visualization_dump_file;
  unsigned long m_visualization_dump_index = 0;
  unsigned long m_visualization_cloud_counter = 0;
  std::vector<DebugSample> m_visualization_aggregate_samples;
  std::vector<VisualizationCircle> m_visualization_circles;
  std::map<int, DebugPadContribution> m_visualization_pad_union;


};

#endif
