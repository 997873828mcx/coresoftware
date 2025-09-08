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

class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;
class PHG4TpcCylinderGeom;
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

  void SetDriftVelocity(double vd) override { drift_velocity = vd; }
  void SetReadoutTime(float t) override { extended_readout_time = t; }
  // otherwise warning of inconsistent overload since only one MapToPadPlane methow is overridden
  using PHG4TpcPadPlane::MapToPadPlane;

  void MapToPadPlane(TpcClusterBuilder &tpc_clustbuilder, TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc * /*hittruthassoc*/, const double x_gem, const double y_gem, const double t_gem, const unsigned int side, PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/ ) override;
  //void MapToPadPlane(TpcClusterBuilder &tpc_clustbuilder, TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc * /*hittruthassoc*/, const double x_gem, const double y_gem, const double t_gem, const unsigned int side, PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/, TH2* h_adc_ref, TH2* h_adc_serf  ) override;

  void SetDefaultParameters() override;
  void UpdateInternalParameters() override;
   

 private:

  //  void populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void SERF_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double rad_gem, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_tbins(const double t, const std::array<double, 2> &cloud_sig_tt, std::vector<int> &adc_tbin, std::vector<double> &adc_tbin_share);

  double check_phi(const unsigned int side, const double phi, const double radius);

  // utility: pick layers whose annulus intersects a radial window around rad
  std::vector<unsigned int> layersInRadialWindow(double rad, double sigma, double nsig) const;

  // utility: find geometry for a given layer
  PHG4TpcCylinderGeom* getGeomForLayer(unsigned int layer) const;

  // utility: determine sector for (x,y) and rotate to a canonical frame
  void rotatePointToSector(double x, double y, unsigned int side, int& sectorFound, double& xNew, double& yNew);

  PHG4TpcCylinderGeomContainer *GeomContainer = nullptr;
  PHG4TpcCylinderGeom *LayerGeom = nullptr;

  double neffelectrons_threshold = std::numeric_limits<double>::signaling_NaN();

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
  static constexpr double _nsigmas = 5;

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
  
std::array<std::vector<PadInfo>,3*16+7> Pads;
bool pointInPolygon( double x, double y,const std::vector<Point>& poly); 
  double integratedDensityOfCircleAndPad(double hitX,double hitY, double sigma, const std::vector<Point>& pad,double gridStep = 0.0);
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



};

#endif
