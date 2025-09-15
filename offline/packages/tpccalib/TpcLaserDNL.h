#ifndef TPCCALIB_TPCLASERDNL_H
#define TPCCALIB_TPCLASERDNL_H
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <string>
#include <memory>
#include <vector>
#include <Rtypes.h>

class SvtxTrackMap;
class TrkrHitSetContainer;
class PHG4TpcCylinderGeomContainer;
class ActsGeometry;
class TrkrClusterContainer;
class TrkrCluster;
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

private:
// nodes
  SvtxTrackMap* m_track_map{nullptr};
  TrkrHitSetContainer* m_hitsets{nullptr};
  PHG4TpcCylinderGeomContainer* m_geom{nullptr};
  ActsGeometry* m_acts{nullptr};
  TrkrClusterContainer* m_clusters{nullptr};

// io
std::string m_outfile{"laser_dnl.root"};
std::unique_ptr<TFile> m_tf;
TTree* m_tt{nullptr};

// cuts
double m_max_dca{0.3}; // cm
double m_max_dz{1.0}; // cm
  double m_pedestal{74.4};
  double m_min_adc{1.0};
  bool m_use_pedestal{true};
  bool m_weight_by_adc{true};
  bool m_use_clusters{false};

// tree vars
int m_event{0};
int m_trkid{0};
unsigned int m_layer{0};
int m_side{0};
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

// debug vectors: store per-hit info for used hits
std::vector<ULong64_t> m_hitkeys;
std::vector<ULong64_t> m_hitsetkeys;
std::vector<unsigned int> m_iphi;
std::vector<unsigned int> m_tbin;

// helpers
static bool cylinder_intersection(double x0,double y0,double z0,
double vx,double vy,double vz,
double R, double& t_out,
double& xi,double& yi,double& zi);
static double wrap_dphi(double d);
};

#endif
