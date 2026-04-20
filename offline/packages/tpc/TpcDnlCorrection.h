#ifndef TPC_TPCDNLCORRECTION_H
#define TPC_TPCDNLCORRECTION_H

/*!
 * \file TpcDnlCorrection.h
 * \brief applies per-cluster TPC DNL phi corrections from a calibration ROOT file
 */

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHG4TpcGeomContainer;
class TF1;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;

class TpcDnlCorrection : public SubsysReco
{
 public:
  explicit TpcDnlCorrection(const std::string& name = "TpcDnlCorrection");
  ~TpcDnlCorrection() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void set_correction_filename(const std::string& value)
  {
    m_correction_filename = value;
  }

  void set_allow_combined_side_fallback(bool value)
  {
    m_allow_combined_side_fallback = value;
  }

 private:
  struct CorrectionEntry
  {
    int source_side = -999;
    double phase_min = std::numeric_limits<double>::quiet_NaN();
    double phase_max = std::numeric_limits<double>::quiet_NaN();
    std::string fit_region;
    std::string func_name;
    std::unique_ptr<TF1> function;
  };

  using CorrectionKey = std::tuple<unsigned int, int, int>;
  using CorrectionMap = std::map<CorrectionKey, std::vector<CorrectionEntry>>;

  struct MatchResult
  {
    const CorrectionEntry* entry = nullptr;
    bool used_combined_side = false;
    bool have_any_for_key = false;
  };

  struct ClusterPhaseInputs
  {
    bool valid = false;
    int npads = 0;
    double r = std::numeric_limits<double>::quiet_NaN();
    double z = std::numeric_limits<double>::quiet_NaN();
    double phi_reco = std::numeric_limits<double>::quiet_NaN();
    double phi_pad_max = std::numeric_limits<double>::quiet_NaN();
    double phi_width = std::numeric_limits<double>::quiet_NaN();
    double phase_reco = std::numeric_limits<double>::quiet_NaN();
  };

  static double wrap_phi(double phi);
  MatchResult find_correction(unsigned int layer, int side, int npads, double phase_reco) const;
  ClusterPhaseInputs build_cluster_phase_inputs(TrkrDefs::cluskey, TrkrDefs::hitsetkey, TrkrCluster*) const;
  bool load_corrections();

  std::string m_correction_filename;
  bool m_allow_combined_side_fallback = true;
  CorrectionMap m_corrections;

  TrkrClusterContainer* m_clusters = nullptr;
  TrkrClusterHitAssoc* m_cluster_hit_assoc = nullptr;
  TrkrHitSetContainer* m_hitsets = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  PHG4TpcGeomContainer* m_geom = nullptr;

  std::uint64_t m_nclusters_total = 0;
  std::uint64_t m_ncorrected = 0;
  std::uint64_t m_ncorrected_combined = 0;
  std::uint64_t m_ninvalid_inputs = 0;
  std::uint64_t m_nmissing_key = 0;
  std::uint64_t m_nphase_gap = 0;
  std::uint64_t m_nsurface_failures = 0;
};

#endif
