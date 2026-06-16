#ifndef KSHORTRECONSTRUCTION_H
#define KSHORTRECONSTRUCTION_H

#include <TTree.h>   

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <Acts/Definitions/Algebra.hpp>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Utilities/Result.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <ActsExamples/EventData/Trajectories.hpp>
#pragma GCC diagnostic pop

#include <Eigen/Dense>

class TFile;
class TH1D;
class TNtuple;

using BoundTrackParam = const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;

class ActsGeometry;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class GlobalVertexMap;

class KshortReconstruction : public SubsysReco
{
 public:
  KshortReconstruction(const std::string& name = "KshortReconstruction");
  virtual ~KshortReconstruction() {}

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void setPtCut(double ptcut) { invariant_pt_cut = ptcut; }
  void setTrackPtCut(double ptcut) { track_pt_cut = ptcut; }
  void setTrackQualityCut(double cut) { _qual_cut = cut; }
  void setPairDCACut(double cut) { pair_dca_cut = cut; }
  void setTrackDCACut(double cut) { track_dca_cut = cut; }
  void setRequireMVTX(bool set) { _require_mvtx = set; }
  void setDecayMass(Float_t decayMassSet) { decaymass = decayMassSet; }  //(muons decaymass = 0.1057) (pions = 0.13957) (electron = 0.000511)
  void set_output_file(const std::string& outputfile) { filepath = outputfile; }
  void save_tracks(bool save = true) { m_save_tracks = save; }
  void setApplyQualityCut(bool apply)
  {
    apply_quality_cut = apply;
  }

  void setApplyDCACut(bool apply)
  {
    apply_dca_cut = apply;
  }

  void setApplyPairDCACut(bool apply)
  {
    apply_pair_dca_cut = apply;
  }

  void setApplyInvariantPtCut(bool apply)
  {
    apply_invariant_pt_cut = apply;
  }
  void setApplyTrackPtCut(bool apply)
  {
    apply_track_pt_cut = apply;
  }
  void setApplyChargeCut(bool apply)
  {
    apply_charge_cut = apply;
  }
  

 private:
  //void fillNtp(SvtxTrack* track1, SvtxTrack* track2, Acts::Vector3 dcavals1, Acts::Vector3 dcavals2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, float invariantPhi, float rapidity, float pseudorapidity, Eigen::Vector3d projected_pos1, Eigen::Vector3d projected_pos2, Eigen::Vector3d projected_mom1, Eigen::Vector3d projected_mom2, Acts::Vector3 pca_rel1_proj, Acts::Vector3 pca_rel2_proj, double pair_dca_proj, unsigned int track1_silicon_cluster_size, unsigned int track2_silicon_cluster_size, unsigned int track1_mvtx_cluster_size, unsigned int track1_mvtx_state_size, unsigned int track1_intt_cluster_size, unsigned int track1_intt_state_size, unsigned int track2_mvtx_cluster_size, unsigned int track2_mvtx_state_size, unsigned int track2_intt_cluster_size, unsigned int track2_intt_state_size, int runNumber, int eventNumber);

  

  // void findPcaTwoTracks(SvtxTrack *track1, SvtxTrack *track2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca);
  void findPcaTwoTracks(const Acts::Vector3& pos1, const Acts::Vector3& pos2, Acts::Vector3 mom1, Acts::Vector3 mom2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca);

  int getNodes(PHCompositeNode* topNode);

  Acts::Vector3 calculateDca(SvtxTrack* track, const Acts::Vector3& momentum, Acts::Vector3 position);

  bool projectTrackToCylinder(SvtxTrack* track, double Radius, Eigen::Vector3d& pos, Eigen::Vector3d& mom);
  bool projectTrackToPoint(SvtxTrack* track, Eigen::Vector3d PCA, Eigen::Vector3d& pos, Eigen::Vector3d& mom);

  Acts::Vector3 getVertex(SvtxTrack* track);
  std::vector<unsigned int> getTrackStates(SvtxTrack* track);

  //TNtuple* ntp_reco_info = nullptr;

  TTree* m_pairTree = nullptr;
  ActsGeometry* _tGeometry = nullptr;
  SvtxTrackMap* m_svtxTrackMap = nullptr;
  SvtxVertexMap* m_vertexMap = nullptr;
  // GlobalVertexMap* m_globalvertexMap = nullptr;

  std::string filepath = "";
  Float_t decaymass = 0.13957;  // pion decay mass
  bool _require_mvtx = true;
  bool apply_quality_cut = true;
  bool apply_dca_cut = true;
  bool apply_pair_dca_cut = true;
  bool apply_invariant_pt_cut = true;
  bool apply_track_pt_cut = true;
  bool apply_charge_cut = true;
  double _qual_cut = 1000.0;
  double pair_dca_cut = 0.15;  // kshort relative cut 500 microns
  double track_dca_cut = 0.01;
  double invariant_pt_cut = 0.1;
  double track_pt_cut = 0.2;
  TFile* fout = nullptr;
  //TH1D* recomass = nullptr;
  //TH2D* h_AP_all = nullptr;
  bool m_save_tracks = false;
  SvtxTrackMap* m_output_trackMap = nullptr;
  std::string m_output_trackMap_node_name = "KshortReconstruction_SvtxTrackMap";

  struct Pair_t
{
  Int_t   run   = 0;
  Int_t   evt   = 0;

 

  Short_t cross1 = -1;
  Short_t cross2 = -1;

  // kinematics
  Float_t px1 = 0, py1 = 0, pz1 = 0;
  Float_t px2 = 0, py2 = 0, pz2 = 0;

  // DCA
  Float_t dca_xy1 = 0, dca_z1 = 0;
  Float_t dca_xy2 = 0, dca_z2 = 0;
  Float_t pairDCA = 0;     // <<< name fixed (was pairDCA in struct, pair_dca in code)

  // A-P
  Float_t alpha = 0, qT = 0;

  // invariant
  //Float_t invMass = 0, invPt = 0, invPhi = 0;
  

  // decay length
  Float_t L = 0;           // (optional – not filled yet)
  Float_t Lproj = 0;       

  // charges & pointing
  Float_t charge1 = 0;
  Float_t charge2 = 0;
  Float_t cosThetaReco = 0;
};
Pair_t m_pair;   // member of KshortReconstruction
};

#endif  // KSHORTRECONSTRUCTION_H
