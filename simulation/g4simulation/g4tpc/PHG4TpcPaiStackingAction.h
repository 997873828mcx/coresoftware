// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_PHG4TPCPAISTACKINGACTION_H
#define G4TPC_PHG4TPCPAISTACKINGACTION_H

#include <g4main/PHG4StackingAction.h>

#include <map>
#include <string>

class G4Track;
class PHCompositeNode;
class PHG4HitContainer;

class PHG4TpcPaiStackingAction : public PHG4StackingAction
{
 public:
  explicit PHG4TpcPaiStackingAction(const std::string &name = "PHG4TpcPaiStackingAction");
  ~PHG4TpcPaiStackingAction() override;

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *) override;
  void PrepareNewEvent() override;
  void SetInterfacePointers(PHCompositeNode *) override;

  void set_cluster_node_name(const std::string &name) { m_ClusterNodeName = name; }
  const std::string &cluster_node_name() const { return m_ClusterNodeName; }
  void set_region_name(const std::string &name) { m_RegionName = name; }
  void set_w_value_eV(const double value) { m_WValueEV = value; }
  void set_min_kinetic_energy_eV(const double value) { m_MinKineticEnergyEV = value; }
  void set_max_cluster_size(const unsigned int value) { m_MaxClusterSize = value; }

 private:
  bool DebugEnabled() const;
  void PrintDiagnostics() const;

  PHG4HitContainer *m_ClusterContainer{nullptr};
  std::string m_ClusterNodeName{"G4HIT_TPC_PAI_CLUSTER"};
  std::string m_RegionName{"REGION_TPCGAS"};
  double m_WValueEV{35.0};
  double m_MinKineticEnergyEV{0.0};
  unsigned int m_MaxClusterSize{100000};
  bool m_DebugEnabled{false};
  unsigned int m_NumClassifiedTracks{0};
  unsigned int m_NumTracksWithoutContainer{0};
  unsigned int m_NumPrimaryTracks{0};
  unsigned int m_NumSecondaryTracks{0};
  unsigned int m_NumSecondaryElectrons{0};
  unsigned int m_NumElectronsWithoutCreator{0};
  unsigned int m_NumRejectedCreator{0};
  unsigned int m_NumElectronsWithoutVolume{0};
  unsigned int m_NumRejectedRegion{0};
  unsigned int m_NumRejectedLowEnergy{0};
  unsigned int m_NumStoredClusters{0};
  unsigned int m_NumKilledElectrons{0};
  std::map<std::string, unsigned int> m_ParticleCounts;
  std::map<std::string, unsigned int> m_CreatorProcessCounts;
  std::map<std::string, unsigned int> m_RegionCounts;
  std::map<std::string, unsigned int> m_PhysicalVolumeCounts;
};

#endif
