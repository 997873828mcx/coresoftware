#include "PHG4TpcPaiStackingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <Geant4/G4Electron.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4VProcess.hh>
#include <Geant4/G4VPhysicalVolume.hh>

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iostream>

PHG4TpcPaiStackingAction::PHG4TpcPaiStackingAction(const std::string &name)
  : PHG4StackingAction(name)
{
  m_DebugEnabled = DebugEnabled();
}

PHG4TpcPaiStackingAction::~PHG4TpcPaiStackingAction()
{
  if (m_DebugEnabled)
  {
    PrintDiagnostics();
  }
}

bool PHG4TpcPaiStackingAction::DebugEnabled() const
{
  const char *debug_env = std::getenv("TPC_PAI_DEBUG");
  return debug_env && debug_env[0] != '\0' && debug_env[0] != '0';
}

void PHG4TpcPaiStackingAction::PrintDiagnostics() const
{
  std::cout << GetName() << " diagnostics:"
            << " classified=" << m_NumClassifiedTracks
            << " no_container=" << m_NumTracksWithoutContainer
            << " primary=" << m_NumPrimaryTracks
            << " secondary=" << m_NumSecondaryTracks
            << " secondary_electrons=" << m_NumSecondaryElectrons
            << " electron_no_creator=" << m_NumElectronsWithoutCreator
            << " rejected_creator=" << m_NumRejectedCreator
            << " electron_no_volume=" << m_NumElectronsWithoutVolume
            << " rejected_region=" << m_NumRejectedRegion
            << " rejected_low_energy=" << m_NumRejectedLowEnergy
            << " stored=" << m_NumStoredClusters
            << " killed=" << m_NumKilledElectrons
            << std::endl;

  auto print_counts = [&](const std::string &label, const std::map<std::string, unsigned int> &counts)
  {
    std::cout << GetName() << " " << label << ":";
    if (counts.empty())
    {
      std::cout << " <none>";
    }
    for (const auto &entry : counts)
    {
      std::cout << " " << entry.first << "=" << entry.second;
    }
    std::cout << std::endl;
  };

  print_counts("particles", m_ParticleCounts);
  print_counts("creator_processes", m_CreatorProcessCounts);
  print_counts("regions", m_RegionCounts);
  print_counts("physical_volumes", m_PhysicalVolumeCounts);
}

void PHG4TpcPaiStackingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  m_ClusterContainer = findNode::getClass<PHG4HitContainer>(topNode, m_ClusterNodeName);
  if (m_ClusterContainer)
  {
    return;
  }

  PHNodeIterator iter(topNode);
  auto *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << GetName() << " - DST node missing, cannot create " << m_ClusterNodeName << std::endl;
    return;
  }

  PHNodeIterator dstiter(dstNode);
  auto *detNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TPC"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("TPC");
    dstNode->addNode(detNode);
  }

  m_ClusterContainer = new PHG4HitContainer(m_ClusterNodeName);
  auto *newNode = new PHIODataNode<PHObject>(m_ClusterContainer, m_ClusterNodeName, "PHObject");
  detNode->addNode(newNode);
}

void PHG4TpcPaiStackingAction::PrepareNewEvent()
{
  if (m_ClusterContainer)
  {
    m_ClusterContainer->Reset();
  }
}

G4ClassificationOfNewTrack PHG4TpcPaiStackingAction::ClassifyNewTrack(const G4Track *track)
{
  if (!track)
  {
    return fUrgent;
  }

  if (m_DebugEnabled)
  {
    ++m_NumClassifiedTracks;
    if (track->GetDefinition())
    {
      ++m_ParticleCounts[track->GetDefinition()->GetParticleName()];
    }
  }

  if (!m_ClusterContainer)
  {
    if (m_DebugEnabled)
    {
      ++m_NumTracksWithoutContainer;
    }
    return fUrgent;
  }
  if (track->GetParentID() <= 0)
  {
    if (m_DebugEnabled)
    {
      ++m_NumPrimaryTracks;
    }
    return fUrgent;
  }
  if (m_DebugEnabled)
  {
    ++m_NumSecondaryTracks;
  }
  if (track->GetDefinition() != G4Electron::ElectronDefinition())
  {
    return fUrgent;
  }
  if (m_DebugEnabled)
  {
    ++m_NumSecondaryElectrons;
  }

  const G4VProcess *creator_process = track->GetCreatorProcess();
  if (!creator_process)
  {
    if (m_DebugEnabled)
    {
      ++m_NumElectronsWithoutCreator;
    }
    return fUrgent;
  }
  const G4String &creator_name = creator_process->GetProcessName();
  if (m_DebugEnabled)
  {
    ++m_CreatorProcessCounts[creator_name];
  }
  if (creator_name != "eIoni" && creator_name != "hIoni" && creator_name != "muIoni")
  {
    if (m_DebugEnabled)
    {
      ++m_NumRejectedCreator;
    }
    return fUrgent;
  }

  const G4VPhysicalVolume *physical_volume = track->GetVolume();
  if (m_DebugEnabled && physical_volume)
  {
    ++m_PhysicalVolumeCounts[physical_volume->GetName()];
  }

  const G4LogicalVolume *logical_volume = physical_volume ? physical_volume->GetLogicalVolume() : nullptr;
  if (!logical_volume)
  {
    logical_volume = track->GetLogicalVolumeAtVertex();
  }
  if (!logical_volume || !logical_volume->GetRegion())
  {
    if (m_DebugEnabled)
    {
      ++m_NumElectronsWithoutVolume;
    }
    return fUrgent;
  }
  if (m_DebugEnabled)
  {
    ++m_RegionCounts[logical_volume->GetRegion()->GetName()];
  }
  const G4String region_name(m_RegionName.c_str());
  const bool is_tpc_gas_region = logical_volume->GetRegion()->GetName() == region_name;
  const bool is_tpc_gas_volume = physical_volume &&
                                 (physical_volume->GetName() == "tpc_gas_north" ||
                                  physical_volume->GetName() == "tpc_gas_south");
  if (!is_tpc_gas_region && !is_tpc_gas_volume)
  {
    if (m_DebugEnabled)
    {
      ++m_NumRejectedRegion;
    }
    return fUrgent;
  }

  const double kinetic_energy_eV = track->GetKineticEnergy() / eV;
  if (kinetic_energy_eV < m_MinKineticEnergyEV)
  {
    if (m_DebugEnabled)
    {
      ++m_NumRejectedLowEnergy;
    }
    return fKill;
  }

  unsigned int cluster_size = 1;
  if (m_WValueEV > 0.0 && m_MaxClusterSize > 0 && std::isfinite(kinetic_energy_eV))
  {
    const double secondary_electrons = std::floor(kinetic_energy_eV / m_WValueEV);
    if (secondary_electrons > 0.0)
    {
      const auto capped = std::min<double>(secondary_electrons, static_cast<double>(m_MaxClusterSize - 1));
      cluster_size += static_cast<unsigned int>(capped);
    }
  }

  const G4ThreeVector &position = track->GetPosition();
  auto *hit = new PHG4Hitv1();
  hit->set_x(0, position.x() / cm);
  hit->set_y(0, position.y() / cm);
  hit->set_z(0, position.z() / cm);
  hit->set_t(0, track->GetGlobalTime() / nanosecond);
  hit->set_x(1, position.x() / cm);
  hit->set_y(1, position.y() / cm);
  hit->set_z(1, position.z() / cm);
  hit->set_t(1, track->GetGlobalTime() / nanosecond);
  hit->set_trkid(track->GetParentID());
  hit->set_edep(track->GetKineticEnergy() / GeV);
  hit->set_eion(track->GetKineticEnergy() / GeV);
  hit->set_index_i(static_cast<int>(cluster_size));
  hit->set_hit_type(1);
  m_ClusterContainer->AddHit(0, hit);

  if (m_DebugEnabled)
  {
    ++m_NumStoredClusters;
    ++m_NumKilledElectrons;
  }

  if (Verbosity() > 100)
  {
    std::cout << GetName() << " stored PAI cluster seed at x/y/z "
              << hit->get_x(0) << "/" << hit->get_y(0) << "/" << hit->get_z(0)
              << " cm, kinetic energy " << kinetic_energy_eV
              << " eV, cluster size " << cluster_size
              << ", parent track " << track->GetParentID() << std::endl;
  }

  return fKill;
}
