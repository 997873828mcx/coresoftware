#include "PHG4TpcDigitizer.h"

#include <cdbobjects/CDBTTree.h>

#include <g4detectors/PHG4CylinderCellDefs.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4detectors/PHG4TpcGeomv1.h>

#include <g4main/PHG4HitDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/TrkrHitv2.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include <unistd.h>

namespace
{
  constexpr unsigned int kFirstTpcLayer = 7;
  constexpr unsigned int kLastTpcLayer = 54;
  constexpr unsigned int kMaskedSignalPad = 0;
  constexpr unsigned int kLiveSignalPad = 1;
  constexpr unsigned int kMaskedEmptyPad = 2;
  constexpr unsigned int kLiveEmptyPad = 3;
  constexpr unsigned int kTimeBin = 0;

  void require(const bool condition, const std::string& message)
  {
    if (!condition)
    {
      std::cerr << "test_PHG4TpcDigitizerMask: " << message << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::string makeTemporaryMapName()
  {
    return "/tmp/test_PHG4TpcDigitizerMask_" +
           std::to_string(static_cast<long long>(getpid())) + ".root";
  }

  void writeDeadChannelMap(const std::string& filename)
  {
    std::remove(filename.c_str());

    CDBTTree output(filename);
    output.SetSingleIntValue("TotalDeadChannels", 2);

    output.SetIntValue(0, "layer", kFirstTpcLayer);
    output.SetIntValue(0, "sector", 0);
    output.SetIntValue(0, "side", 0);
    output.SetIntValue(0, "pad", kMaskedSignalPad);

    output.SetIntValue(1, "layer", kFirstTpcLayer);
    output.SetIntValue(1, "sector", 0);
    output.SetIntValue(1, "side", 0);
    output.SetIntValue(1, "pad", kMaskedEmptyPad);

    output.Commit();
    output.CommitSingle();
    output.WriteCDBTTree();
  }

  PHG4TpcGeomContainer* addGeometry(PHCompositeNode* topNode)
  {
    auto* runNode = new PHCompositeNode("RUN");
    auto* geometryNode = new PHCompositeNode("RECO_TRACKING_GEOMETRY");
    topNode->addNode(runNode);
    runNode->addNode(geometryNode);

    auto* geometry = new PHG4TpcGeomContainer;
    geometryNode->addNode(
        new PHIODataNode<PHObject>(
            geometry, "TPCGEOMCONTAINER", "PHObject"));

    for (unsigned int layer = kFirstTpcLayer;
         layer <= kLastTpcLayer;
         ++layer)
    {
      auto* layerGeometry = new PHG4TpcGeomv1;
      layerGeometry->set_binning(PHG4CylinderCellDefs::sizebinning);
      layerGeometry->set_radius(30.0 + static_cast<double>(layer));
      layerGeometry->set_thickness(1.0);
      layerGeometry->set_zbins(1);
      layerGeometry->set_zmin(0.0);
      layerGeometry->set_zstep(53.326184);
      layerGeometry->set_adc_clock(53.326184);

      // The target layer uses the real inner-region pad count so that pads
      // 0--93 map to sector zero exactly as they do in production. Other
      // layers need only one pad for this focused digitizer test.
      const int phiBins = layer == kFirstTpcLayer ? 1128 : 1;
      layerGeometry->set_phibins(phiBins);
      layerGeometry->set_phistep(6.283185307179586 /
                                 static_cast<double>(phiBins));
      layerGeometry->set_drift_velocity_sim(1.0);

      const int addResult = geometry->AddLayerCellGeom(
          static_cast<int>(layer), layerGeometry);
      require(addResult == 0, "could not add synthetic TPC geometry layer");
    }

    return geometry;
  }

  TrkrHitSetContainerv1* addHitNodes(
      PHCompositeNode* topNode,
      TrkrHitTruthAssocv1*& truthAssociations)
  {
    auto* dstNode = new PHCompositeNode("DST");
    auto* trackerNode = new PHCompositeNode("TRKR");
    topNode->addNode(dstNode);
    dstNode->addNode(trackerNode);

    auto* hitsets = new TrkrHitSetContainerv1;
    trackerNode->addNode(
        new PHIODataNode<PHObject>(hitsets, "TRKR_HITSET", "PHObject"));

    truthAssociations = new TrkrHitTruthAssocv1;
    trackerNode->addNode(
        new PHIODataNode<PHObject>(
            truthAssociations, "TRKR_HITTRUTHASSOC", "PHObject"));

    return hitsets;
  }

  void addSignalHit(
      TrkrHitSet* hitset,
      const TrkrDefs::hitkey hitkey)
  {
    auto* hit = new TrkrHitv2;
    hit->addEnergy(100000.0);
    hitset->addHitSpecificKey(hitkey, hit);
  }

  std::size_t associationCount(
      const TrkrHitTruthAssoc* truthAssociations,
      const TrkrDefs::hitsetkey hitsetkey,
      const TrkrDefs::hitkey hitkey)
  {
    TrkrHitTruthAssoc::MMap associations;
    truthAssociations->getG4Hits(hitsetkey, hitkey, associations);
    return associations.size();
  }
}  // namespace

int main()
{
  const std::string mapFilename = makeTemporaryMapName();
  writeDeadChannelMap(mapFilename);

  auto topNode = std::make_unique<PHCompositeNode>("TOP");
  PHG4TpcGeomContainer* geometry = addGeometry(topNode.get());
  require(geometry, "could not create synthetic TPC geometry");

  TrkrHitTruthAssocv1* truthAssociations = nullptr;
  TrkrHitSetContainerv1* hitsets =
      addHitNodes(topNode.get(), truthAssociations);
  require(hitsets, "could not create TRKR_HITSET node");
  require(truthAssociations, "could not create TRKR_HITTRUTHASSOC node");

  const auto hitsetkey =
      TpcDefs::genHitSetKey(kFirstTpcLayer, 0, 0);
  const auto maskedSignalKey =
      TpcDefs::genHitKey(kMaskedSignalPad, kTimeBin);
  const auto liveSignalKey =
      TpcDefs::genHitKey(kLiveSignalPad, kTimeBin);
  const auto maskedEmptyKey =
      TpcDefs::genHitKey(kMaskedEmptyPad, kTimeBin);
  const auto liveEmptyKey =
      TpcDefs::genHitKey(kLiveEmptyPad, kTimeBin);

  TrkrHitSet* inputHitset =
      hitsets->findOrAddHitSet(hitsetkey)->second;
  addSignalHit(inputHitset, maskedSignalKey);
  addSignalHit(inputHitset, liveSignalKey);

  // Give the masked signal two contributors to verify that cleanup removes
  // every association rather than only the first one.
  truthAssociations->addAssoc(
      hitsetkey, maskedSignalKey,
      static_cast<PHG4HitDefs::keytype>(101));
  truthAssociations->addAssoc(
      hitsetkey, maskedSignalKey,
      static_cast<PHG4HitDefs::keytype>(102));
  truthAssociations->addAssoc(
      hitsetkey, liveSignalKey,
      static_cast<PHG4HitDefs::keytype>(201));
  require(
      associationCount(truthAssociations, hitsetkey, maskedSignalKey) == 2,
      "masked input hit does not have its two test associations");
  require(
      associationCount(truthAssociations, hitsetkey, liveSignalKey) == 1,
      "live input hit does not have its test association");

  PHG4TpcDigitizer digitizer("PHG4TpcDigitizerMaskTest");
  digitizer.SetMaskChannelsFromFile();
  digitizer.SetDeadChannelMapName(mapFilename);

  // Keep the empty-bin electronics path enabled, but make it deterministic:
  // every live empty bin contains the pedestal and therefore exceeds a
  // threshold placed one electron below it. A masked empty channel must still
  // produce no hit at all.
  digitizer.SetENC(0.0F);
  digitizer.SetADCThreshold(-1.0F);

  const int initResult = digitizer.InitRun(topNode.get());
  require(
      initResult == Fun4AllReturnCodes::EVENT_OK,
      "digitizer InitRun failed");
  const int eventResult = digitizer.process_event(topNode.get());
  require(
      eventResult == Fun4AllReturnCodes::EVENT_OK,
      "digitizer process_event failed");

  TrkrHitSet* outputHitset = hitsets->findHitSet(hitsetkey);
  require(outputHitset, "target output hitset is missing");

  require(
      !outputHitset->getHit(maskedSignalKey),
      "masked signal channel survived digitization");
  require(
      !outputHitset->getHit(maskedEmptyKey),
      "masked empty channel acquired a digitized hit");

  const TrkrHit* liveSignal = outputHitset->getHit(liveSignalKey);
  const TrkrHit* liveEmpty = outputHitset->getHit(liveEmptyKey);
  require(liveSignal, "live signal channel was removed");
  require(liveEmpty, "live empty channel was not digitized");
  require(liveSignal->getAdc() > 0, "live signal channel has zero ADC");
  require(liveEmpty->getAdc() > 0, "live empty channel has zero ADC");

  require(
      associationCount(truthAssociations, hitsetkey, maskedSignalKey) == 0,
      "masked signal truth associations survived cleanup");
  require(
      associationCount(truthAssociations, hitsetkey, liveSignalKey) == 1,
      "live signal truth association was removed");

  topNode.reset();
  if (std::remove(mapFilename.c_str()) != 0)
  {
    std::cerr << "Could not remove temporary map " << mapFilename
              << std::endl;
    return 1;
  }

  std::cout << "PHG4TpcDigitizer channel-mask regression test passed"
            << std::endl;
  return 0;
}
