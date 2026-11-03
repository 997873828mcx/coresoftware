#include "TpcHitsetPedestalCorrection.h"
#include "Tpc_ModuleTrack.h"
#include "Tpc_ModuleTrackContainer.h"
#include "Tpc_ModuleTrackReco.h"
#include "Tpc_PolyClusterizer.h"
#include "Tpc_AssembledTrackReco.h"
#include <tpc/TpcAdcBaselineCorrection.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/TrkrHitv2.h>

#include <TFile.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace
{
  void require(bool ok, const std::string& message)
  {
    if (!ok)
    {
      std::cerr << message << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  struct Event
  {
    std::unique_ptr<PHCompositeNode> top {std::make_unique<PHCompositeNode>("TOP")};
    PHCompositeNode* dst {new PHCompositeNode("DST")};
    TrkrHitSetContainerv1* hitsets {new TrkrHitSetContainerv1};
    TrkrHitTruthAssocv1* truth {new TrkrHitTruthAssocv1};
    TrkrDefs::hitsetkey tpcKey {TpcDefs::genHitSetKey(7, 0, 0)};
    TrkrDefs::hitsetkey otherKey {TrkrDefs::genHitSetKey(TrkrDefs::TrkrId::mvtxId, 0)};

    Event()
    {
      top->addNode(dst);
      dst->addNode(new PHIODataNode<PHObject>(hitsets, "TRKR_HITSET", "PHObject"));
      dst->addNode(new PHIODataNode<PHObject>(truth, "TRKR_HITTRUTHASSOC", "PHObject"));
      for (unsigned int i = 0; i < 4; ++i)
      {
        auto* hit = new TrkrHitv2;
        hit->setAdc(i == 0 ? 59 : i == 1 ? 60 : i == 2 ? 61 : 100);
        hitsets->findOrAddHitSet(tpcKey)->second->addHitSpecificKey(i, hit);
        truth->addAssoc(tpcKey, i, 1000 + i);
      }
      auto* otherHit = new TrkrHitv2;
      otherHit->setAdc(123);
      hitsets->findOrAddHitSet(otherKey)->second->addHitSpecificKey(0, otherHit);
    }

    unsigned int adc(unsigned int i) const { return hitsets->findHitSet(tpcKey)->getHit(i)->getAdc(); }
    PdbParameterMap* metadata() const { return findNode::getClass<PdbParameterMap>(top.get(), "TPC_ADC_BASELINE"); }
    PdbParameterMap* addMetadata(double pedestal = 60.0, int state = 0)
    {
      auto* map = new PdbParameterMap;
      map->set_int_param("schema_version", 1);
      map->set_int_param("state", state);
      map->set_double_param("pedestal_adc", pedestal);
      map->set_double_param("removed_pedestal_adc", state == 0 ? 0.0 : pedestal);
      dst->addNode(new PHIODataNode<PHObject>(map, "TPC_ADC_BASELINE", "PHObject"));
      return map;
    }
  };
}

int main()
{
  constexpr int ok = Fun4AllReturnCodes::EVENT_OK;
  constexpr int bad = Fun4AllReturnCodes::ABORTRUN;
  Event event;
  event.addMetadata();
  TpcHitsetPedestalCorrection correction;
  correction.setInputPedestalADC(60);
  require(correction.InitRun(event.top.get()) == ok, "InitRun failed");
  require(correction.process_event(event.top.get()) == ok, "raw conversion failed");
  require(event.adc(0) == 0 && event.adc(1) == 0 && event.adc(2) == 1 && event.adc(3) == 40, "60 ADC arithmetic/underflow failed");
  require(event.hitsets->findHitSet(event.tpcKey)->size() == 4, "zero ADC keys were removed");
  require(event.hitsets->findHitSet(event.otherKey)->getHit(0)->getAdc() == 123, "non-TPC hit changed");
  TrkrHitTruthAssoc::MMap associations;
  event.truth->getG4Hits(event.tpcKey, 0, associations);
  require(associations.size() == 1, "zero-weight truth association changed");
  require(correction.getEventsCorrected() == 1 && correction.getHitsCorrected() == 4 && correction.getHitsClamped() == 1 && correction.getHitsZero() == 2, "counters failed");
  require(correction.getInputPedestalADC() == 60 && correction.getEffectiveIntegerShiftADC() == 60, "resolved pedestal getters failed");
  require(correction.process_event(event.top.get()) == ok && event.adc(3) == 40 && correction.getEventsAlreadySubtracted() == 1, "second invocation subtracted twice");

  TpcHitsetPedestalCorrection wrongPedestal;
  wrongPedestal.setInputPedestalADC(74.4727272727);
  require(wrongPedestal.process_event(event.top.get()) == bad && event.adc(3) == 40, "mismatched already-subtracted metadata was accepted");

  Event undeclared;
  TpcHitsetPedestalCorrection strict;
  require(strict.process_event(undeclared.top.get()) == bad && undeclared.adc(3) == 100, "undeclared legacy input was changed");

  Event legacy;
  TpcHitsetPedestalCorrection old;
  old.setInputPedestalADC(74.4727272727);
  require(old.process_event(legacy.top.get()) == ok && legacy.adc(3) == 26 && old.getEffectiveIntegerShiftADC() == 74, "legacy fractional pedestal quantization failed");
  legacy.metadata()->Reset();
  legacy.hitsets->findHitSet(legacy.tpcKey)->getHit(3)->setAdc(100);
  require(old.process_event(legacy.top.get()) == ok && legacy.adc(3) == 26 && old.getEventsCorrected() == 2, "event metadata reset did not restore legacy declaration");

  Event data;
  TpcHitsetPedestalCorrection already;
  already.setInputAlreadySubtracted(true);
  require(already.process_event(data.top.get()) == ok && data.adc(3) == 100 && already.getEventsCorrected() == 0, "explicit real-data input changed");
  Event rawClaimedAsData;
  rawClaimedAsData.addMetadata();
  require(already.process_event(rawClaimedAsData.top.get()) == bad, "raw metadata was accepted as subtracted data");
  Event malformed;
  malformed.addMetadata()->set_int_param("schema_version", 99);
  require(strict.process_event(malformed.top.get()) == bad, "invalid metadata schema accepted");
  Event partial;
  partial.addMetadata()->Reset();
  partial.metadata()->set_int_param("state", 1);
  require(correction.process_event(partial.top.get()) == bad && partial.adc(3) == 100, "partial metadata was treated as a legacy event");
  bool rejected = false;
  try { strict.setInputPedestalADC(std::numeric_limits<double>::quiet_NaN()); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "non-finite declared pedestal accepted");

  // A custom macro need not register the correction subsystem. The first
  // consumer normalizes metadata-marked ADC, then overrides its historical
  // internal pedestal with zero. Six 40-count signals must remain findable.
  Event custom;
  custom.addMetadata();
  for (unsigned int layer = 7; layer < 13; ++layer)
  {
    auto* hit = new TrkrHitv2;
    hit->setAdc(100);
    custom.hitsets->findOrAddHitSet(TpcDefs::genHitSetKey(layer, 0, 0))->second->addHitSpecificKey(TpcDefs::genHitKey(10, 10), hit);
  }
  Tpc_ModuleTrackReco finder;
  finder.setMaxThreads(1);
  finder.setPedestal(60); // must not subtract this again after normalization
  require(finder.InitRun(custom.top.get()) == ok, "custom macro finder InitRun failed");
  require(finder.process_event(custom.top.get()) == ok && custom.adc(3) == 40, "finder did not auto-normalize raw ADC");
  auto* moduleTracks = findNode::getClass<Tpc_ModuleTrackContainer>(custom.top.get(), "TPC_MODULETRACKS");
  require(moduleTracks && moduleTracks->size() > 0 && moduleTracks->get_track(0)->get_nrawhits() >= 6, "finder subtracted its internal pedestal a second time");
  require(finder.process_event(custom.top.get()) == ok && custom.adc(3) == 40, "finder normalized ADC twice");
  TpcHitsetPedestalCorrection afterFinder;
  afterFinder.setInputPedestalADC(60);
  require(afterFinder.process_event(custom.top.get()) == ok && afterFinder.getEventsCorrected() == 0 && afterFinder.getEventsAlreadySubtracted() == 1, "explicit correction after automatic consumer was not idempotent");

  Event polyOnly;
  polyOnly.addMetadata();
  Tpc_PolyClusterizer poly;
  // Empty reconstruction inputs exercise the ADC consumer hook without
  // invoking Garfield or introducing a different trajectory calculation.
  require(poly.process_event(polyOnly.top.get()) == ok && polyOnly.adc(3) == 40, "polyclusterizer did not auto-normalize raw ADC");
  require(poly.process_event(polyOnly.top.get()) == ok && polyOnly.adc(3) == 40, "polyclusterizer normalized ADC twice");
  Event assemblyOnly;
  assemblyOnly.addMetadata();
  Tpc_AssembledTrackReco assembler;
  require(assembler.process_event(assemblyOnly.top.get()) == ok && assemblyOnly.adc(3) == 40, "assembler did not auto-normalize raw ADC");
  require(assembler.process_event(assemblyOnly.top.get()) == ok && assemblyOnly.adc(3) == 40, "assembler normalized ADC twice");
  Event noMetadata;
  const auto noGuess = TpcAdcBaseline::normalizeIfMetadata(noMetadata.top.get());
  require(noGuess.ok && !noGuess.has_metadata && noMetadata.adc(3) == 100 && !noMetadata.metadata(), "automatic helper guessed an unmarked legacy/data pedestal");
  require(finder.process_event(partial.top.get()) == bad, "consumer ignored partial metadata");

  const std::string filename = "/tmp/test_TpcHitsetPedestalCorrection_" + std::to_string(getpid()) + ".root";
  {
    TFile output(filename.c_str(), "RECREATE");
    output.WriteObject(event.metadata(), "TPC_ADC_BASELINE");
  }
  {
    TFile input(filename.c_str(), "READ");
    require(!input.IsZombie() && !input.TestBit(TFile::kRecovered), "metadata ROOT file is unreadable/recovered");
    auto* saved = input.Get<PdbParameterMap>("TPC_ADC_BASELINE");
    require(saved, "metadata did not persist");
    Event reloaded;
    reloaded.dst->addNode(new PHIODataNode<PHObject>(new PdbParameterMap(*saved), "TPC_ADC_BASELINE", "PHObject"));
    require(correction.process_event(reloaded.top.get()) == ok && reloaded.adc(3) == 100, "persisted marker did not prevent subtraction");
  }
  std::remove(filename.c_str());
  std::cout << "TPC ADC pedestal contract and subtraction tests passed" << std::endl;
  return 0;
}
