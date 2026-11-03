#include "TpcClusterizer.h"
#include "TpcSimpleClusterizer.h"
#include "Tpc3DClusterizer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition)
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
    TrkrDefs::hitsetkey key {TpcDefs::genHitSetKey(7, 0, 0)};

    Event()
    {
      top->addNode(dst);
      dst->addNode(new PHIODataNode<PHObject>(hitsets, "TRKR_HITSET", "PHObject"));
      const unsigned int values[] = {59, 60, 61, 100};
      for (unsigned int index = 0; index < 4; ++index)
      {
        auto* hit = new TrkrHitv2;
        hit->setAdc(values[index]);
        hitsets->findOrAddHitSet(key)->second->addHitSpecificKey(index, hit);
      }
    }

    unsigned int adc(unsigned int index) const { return hitsets->findHitSet(key)->getHit(index)->getAdc(); }
    PdbParameterMap* metadata() const { return findNode::getClass<PdbParameterMap>(top.get(), "TPC_ADC_BASELINE"); }
    void addMetadata()
    {
      auto* baseline = new PdbParameterMap;
      baseline->set_int_param("schema_version", 1);
      baseline->set_int_param("state", 0);
      baseline->set_double_param("pedestal_adc", 60.0);
      baseline->set_double_param("removed_pedestal_adc", 0.0);
      dst->addNode(new PHIODataNode<PHObject>(baseline, "TPC_ADC_BASELINE", "PHObject"));
    }
  };

  template <class Consumer>
  void exerciseConsumer(Consumer& consumer)
  {
    // Intentionally omit tracking geometry and output containers. Exercise the
    // real consumer's event entry point, then stop before any trajectory or
    // centroid calculation; the expected ABORTRUN is this fixture's missing
    // reconstruction inputs, not a successful reconstruction claim.
    Event raw;
    raw.addMetadata();
    const int stopped = Fun4AllReturnCodes::ABORTRUN;
    require(consumer.process_event(raw.top.get()) == stopped, consumer.Name() + ": fixture did not stop before geometry");
    require(raw.adc(0) == 0 && raw.adc(1) == 0 && raw.adc(2) == 1 && raw.adc(3) == 40,
            consumer.Name() + ": raw 60-ADC event was not automatically normalized");
    require(raw.hitsets->findHitSet(raw.key)->size() == 4, consumer.Name() + ": zero ADC hit keys were removed");
    require(consumer.process_event(raw.top.get()) == stopped && raw.adc(3) == 40,
            consumer.Name() + ": repeated consumer subtracted twice");
    Event unmarked;
    require(consumer.process_event(unmarked.top.get()) == stopped && unmarked.adc(3) == 100 && !unmarked.metadata(),
            consumer.Name() + ": unmarked legacy/data input was guessed or changed");
    Event malformed;
    malformed.addMetadata();
    malformed.metadata()->set_int_param("schema_version", 99);
    require(consumer.process_event(malformed.top.get()) == stopped && malformed.adc(3) == 100,
            consumer.Name() + ": malformed metadata was consumed or partially corrected");
  }
}

int main()
{
  TpcClusterizer standard;
  exerciseConsumer(standard);
  Event marked;
  marked.addMetadata();
  standard.set_pedestal(60.0);
  standard.process_event(marked.top.get());
  require(marked.adc(3) == 40, "standard consumer ignored baseline metadata with configured offset");
  Event data;
  standard.set_rawdata_reco();
  standard.process_event(data.top.get());
  require(data.adc(3) == 100,
          "standard data-unpacker convention changed");
  Event rawObject;
  rawObject.addMetadata();
  standard.set_read_raw(true);
  standard.process_event(rawObject.top.get());
  require(rawObject.adc(3) == 100, "RawHitSet consumer changed an unrelated TRKR_HITSET");

  TpcSimpleClusterizer simple;
  exerciseConsumer(simple);
  Event simpleMarked;
  simpleMarked.addMetadata();
  simple.process_event(simpleMarked.top.get());
  require(simpleMarked.adc(3) == 40, "simple consumer ignored baseline metadata");

  Tpc3DClusterizer threeD;
  exerciseConsumer(threeD);
  std::cout << "TPC clusterizer ADC consumer integration PASS" << std::endl;
  return EXIT_SUCCESS;
}
