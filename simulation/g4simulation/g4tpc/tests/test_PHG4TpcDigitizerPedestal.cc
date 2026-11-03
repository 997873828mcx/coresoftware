#include "PHG4TpcDigitizer.h"
#include "TpcSampaResponse.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderCellDefs.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4detectors/PHG4TpcGeomv1.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <phparameter/PHParameters.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/TrkrHitv2.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

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

  void checkDigitizer(double pedestal, bool noNoise, bool explicitSetting = true)
  {
    auto top = std::make_unique<PHCompositeNode>("TOP");
    auto* run = new PHCompositeNode("RUN");
    auto* dst = new PHCompositeNode("DST");
    top->addNode(run);
    top->addNode(dst);
    auto* geometry = new PHG4TpcGeomContainer;
    run->addNode(new PHIODataNode<PHObject>(geometry, "TPCGEOMCONTAINER", "PHObject"));
    for (int layer = 7; layer < 55; ++layer)
    {
      auto* geom = new PHG4TpcGeomv1;
      geom->set_binning(PHG4CylinderCellDefs::sizebinning);
      geom->set_radius(30.0 + layer);
      geom->set_thickness(1.0);
      geom->set_zbins(8);
      geom->set_zmin(0.0);
      geom->set_zstep(53.326184);
      geom->set_adc_clock(53.326184);
      geom->set_phibins(12);
      geom->set_phistep(6.283185307179586 / 12.0);
      geom->set_drift_velocity_sim(1.0);
      require(geometry->AddLayerCellGeom(layer, geom) == 0, "geometry creation failed");
    }
    auto* hitsets = new TrkrHitSetContainerv1;
    dst->addNode(new PHIODataNode<PHObject>(hitsets, "TRKR_HITSET", "PHObject"));
    auto* truth = new TrkrHitTruthAssocv1;
    dst->addNode(new PHIODataNode<PHObject>(truth, "TRKR_HITTRUTHASSOC", "PHObject"));
    const auto hsk = TpcDefs::genHitSetKey(7, 0, 0);
    auto* signal = new TrkrHitv2;
    signal->addEnergy(100000.0);
    hitsets->findOrAddHitSet(hsk)->second->addHitSpecificKey(TpcDefs::genHitKey(0, 0), signal);
    truth->addAssoc(hsk, TpcDefs::genHitKey(0, 0), 123);
    PHG4TpcDigitizer digitizer;
    require(digitizer.GetPedestalADC() == 60.0, "source default is not 60 ADC");
    if (explicitSetting) digitizer.SetPedestalADC(pedestal);
    digitizer.SetENC(0);
    digitizer.set_skip_noise_flag(true);
    digitizer.set_no_noise_flag(noNoise);
    require(digitizer.InitRun(top.get()) == Fun4AllReturnCodes::EVENT_OK, "digitizer initialization failed");
    require(std::abs(digitizer.GetSignalConversionGain() - 20.0 * 1.60e-4) < 1.e-9,
            "analog-peak conversion contains an extra response scale");
    require(digitizer.GetSignalConversionGain() == digitizer.GetNoiseConversionGain(),
            "signal/noise physical electron conversion differs");
    require(digitizer.process_event(top.get()) == Fun4AllReturnCodes::EVENT_OK, "digitizer event failed");
    auto* output = hitsets->findHitSet(hsk);
    require(output && output->size() == (noNoise ? 1U : 5U), "unexpected five-bin readout after zero-ADC cleanup");
    for (unsigned int bin = 1; bin < 5; ++bin)
    {
      auto* tail = output->getHit(TpcDefs::genHitKey(0, bin));
      require(noNoise ? !tail : tail && tail->getAdc() == static_cast<unsigned int>(pedestal), "following-bin pedestal is wrong");
    }
    require(output->getHit(TpcDefs::genHitKey(0, 0))->getAdc() == (noNoise ? 148U : 148U + static_cast<unsigned int>(pedestal)), "signal amplitude/gain changed");
    auto* metadata = findNode::getClass<PdbParameterMap>(top.get(), "TPC_ADC_BASELINE");
    require(metadata, "digitizer did not publish baseline metadata");
    const auto doubles = metadata->get_dparam_iters();
    bool found = false;
    for (auto value = doubles.first; value != doubles.second; ++value)
    {
      if (value->first == "pedestal_adc")
      {
        found = value->second == (noNoise ? 0.0 : pedestal);
      }
    }
    require(found, "digitizer baseline metadata differs from actual setting");
  }

  void checkResponseMetadata()
  {
    for (int mode = 0; mode < 4; ++mode)
    {
      auto top = std::make_unique<PHCompositeNode>("TOP");
      top->addNode(new PHCompositeNode("DST"));
      auto* par = new PHCompositeNode("PAR");
      top->addNode(par);
      auto* parameters = new PHParameters("TPC");
      if (mode != 0)
      {
        parameters->set_string_param("tpc_sampa_response_model",
                                    mode == 1 ? "legacy_unit_area_shares" : TpcSampaResponse::Mode);
        parameters->set_int_param("tpc_sampa_response_version", mode == 2 ? 999 : TpcSampaResponse::Version);
      }
      par->addNode(new PHDataNode<PHParameters>(parameters, "G4TPCPADPLANEPAR"));
      PHG4TpcDigitizer digitizer;
      const int result = digitizer.InitRun(top.get());
      require(result == (mode == 3 ? Fun4AllReturnCodes::EVENT_OK : Fun4AllReturnCodes::ABORTRUN),
              "legacy or mismatched pad-readout representation was accepted");
    }
  }
}

int main()
{
  checkDigitizer(60.0, false, false);
  checkDigitizer(100.0, false);
  checkDigitizer(60.0, true);
  checkResponseMetadata();
  bool rejected = false;
  try { PHG4TpcDigitizer().SetPedestalADC(-1); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "negative pedestal setting accepted");
  std::cout << "TPC digitizer explicit ADC pedestal and following-bin tests passed" << std::endl;
  return 0;
}
