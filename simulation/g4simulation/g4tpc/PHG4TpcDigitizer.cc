#include "PHG4TpcDigitizer.h"
#include "TpcSampaResponse.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrHitv2.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco
#
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>  // for exit
#include <iostream>
#include <limits>
#include <memory>  // for allocator_tra...
#include <stdexcept>

PHG4TpcDigitizer::PHG4TpcDigitizer(const std::string &name)
  : SubsysReco(name)
  , RandomGenerator(gsl_rng_alloc(gsl_rng_mt19937))  // will be assigned in PHG4TpcDigitizer::InitRun
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  std::cout << Name() << " random seed: " << seed << std::endl;

  gsl_rng_set(RandomGenerator, seed);

  if (Verbosity() > 0)
  {
    std::cout << "Creating PHG4TpcDigitizer with name = " << name << std::endl;
  }
}

PHG4TpcDigitizer::~PHG4TpcDigitizer()
{
  gsl_rng_free(RandomGenerator);
}

void PHG4TpcDigitizer::SetPedestalADC(const double pedestal)
{
  if (!std::isfinite(pedestal) || pedestal < 0.0 || pedestal > 1023.0)
  {
    throw std::invalid_argument("PHG4TpcDigitizer pedestal must be in [0,1023] ADC");
  }
  m_pedestalADC = pedestal;
}

double PHG4TpcDigitizer::GetPedestalADC() const
{
  return std::isfinite(m_pedestalADC)
             ? m_pedestalADC
             : static_cast<double>(Pedestal) * ChargeToPeakVolts * 1.60e-04 * 1024.0 / 2200.0;
}

int PHG4TpcDigitizer::InitRun(PHCompositeNode *topNode)
{
  // A coupled pad-readout/digitizer chain must agree on its pre-ADC units.
  // Standalone users may provide equivalent-electron amplitudes directly
  // without a pad-plane parameter node. Legacy unit-area charge shares need
  // their original digitizer and must not silently enter this conversion.
  const auto* readout = findNode::getClass<PHParameters>(topNode, "G4TPCPADPLANEPAR");
  if (readout &&
      (!readout->exist_string_param("tpc_sampa_response_model") ||
       !readout->exist_int_param("tpc_sampa_response_version") ||
       readout->get_string_param("tpc_sampa_response_model") != TpcSampaResponse::Mode ||
       readout->get_int_param("tpc_sampa_response_version") != TpcSampaResponse::Version))
  {
    std::cerr << PHWHERE << "incompatible pad-readout SAMPA response representation" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // The pad readout already applies the peak-normalized SAMPA response averaged
  // over each ADC clock bin. Its stored value is an equivalent-electron signal
  // amplitude, not a unit-area share of the total GEM charge. Only the physical
  // charge-to-analog-peak gain is needed here; an extra empirical multiplier
  // would count the response normalization twice.
  ADCSignalConversionGain = ChargeToPeakVolts * 1.60e-04;  // mV/fC * fC/electron
  // The noise is by definition the RMS noise width voltage divided by ChargeToPeakVolts
  ADCNoiseConversionGain = ChargeToPeakVolts * 1.60e-04;  // 20 (or 30) mV/fC * fC/electron

  const float threshold_electrons = ADCThreshold + (no_noise ? 0.0F : Pedestal);
  ADCThreshold_mV = threshold_electrons * ADCNoiseConversionGain;
  if (std::isfinite(m_pedestalADC))
  {
    ADCThreshold_mV = ADCThreshold * ADCNoiseConversionGain +
                      (no_noise ? 0.0 : m_pedestalADC * 2200.0 / 1024.0);
  }

  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);
  if (!findNode::getClass<PdbParameterMap>(topNode, "TPC_ADC_BASELINE"))
  {
    dstNode->addNode(new PHIODataNode<PHObject>(
        new PdbParameterMap, "TPC_ADC_BASELINE", "PHObject"));
  }

  if (m_maskDeadChannels &&
      !m_deadChannelMap.load(
          m_deadChannelMapName,
          "TotalDeadChannels",
          m_maskFromFile,
          Verbosity()))
  {
    std::cerr << PHWHERE << "failed to load requested TPC dead-channel map"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (m_maskHotChannels &&
      !m_hotChannelMap.load(
          m_hotChannelMapName,
          "TotalHotChannels",
          m_maskFromFile,
          Verbosity()))
  {
    std::cerr << PHWHERE << "failed to load requested TPC hot-channel map"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CalculateCylinderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    std::cout << "====================== PHG4TpcDigitizer::InitRun() =====================" << std::endl;
    for (auto &tpiter : _max_adc)
    {
      std::cout << " Max ADC in Layer #" << tpiter.first << " = " << tpiter.second << std::endl;
    }
    for (auto &tpiter : _energy_scale)
    {
      std::cout << " Energy per ADC in Layer #" << tpiter.first << " = " << 1.0e6 * tpiter.second << " keV" << std::endl;
    }
    std::cout << "===========================================================================" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcDigitizer::process_event(PHCompositeNode *topNode)
{
  // Event-level metadata travels with TRKR_HITSET. It is reset and filled for
  // every event, including empty events; it must be saved by the DST writer.
  auto* baseline = findNode::getClass<PdbParameterMap>(topNode, "TPC_ADC_BASELINE");
  if (!baseline)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  baseline->Reset();
  baseline->set_int_param("schema_version", 1);
  baseline->set_int_param("state", 0);
  baseline->set_double_param("pedestal_adc", no_noise ? 0.0 : GetPedestalADC());
  baseline->set_double_param("removed_pedestal_adc", 0.0);
  baseline->set_string_param("source", "PHG4TpcDigitizer");
  DigitizeCylinderCells(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcDigitizer::CalculateCylinderCellADCScale(PHCompositeNode *topNode)
{
  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

  PHG4TpcGeomContainer *geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");

  if (!geom_container)
  {
    return;
  }

  PHG4TpcGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4TpcGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_thickness();
    float pitch = (layeriter->second)->get_phistep() * (layeriter->second)->get_radius();
    float length = (layeriter->second)->get_zstep() * (layeriter->second)->get_drift_velocity_sim();

    float minpath = pitch;
    minpath = std::min(length, minpath);
    minpath = std::min(thickness, minpath);
    float mip_e = 0.003876 * minpath;

    if (!_max_adc.contains(layer))
    {
      // cppcheck-suppress stlFindInsert
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }

  return;
}

void PHG4TpcDigitizer::DigitizeCylinderCells(PHCompositeNode *topNode)
{
  unsigned int print_layer = 18;  // to print diagnostic output for layer 47

  // Digitizes equivalent-electron amplitudes from PHG4TpcPadPlaneReadout.
  // Each time-bin value includes the peak-normalized SAMPA bin-average response;
  // summing those amplitudes is not a count of physical GEM electrons.

  // NOTES:
  // Modified by ADF June 2018 to do the following:
  //   Add noise to cells before digitizing
  //   Digitize the first adc time bin to exceed the threshold, and the 4 bins after that
  //   If the adc value is still above the threshold after 5 bins, repeat for the next 5 bins

  // For pad charge Q and unit-peak analog response g, the pad readout stores
  // E_i = Q * integral_bin(g(t-arrival) dt) / ADC_clock as equivalent electrons.
  // V_i [mV] = E_i * 1.60e-4 [fC/electron] * ChargeToPeakVolts [mV/fC].
  // The stored pulse sum is Q * response_area / ADC_clock before truncation;
  // it is not normalized to Q, and a bin's voltage need not equal the analog
  // peak. Clock-phase dependence and loss outside the readout window remain.
  // ADC = V_i * 1024 / 2200, truncated and saturated to the 10-bit ADC range.
  //
  // The noise/pedestal contract is independent of the signal pulse shape:
  // sigma_noise [mV] = ENC * 1.60e-4 * ChargeToPeakVolts. The configured
  // pedestal is added before zero suppression (unless no_noise is selected).
  // Neither ENC, the physical gain, nor the threshold receives a pulse-area
  // correction. At 20 mV/fC, 670-electron ENC corresponds to about one ADC count.

  //----------
  // Get Nodes
  //----------

  PHG4TpcGeomContainer *geom_container =
      findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TPCGEOMCONTAINER" << std::endl;
  }

  // Get the TrkrHitSetContainer node
  TrkrHitSetContainer *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkrhitsetcontainer)
  {
    std::cout << "Could not locate TRKR_HITSET node, quit! " << std::endl;
    exit(1);
  }

  TrkrHitTruthAssoc *hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    std::cout << PHWHERE << " Could not locate TRKR_HITTRUTHASSOC node, quit! " << std::endl;
    exit(1);
  }

  //-------------
  // Digitization
  //-------------

  // Loop over all TPC layers
  for (unsigned int layer = TpcMinLayer; layer < TpcMinLayer + TpcNLayers; ++layer)
  {
    // we need the geometry object for this layer
    PHG4TpcGeom *layergeom = geom_container->GetLayerCellGeom(layer);
    if (!layergeom)
    {
      exit(1);
    }

    int nphibins = layergeom->get_phibins();
    if (Verbosity() > 1)
    {
      std::cout << "    nphibins " << nphibins << std::endl;
    }

    for (unsigned int side = 0; side < 2; ++side)
    {
      if (Verbosity() > 1)
      {
        std::cout << "TPC layer " << layer << " side " << side << std::endl;
      }

      // for this layer and side, use a vector of a vector of cells for each phibin
      if (phi_sorted_hits.size() != static_cast<std::size_t>(nphibins))
      {
        phi_sorted_hits.clear();
        phi_sorted_hits.resize(nphibins);
      }
      else
      {
        for (auto &hits : phi_sorted_hits)
        {
          hits.clear();
        }
      }

      // Loop over all hitsets containing signals for this layer and add them to phi_sorted_hits for their phibin
      TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId, layer);
      for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
           hitset_iter != hitset_range.second;
           ++hitset_iter)
      {
        // we have an iterator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
        // get the hitset key
        TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
        unsigned int this_side = TpcDefs::getSide(hitsetkey);
        // skip this hitset if it is not on this side
        if (this_side != side)
        {
          continue;
        }

        if (Verbosity() > 2)
        {
          if (layer == print_layer)
          {
            std::cout << "new: PHG4TpcDigitizer:  processing signal hits for layer " << layer
                      << " hitsetkey " << hitsetkey << " side " << side << std::endl;
          }
        }

        // get all of the hits from this hitset
        TrkrHitSet *hitset = hitset_iter->second;
        TrkrHitSet::ConstRange hit_range = hitset->getHits();
        for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
             hit_iter != hit_range.second;
             ++hit_iter)
        {
          // Fill the vector of signal hits for each phibin
          unsigned int phibin = TpcDefs::getPad(hit_iter->first);
          phi_sorted_hits[phibin].push_back(hit_iter);
        }
        // For this hitset we now have the signal hits sorted into vectors for each phi
      }

      // Process one phi bin at a time
      if (Verbosity() > 1)
      {
        std::cout << "    phi_sorted_hits size " << phi_sorted_hits.size() << " ntbins " << layergeom->get_zbins() << std::endl;
      }

      for (unsigned int iphi = 0; iphi < phi_sorted_hits.size(); iphi++)
      {
        const unsigned int sector = 12U * iphi / nphibins;
        const TrkrDefs::hitsetkey channel_hitsetkey =
            TpcDefs::genHitSetKey(layer, sector, side);

        // A masked electronics channel has neither signal nor electronic
        // noise. Skip the complete time stream before sampling any noise.
        if (isChannelMasked(channel_hitsetkey, iphi))
        {
          continue;
        }

        int ntbins = layergeom->get_zbins();
        signal_hit_by_tbin.assign(ntbins, nullptr);

        // add a signal hit from phi_sorted_hits for each t bin that has one
        for (unsigned int it = 0; it < phi_sorted_hits[iphi].size(); it++)
        {
          int tbin = TpcDefs::getTBin(phi_sorted_hits[iphi][it]->first);
          if (!signal_hit_by_tbin[tbin])
          {
            signal_hit_by_tbin[tbin] = phi_sorted_hits[iphi][it]->second;
          }

          if (Verbosity() > 2)
          {
            if (layer == print_layer)
            {
              TrkrDefs::hitkey hitkey = phi_sorted_hits[iphi][it]->first;
              std::cout << "iphi " << iphi << " adding existing signal hit for layer " << layer
                        << " side " << side
                        << " tbin " << tbin << "  hitkey " << hitkey
                        << " pad " << TpcDefs::getPad(hitkey)
                        << " t bin " << TpcDefs::getTBin(hitkey)
                        << "  energy " << (phi_sorted_hits[iphi][it]->second)->getEnergy()
                        << std::endl;
            }
          }
        }

        // initialize entries to zero for each t bin
        adc_input.assign(ntbins, 0.0);

        // Now for this phibin we process all bins ordered by t into hits with noise
        //======================================================
        // For this step we take the edep value and convert it to mV at the ADC input
        // See comments above for how to do this for signal and noise

        for (int it = 0; it < ntbins; it++)
        {
          TrkrHit *signal_hit = signal_hit_by_tbin[it];
          if (signal_hit)
          {
            // This tbin has a hit, add noise
            float signal_with_noise = add_noise_to_bin(signal_hit->getEnergy());
            adc_input[it] = signal_with_noise;

            if (Verbosity() > 2)
            {
              if (layer == print_layer)
              {
                std::cout << "existing signal hit: layer " << layer << " iphi " << iphi << " it " << it
                          << " edep " << signal_hit->getEnergy()
                          << " adc gain " << ADCSignalConversionGain
                          << " signal with noise " << signal_with_noise
                          << " adc_input " << adc_input[it] << std::endl;
              }
            }
          }
          else
          {
            if (!skip_noise)
            {
              // This t bin does not have a filled cell, add noise
              float noise = add_noise_to_bin(0.0);
              adc_input[it] = noise;

              if (Verbosity() > 2)
              {
                if (layer == print_layer)
                {
                  std::cout << "noise hit: layer " << layer << " side " << side << " iphi " << iphi << " it " << it
                            << " adc gain " << ADCSignalConversionGain
                            << " noise " << noise
                            << " adc_input " << adc_input[it] << std::endl;
                }
              }
            }
          }
        }

        // Now we can digitize the entire stream of t bins for this phi bin
        int binpointer = 0;

        // Since we now store the local z of the hit as time of arrival at the readout plane,
        // there is no difference between north and south
        // The first to arrive is always bin 0

        for (int it = 0; it < ntbins; it++)
        {
          if (it < binpointer)
          {
            continue;
          }

          // optionally do not trigger on bins with no signal
          if (!signal_hit_by_tbin[it] && skip_noise)
          {
            binpointer++;
            continue;
          }

          // convert threshold in "equivalent electrons" to mV
          if (adc_input[it] > ADCThreshold_mV)
          {
            // digitize this bin and the following 4 bins

            if (Verbosity() > 2)
            {
              if (layer == print_layer)
              {
                std::cout << std::endl
                          << "Hit above threshold of "
                          << ADCThreshold * ADCNoiseConversionGain << " for phibin " << iphi
                          << " it " << it << " with adc_input " << adc_input[it]
                          << " digitize this and 4 following bins: " << std::endl;
              }
            }

            for (int itup = 0; itup < 5; itup++)
            {
              if (it + itup < ntbins && it + itup >= 0)  // stay within the bin limits
              {
                float input = 0;
                TrkrHit *signal_hit = signal_hit_by_tbin[it + itup];
                if (!signal_hit && skip_noise)
                {
                  input = add_noise_to_bin(0.0);  // no noise added to this bin previously because skip_noise is true
                }
                else
                {
                  input = adc_input[it + itup];
                }
                // input voltage x 1024 channels over 2200 mV max range
                unsigned int adc_output = (unsigned int) (input * 1024.0 / 2200.0);

                if (input < 0)
                {
                  adc_output = 0;
                }
                adc_output = std::min<unsigned int>(adc_output, 1023);

                // Get the hitkey
                TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, it + itup);
                TrkrHit *hit = nullptr;

                if (Verbosity() > 2)
                {
                  if (layer == print_layer)
                  {
                    std::cout << "    Digitizing:  iphi " << iphi << "  it+itup " << it + itup
                              << " adc_hitid " << (signal_hit ? hitkey : 0)
                              << " is_populated " << (signal_hit ? 1 : 2)
                              << "  adc_input " << adc_input[it + itup]
                              << " ADCThreshold " << ADCThreshold * ADCNoiseConversionGain
                              << " adc_output " << adc_output
                              << " hitkey " << hitkey
                              << " side " << side
                              << "  binpointer " << binpointer
                              << std::endl;
                  }
                }

                if (signal_hit)
                {
                  // this is a signal hit, it already exists
                  hit = signal_hit;
                }
                else
                {
                  // Hit does not exist yet, have to make one
                  // we need the hitset key, requires (layer, sector, side)
                  auto hitset_iter =
                      trkrhitsetcontainer->findOrAddHitSet(channel_hitsetkey);

                  hit = new TrkrHitv2();
                  hitset_iter->second->addHitSpecificKey(hitkey, hit);

                  if (Verbosity() > 2)
                  {
                    if (layer == print_layer)
                    {
                      std::cout << "      adding noise TrkrHit for iphi " << iphi
                                << " tbin " << it + itup
                                << " side " << side
                                << " created new hit with hitkey " << hitkey
                                << " energy " << adc_input[it + itup] << " adc " << adc_output
                                << "  binpointer " << binpointer
                                << std::endl;
                    }
                  }
                }

                hit->setAdc(adc_output);

              }  // end boundary check
              binpointer++;  // skip this bin in future
            }  // end itup loop

          }  //  adc threshold if
          else
          {
            // set adc value to zero if there is a hit
            TrkrHit *hit = signal_hit_by_tbin[it];
            if (hit)
            {
              hit->setAdc(0);
            }
            // bin below threshold, move on
            binpointer++;
          }  // end adc threshold if/else
        }  // end time bin loop
      }  // end phibins loop
    }  // end loop over sides
  }  // end loop over TPC layers

  //======================================================
  if (Verbosity() > 5)
  {
    std::cout << "From PHG4TpcDigitizer: hitsetcontainer dump at end before cleaning:" << std::endl;
  }
  std::vector<TrkrDefs::hitkey> delete_hitkey_list;

  // Clean up undigitized hits - we want all hitsets for the Tpc
  // This loop is pretty efficient because the remove methods all take a specified hitset as input
  TrkrHitSetContainer::ConstRange hitset_range_now = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range_now.first;
       hitset_iter != hitset_range_now.second;
       ++hitset_iter)
  {
    // we have an iterator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
    TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    const int sector = TpcDefs::getSectorId(hitsetkey);
    const int side = TpcDefs::getSide(hitsetkey);
    if (Verbosity() > 5)
    {
      std::cout << "PHG4TpcDigitizer: hitset with key: " << hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << std::endl;
    }

    // get all of the hits from this hitset
    TrkrHitSet *hitset = hitset_iter->second;
    delete_hitkey_list.clear();
    TrkrHitSet::ConstRange hit_range = hitset->getHits();
    for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
         hit_iter != hit_range.second;
         ++hit_iter)
    {
      TrkrDefs::hitkey hitkey = hit_iter->first;
      TrkrHit *tpchit = hit_iter->second;
      const bool masked =
          isChannelMasked(hitsetkey, TpcDefs::getPad(hitkey));

      if (Verbosity() > 5)
      {
        std::cout << "    layer " << layer << "  hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey)
                  << " t bin " << TpcDefs::getTBin(hitkey)
                  << " adc " << tpchit->getAdc() << std::endl;
      }

      if (tpchit->getAdc() == 0 || masked)
      {
        if (Verbosity() > 20)
        {
          std::cout << "                       --   this hit not digitized - delete it" << std::endl;
        }
        // screws up the iterator to delete it here, store the hitkey for later deletion
        delete_hitkey_list.push_back(hitkey);
      }
    }

    // delete all undigitized hits
    for (auto &hitkey : delete_hitkey_list)
    {
      if (isChannelMasked(hitsetkey, TpcDefs::getPad(hitkey)))
      {
        // A hit may have several contributing G4 hits. The legacy v1
        // implementation removes one association per call, so repeat once
        // for every association reported before deleting the digitized hit.
        TrkrHitTruthAssoc::MMap associations;
        hittruthassoc->getG4Hits(hitsetkey, hitkey, associations);
        for (std::size_t index = 0; index < associations.size(); ++index)
        {
          hittruthassoc->removeAssoc(hitsetkey, hitkey);
        }
      }
      hitset->removeHit(hitkey);
      if (Verbosity() > 20)
      {
        if (layer == print_layer)
        {
          std::cout << "removed hit with hitsetkey " << hitsetkey
                    << " and hitkey " << hitkey << std::endl;
        }
      }

      // should also delete all entries with this hitkey from the TrkrHitTruthAssoc map
      // hittruthassoc->removeAssoc(hitsetkey, hitkey);   // Slow! Commented out by ADF 9/6/2022
    }
  }

  // Final hitset dump
  if (Verbosity() > 5)
  {
    std::cout << "From PHG4TpcDigitizer: hitsetcontainer dump at end after cleaning:" << std::endl;
  }
  // We want all hitsets for the Tpc
  TrkrHitSetContainer::ConstRange hitset_range_final = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range_final.first;
       hitset_iter != hitset_range_final.second;
       ++hitset_iter)
  {
    // we have an itrator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
    TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    if (layer != print_layer)
    {
      continue;
    }
    const int sector = TpcDefs::getSectorId(hitsetkey);
    const int side = TpcDefs::getSide(hitsetkey);
    if (Verbosity() > 5 && layer == print_layer)
    {
      std::cout << "PHG4TpcDigitizer: hitset with key: " << hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << std::endl;
    }

    // get all of the hits from this hitset
    TrkrHitSet *hitset = hitset_iter->second;
    TrkrHitSet::ConstRange hit_range = hitset->getHits();
    for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
         hit_iter != hit_range.second;
         ++hit_iter)
    {
      TrkrDefs::hitkey hitkey = hit_iter->first;
      TrkrHit *tpchit = hit_iter->second;
      if (Verbosity() > 5)
      {
        std::cout << "      LAYER " << layer << " hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) << " t bin " << TpcDefs::getTBin(hitkey)
                  << " adc " << tpchit->getAdc() << std::endl;
      }

      if (tpchit->getAdc() == 0)
      {
        std::cout << "   Oops!                    --   this hit not digitized and not deleted!" << std::endl;
      }
    }
  }

  // hittruthassoc->identify();

  return;
}

bool PHG4TpcDigitizer::isChannelMasked(
    const TrkrDefs::hitsetkey hitsetkey,
    const unsigned int pad) const
{
  return (m_maskDeadChannels && m_deadChannelMap.contains(hitsetkey, pad)) ||
         (m_maskHotChannels && m_hotChannelMap.contains(hitsetkey, pad));
}

float PHG4TpcDigitizer::add_noise_to_bin(float signal)
{
  // add noise to the signal and return adc input voltage
  float adc_input_voltage = signal * ADCSignalConversionGain;  // mV, see comments above
  if (!no_noise)
  {
    // Keep the historical arithmetic unchanged unless an ADC pedestal was
    // explicitly selected. Setting 60 ADC uses the exact voltage 128.90625 mV,
    // avoiding an electrons-to-voltage round trip that can truncate 60 to 59.
    const float noise_voltage = std::isfinite(m_pedestalADC)
                                    ? m_pedestalADC * 2200.0 / 1024.0 + added_noise() * ADCNoiseConversionGain
                                    : (Pedestal + added_noise()) * ADCNoiseConversionGain;
    adc_input_voltage += noise_voltage;
  }

  return adc_input_voltage;
}

float PHG4TpcDigitizer::added_noise()
{
  float noise = gsl_ran_gaussian(RandomGenerator, TpcEnc);

  return noise;
}
