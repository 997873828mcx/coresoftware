#include "TpcHitsetPedestalCorrection.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <tpc/TpcAdcBaselineCorrection.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

TpcHitsetPedestalCorrection::TpcHitsetPedestalCorrection(const std::string& name)
  : SubsysReco(name)
{
}

void TpcHitsetPedestalCorrection::setInputPedestalADC(const double pedestal)
{
  if (!std::isfinite(pedestal) || pedestal < 0.0 || pedestal > 1023.0)
  {
    throw std::invalid_argument("TPC input pedestal must be in [0,1023] ADC");
  }
  m_declaredPedestalADC = pedestal;
}

double TpcHitsetPedestalCorrection::getEffectiveIntegerShiftADC() const
{
  return std::isfinite(m_lastPedestalADC) ? std::ceil(m_lastPedestalADC - 0.5) : m_lastPedestalADC;
}

int TpcHitsetPedestalCorrection::InitRun(PHCompositeNode* topNode)
{
  if (!findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET"))
  {
    std::cerr << Name() << ": TRKR_HITSET is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcHitsetPedestalCorrection::process_event(PHCompositeNode* topNode)
{
  TpcAdcBaseline::Options options;
  options.declared_pedestal_adc = m_declaredPedestalADC;
  options.input_already_subtracted = m_inputAlreadySubtracted;
  options.correction_name = Name();
  const auto result = TpcAdcBaseline::normalize(topNode, options);
  if (!result.ok)
  {
    std::cerr << Name() << ": " << result.error << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_lastPedestalADC = result.pedestal_adc;
  m_eventsCorrected += result.corrected;
  m_eventsAlreadySubtracted += result.already_subtracted;
  m_hitsCorrected += result.hits_corrected;
  m_hitsClamped += result.hits_clamped;
  m_hitsZero += result.hits_zero;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcHitsetPedestalCorrection::End(PHCompositeNode*)
{
  std::cout << Name() << ": corrected_events=" << m_eventsCorrected
            << " already_subtracted_events=" << m_eventsAlreadySubtracted
            << " corrected_hits=" << m_hitsCorrected
            << " clamped_hits=" << m_hitsClamped
            << " zero_hits=" << m_hitsZero
            << " pedestal_adc=" << m_lastPedestalADC << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
