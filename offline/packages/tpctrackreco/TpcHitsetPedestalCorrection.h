#ifndef TPCTRACKRECO_TPCHITSETPEDESTALCORRECTION_H
#define TPCTRACKRECO_TPCHITSETPEDESTALCORRECTION_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <limits>
#include <string>

class PHCompositeNode;

/** Convert raw simulated TRKR_HITSET ADC to pedestal-subtracted ADC once.
 *
 * Reads the event-level TPC_ADC_BASELINE PdbParameterMap. Legacy inputs without
 * that node require an explicit setInputPedestalADC declaration. Already
 * subtracted data without metadata require setInputAlreadySubtracted(true).
 * Neither source mode is inferred from hit values. Zero ADC hits and all keys
 * and truth associations are retained. Unsigned ADC storage requires clamping
 * at zero and nearest-integer rounding after subtracting a fractional pedestal.
 */
class TpcHitsetPedestalCorrection : public SubsysReco
{
 public:
  explicit TpcHitsetPedestalCorrection(const std::string& name = "TpcHitsetPedestalCorrection");
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setInputPedestalADC(double pedestal);
  void setInputAlreadySubtracted(bool value) { m_inputAlreadySubtracted = value; }

  double getInputPedestalADC() const { return m_lastPedestalADC; }
  double getEffectiveIntegerShiftADC() const;
  std::uint64_t getEventsCorrected() const { return m_eventsCorrected; }
  std::uint64_t getEventsAlreadySubtracted() const { return m_eventsAlreadySubtracted; }
  std::uint64_t getHitsCorrected() const { return m_hitsCorrected; }
  std::uint64_t getHitsClamped() const { return m_hitsClamped; }
  std::uint64_t getHitsZero() const { return m_hitsZero; }

 private:
  double m_declaredPedestalADC {std::numeric_limits<double>::quiet_NaN()};
  double m_lastPedestalADC {std::numeric_limits<double>::quiet_NaN()};
  bool m_inputAlreadySubtracted {false};
  std::uint64_t m_eventsCorrected {0};
  std::uint64_t m_eventsAlreadySubtracted {0};
  std::uint64_t m_hitsCorrected {0};
  std::uint64_t m_hitsClamped {0};
  std::uint64_t m_hitsZero {0};
};

#endif
