#ifndef G4TPC_TPCSAMPARESPONSE_H
#define G4TPC_TPCSAMPARESPONSE_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

// The electronics model averages the analog voltage over each ADC clock bin;
// it is not an instantaneous sample at the bin center.  All pulse weights are
// normalized to a unit ANALOG peak, so the digitizer needs only the physical
// charge-to-peak-voltage gain.  No clock-phase-dependent peak correction is used.
namespace TpcSampaResponse
{
  inline constexpr double DefaultPeakingTimeNs = 55.0;
  inline constexpr const char* Mode = "peak_normalized_bin_average_v1";
  inline constexpr int Version = 1;
  inline constexpr int ClockBins = 8;

  inline void validate(double peakingTimeNs, double clockNs)
  {
    if (!std::isfinite(peakingTimeNs) || peakingTimeNs <= 0.0 ||
        !std::isfinite(clockNs) || clockNs <= 0.0)
    {
      throw std::invalid_argument("SAMPA peaking time and ADC clock must be finite and positive");
    }
  }

  // Integral of g(u)=exp(4)*(u/T)^4*exp(-4*u/T), u>=0; g(T)=1.
  inline double areaNs(double peakingTimeNs)
  {
    validate(peakingTimeNs, 1.0);
    return 3.0 * std::exp(4.0) * peakingTimeNs / 128.0;
  }

  inline double binIntegralToPeakScale(double peakingTimeNs, double clockNs)
  {
    validate(peakingTimeNs, clockNs);
    return areaNs(peakingTimeNs) / clockNs;
  }

  inline double analogResponse(double delayNs, double peakingTimeNs)
  {
    validate(peakingTimeNs, 1.0);
    if (!std::isfinite(delayNs))
    {
      throw std::invalid_argument("SAMPA delay must be finite");
    }
    if (delayNs <= 0.0) { return 0.0; }
    const double ratio = delayNs / peakingTimeNs;
    if (!std::isfinite(ratio)) { return 0.0; }
    // The logarithmic form avoids inf*0 in the very late tail.
    return std::exp(4.0 + 4.0 * std::log(ratio) - 4.0 * ratio);
  }

  namespace detail
  {
    // Survival probability of a Gamma(shape=5,scale=1) distribution.
    inline double survival(double x)
    {
      if (x <= 0.0) { return 1.0; }
      if (x >= 1000.0) { return 0.0; }
      return std::exp(-x) * (1.0 + x * (1.0 + x * (0.5 + x * (1.0 / 6.0 + x / 24.0))));
    }

    inline double cumulative(double x)
    {
      if (x <= 0.0) { return 0.0; }
      if (x >= 1.0) { return 1.0 - survival(x); }
      // A positive series avoids cancellation of the first four moments when
      // an arrival occurs immediately before a clock-bin boundary.
      double term = std::pow(x, 5) / 120.0;
      double sum = term;
      for (int order = 6; order < 128; ++order)
      {
        term *= x / order;
        const double next = sum + term;
        if (next == sum) { break; }
        sum = next;
      }
      return std::exp(-x) * sum;
    }
  }  // namespace detail

  inline double binAverage(double lowerDelayNs, double upperDelayNs,
                           double peakingTimeNs, double clockNs)
  {
    validate(peakingTimeNs, clockNs);
    if (!std::isfinite(lowerDelayNs) || !std::isfinite(upperDelayNs) ||
        upperDelayNs < lowerDelayNs ||
        upperDelayNs - lowerDelayNs > clockNs * (1.0 + 1.e-10))
    {
      throw std::invalid_argument("invalid SAMPA ADC integration interval");
    }
    const double lo = std::max(0.0, lowerDelayNs) * 4.0 / peakingTimeNs;
    const double hi = std::max(0.0, upperDelayNs) * 4.0 / peakingTimeNs;
    const double fraction = lo >= 1.0
                                ? detail::survival(lo) - detail::survival(hi)
                                : detail::cumulative(hi) - detail::cumulative(lo);
    return binIntegralToPeakScale(peakingTimeNs, clockNs) * std::max(0.0, fraction);
  }
}  // namespace TpcSampaResponse

#endif
