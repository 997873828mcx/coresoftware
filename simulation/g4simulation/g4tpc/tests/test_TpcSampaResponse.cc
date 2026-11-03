#include "TpcSampaResponse.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }

  // Independent Simpson integration of the peak-normalized analog function;
  // no production cumulative probability or bin-weight helper is used here.
  double quadrature(double lo, double hi, double peak, int intervals = 2048)
  {
    lo = std::max(0.0, lo);
    hi = std::max(0.0, hi);
    const double step = (hi - lo) / intervals;
    double sum = 0.0;
    for (int index = 0; index <= intervals; ++index)
    {
      const double ratio = (lo + index * step) / peak;
      const double value = std::exp(4.0 - 4.0 * ratio) * std::pow(ratio, 4);
      sum += (index == 0 || index == intervals ? 1.0 : index % 2 == 0 ? 2.0 : 4.0) * value;
    }
    return sum * step / 3.0;
  }

  template<class Function>
  void rejects(Function function)
  {
    bool rejected = false;
    try { function(); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "invalid SAMPA input was accepted");
  }
}

int main()
{
  try
  {
    require(TpcSampaResponse::DefaultPeakingTimeNs == 55.0, "default is not SAMPA v5 55 ns");
    require(std::abs(TpcSampaResponse::binIntegralToPeakScale(55.0, 53.326184) - 1.319809941343262) < 1.e-13,
            "55 ns response-area/clock normalization differs");
    require(std::abs(TpcSampaResponse::binIntegralToPeakScale(80.0, 53.326184) - 1.919723551044745) < 1.e-13,
            "80 ns response-area/clock normalization differs");
    for (const double peak : {55.0, 80.0})
    {
      require(TpcSampaResponse::analogResponse(-1.0, peak) == 0.0, "response is not causal");
      require(std::abs(TpcSampaResponse::analogResponse(peak, peak) - 1.0) < 1.e-14, "analog peak is not unity");
      const double area = TpcSampaResponse::areaNs(peak);
      require(std::abs(quadrature(0.0, 16.0 * peak, peak, 32768) - area) < 1.e-9,
              "full analog response area disagrees with independent quadrature");
      for (const double clock : {50.0, 53.326184})
      {
        double minimum_peak = 1.0;
        double maximum_peak = 0.0;
        for (const double phase : {1.e-6, .0625, .25, .5, .9375, .999999})
        {
          double sum = 0.0;
          double discrete_peak = 0.0;
          for (int bin = 0; bin < TpcSampaResponse::ClockBins; ++bin)
          {
            const double lo = (bin - phase) * clock;
            const double hi = (bin + 1.0 - phase) * clock;
            const double weight = TpcSampaResponse::binAverage(lo, hi, peak, clock);
            require(weight >= 0.0 && weight <= 1.0, "invalid bin-average amplitude");
            require(std::abs(weight - quadrature(lo, hi, peak) / clock) < 2.e-11,
                    "bin-average amplitude differs from independent quadrature");
            sum += weight;
            discrete_peak = std::max(discrete_peak, weight);
          }
          require(sum * clock < area, "eight-clock tail was renormalized into retained bins");
          require(std::abs(sum * clock - quadrature(0.0, (8.0 - phase) * clock, peak, 32768)) < 1.e-9,
                  "retained pulse area differs from independent quadrature");
          minimum_peak = std::min(minimum_peak, discrete_peak);
          maximum_peak = std::max(maximum_peak, discrete_peak);
        }
        require(maximum_peak - minimum_peak > .05, "clock-phase dependence was artificially removed");
        require(maximum_peak < 1.0, "discrete bin-average peak was forced to the analog peak");
      }
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    rejects([] { TpcSampaResponse::areaNs(0.0); });
    rejects([] { TpcSampaResponse::binIntegralToPeakScale(55.0, -1.0); });
    rejects([nan] { TpcSampaResponse::binIntegralToPeakScale(nan, 53.0); });
    rejects([nan] { TpcSampaResponse::analogResponse(nan, 55.0); });
    rejects([] { TpcSampaResponse::binAverage(10.0, 0.0, 55.0, 53.0); });
    rejects([] { TpcSampaResponse::binAverage(0.0, 60.0, 55.0, 53.0); });
    std::cout << "PASS: peak-normalized SAMPA bin averages, exact areas, clock phases and invalid inputs\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
