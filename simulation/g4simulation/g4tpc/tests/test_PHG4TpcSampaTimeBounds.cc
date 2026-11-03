#include "PHG4TpcPadPlaneReadout.h"

#include <g4detectors/PHG4CylinderCellDefs.h>
#include <g4detectors/PHG4TpcGeomv2.h>
#include <phool/recoConsts.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Existing production test access, without adding a public geometry setter.
struct PHG4TpcPadPlaneReadoutTestAccess
{
  static std::pair<std::vector<int>, std::vector<double>> response(
      PHG4TpcPadPlaneReadout &readout, PHG4TpcGeomv2 &geometry, double arrival)
  {
    readout.LayerGeom = &geometry;
    std::pair<std::vector<int>, std::vector<double>> result;
    readout.sampaTimeDistribution(arrival, result.first, result.second);
    return result;
  }
  static void peakingTime(PHG4TpcPadPlaneReadout& readout, double value)
  {
    readout.Ts = value;
  }
};

namespace
{
  void require(bool condition, const char *message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }

  double binAverageReference(double lo, double hi, double peak, double clock)
  {
    lo = std::max(0.0, lo);
    hi = std::max(0.0, hi);
    constexpr int intervals = 256;
    const double step = (hi - lo) / intervals;
    double sum = 0.0;
    for (int index = 0; index <= intervals; ++index)
    {
      const double ratio = (lo + index * step) / peak;
      const double value = std::exp(4.0 - 4.0 * ratio) * std::pow(ratio, 4);
      sum += (index == 0 || index == intervals ? 1.0 : index % 2 == 0 ? 2.0 : 4.0) * value;
    }
    return sum * step / (3.0 * clock);
  }
}

int main(int argc, char **argv)
{
  try
  {
    require(argc == 1 || argc == 3, "usage: test_PHG4TpcSampaTimeBounds [--dump-normal PATH]");
    const bool normal_only = argc == 3;
    require(!normal_only || std::string(argv[1]) == "--dump-normal", "unknown test option");
    std::ofstream dump;
    if (normal_only) { dump.open(argv[2], std::ios::binary); require(bool(dump), "cannot open normal response dump"); }
    recoConsts::instance()->set_IntFlag("RANDOMSEED", 202610171);
    PHG4TpcPadPlaneReadout readout("SampaTimeBounds");
    std::size_t cases = 0;
    require(readout.GetSampaPeakingTimeNs() == 55.0, "readout default is not 55 ns");
    for (const double peak : {55.0, 80.0})
    {
      PHG4TpcPadPlaneReadoutTestAccess::peakingTime(readout, peak);
      for (const int bins : {512, 528})
      {
        for (const double clock : {50.0, 53.326184})
        {
          for (const double origin : {0.0, -125.0})
          {
            PHG4TpcGeomv2 clipped, extended;
            clipped.set_binning(PHG4CylinderCellDefs::sizebinning);
            extended.set_binning(PHG4CylinderCellDefs::sizebinning);
            clipped.set_zbins(bins); extended.set_zbins(bins + 16);
            clipped.set_zmin(origin); extended.set_zmin(origin);
            clipped.set_zstep(clock); extended.set_zstep(clock);
            std::vector<int> starts{0, 1, 200, bins - 9, bins - 8};
            if (!normal_only)
            {
              for (int first = bins - 7; first < bins; ++first) { starts.push_back(first); }
            }
            for (const int start : starts)
            {
              for (const double phase : {.0625, .25, .5, .9375})
              {
                const double arrival = origin + (start + phase) * clock;
                const auto full = PHG4TpcPadPlaneReadoutTestAccess::response(readout, extended, arrival);
                const auto result = PHG4TpcPadPlaneReadoutTestAccess::response(readout, clipped, arrival);
                const int first = clipped.get_zbin(arrival);
                require(result.first.size() == static_cast<std::size_t>(std::min(8, bins - first)), "incorrect clipped pulse length");
                require(result.first.size() == result.second.size(), "bin/share dimensions differ");
                double total = 0.0;
                for (std::size_t index = 0; index < result.first.size(); ++index)
                {
                  require(result.first[index] >= 0 && result.first[index] < bins, "out-of-range ADC bin emitted");
                  require(result.first[index] == first + static_cast<int>(index), "pulse bin order changed");
                  require(result.first[index] == full.first[index], "ADC clipping moved a retained bin");
                  require(std::memcmp(&result.second[index], &full.second[index], sizeof(double)) == 0,
                          "ADC clipping changed a retained pulse share");
                  require(std::isfinite(result.second[index]) && result.second[index] >= 0.0, "invalid pulse share");
                  const double center = clipped.get_zcenter(result.first[index]);
                  require(std::abs(result.second[index] - binAverageReference(
                      center - clock / 2.0 - arrival, center + clock / 2.0 - arrival, peak, clock)) < 2.e-8,
                      "production readout differs from independent voltage-average integral");
                  total += result.second[index];
                  if (normal_only)
                  {
                    dump.write(reinterpret_cast<const char *>(&result.first[index]), sizeof(int));
                    dump.write(reinterpret_cast<const char *>(&result.second[index]), sizeof(double));
                  }
                }
                require(total > 0.0, "pulse lost all integrated response");
                ++cases;
              }
            }
          }
        }
      }
    }
    if (normal_only) { dump.close(); require(bool(dump), "normal response dump failed"); }
    std::cout << "PASS: " << cases << " SAMPA cases; all ADC bins valid and clipping preserves retained shares exactly\n";
    return 0;
  }
  catch (const std::exception &error)
  {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
