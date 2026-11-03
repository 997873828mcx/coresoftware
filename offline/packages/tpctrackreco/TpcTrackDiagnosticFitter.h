#ifndef TPCTRACKRECO_TPCTRACKDIAGNOSTICFITTER_H
#define TPCTRACKRECO_TPCTRACKDIAGNOSTICFITTER_H

#include "TpcTrackKalmanFitter.h"

#include <cstdint>
#include <utility>
#include <vector>

struct TpcTrackDiagnosticMeasurement
{
  TpcTrackPoint point;
  std::uint64_t cluster_id{0};
  // Geometry radius, never the held-out measurement's measured radius.
  double layer_radius_cm{0.0};
  TpcTrackVec3 cylinder_origin;
  TpcTrackVec3 cylinder_axis{0.0, 0.0, 1.0};
  std::vector<std::pair<std::uint32_t, std::uint32_t>> raw_hits;
};

struct TpcTrackLeaveOneOutResult
{
  bool success{false};
  std::string message;
  std::size_t excluded_measurements{0};
  std::size_t excluded_shared_raw_hit_measurements{0};
  int independent_charge{0};
  bool independent_charge_ambiguous{false};
  std::vector<std::uint64_t> retained_cluster_ids;
  TpcKalmanResult fit;
  TpcKalmanSurfaceResult prediction;
};

class TpcTrackDiagnosticFitter
{
 public:
  // Same numerical recipe as configure_yuri_v0_kalman_tuning(), with every
  // scalar set explicitly. Caller supplies the mapped PHField pointer.
  static TpcKalmanConfig yuri_config(const PHField *field, double bfield_t = 1.4);
  static void assign_yuri_measurement(TpcTrackPoint &point);

  static bool fit(const std::vector<TpcTrackDiagnosticMeasurement> &measurements,
                  int charge, const TpcKalmanConfig &config,
                  TpcKalmanResult &result);

  // The tested point and all measurements sharing any of its raw-hit keys are
  // physically absent from the input passed to fit(), including ordering,
  // helix seeding, path construction, robust decisions and every iteration.
  // The supplied pattern charge is never used by the held-out fit: repeat the
  // production poly fit on remaining points and derive its charge independently.
  static bool leave_one_out(
      const std::vector<TpcTrackDiagnosticMeasurement> &measurements,
      std::size_t held_out, int charge, const TpcKalmanConfig &config,
      TpcTrackLeaveOneOutResult &result, bool fit_allowed = true);
};

#endif
