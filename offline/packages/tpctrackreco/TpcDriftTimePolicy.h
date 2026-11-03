#ifndef TPCTRACKRECO_TPCDRIFTTIMEPOLICY_H
#define TPCTRACKRECO_TPCDRIFTTIMEPOLICY_H

#include <algorithm>
#include <cmath>

// No global time shift. last_stored_time_ns is (N-1)*stored_step_ns, not
// the time of a discarded terminal Euler step or the true membrane crossing.
// Crossing-finder legacy edge behavior is retained independently of the
// explicitly enabled simulation-only CM response-tail extension.
struct TpcDriftTimeDecision
{
  bool accepted{false};
  bool central_membrane_clamped{false};
  double time_ns{0.0};
};

inline TpcDriftTimeDecision resolve_tpc_drift_time(
    const double requested_time_ns, const double last_stored_time_ns,
    const double legacy_edge_tolerance_ns, const double cm_tail_tolerance_ns,
    const bool terminal_status_is_central_membrane)
{
  TpcDriftTimeDecision result;
  if (!std::isfinite(requested_time_ns) || !std::isfinite(last_stored_time_ns) ||
      !std::isfinite(legacy_edge_tolerance_ns) || !std::isfinite(cm_tail_tolerance_ns) ||
      requested_time_ns<0.0 || last_stored_time_ns<0.0 ||
      legacy_edge_tolerance_ns<0.0 || cm_tail_tolerance_ns<0.0)
  {
    return result;
  }
  if (requested_time_ns<=last_stored_time_ns)
  {
    result.accepted=true;
    result.time_ns=requested_time_ns;
    return result;
  }
  result.central_membrane_clamped=terminal_status_is_central_membrane &&
      cm_tail_tolerance_ns>0.0 && requested_time_ns<=last_stored_time_ns+cm_tail_tolerance_ns;
  result.accepted=requested_time_ns<=last_stored_time_ns+legacy_edge_tolerance_ns ||
                  result.central_membrane_clamped;
  if (result.accepted) { result.time_ns=last_stored_time_ns; }
  return result;
}

#endif
