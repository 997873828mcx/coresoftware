#include "TpcTrackDiagnosticFitter.h"
#include "Tpc_FittingTools.h"

#include <algorithm>
#include <cmath>
#include <set>

TpcKalmanConfig TpcTrackDiagnosticFitter::yuri_config(const PHField *field,
                                                    const double bfield_t)
{
  TpcKalmanConfig config;
  config.bfield_t = bfield_t;
  config.magnetic_field = field;
  config.analytic_uniform_propagation = false;
  config.rkn_max_step_cm = 10.0;
  config.rkn_step_tolerance = 0.0;
  config.rkn_max_step_trials = 1;
  config.rkn_max_total_steps = 2000;
  config.rkn_variable_field_within_step = false;
  config.rkn_fast_field_jacobian = true;
  config.rkn_fast_field_pca = true;
  config.rkn_field_pca_refine_iterations = 1;
  config.point_order = TpcTrackPointOrder::RadiusDescending;
  config.meas_sigma_rphi_cm = 0.08;
  config.meas_sigma_r_cm = 0.30;
  config.meas_sigma_z_cm = 0.23;
  config.min_measurement_sigma_cm = 1.0e-6;
  config.initial_sigma_pos_cm = 0.05;
  config.initial_sigma_phi = 0.010;
  config.initial_sigma_qop_t = 0.20;
  config.use_relative_initial_sigma_qop_t = true;
  config.initial_sigma_qop_t_relative = 0.05;
  config.initial_sigma_qop_t_floor = 5.0e-4;
  config.initial_sigma_tanl = 0.020;
  config.process_sigma_pos_cm = 0.0;
  config.process_sigma_phi = 0.0;
  config.process_sigma_qop_t = 0.0;
  config.process_sigma_tanl = 0.0;
  config.collect_innovation_components = true;
  config.robust_measurement_updates = true;
  config.robust_chi2_soft = 16.0;
  config.robust_chi2_hard = 100.0;
  config.fit_iterations = 1;
  config.material_x0_per_cm = 0.0;
  config.multiple_scattering_scale = 0.0;
  config.energy_loss_gev_per_cm = 0.0;
  config.energy_loss_sigma_fraction = 0.0;
  config.min_pt_gev = 0.05;
  return config;
}

void TpcTrackDiagnosticFitter::assign_yuri_measurement(TpcTrackPoint &point)
{
  const int layer = point.layer;
  const bool transition = layer == 22 || layer == 23 || layer == 38 || layer == 39;
  point.use_in_kalman = !transition;
  if (transition)
  {
    point.sigma_rphi_cm = 0.200;
    point.sigma_r_cm = 0.500;
    point.sigma_z_cm = 0.200;
  }
  else if (layer >= 7 && layer <= 54)
  {
    point.sigma_rphi_cm = layer <= 22 ? 0.040 : 0.025;
    point.sigma_r_cm = layer >= 23 && layer <= 38 ? 0.050 : 0.200;
    point.sigma_z_cm = 0.100;
  }
  else
  {
    point.sigma_rphi_cm = 0.08;
    point.sigma_r_cm = 0.30;
    point.sigma_z_cm = 0.23;
  }
}

bool TpcTrackDiagnosticFitter::fit(
    const std::vector<TpcTrackDiagnosticMeasurement> &measurements,
    const int charge, const TpcKalmanConfig &config, TpcKalmanResult &result)
{
  std::vector<TpcTrackPoint> points;
  points.reserve(measurements.size());
  for (const auto &measurement : measurements)
  {
    points.push_back(measurement.point);
  }
  return TpcTrackKalmanFitter::fit(points, charge, config, result);
}

bool TpcTrackDiagnosticFitter::leave_one_out(
    const std::vector<TpcTrackDiagnosticMeasurement> &measurements,
    const std::size_t held_out, const int /*pattern_charge*/, const TpcKalmanConfig &config,
    TpcTrackLeaveOneOutResult &result, const bool fit_allowed)
{
  result = {};
  if (held_out >= measurements.size())
  {
    result.message = "invalid held-out measurement index";
    return false;
  }
  const auto &held = measurements[held_out];
  const std::set<std::pair<std::uint32_t, std::uint32_t>> held_hits(
      held.raw_hits.begin(), held.raw_hits.end());
  std::vector<TpcTrackDiagnosticMeasurement> remaining;
  remaining.reserve(measurements.size() - 1);
  for (std::size_t index = 0; index < measurements.size(); ++index)
  {
    const auto &measurement = measurements[index];
    const bool overlapping = std::any_of(measurement.raw_hits.begin(), measurement.raw_hits.end(),
        [&](const auto &hit) { return held_hits.count(hit) != 0; });
    if (index == held_out || overlapping)
    {
      ++result.excluded_measurements;
      if (index != held_out)
      {
        ++result.excluded_shared_raw_hit_measurements;
      }
      continue;
    }
    remaining.push_back(measurement);
    result.retained_cluster_ids.push_back(measurement.cluster_id);
  }
  if (!fit_allowed)
  {
    result.message="input polytrack unavailable; exclusion metadata only";
    return false;
  }
  // Charge participates in the physical endpoint/orientation choice even
  // though it does not enter the raw Kalman updates. Recompute the exact
  // PolyTrackReco charge rule from remaining points, preserving their pattern
  // association order. No full-fit charge or state can leak into prediction.
  std::vector<Tpc_FittingTools::Point> charge_points;
  for (const auto &measurement:remaining)
  {
    Tpc_FittingTools::Point point;
    point.x=measurement.point.position.x;point.y=measurement.point.position.y;point.z=measurement.point.position.z;
    charge_points.push_back(point);
  }
  Tpc_FittingTools::FitResult charge_fit;
  if (!Tpc_FittingTools::fit(charge_points,charge_fit) ||
      !std::isfinite(charge_fit.curvature) || charge_fit.curvature==0.0)
  {
    result.message="independent remaining-point poly charge fit failed";
    return false;
  }
  result.independent_charge=charge_fit.curvature>=0.0?-1:1;
  bool radial_rise=false,radial_fall=false;
  for (std::size_t index=1;index<charge_points.size();++index)
  {
    const auto &a=charge_points[index-1];const auto &b=charge_points[index];
    const double dr=std::hypot(b.x,b.y)-std::hypot(a.x,a.y);
    radial_rise=radial_rise || dr>0.1;radial_fall=radial_fall || dr<-.1;
  }
  const auto &first=charge_points.front();const auto &middle=charge_points[charge_points.size()/2];const auto &last=charge_points.back();
  const double ax=middle.x-first.x,ay=middle.y-first.y,bx=last.x-middle.x,by=last.y-middle.y;
  const double cross_scale=std::hypot(ax,ay)*std::hypot(bx,by);
  result.independent_charge_ambiguous=(radial_rise&&radial_fall) ||
      cross_scale==0.0 || std::abs(ax*by-ay*bx)<1e-8*cross_scale;
  if (!fit(remaining, result.independent_charge, config, result.fit))
  {
    result.message = "held-out fit: " + result.fit.message;
    return false;
  }
  std::array<double, 6> state{};
  std::array<double, 36> covariance{};
  if (!TpcTrackKalmanFitter::physical_endpoint(result.fit, state, covariance))
  {
    result.message = "held-out fit has no finite physical endpoint";
    return false;
  }
  // Select by geometry radius and the independently fitted endpoint. No
  // tested measurement coordinate or full-fit state enters branch selection.
  const double dx = state[0] - held.cylinder_origin.x;
  const double dy = state[1] - held.cylinder_origin.y;
  const double dz = state[2] - held.cylinder_origin.z;
  const double axis_norm = std::sqrt(held.cylinder_axis.x * held.cylinder_axis.x +
                                     held.cylinder_axis.y * held.cylinder_axis.y +
                                     held.cylinder_axis.z * held.cylinder_axis.z);
  const double along = (dx * held.cylinder_axis.x + dy * held.cylinder_axis.y +
                        dz * held.cylinder_axis.z) / axis_norm;
  const double endpoint_radius = std::sqrt(std::max(0.0, dx * dx + dy * dy + dz * dz - along * along));
  const double hint = held.layer_radius_cm - endpoint_radius;
  if (!TpcTrackKalmanFitter::predict_to_cylinder(
          state, covariance, held.layer_radius_cm, hint, -80.0, 120.0,
          config, result.prediction, result.fit.mass_gev, held.cylinder_origin, held.cylinder_axis))
  {
    result.message = "held-out prediction: " + result.prediction.message;
    return false;
  }
  result.success = true;
  result.message = "ok";
  return true;
}
