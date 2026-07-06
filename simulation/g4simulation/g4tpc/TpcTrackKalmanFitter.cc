#include "TpcTrackKalmanFitter.h"

#include "TpcTrackHelixFitter.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  constexpr double kPi = 3.14159265358979323846;

  template <class T>
  constexpr T square(const T &value)
  {
    return value * value;
  }

  double normalize_phi(const double phi)
  {
    return std::atan2(std::sin(phi), std::cos(phi));
  }

  using StateVector = Eigen::Matrix<double, 6, 1>;
  using StateMatrix = Eigen::Matrix<double, 6, 6>;
  using MeasurementVector = Eigen::Matrix<double, 3, 1>;
  using MeasurementMatrix = Eigen::Matrix<double, 3, 6>;
  using MeasurementCov = Eigen::Matrix<double, 3, 3>;

  StateVector to_eigen(const std::array<double, 6> &state)
  {
    StateVector output;
    for (int i = 0; i < 6; ++i)
    {
      output(i) = state[static_cast<std::size_t>(i)];
    }
    return output;
  }

  std::array<double, 6> to_array(const StateVector &state)
  {
    std::array<double, 6> output{};
    for (int i = 0; i < 6; ++i)
    {
      output[static_cast<std::size_t>(i)] = state(i);
    }
    return output;
  }

  std::array<double, 36> to_array(const StateMatrix &matrix)
  {
    std::array<double, 36> output{};
    for (int row = 0; row < 6; ++row)
    {
      for (int col = 0; col < 6; ++col)
      {
        output[static_cast<std::size_t>(6 * row + col)] = matrix(row, col);
      }
    }
    return output;
  }

  StateVector residual(const StateVector &lhs, const StateVector &rhs)
  {
    StateVector diff = lhs - rhs;
    diff(TpcTrackKalmanFitter::Phi) = normalize_phi(diff(TpcTrackKalmanFitter::Phi));
    return diff;
  }

  double omega_from_state(const StateVector &state, const double bfield_t)
  {
    return 0.003 * bfield_t * state(TpcTrackKalmanFitter::QOverPt);
  }

  StateVector apply_mean_energy_loss(const StateVector &state,
                                     const double ds_cm,
                                     const TpcKalmanConfig &config,
                                     const double mass_gev)
  {
    if (config.energy_loss_gev_per_cm <= 0.0 || ds_cm == 0.0)
    {
      return state;
    }

    StateVector output = state;
    const double qop_t = output(TpcTrackKalmanFitter::QOverPt);
    if (std::abs(qop_t) < 1.0e-12)
    {
      return output;
    }

    const double tanl = output(TpcTrackKalmanFitter::TanLambda);
    const double path3d_cm = std::abs(ds_cm) * std::sqrt(1.0 + tanl * tanl);
    if (path3d_cm <= 0.0)
    {
      return output;
    }

    const double pt = 1.0 / std::abs(qop_t);
    const double momentum = pt * std::sqrt(1.0 + tanl * tanl);
    const double energy = std::sqrt(momentum * momentum + mass_gev * mass_gev);
    const double signed_loss = std::copysign(config.energy_loss_gev_per_cm * path3d_cm, ds_cm);
    const double new_energy = std::max(mass_gev + 1.0e-9, energy - signed_loss);
    const double new_momentum = std::sqrt(std::max(0.0, new_energy * new_energy - mass_gev * mass_gev));
    if (new_momentum <= 0.0)
    {
      return output;
    }

    const double min_pt = std::max(config.min_pt_gev, 1.0e-6);
    const double new_pt = std::max(min_pt, new_momentum / std::sqrt(1.0 + tanl * tanl));
    output(TpcTrackKalmanFitter::QOverPt) = std::copysign(1.0 / new_pt, qop_t);
    return output;
  }

  StateVector propagate_eigen(const StateVector &state,
                              const double ds_cm,
                              const TpcKalmanConfig &config,
                              const double mass_gev)
  {
    StateVector output = state;
    const double phi = state(TpcTrackKalmanFitter::Phi);
    const double omega = omega_from_state(state, config.bfield_t);

    if (std::abs(omega) < 1.0e-10)
    {
      output(TpcTrackKalmanFitter::X) += ds_cm * std::cos(phi);
      output(TpcTrackKalmanFitter::Y) += ds_cm * std::sin(phi);
    }
    else
    {
      const double phi2 = phi + omega * ds_cm;
      output(TpcTrackKalmanFitter::X) += (std::sin(phi2) - std::sin(phi)) / omega;
      output(TpcTrackKalmanFitter::Y) += -(std::cos(phi2) - std::cos(phi)) / omega;
      output(TpcTrackKalmanFitter::Phi) = phi2;
    }

    output(TpcTrackKalmanFitter::Z) += state(TpcTrackKalmanFitter::TanLambda) * ds_cm;
    output(TpcTrackKalmanFitter::Phi) = normalize_phi(output(TpcTrackKalmanFitter::Phi));
    return apply_mean_energy_loss(output, ds_cm, config, mass_gev);
  }

  StateMatrix transport_jacobian(const StateVector &state,
                                 const double ds_cm,
                                 const TpcKalmanConfig &config,
                                 const double mass_gev)
  {
    StateMatrix jac = StateMatrix::Identity();
    const double scales[6] = {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-6};
    for (int col = 0; col < 6; ++col)
    {
      const double step = scales[col] * std::max(1.0, std::abs(state(col)));
      StateVector plus = state;
      StateVector minus = state;
      plus(col) += step;
      minus(col) -= step;
      if (col == TpcTrackKalmanFitter::Phi)
      {
        plus(col) = normalize_phi(plus(col));
        minus(col) = normalize_phi(minus(col));
      }

      const StateVector f_plus = propagate_eigen(plus, ds_cm, config, mass_gev);
      const StateVector f_minus = propagate_eigen(minus, ds_cm, config, mass_gev);
      jac.col(col) = residual(f_plus, f_minus) / (2.0 * step);
    }
    return jac;
  }

  double multiple_scattering_theta0(const StateVector &state,
                                    const double ds_cm,
                                    const TpcKalmanConfig &config,
                                    const double mass_gev)
  {
    if (config.material_x0_per_cm <= 0.0 || ds_cm == 0.0)
    {
      return 0.0;
    }

    const double qop_t = state(TpcTrackKalmanFitter::QOverPt);
    if (std::abs(qop_t) < 1.0e-12)
    {
      return 0.0;
    }

    const double tanl = state(TpcTrackKalmanFitter::TanLambda);
    const double path3d_cm = std::abs(ds_cm) * std::sqrt(1.0 + tanl * tanl);
    const double x_over_x0 = config.material_x0_per_cm * path3d_cm;
    if (x_over_x0 <= 0.0)
    {
      return 0.0;
    }

    const double pt = 1.0 / std::abs(qop_t);
    const double momentum = pt * std::sqrt(1.0 + tanl * tanl);
    const double energy = std::sqrt(momentum * momentum + mass_gev * mass_gev);
    const double beta = (energy > 0.0) ? momentum / energy : 0.0;
    if (beta <= 0.0 || momentum <= 0.0)
    {
      return 0.0;
    }

    const double log_term = 1.0 + 0.038 * std::log(x_over_x0);
    return config.multiple_scattering_scale * 0.0136 / (beta * momentum) *
           std::sqrt(x_over_x0) * log_term;
  }

  StateMatrix process_noise(const StateVector &state,
                            const double ds_cm,
                            const TpcKalmanConfig &config,
                            const double mass_gev)
  {
    const double scale = std::max(1.0, std::abs(ds_cm));
    StateMatrix noise = StateMatrix::Zero();
    noise(TpcTrackKalmanFitter::X, TpcTrackKalmanFitter::X) = square(config.process_sigma_pos_cm * scale);
    noise(TpcTrackKalmanFitter::Y, TpcTrackKalmanFitter::Y) = square(config.process_sigma_pos_cm * scale);
    noise(TpcTrackKalmanFitter::Z, TpcTrackKalmanFitter::Z) = square(config.process_sigma_pos_cm * scale);
    noise(TpcTrackKalmanFitter::Phi, TpcTrackKalmanFitter::Phi) = square(config.process_sigma_phi * scale);
    noise(TpcTrackKalmanFitter::QOverPt, TpcTrackKalmanFitter::QOverPt) = square(config.process_sigma_qop_t * scale);
    noise(TpcTrackKalmanFitter::TanLambda, TpcTrackKalmanFitter::TanLambda) = square(config.process_sigma_tanl * scale);

    const double theta0 = multiple_scattering_theta0(state, ds_cm, config, mass_gev);
    if (theta0 > 0.0)
    {
      noise(TpcTrackKalmanFitter::Phi, TpcTrackKalmanFitter::Phi) += square(theta0);
      noise(TpcTrackKalmanFitter::TanLambda, TpcTrackKalmanFitter::TanLambda) +=
          square(theta0 * std::sqrt(1.0 + square(state(TpcTrackKalmanFitter::TanLambda))));
    }

    if (config.energy_loss_sigma_fraction > 0.0 && config.energy_loss_gev_per_cm > 0.0)
    {
      const double tanl = state(TpcTrackKalmanFitter::TanLambda);
      const double path3d_cm = std::abs(ds_cm) * std::sqrt(1.0 + tanl * tanl);
      const double loss_sigma = config.energy_loss_sigma_fraction *
                                config.energy_loss_gev_per_cm * path3d_cm;
      const double qop_t = state(TpcTrackKalmanFitter::QOverPt);
      if (std::abs(qop_t) > 1.0e-12)
      {
        const double pt = 1.0 / std::abs(qop_t);
        const double momentum = pt * std::sqrt(1.0 + tanl * tanl);
        if (momentum > 0.0)
        {
          noise(TpcTrackKalmanFitter::QOverPt, TpcTrackKalmanFitter::QOverPt) +=
              square(qop_t * loss_sigma / momentum);
        }
      }
    }

    return noise;
  }

  MeasurementCov measurement_covariance(const TpcTrackPoint &point,
                                        const double var_rphi,
                                        const double var_r,
                                        const double var_z)
  {
    const double radius = std::hypot(point.position.x, point.position.y);
    const double cos_phi = (radius > 0.0) ? point.position.x / radius : 1.0;
    const double sin_phi = (radius > 0.0) ? point.position.y / radius : 0.0;

    MeasurementCov cov = MeasurementCov::Zero();
    cov(0, 0) = var_r * square(cos_phi) + var_rphi * square(sin_phi);
    cov(1, 1) = var_r * square(sin_phi) + var_rphi * square(cos_phi);
    cov(0, 1) = (var_r - var_rphi) * sin_phi * cos_phi;
    cov(1, 0) = cov(0, 1);
    cov(2, 2) = var_z;
    return cov;
  }

  bool make_seed(const std::vector<TpcTrackPoint> &points,
                 const TpcKalmanConfig &config,
                 TpcTrackHelix &seed,
                 std::vector<double> &theta_values,
                 std::vector<double> &path_s)
  {
    if (!TpcTrackHelixFitter::fit(points, 0, config.bfield_t, seed))
    {
      return false;
    }

    theta_values.clear();
    theta_values.reserve(points.size());
    for (const auto &point : points)
    {
      double theta = std::atan2(point.position.y - seed.cy, point.position.x - seed.cx);
      if (!theta_values.empty())
      {
        while (theta - theta_values.back() > kPi)
        {
          theta -= 2.0 * kPi;
        }
        while (theta - theta_values.back() < -kPi)
        {
          theta += 2.0 * kPi;
        }
      }
      theta_values.push_back(theta);
    }

    if (theta_values.size() < 2)
    {
      return false;
    }

    path_s.assign(theta_values.size(), 0.0);
    for (std::size_t i = 1; i < theta_values.size(); ++i)
    {
      path_s[i] = path_s[i - 1] + seed.radius * std::abs(theta_values[i] - theta_values[i - 1]);
    }
    return true;
  }

  bool initial_state(const std::vector<TpcTrackPoint> &points,
                     const TpcTrackHelix &seed,
                     const std::vector<double> &theta_values,
                     const TpcKalmanConfig &config,
                     StateVector &state)
  {
    if (points.empty() || theta_values.empty())
    {
      return false;
    }

    const double direction = seed.direction;
    const double theta0 = theta_values.front();
    const double denom = 0.003 * config.bfield_t * seed.radius;
    if (std::abs(denom) <= 0.0 || !std::isfinite(denom))
    {
      return false;
    }

    state.setZero();
    state(TpcTrackKalmanFitter::X) = points.front().position.x;
    state(TpcTrackKalmanFitter::Y) = points.front().position.y;
    state(TpcTrackKalmanFitter::Z) = points.front().position.z;
    state(TpcTrackKalmanFitter::Phi) = normalize_phi(std::atan2(direction * std::cos(theta0),
                                                               direction * -std::sin(theta0)));
    state(TpcTrackKalmanFitter::QOverPt) = direction / denom;
    state(TpcTrackKalmanFitter::TanLambda) = direction * seed.pitch / seed.radius;

    const double max_abs_qop_t = 1.0 / std::max(config.min_pt_gev, 1.0e-6);
    if (std::abs(state(TpcTrackKalmanFitter::QOverPt)) > max_abs_qop_t)
    {
      state(TpcTrackKalmanFitter::QOverPt) = std::copysign(max_abs_qop_t,
                                                          state(TpcTrackKalmanFitter::QOverPt));
    }
    return true;
  }
}  // namespace

bool TpcTrackKalmanFitter::fit(const std::vector<TpcTrackPoint> &input_points,
                               const int charge,
                               const TpcKalmanConfig &config,
                               TpcKalmanResult &result,
                               const double mass_gev)
{
  result = TpcKalmanResult{};
  result.charge = charge;
  result.bfield_t = config.bfield_t;
  result.mass_gev = mass_gev;

  if (input_points.size() < 5)
  {
    result.message = "need at least five TPC points";
    return false;
  }

  std::vector<TpcTrackPoint> points = input_points;
  TpcTrackHelixFitter::order_points(points, config.point_order);

  std::vector<double> theta_values;
  if (!make_seed(points, config, result.seed, theta_values, result.path_s))
  {
    result.message = "helix seed failed";
    return false;
  }

  StateVector state;
  if (!initial_state(points, result.seed, theta_values, config, state))
  {
    result.message = "initial state failed";
    return false;
  }

  const std::size_t npoints = points.size();
  std::vector<StateVector> states_filtered(npoints);
  std::vector<StateVector> states_predicted(npoints);
  std::vector<StateMatrix> covs_filtered(npoints);
  std::vector<StateMatrix> covs_predicted(npoints);
  std::vector<StateMatrix> transport(npoints);

  const double min_meas = std::max(config.min_measurement_sigma_cm, 1.0e-12);
  const double sigma_rphi = std::max(config.meas_sigma_rphi_cm, min_meas);
  const double sigma_r = std::max(config.meas_sigma_r_cm, min_meas);
  const double sigma_z = std::max(config.meas_sigma_z_cm, min_meas);
  const double var_rphi = square(sigma_rphi);
  const double var_r = square(sigma_r);
  const double var_z = square(sigma_z);

  MeasurementMatrix hmat = MeasurementMatrix::Zero();
  hmat(0, X) = 1.0;
  hmat(1, Y) = 1.0;
  hmat(2, Z) = 1.0;

  StateMatrix cov = StateMatrix::Zero();
  cov(X, X) = square(config.initial_sigma_pos_cm);
  cov(Y, Y) = square(config.initial_sigma_pos_cm);
  cov(Z, Z) = square(config.initial_sigma_pos_cm);
  cov(Phi, Phi) = square(config.initial_sigma_phi);
  cov(QOverPt, QOverPt) = square(config.initial_sigma_qop_t);
  cov(TanLambda, TanLambda) = square(config.initial_sigma_tanl);

  double chi2 = 0.0;
  int ndof = 0;
  const StateMatrix eye = StateMatrix::Identity();

  for (std::size_t index = 0; index < npoints; ++index)
  {
    StateVector pred_state = state;
    StateMatrix pred_cov = cov;
    StateMatrix fmat = eye;
    if (index > 0)
    {
      const double ds = result.path_s[index] - result.path_s[index - 1];
      fmat = transport_jacobian(state, ds, config, mass_gev);
      pred_state = propagate_eigen(state, ds, config, mass_gev);
      pred_cov = fmat * cov * fmat.transpose() + process_noise(state, ds, config, mass_gev);
      pred_cov = 0.5 * (pred_cov + pred_cov.transpose()).eval();
    }

    MeasurementVector measurement;
    measurement << points[index].position.x, points[index].position.y, points[index].position.z;
    const MeasurementCov meas_cov = measurement_covariance(points[index], var_rphi, var_r, var_z);
    const MeasurementVector meas_residual = measurement - hmat * pred_state;
    const MeasurementCov innovation = hmat * pred_cov * hmat.transpose() + meas_cov;
    const MeasurementCov innovation_inv =
        innovation.completeOrthogonalDecomposition().solve(MeasurementCov::Identity());
    const Eigen::Matrix<double, 6, 3> gain = pred_cov * hmat.transpose() * innovation_inv;
    state = pred_state + gain * meas_residual;
    state(Phi) = normalize_phi(state(Phi));
    cov = (eye - gain * hmat) * pred_cov * (eye - gain * hmat).transpose() +
          gain * meas_cov * gain.transpose();
    cov = 0.5 * (cov + cov.transpose()).eval();

    chi2 += (meas_residual.transpose() * innovation_inv * meas_residual)(0, 0);
    ndof += 3;

    states_predicted[index] = pred_state;
    covs_predicted[index] = pred_cov;
    transport[index] = fmat;
    states_filtered[index] = state;
    covs_filtered[index] = cov;
  }

  std::vector<StateVector> states_smoothed = states_filtered;
  std::vector<StateMatrix> covs_smoothed = covs_filtered;
  for (int index = static_cast<int>(npoints) - 2; index >= 0; --index)
  {
    const StateMatrix pred_inv = covs_predicted[static_cast<std::size_t>(index + 1)]
                                     .completeOrthogonalDecomposition()
                                     .solve(StateMatrix::Identity());
    const StateMatrix smoother_gain =
        covs_filtered[static_cast<std::size_t>(index)] *
        transport[static_cast<std::size_t>(index + 1)].transpose() *
        pred_inv;
    const StateVector smooth_residual = residual(states_smoothed[static_cast<std::size_t>(index + 1)],
                                                 states_predicted[static_cast<std::size_t>(index + 1)]);
    states_smoothed[static_cast<std::size_t>(index)] =
        states_filtered[static_cast<std::size_t>(index)] + smoother_gain * smooth_residual;
    states_smoothed[static_cast<std::size_t>(index)](Phi) =
        normalize_phi(states_smoothed[static_cast<std::size_t>(index)](Phi));
    covs_smoothed[static_cast<std::size_t>(index)] =
        covs_filtered[static_cast<std::size_t>(index)] +
        smoother_gain *
            (covs_smoothed[static_cast<std::size_t>(index + 1)] -
             covs_predicted[static_cast<std::size_t>(index + 1)]) *
            smoother_gain.transpose();
    covs_smoothed[static_cast<std::size_t>(index)] =
        0.5 * (covs_smoothed[static_cast<std::size_t>(index)] +
               covs_smoothed[static_cast<std::size_t>(index)].transpose())
                  .eval();
  }

  result.states_filtered.reserve(npoints);
  result.covs_filtered.reserve(npoints);
  result.states_smoothed.reserve(npoints);
  result.covs_smoothed.reserve(npoints);
  for (std::size_t index = 0; index < npoints; ++index)
  {
    result.states_filtered.push_back(to_array(states_filtered[index]));
    result.covs_filtered.push_back(to_array(covs_filtered[index]));
    result.states_smoothed.push_back(to_array(states_smoothed[index]));
    result.covs_smoothed.push_back(to_array(covs_smoothed[index]));
  }

  result.chi2 = chi2;
  result.ndof = ndof - StateDim;
  result.success = true;
  result.message = "ok";
  return true;
}

TpcTrackVec3 TpcTrackKalmanFitter::state_position(const std::array<double, StateDim> &state)
{
  return {state[X], state[Y], state[Z]};
}

TpcTrackVec3 TpcTrackKalmanFitter::state_momentum(const std::array<double, StateDim> &state)
{
  const double qop_t = state[QOverPt];
  const double pt = (std::abs(qop_t) < 1.0e-12) ? 1.0e12 : 1.0 / std::abs(qop_t);
  return {
      pt * std::cos(state[Phi]),
      pt * std::sin(state[Phi]),
      pt * state[TanLambda]};
}

TpcTrackVec3 TpcTrackKalmanFitter::state_tangent(const std::array<double, StateDim> &state)
{
  return {
      std::cos(state[Phi]),
      std::sin(state[Phi]),
      state[TanLambda]};
}

std::array<double, TpcTrackKalmanFitter::StateDim> TpcTrackKalmanFitter::propagate_state(
    const std::array<double, StateDim> &state,
    const double ds_cm,
    const TpcKalmanConfig &config,
    const double mass_gev)
{
  return to_array(propagate_eigen(to_eigen(state), ds_cm, config, mass_gev));
}

std::pair<double, double> TpcTrackKalmanFitter::dca_to_vertex(const TpcKalmanResult &fit,
                                                             const TpcTrackVec3 &vertex)
{
  if (!fit.success || fit.states_smoothed.empty())
  {
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

  const auto &state_array = fit.states_smoothed.front();
  const StateVector state = to_eigen(state_array);
  const double omega = omega_from_state(state, fit.bfield_t);
  if (std::abs(omega) < 1.0e-10)
  {
    return TpcTrackHelixFitter::line_dca_to_vertex(state_position(state_array),
                                                   state_momentum(state_array),
                                                   vertex);
  }

  const double radius = 1.0 / omega;
  const double center_x = state(X) - std::sin(state(Phi)) / omega;
  const double center_y = state(Y) + std::cos(state(Phi)) / omega;
  const double vx = vertex.x - center_x;
  const double vy = vertex.y - center_y;
  const double distance_to_center = std::sqrt(vx * vx + vy * vy);
  double theta_closest = std::atan2(state(Y) - center_y, state(X) - center_x);
  double dca_xy = std::abs(radius);
  if (distance_to_center > 0.0)
  {
    const double closest_x = center_x + std::abs(radius) * vx / distance_to_center;
    const double closest_y = center_y + std::abs(radius) * vy / distance_to_center;
    const double theta0 = std::atan2(state(Y) - center_y, state(X) - center_x);
    const double theta_raw = std::atan2(closest_y - center_y, closest_x - center_x);
    theta_closest = theta0 + normalize_phi(theta_raw - theta0);
    dca_xy = std::abs(distance_to_center - std::abs(radius));
  }

  const double theta0 = std::atan2(state(Y) - center_y, state(X) - center_x);
  const double dtheta = normalize_phi(theta_closest - theta0);
  const double s_cm = dtheta / omega;
  TpcKalmanConfig config;
  config.bfield_t = fit.bfield_t;
  const auto closest = propagate_state(state_array, s_cm, config, fit.mass_gev);
  return {dca_xy, std::abs(closest[Z] - vertex.z)};
}
