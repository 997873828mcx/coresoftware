#include "TpcTrackKalmanFitter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  constexpr double kBFieldT = 1.4;
  constexpr double kPionMassGeV = 0.13957039;

  using State = std::array<double, TpcTrackKalmanFitter::StateDim>;

  double square(const double value)
  {
    return value * value;
  }

  double distance(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
  {
    return std::sqrt(square(lhs.x - rhs.x) +
                     square(lhs.y - rhs.y) +
                     square(lhs.z - rhs.z));
  }

  double norm(const TpcTrackVec3 &value)
  {
    return std::sqrt(square(value.x) + square(value.y) + square(value.z));
  }

  double cosine(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
  {
    const double denom = norm(lhs) * norm(rhs);
    if (denom <= 0.0)
    {
      return -1.0;
    }
    return (lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z) / denom;
  }

  void require(const bool condition, const std::string &message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  std::string charge_label(const int charge)
  {
    return charge > 0 ? "+1" : "-1";
  }

  TpcKalmanConfig make_config(const TpcTrackPointOrder order,
                              const bool analytic_uniform)
  {
    TpcKalmanConfig config;
    config.bfield_t = kBFieldT;
    config.analytic_uniform_propagation = analytic_uniform;
    config.point_order = order;
    config.meas_sigma_rphi_cm = 0.03;
    config.meas_sigma_r_cm = 0.30;
    config.meas_sigma_z_cm = 0.05;
    config.initial_sigma_pos_cm = 0.2;
    config.initial_sigma_phi = 0.3;
    config.initial_sigma_qop_t = 0.5;
    config.initial_sigma_tanl = 0.3;
    config.process_sigma_pos_cm = 1.0e-5;
    config.process_sigma_phi = 1.0e-6;
    config.process_sigma_qop_t = 1.0e-7;
    config.process_sigma_tanl = 1.0e-7;
    config.collect_innovation_components = true;
    return config;
  }

  State make_initial_state(const int charge, const double pt_gev, const double tan_lambda)
  {
    State state{};
    state[TpcTrackKalmanFitter::X] = 32.0;
    state[TpcTrackKalmanFitter::Y] = -4.0;
    state[TpcTrackKalmanFitter::Z] = 3.0;
    state[TpcTrackKalmanFitter::Phi] = 0.35;
    state[TpcTrackKalmanFitter::QOverPt] =
        -static_cast<double>(charge) / pt_gev;
    state[TpcTrackKalmanFitter::TanLambda] = tan_lambda;
    return state;
  }

  TpcTrackPoint make_point(const State &state,
                           const int index,
                           const double path_cm)
  {
    TpcTrackPoint point;
    point.track_id = 11;
    point.shower_id = 7;
    point.layer = index % 48;
    point.position = TpcTrackKalmanFitter::state_position(state);
    point.momentum = TpcTrackKalmanFitter::state_momentum(state);
    point.t = path_cm;
    point.path = path_cm;
    return point;
  }

  std::vector<TpcTrackPoint> make_points(const State &initial_state,
                                         const TpcKalmanConfig &config,
                                         const int npoints,
                                         const double step_cm)
  {
    std::vector<TpcTrackPoint> points;
    points.reserve(static_cast<std::size_t>(npoints));
    for (int i = 0; i < npoints; ++i)
    {
      const double s_cm = step_cm * static_cast<double>(i);
      const State state = TpcTrackKalmanFitter::propagate_state(
          initial_state, s_cm, config, kPionMassGeV);
      points.push_back(make_point(state, i, s_cm));
    }
    return points;
  }

  void check_innovation_diagnostics(const TpcKalmanResult &result,
                                    const std::size_t npoints,
                                    const std::string &label)
  {
    const auto require_size = [&](const std::size_t size, const std::string &name)
    {
      require(size == npoints, label + ": wrong size for " + name);
    };

    require_size(result.measurement_chi2.size(), "measurement_chi2");
    require_size(result.measurement_used.size(), "measurement_used");
    require_size(result.measurement_in_seed.size(), "measurement_in_seed");
    require_size(result.innovation_residual_r.size(), "innovation_residual_r");
    require_size(result.innovation_residual_rphi.size(), "innovation_residual_rphi");
    require_size(result.innovation_residual_z.size(), "innovation_residual_z");
    require_size(result.prediction_sigma_r.size(), "prediction_sigma_r");
    require_size(result.prediction_sigma_rphi.size(), "prediction_sigma_rphi");
    require_size(result.prediction_sigma_z.size(), "prediction_sigma_z");
    require_size(result.innovation_sigma_r.size(), "innovation_sigma_r");
    require_size(result.innovation_sigma_rphi.size(), "innovation_sigma_rphi");
    require_size(result.innovation_sigma_z.size(), "innovation_sigma_z");
    require_size(result.innovation_rho_r_rphi.size(), "innovation_rho_r_rphi");
    require_size(result.innovation_rho_r_z.size(), "innovation_rho_r_z");
    require_size(result.innovation_rho_rphi_z.size(), "innovation_rho_rphi_z");
    require_size(result.innovation_whitened_0.size(), "innovation_whitened_0");
    require_size(result.innovation_whitened_1.size(), "innovation_whitened_1");
    require_size(result.innovation_whitened_2.size(), "innovation_whitened_2");

    for (std::size_t index = 0; index < npoints; ++index)
    {
      require(result.measurement_used[index] == 1U,
              label + ": synthetic measurement was not accepted");
      require(result.measurement_in_seed[index] == 1U,
              label + ": global seed marker should include every measurement");
      require(result.innovation_sigma_r[index] > 0.0 &&
                  result.innovation_sigma_rphi[index] > 0.0 &&
                  result.innovation_sigma_z[index] > 0.0,
              label + ": innovation sigma is not positive");
      require(std::abs(result.innovation_rho_r_rphi[index]) <= 1.0 &&
                  std::abs(result.innovation_rho_r_z[index]) <= 1.0 &&
                  std::abs(result.innovation_rho_rphi_z[index]) <= 1.0,
              label + ": innovation correlation is outside [-1,1]");

      const double whitened_chi2 =
          square(result.innovation_whitened_0[index]) +
          square(result.innovation_whitened_1[index]) +
          square(result.innovation_whitened_2[index]);
      const double tolerance = 1.0e-7 * std::max(1.0, result.measurement_chi2[index]);
      require(std::isfinite(whitened_chi2) &&
                  std::abs(whitened_chi2 - result.measurement_chi2[index]) <= tolerance,
              label + ": whitened components do not reproduce incremental chi2");
    }
  }

  TpcKalmanResult fit_or_throw(const std::vector<TpcTrackPoint> &points,
                               const int charge,
                               const TpcKalmanConfig &config,
                               const std::string &label)
  {
    TpcKalmanResult result;
    if (!TpcTrackKalmanFitter::fit(points, charge, config, result, kPionMassGeV))
    {
      std::ostringstream message;
      message << label << ": fit failed: " << result.message;
      throw std::runtime_error(message.str());
    }
    require(!result.states_smoothed.empty(), label + ": no smoothed states");
    check_innovation_diagnostics(result, points.size(), label);
    return result;
  }

  void check_physical_start(const TpcKalmanResult &fit,
                            const State &truth_initial_state,
                            const TpcTrackVec3 &truth_vertex,
                            const int charge,
                            const std::string &label)
  {
    const State selected = TpcTrackKalmanFitter::propagation_state(fit, truth_vertex);
    const TpcTrackVec3 selected_pos = TpcTrackKalmanFitter::state_position(selected);
    const TpcTrackVec3 selected_mom = TpcTrackKalmanFitter::state_momentum(selected);
    const TpcTrackVec3 truth_pos = TpcTrackKalmanFitter::state_position(truth_initial_state);
    const TpcTrackVec3 truth_mom = TpcTrackKalmanFitter::state_momentum(truth_initial_state);

    const double pos_error_cm = distance(selected_pos, truth_pos);
    const double mom_cos = cosine(selected_mom, truth_mom);
    const double physical_qop_sign = -static_cast<double>(charge);
    const double qop_t = selected[TpcTrackKalmanFitter::QOverPt];

    std::ostringstream detail;
    detail << label
           << " pos_error_cm=" << pos_error_cm
           << " mom_cos=" << mom_cos
           << " qop_t=" << qop_t;

    require(pos_error_cm < 5.0, detail.str() + " selected wrong endpoint");
    require(mom_cos > 0.98, detail.str() + " selected wrong momentum direction");
    require(qop_t * physical_qop_sign > 0.0,
            detail.str() + " selected state has wrong physical q/pT sign");

    const auto dca = TpcTrackKalmanFitter::dca_to_vertex(fit, truth_vertex);
    std::ostringstream dca_detail;
    dca_detail << label
               << " dca_xy_cm=" << dca.first
               << " dca_z_cm=" << dca.second;
    require(std::isfinite(dca.first) && dca.first < 0.25,
            dca_detail.str() + " transverse DCA should stay on the generated branch");
    require(std::isfinite(dca.second) && dca.second < 0.50,
            dca_detail.str() + " DCA z should stay on the generated branch");

    std::cout << "[ok] " << detail.str()
              << " dca_xy_cm=" << dca.first
              << " dca_z_cm=" << dca.second << '\n';
  }

  void run_looper_case(const int charge, const bool reverse_input,
                       const bool analytic_uniform)
  {
    const TpcKalmanConfig config =
        make_config(TpcTrackPointOrder::Auto, analytic_uniform);
    const State truth_initial_state = make_initial_state(charge, 0.20, 0.12);
    const State truth_vertex_state = TpcTrackKalmanFitter::propagate_state(
        truth_initial_state, -18.0, config, kPionMassGeV);
    const TpcTrackVec3 truth_vertex =
        TpcTrackKalmanFitter::state_position(truth_vertex_state);

    std::vector<TpcTrackPoint> points =
        make_points(truth_initial_state, config, 96, 8.0);
    if (reverse_input)
    {
      std::reverse(points.begin(), points.end());
    }

    const std::string label = std::string(analytic_uniform ? "analytic " : "rk ") +
                              "looper charge=" +
                              charge_label(charge) +
                              (reverse_input ? " reversed-input" : " physical-input");
    const TpcKalmanResult result = fit_or_throw(points, charge, config, label);
    check_physical_start(result, truth_initial_state, truth_vertex, charge, label);
  }

  void run_high_pt_case(const int charge, const bool reverse_input,
                        const bool analytic_uniform)
  {
    const TpcKalmanConfig config =
        make_config(TpcTrackPointOrder::Auto, analytic_uniform);
    const State truth_initial_state = make_initial_state(charge, 8.0, -0.05);
    const State truth_vertex_state = TpcTrackKalmanFitter::propagate_state(
        truth_initial_state, -12.0, config, kPionMassGeV);
    const TpcTrackVec3 truth_vertex =
        TpcTrackKalmanFitter::state_position(truth_vertex_state);

    std::vector<TpcTrackPoint> points =
        make_points(truth_initial_state, config, 24, 5.0);
    if (reverse_input)
    {
      std::reverse(points.begin(), points.end());
    }

    const std::string label = std::string(analytic_uniform ? "analytic " : "rk ") +
                              "high-pt charge=" +
                              charge_label(charge) +
                              (reverse_input ? " reversed-input" : " physical-input");
    const TpcKalmanResult result = fit_or_throw(points, charge, config, label);
    check_physical_start(result, truth_initial_state, truth_vertex, charge, label);
  }
}  // namespace

int main()
{
  try
  {
    for (const bool analytic_uniform : {false, true})
    {
      for (const int charge : {-1, 1})
      {
        run_looper_case(charge, false, analytic_uniform);
        run_looper_case(charge, true, analytic_uniform);
        run_high_pt_case(charge, false, analytic_uniform);
        run_high_pt_case(charge, true, analytic_uniform);
      }
    }
  }
  catch (const std::exception &error)
  {
    std::cerr << "[fail] " << error.what() << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "All TPC Kalman looper endpoint tests passed.\n";
  return EXIT_SUCCESS;
}
