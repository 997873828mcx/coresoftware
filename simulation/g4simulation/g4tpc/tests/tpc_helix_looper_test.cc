#include "TpcTrackHelixFitter.h"

#include <algorithm>
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

  TpcTrackHelix make_truth_helix(const int charge)
  {
    const double pt = 0.20;
    const double phi = 0.35;
    const double tan_lambda = 0.12;
    const TpcTrackVec3 position{32.0, -4.0, 3.0};
    const TpcTrackVec3 momentum{
        pt * std::cos(phi),
        pt * std::sin(phi),
        pt * tan_lambda};

    TpcTrackHelix helix;
    require(TpcTrackHelixFitter::from_state(position, momentum, charge, kBFieldT, helix),
            "failed to build truth helix");
    return helix;
  }

  std::vector<TpcTrackPoint> make_points(const TpcTrackHelix &helix)
  {
    std::vector<TpcTrackPoint> points;
    constexpr int npoints = 96;
    constexpr double step_theta = 0.12;
    points.reserve(npoints);
    for (int i = 0; i < npoints; ++i)
    {
      const double theta = helix.theta_first +
                           helix.direction * step_theta * static_cast<double>(i);
      TpcTrackPoint point;
      point.track_id = 13;
      point.shower_id = 9;
      point.layer = i % 48;
      point.position = TpcTrackHelixFitter::point(helix, theta);
      point.momentum = TpcTrackHelixFitter::momentum(helix, theta);
      point.t = static_cast<double>(i);
      point.path = static_cast<double>(i);
      points.push_back(point);
    }
    return points;
  }

  TpcTrackHelix fit_or_throw(std::vector<TpcTrackPoint> points,
                             const int charge,
                             const std::string &label)
  {
    TpcTrackHelixFitter::order_points(points, TpcTrackPointOrder::Auto);

    TpcTrackHelix fitted;
    if (!TpcTrackHelixFitter::fit(points, 0, kBFieldT, fitted))
    {
      throw std::runtime_error(label + ": helix fit failed");
    }
    if (!TpcTrackHelixFitter::orient_to_charge(fitted, charge))
    {
      throw std::runtime_error(label + ": charge orientation failed");
    }
    return fitted;
  }

  void run_case(const int charge, const bool reverse_input)
  {
    const TpcTrackHelix truth = make_truth_helix(charge);
    std::vector<TpcTrackPoint> points = make_points(truth);
    if (reverse_input)
    {
      std::reverse(points.begin(), points.end());
    }

    const std::string label = std::string("helix looper charge=") +
                              charge_label(charge) +
                              (reverse_input ? " reversed-input" : " physical-input");
    const TpcTrackHelix fitted = fit_or_throw(points, charge, label);

    const TpcTrackVec3 selected_pos =
        TpcTrackHelixFitter::point(fitted, fitted.theta_first);
    const TpcTrackVec3 selected_mom =
        TpcTrackHelixFitter::momentum(fitted, fitted.theta_first);
    const TpcTrackVec3 truth_pos =
        TpcTrackHelixFitter::point(truth, truth.theta_first);
    const TpcTrackVec3 truth_mom =
        TpcTrackHelixFitter::momentum(truth, truth.theta_first);

    const double pos_error_cm = distance(selected_pos, truth_pos);
    const double mom_cos = cosine(selected_mom, truth_mom);
    const double expected_direction = -static_cast<double>(charge);

    std::ostringstream detail;
    detail << label
           << " pos_error_cm=" << pos_error_cm
           << " mom_cos=" << mom_cos
           << " direction=" << fitted.direction;

    require(pos_error_cm < 1.0e-6, detail.str() + " selected wrong endpoint");
    require(mom_cos > 0.999999, detail.str() + " selected wrong momentum direction");
    require(fitted.direction * expected_direction > 0.0,
            detail.str() + " fitted direction has wrong charge convention");

    const double vertex_theta = truth.theta_first - truth.direction * 0.25;
    const TpcTrackVec3 vertex = TpcTrackHelixFitter::point(truth, vertex_theta);
    const auto dca = TpcTrackHelixFitter::helix_dca_to_vertex(fitted, vertex);

    std::ostringstream dca_detail;
    dca_detail << label
               << " dca_xy_cm=" << dca.first
               << " dca_z_cm=" << dca.second;
    require(std::isfinite(dca.first) && dca.first < 1.0e-6,
            dca_detail.str() + " transverse DCA should stay on the generated branch");
    require(std::isfinite(dca.second) && dca.second < 1.0e-6,
            dca_detail.str() + " DCA z should stay on the generated branch");

    std::cout << "[ok] " << detail.str()
              << " dca_xy_cm=" << dca.first
              << " dca_z_cm=" << dca.second << '\n';
  }
}  // namespace

int main()
{
  try
  {
    for (const int charge : {-1, 1})
    {
      run_case(charge, false);
      run_case(charge, true);
    }
  }
  catch (const std::exception &error)
  {
    std::cerr << "[fail] " << error.what() << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "All TPC helix looper endpoint tests passed.\n";
  return EXIT_SUCCESS;
}
