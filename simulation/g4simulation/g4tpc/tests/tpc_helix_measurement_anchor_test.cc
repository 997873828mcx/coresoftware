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
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kBFieldT = 1.4;

  void require(const bool condition, const std::string &message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  double distance(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
  {
    return TpcTrackHelixFitter::distance(lhs, rhs);
  }

  TpcTrackHelix make_perigee_helix(const int charge,
                                   const double pt_gev,
                                   const double tan_lambda)
  {
    const double direction = -static_cast<double>(charge);
    const double radius = pt_gev / (0.003 * kBFieldT);
    const double theta_perigee = kPi;

    TpcTrackHelix helix;
    helix.cx = radius + 8.0;
    helix.cy = 0.0;
    helix.radius = radius;
    helix.pitch = direction * radius * tan_lambda;
    helix.z0 = -helix.pitch * theta_perigee;
    helix.theta_first = theta_perigee;
    helix.theta_last = theta_perigee;
    helix.theta_min = theta_perigee;
    helix.theta_max = theta_perigee;
    helix.direction = direction;
    helix.bfield_t = kBFieldT;
    return helix;
  }

  std::vector<TpcTrackPoint> make_measurements(const TpcTrackHelix &helix,
                                               const double first_path_cm,
                                               const int npoints,
                                               const double step_cm)
  {
    std::vector<TpcTrackPoint> points;
    points.reserve(static_cast<std::size_t>(npoints));
    for (int index = 0; index < npoints; ++index)
    {
      const double path_cm = first_path_cm + step_cm * static_cast<double>(index);
      const double theta = helix.theta_first +
                           helix.direction * path_cm / helix.radius;
      TpcTrackPoint point;
      point.track_id = 17;
      point.layer = index;
      point.position = TpcTrackHelixFitter::point(helix, theta);
      point.momentum = TpcTrackHelixFitter::momentum(helix, theta);
      point.path = path_cm;
      points.push_back(point);
    }
    return points;
  }

  bool theta_in_range(const double theta,
                      const TpcTrackHelixSearchRange &range)
  {
    return theta >= range.theta_min - 1.0e-10 &&
           theta <= range.theta_max + 1.0e-10;
  }

  void check_anchor_case(const int charge,
                         const bool reverse_input,
                         const double first_path_cm,
                         const char *label)
  {
    const TpcTrackHelix helix = make_perigee_helix(charge, 0.20, 0.12);
    const double circumference_cm = 2.0 * kPi * helix.radius;
    std::vector<TpcTrackPoint> points =
        make_measurements(helix, first_path_cm, 12, 4.0);
    if (reverse_input)
    {
      std::reverse(points.begin(), points.end());
    }

    TpcTrackHelixSearchRange range;
    const bool ok = TpcTrackHelixFitter::measurement_anchored_search_range(
        helix, points, 80.0, 5.0, range);

    std::ostringstream detail;
    detail << "charge=" << charge
           << " reverse=" << reverse_input
           << " branch=" << label
           << " anchor_path_cm=" << range.anchor_path_cm
           << " expected_cm=" << first_path_cm
           << " residual_cm=" << range.anchor_residual_cm;

    require(ok && range.valid, detail.str() + " failed to build range");
    require(std::abs(range.anchor_path_cm - first_path_cm) < 1.0e-6,
            detail.str() + " selected the wrong measured branch");
    require(range.anchor_residual_cm < 1.0e-8,
            detail.str() + " anchor does not close on its measurement");

    const double upstream_vertex_path_cm = first_path_cm - 60.0;
    const double vertex_theta = helix.theta_first +
                                helix.direction * upstream_vertex_path_cm / helix.radius;
    require(theta_in_range(vertex_theta, range),
            detail.str() + " omitted an upstream production point");

    const double range_span_cm =
        (range.theta_max - range.theta_min) * helix.radius;
    require(range_span_cm <= 0.95 * circumference_cm + 1.0e-9,
            detail.str() + " search exceeds the branch-safe turn cap");
    std::cout << "[ok] " << detail.str()
              << " range_span_cm=" << range_span_cm << '\n';
  }

  TpcTrackHelix make_vertex_helix(const TpcTrackVec3 &vertex,
                                  const TpcTrackVec3 &momentum,
                                  const int charge)
  {
    TpcTrackHelix helix;
    require(TpcTrackHelixFitter::from_state(
                vertex, momentum, charge, kBFieldT, helix),
            "failed to construct pair helix");
    return helix;
  }

  void check_pair_pca()
  {
    const TpcTrackVec3 vertex{6.0, -2.0, 4.0};
    const TpcTrackHelix helix1 = make_vertex_helix(
        vertex, {0.90, 0.35, 0.12}, 1);
    const TpcTrackHelix helix2 = make_vertex_helix(
        vertex, {-0.55, 0.80, -0.08}, -1);

    const auto points1 = make_measurements(helix1, 42.0, 14, 3.0);
    const auto points2 = make_measurements(helix2, 38.0, 14, 3.0);
    TpcTrackHelixSearchRange range1;
    TpcTrackHelixSearchRange range2;
    require(TpcTrackHelixFitter::measurement_anchored_search_range(
                helix1, points1, 70.0, 3.0, range1),
            "failed first pair range");
    require(TpcTrackHelixFitter::measurement_anchored_search_range(
                helix2, points2, 70.0, 3.0, range2),
            "failed second pair range");

    const auto candidates = TpcTrackHelixFitter::pca_candidates_in_ranges(
        helix1, helix2, range1, range2, 64, 32);
    require(!candidates.empty(), "measurement-anchored PCA returned no candidates");
    const auto &best = candidates.front();
    const TpcTrackVec3 midpoint = TpcTrackHelixFitter::scale(
        TpcTrackHelixFitter::add(best.pca1, best.pca2), 0.5);

    std::ostringstream detail;
    detail << "pair_dca_cm=" << best.dca
           << " vertex_error_cm=" << distance(midpoint, vertex);
    require(best.dca < 2.0e-4, detail.str() + " pair DCA did not close");
    require(distance(midpoint, vertex) < 2.0e-4,
            detail.str() + " selected the wrong intersection");
    std::cout << "[ok] " << detail.str() << '\n';
  }

  void check_low_pt_turn_cap()
  {
    const TpcTrackHelix helix = make_perigee_helix(1, 0.03, 0.10);
    const auto points = make_measurements(helix, 10.0, 8, 2.0);
    TpcTrackHelixSearchRange range;
    require(TpcTrackHelixFitter::measurement_anchored_search_range(
                helix, points, 80.0, 5.0, range),
            "failed to build low-pT search range");

    const double circumference_cm = 2.0 * kPi * helix.radius;
    const double range_span_cm =
        (range.theta_max - range.theta_min) * helix.radius;
    require(range_span_cm <= 0.95 * circumference_cm + 1.0e-9,
            "low-pT search exceeded the branch-safe turn cap");
    require(range_span_cm >= 0.94 * circumference_cm,
            "low-pT search cap was unexpectedly restrictive");
    require(std::abs(range.downstream_cm - 5.0) < 1.0e-9,
            "low-pT search did not preserve downstream tolerance");
    std::cout << "[ok] low-pt turn cap span_cm=" << range_span_cm
              << " circumference_cm=" << circumference_cm << '\n';
  }
}  // namespace

int main()
{
  try
  {
    for (const int charge : {-1, 1})
    {
      for (const bool reverse_input : {false, true})
      {
        const TpcTrackHelix helix = make_perigee_helix(charge, 0.20, 0.12);
        const double circumference_cm = 2.0 * kPi * helix.radius;
        check_anchor_case(charge, reverse_input, 45.0, "after-perigee");
        check_anchor_case(charge, reverse_input, -45.0, "before-perigee");
        check_anchor_case(charge, reverse_input,
                          circumference_cm + 45.0, "later-turn");
        check_anchor_case(charge, reverse_input,
                          -circumference_cm + 45.0, "earlier-turn");
      }
    }
    check_low_pt_turn_cap();
    check_pair_pca();
  }
  catch (const std::exception &error)
  {
    std::cerr << "[fail] " << error.what() << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "All measurement-anchored helix search tests passed.\n";
  return EXIT_SUCCESS;
}
