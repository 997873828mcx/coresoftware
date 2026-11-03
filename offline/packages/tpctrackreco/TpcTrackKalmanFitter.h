// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCTRACKKALMANFITTER_H
#define TPCTRACKRECO_TPCTRACKKALMANFITTER_H

#include "TpcTrackFit.h"

#include <array>
#include <utility>
#include <vector>

// Diagnostic prediction only: no measurement or vertex constraint is applied.
// Covariance is row-major in (x,y,z,phi,q/pT,tan(lambda)). The Jacobian includes
// the state dependence of the intersection path, rather than a fixed-path
// approximation. Path is the fitter's transverse arc length, in cm.
struct TpcKalmanSurfaceResult
{
  bool success{false};
  std::string message;
  std::array<double, 6> state{};
  std::array<double, 36> covariance{};
  std::array<double, 36> jacobian{};
  double path_cm{0.0};
  unsigned int intersections{0};
};

class TpcTrackKalmanFitter
{
 public:
  enum StateIndex
  {
    X = 0,
    Y = 1,
    Z = 2,
    Phi = 3,
    QOverPt = 4,
    TanLambda = 5,
    StateDim = 6
  };

  static bool fit(const std::vector<TpcTrackPoint> &points,
                  int charge,
                  const TpcKalmanConfig &config,
                  TpcKalmanResult &result,
                  double mass_gev = 0.13957039);

  static TpcTrackVec3 state_position(const std::array<double, StateDim> &state);
  static TpcTrackVec3 state_momentum(const std::array<double, StateDim> &state);
  static TpcTrackVec3 state_tangent(const std::array<double, StateDim> &state);
  static std::array<double, StateDim> propagation_state(const TpcKalmanResult &fit,
                                                        const TpcTrackVec3 &reference_vertex);
  static std::array<double, StateDim> propagate_state(const std::array<double, StateDim> &state,
                                                      double ds_cm,
                                                      const TpcKalmanConfig &config,
                                                      double mass_gev = 0.13957039);
  static std::pair<double, double> dca_to_vertex(const TpcKalmanResult &fit,
                                                 const TpcTrackVec3 &vertex,
                                                 const TpcKalmanConfig *config = nullptr);

  static bool physical_endpoint(const TpcKalmanResult &fit,
                                std::array<double, StateDim> &state,
                                std::array<double, StateDim * StateDim> &covariance);

  // Cylinder roots are enumerated in [path_min_cm,path_max_cm]. Select the
  // root nearest path_hint_cm; ties select the smaller path. The hint and
  // search interval must be chosen without the tested cluster's coordinates.
  static bool predict_to_cylinder(const std::array<double, StateDim> &state,
                                  const std::array<double, StateDim * StateDim> &covariance,
                                  double radius_cm,
                                  double path_hint_cm,
                                  double path_min_cm,
                                  double path_max_cm,
                                  const TpcKalmanConfig &config,
                                  TpcKalmanSurfaceResult &result,
                                  double mass_gev = 0.13957039,
                                  const TpcTrackVec3 &cylinder_origin = {},
                                  const TpcTrackVec3 &cylinder_axis = {0.0, 0.0, 1.0});

  // Find a local transverse minimum near the uniform-field PCA seed and
  // refine using the actual configured propagator. No vertex update is made.
  static bool predict_to_beamline_pca(const TpcKalmanResult &fit,
                                      const TpcTrackVec3 &beamline,
                                      const TpcKalmanConfig &config,
                                      TpcKalmanSurfaceResult &result,
                                      double max_path_cm = 80.0);

  // Propagate a smoothed fit to its transverse PCA to a vertex and apply an
  // optional XY/XYZ Gaussian vertex measurement without modifying the fit.
  static bool constrain_to_vertex(const TpcKalmanResult &fit,
                                  const TpcTrackVec3 &vertex,
                                  const TpcKalmanConfig &config,
                                  double sigma_xy_cm,
                                  double sigma_z_cm,
                                  bool constrain_z,
                                  TpcKalmanVertexConstraintResult &result);
};

#endif
