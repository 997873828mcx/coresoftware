#include "TpcTrackDiagnosticFitter.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <phfield/PHField.h>
#include <CLHEP/Units/SystemOfUnits.h>

namespace
{
  void require(const bool condition, const char *message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  std::vector<TpcTrackDiagnosticMeasurement> measurements(const double pt, const int charge,
                                                         const double phi)
  {
    auto config = TpcTrackDiagnosticFitter::yuri_config(nullptr);
    config.analytic_uniform_propagation = true;
    const std::array<double, 6> truth{{0.0, 0.0, 2.0, phi, -charge / pt, 0.25}};
    std::vector<TpcTrackDiagnosticMeasurement> output;
    for (std::size_t index = 0; index < 30; ++index)
    {
      const auto state = TpcTrackKalmanFitter::propagate_state(truth, 32.0 + 1.5 * index, config);
      TpcTrackDiagnosticMeasurement measurement;
      measurement.cluster_id = 1000 + index * 7;
      measurement.point.layer = 7 + static_cast<int>(index);
      measurement.point.original_index = index;
      measurement.point.position = TpcTrackKalmanFitter::state_position(state);
      measurement.layer_radius_cm = std::hypot(state[0], state[1]);
      measurement.raw_hits.emplace_back(123U, static_cast<unsigned int>(index));
      TpcTrackDiagnosticFitter::assign_yuri_measurement(measurement.point);
      output.push_back(measurement);
    }
    return output;
  }

  void check_track(const double pt, const int charge, const double phi)
  {
    auto input = measurements(pt, charge, phi);
    auto config = TpcTrackDiagnosticFitter::yuri_config(nullptr);
    TpcKalmanResult nominal;
    require(TpcTrackDiagnosticFitter::fit(input, charge, config, nominal), "synthetic nominal fit failed");
    auto quiet_config = config;
    quiet_config.collect_innovation_components = false;
    TpcKalmanResult quiet;
    require(TpcTrackDiagnosticFitter::fit(input, charge, quiet_config, quiet), "quiet nominal fit failed");
    require(quiet.states_smoothed == nominal.states_smoothed &&
                quiet.covs_smoothed == nominal.covs_smoothed &&
                quiet.measurement_used == nominal.measurement_used,
            "diagnostic collection changed nominal fit");

    TpcKalmanSurfaceResult pca;
    require(TpcTrackKalmanFitter::predict_to_beamline_pca(nominal, {}, config, pca),
            "beamline prediction failed");
    require(std::hypot(pca.state[0], pca.state[1]) < 0.02, "beamline PCA position incorrect");
    require(std::abs(1.0 / std::abs(pca.state[4]) / pt - 1.0) < 0.005,
            "beamline PCA momentum incorrect");

    // Both cluster positions are absent from all LOO seeding/path/update
    // inputs. Moving them by meters must leave the retained fit and predicted
    // state/covariance bit-for-bit identical.
    input[11].raw_hits.push_back(input[10].raw_hits.front());
    TpcTrackLeaveOneOutResult first;
    require(TpcTrackDiagnosticFitter::leave_one_out(input, 10, charge, config, first),
            "held-out fit failed");
    require(first.excluded_measurements == 2 && first.excluded_shared_raw_hit_measurements == 1,
            "shared raw hit was not removed");
    TpcTrackLeaveOneOutResult flipped_charge;
    require(TpcTrackDiagnosticFitter::leave_one_out(input,10,-charge,config,flipped_charge),
            "independent charge test failed");
    require(first.independent_charge==flipped_charge.independent_charge &&
                first.fit.states_smoothed==flipped_charge.fit.states_smoothed &&
                first.fit.covs_smoothed==flipped_charge.fit.covs_smoothed &&
                first.prediction.state==flipped_charge.prediction.state &&
                first.prediction.covariance==flipped_charge.prediction.covariance,
            "full-track charge leaked into held-out orientation or prediction");
    input[10].point.position = {900.0, -700.0, 800.0};
    input[11].point.position = {-800.0, 1000.0, -1000.0};
    TpcTrackLeaveOneOutResult moved;
    require(TpcTrackDiagnosticFitter::leave_one_out(input, 10, charge, config, moved),
            "moved held-out fit failed");
    require(first.fit.seed.radius == moved.fit.seed.radius &&
                first.fit.path_s == moved.fit.path_s &&
                first.fit.states_smoothed == moved.fit.states_smoothed &&
                first.fit.covs_smoothed == moved.fit.covs_smoothed &&
                first.fit.measurement_used == moved.fit.measurement_used &&
                first.prediction.state == moved.prediction.state &&
                first.prediction.covariance == moved.prediction.covariance,
            "held-out coordinates leaked into prediction");
    const double radius = std::hypot(first.prediction.state[0], first.prediction.state[1]);
    require(std::abs(radius - input[10].layer_radius_cm) < 1.0e-6,
            "prediction missed geometry cylinder");
    // Surface covariance has no radial variance at a fixed cylinder.
    const double nx = first.prediction.state[0] / radius;
    const double ny = first.prediction.state[1] / radius;
    const auto &cov = first.prediction.covariance;
    const double radial_variance = nx * nx * cov[0] + nx * ny * (cov[1] + cov[6]) + ny * ny * cov[7];
    require(std::abs(radial_variance) < 1.0e-8, "surface covariance omitted intersection dependence");
  }

  struct GradientField : PHField
  {
    void GetFieldValue(const double position[4],double *field) const override
    {
      field[0]=0.0002*(position[2]/CLHEP::cm)*CLHEP::tesla;
      field[1]=0.0001*(position[1]/CLHEP::cm)*CLHEP::tesla;
      field[2]=(1.4+0.0005*(position[0]/CLHEP::cm))*CLHEP::tesla;
    }
  };

  void check_surfaces()
  {
    auto config=TpcTrackDiagnosticFitter::yuri_config(nullptr);
    config.analytic_uniform_propagation=true;
    const std::array<double,6> straight{{30,0,0,0,0,.2}};
    std::array<double,36> covariance{};
    covariance[0]=.01;covariance[7]=.01;covariance[14]=.02;covariance[21]=1e-4;covariance[28]=1e-10;covariance[35]=1e-3;
    TpcKalmanSurfaceResult prediction;
    require(TpcTrackKalmanFitter::predict_to_cylinder(straight,covariance,40,10,0,20,config,prediction),
            "straight-track cylinder failed");
    require(std::abs(prediction.state[0]-40)<1e-7 && std::abs(prediction.state[2]-2)<1e-7,
            "straight-track reference position incorrect");
    // For a straight trajectory, z=z0+(R-x0)*tan(lambda). Check the
    // analytically propagated surface covariance, including dx0 dependence.
    const double expected_z_variance=.02+.2*.2*.01+100*.001;
    require(std::abs(prediction.covariance[14]-expected_z_variance)<1e-5,
            "known straight-track covariance projection incorrect");

    const std::array<double,6> initial{{0,0,2,.7,-1/.8,.25}};
    const auto known=TpcTrackKalmanFitter::propagate_state(initial,50,config);
    const TpcTrackVec3 origin{.3,-.2,1},axis{.01,-.02,1};
    const double dx=known[0]-origin.x,dy=known[1]-origin.y,dz=known[2]-origin.z;
    const double along=(dx*axis.x+dy*axis.y+dz*axis.z)/std::sqrt(axis.x*axis.x+axis.y*axis.y+axis.z*axis.z);
    const double radius=std::sqrt(dx*dx+dy*dy+dz*dz-along*along);
    require(TpcTrackKalmanFitter::predict_to_cylinder(initial,covariance,radius,50,40,60,config,prediction,
              .13957039,origin,axis),"translated/tilted cylinder prediction failed");
    require(std::abs(prediction.path_cm-50)<1e-6,"translated/tilted cylinder reference incorrect");

    GradientField field;
    config.magnetic_field=&field;config.analytic_uniform_propagation=false;
    config.rkn_max_step_cm=.5;config.rkn_variable_field_within_step=true;
    TpcKalmanResult fit;fit.success=true;fit.charge=1;fit.bfield_t=1.4;fit.magnetic_field=&field;
    fit.states_smoothed={TpcTrackKalmanFitter::propagate_state(initial,35,config)};
    fit.covs_smoothed={covariance};
    require(TpcTrackKalmanFitter::predict_to_beamline_pca(fit,{},config,prediction),"mapped-field PCA failed");
    const double stationarity=prediction.state[0]*std::cos(prediction.state[3])+prediction.state[1]*std::sin(prediction.state[3]);
    require(std::abs(stationarity)<1e-6 && std::hypot(prediction.state[0],prediction.state[1])<.001,
            "mapped-field PCA failed actual propagated stationarity");
  }
}  // namespace

int main()
{
  try
  {
    for (const double pt : {0.2, 0.8, 10.0})
    {
      for (const int charge : {-1, 1})
      {
        check_track(pt, charge, 0.7);
        check_track(pt, charge, 3.14159265358979323846 - 1.0e-5);
      }
    }
    check_surfaces();
    auto short_input = measurements(1.0, 1, 0.0);
    short_input.resize(5);
    TpcTrackLeaveOneOutResult failure;
    require(!TpcTrackDiagnosticFitter::leave_one_out(short_input, 0, 1,
              TpcTrackDiagnosticFitter::yuri_config(nullptr), failure),
            "fewer than five independent measurements unexpectedly fitted");
    require(!failure.message.empty() && failure.excluded_measurements == 1,
            "failed refit lost its reason/exclusion count");
    std::cout << "PASS: 12 helices, both charges, low/high pT, phi wrapping, independent LOO, "
                 "shared raw hits, straight-track analytic covariance, translated/tilted cylinder, mapped-field PCA, "
                 "and diagnostic noninterference\n";
    return 0;
  }
  catch (const std::exception &error)
  {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
