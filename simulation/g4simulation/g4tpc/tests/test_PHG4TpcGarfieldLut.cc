#include "../PHG4TpcDistortion.h"

#ifndef PHG4TPCDISTORTION_GARFIELD_COMPACT_LUT
#error Compact LUT adapters require an explicit native capability marker
#endif
#ifndef PHG4TPCDISTORTION_GARFIELD_COMPACT_LIMITS
#error Compact capacity adapters require an explicit native capability marker
#endif

#include <phgarfield/PHGarfieldForwardMap.h>
#include <phgarfield/PHGarfieldAdaptiveMap.h>
#include <phgarfield/PHGarfieldCompactMap.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TH3I.h>
#include <TNamed.h>
#include <TKey.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TSeqCollection.h>

#include <unistd.h>
#include <sys/wait.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace
{
  constexpr int PhiNodes = 3;
  constexpr int RadialNodes = 2;
  constexpr int LongitudinalNodes = 2;
  constexpr double InnerRadiusCm = 20.0;
  constexpr double OuterRadiusCm = 30.0;
  constexpr double CentralMembraneHalfWidthCm = 1.0;
  constexpr double ReadoutAbsZCm = 100.0;
  constexpr double TwoPi = 2.0 * std::numbers::pi;

  struct TemporaryFile
  {
    ~TemporaryFile()
    {
      if (!filename.empty())
      {
        std::remove(filename.c_str());
      }
    }

    std::string filename;
  };

  std::uint64_t fnv1a64(const std::string& text)
  {
    std::uint64_t value = 14695981039346656037ULL;
    for (const unsigned char character : text)
    {
      value ^= character;
      value *= 1099511628211ULL;
    }
    return value;
  }

  std::string hexadecimal(const std::uint64_t value)
  {
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(16) << value;
    return output.str();
  }

  template <class T>
  void write_parameter(const char* name, const T value)
  {
    TParameter<T> parameter(name, value);
    parameter.Write();
  }

  void write_named(const char* name, const std::string& value)
  {
    TNamed named(name, value.c_str());
    named.Write();
  }

  template <class Histogram>
  void fill_guards(Histogram& histogram)
  {
    const int nx = histogram.GetNbinsX();
    const int ny = histogram.GetNbinsY();
    const int nz = histogram.GetNbinsZ();
    for (int ix = 1; ix <= nx; ++ix)
    {
      const int source_x =
          ix == 1 ? nx - 1 : (ix == nx ? 2 : ix);
      for (int iy = 1; iy <= ny; ++iy)
      {
        const int source_y =
            iy == 1 ? 2 : (iy == ny ? ny - 1 : iy);
        for (int iz = 1; iz <= nz; ++iz)
        {
          if (ix != 1 && ix != nx &&
              iy != 1 && iy != ny &&
              iz != 1 && iz != nz)
          {
            continue;
          }
          const int source_z =
              iz == 1 ? 2 : (iz == nz ? nz - 1 : iz);
          histogram.SetBinContent(
              ix,
              iy,
              iz,
              histogram.GetBinContent(source_x, source_y, source_z));
        }
      }
    }
  }

  bool reserve_temporary_filename(TemporaryFile& temporary)
  {
    std::string pattern = "/tmp/PHG4TpcGarfieldLutTest.XXXXXX.root";
    const int descriptor = mkstemps(pattern.data(), 5);
    if (descriptor < 0)
    {
      return false;
    }
    close(descriptor);
    temporary.filename = pattern;
    if (std::remove(pattern.c_str()) != 0)
    {
      return false;
    }
    return true;
  }

  bool write_synthetic_map(const std::string& filename, const bool extended = false, const bool numerical = false,
                           const PHGarfieldAdaptiveMap::Data* adaptive = nullptr)
  {
    std::unique_ptr<TFile> output(TFile::Open(filename.c_str(), "RECREATE"));
    if (!output || output->IsZombie())
    {
      return false;
    }

    TDirectory* transport = output->mkdir("GarfieldFullTransport");
    if (!transport)
    {
      return false;
    }
    transport->cd();

    const double phi_spacing = TwoPi / static_cast<double>(PhiNodes);
    const double radial_spacing =
        (OuterRadiusCm - InnerRadiusCm) /
        static_cast<double>(RadialNodes - 1);
    const std::array<double, 2> z_min = {
        -ReadoutAbsZCm,
        CentralMembraneHalfWidthCm};
    const std::array<double, 2> z_max = {
        -CentralMembraneHalfWidthCm,
        ReadoutAbsZCm};
    const std::array<std::string, 2> suffix = {"negz", "posz"};
    double maximum_drift_time_ns = 0.0;

    for (unsigned int side = 0; side < 2; ++side)
    {
      const double z_spacing =
          (z_max[side] - z_min[side]) /
          static_cast<double>(LongitudinalNodes - 1);
      const auto make_node_histogram =
          [&](const std::string& name)
      {
        auto histogram = std::make_unique<TH3D>(
            name.c_str(),
            name.c_str(),
            PhiNodes + 2,
            -1.5 * phi_spacing,
            TwoPi + 0.5 * phi_spacing,
            RadialNodes + 2,
            InnerRadiusCm - 1.5 * radial_spacing,
            OuterRadiusCm + 1.5 * radial_spacing,
            LongitudinalNodes + 2,
            z_min[side] - 1.5 * z_spacing,
            z_max[side] + 1.5 * z_spacing);
        histogram->SetDirectory(nullptr);
        return histogram;
      };

      auto delta_r =
          make_node_histogram("hDeltaR_" + suffix[side]);
      auto delta_rphi =
          make_node_histogram("hDeltaRPhi_" + suffix[side]);
      auto drift_time =
          make_node_histogram("hDriftTime_" + suffix[side]);
      auto reach = std::make_unique<TH3I>(
          ("hReachesReadoutNode_" + suffix[side]).c_str(),
          "synthetic node reach",
          PhiNodes + 2,
          -1.5 * phi_spacing,
          TwoPi + 0.5 * phi_spacing,
          RadialNodes + 2,
          InnerRadiusCm - 1.5 * radial_spacing,
          OuterRadiusCm + 1.5 * radial_spacing,
          LongitudinalNodes + 2,
          z_min[side] - 1.5 * z_spacing,
          z_max[side] + 1.5 * z_spacing);
      reach->SetDirectory(nullptr);
      auto cell_class = std::make_unique<TH3I>(
          ("hReachCellClass_" + suffix[side]).c_str(),
          "synthetic cell reach class",
          PhiNodes,
          0.0,
          TwoPi,
          RadialNodes - 1,
          InnerRadiusCm,
          OuterRadiusCm,
          LongitudinalNodes - 1,
          z_min[side],
          z_max[side]);
      cell_class->SetDirectory(nullptr);

      for (int iphi = 0; iphi < PhiNodes; ++iphi)
      {
        const double phi = static_cast<double>(iphi) * phi_spacing;
        for (int ir = 0; ir < RadialNodes; ++ir)
        {
          for (int iz = 0; iz < LongitudinalNodes; ++iz)
          {
            const bool readout_node =
                (side == 0 && iz == 0) ||
                (side == 1 && iz + 1 == LongitudinalNodes);
            const double radial_displacement =
                readout_node ? 0.0 : (ir == 0 ? 0.02 : -0.02);
            const double transverse_displacement =
                readout_node ? 0.0 : 0.03 * std::sin(phi);
            const double time_ns =
                readout_node
                    ? 0.0
                    : 100.0 + 2.0 * static_cast<double>(iphi) +
                          static_cast<double>(ir);
            const int xbin = iphi + 2;
            const int ybin = ir + 2;
            const int zbin = iz + 2;
            delta_r->SetBinContent(
                xbin, ybin, zbin, radial_displacement);
            delta_rphi->SetBinContent(
                xbin, ybin, zbin, transverse_displacement);
            drift_time->SetBinContent(xbin, ybin, zbin, time_ns);
            reach->SetBinContent(xbin, ybin, zbin, 1);
            maximum_drift_time_ns =
                std::max(maximum_drift_time_ns, time_ns);
          }
        }
      }

      fill_guards(*delta_r);
      fill_guards(*delta_rphi);
      fill_guards(*drift_time);
      fill_guards(*reach);
      for (int iphi = 1; iphi <= PhiNodes; ++iphi)
      {
        cell_class->SetBinContent(iphi, 1, 1, 2);
      }

      if (numerical)
      {
        auto transverse = std::unique_ptr<TH3D>(new TH3D(("hInterpolationTransverseError_"+suffix[side]).c_str(), "probed error", PhiNodes, 0., TwoPi, RadialNodes-1, InnerRadiusCm, OuterRadiusCm, LongitudinalNodes-1, z_min[side], z_max[side]));
        auto time_error = std::unique_ptr<TH3D>(new TH3D(("hInterpolationTimeEquivalentError_"+suffix[side]).c_str(), "probed error", PhiNodes, 0., TwoPi, RadialNodes-1, InnerRadiusCm, OuterRadiusCm, LongitudinalNodes-1, z_min[side], z_max[side]));
        cell_class->SetBinContent(1,1,1,3);
        transverse->SetBinContent(1,1,1,10.);
        transverse->Write(); time_error->Write();
      }
      if (delta_r->Write() <= 0 ||
          delta_rphi->Write() <= 0 ||
          drift_time->Write() <= 0 ||
          reach->Write() <= 0 ||
          cell_class->Write() <= 0)
      {
        return false;
      }
    }

    if (adaptive && (!extended || !PHGarfieldAdaptiveMap::Write(
            transport->mkdir("Adaptive"), *adaptive, std::cerr)))
    { return false; }

    TDirectory* metadata = transport->mkdir("Metadata");
    if (!metadata)
    {
      return false;
    }
    metadata->cd();

    PHGarfieldForwardMap::Config fingerprint_config;
    fingerprint_config.phi_nodes = PhiNodes;
    fingerprint_config.radial_nodes = RadialNodes;
    fingerprint_config.longitudinal_nodes = LongitudinalNodes;
    fingerprint_config.inner_radius_cm = InnerRadiusCm;
    fingerprint_config.outer_radius_cm = OuterRadiusCm;
    fingerprint_config.central_membrane_half_width_cm =
        CentralMembraneHalfWidthCm;
    fingerprint_config.readout_abs_z_cm = ReadoutAbsZCm;
    fingerprint_config.integration_step_ns = 12.5;
    fingerprint_config.nominal_phg4_drift_velocity_cm_per_ns =
        0.00755;
    fingerprint_config.maximum_nominal_time_mismatch = 0.02;
    fingerprint_config.cdb_global_tag = "SYNTHETIC";
    fingerprint_config.cdb_timestamp = 1;
    fingerprint_config.probe_cell_centers = true;
    fingerprint_config.configuration_label =
        "synthetic reader fixture";

    PHGarfieldForwardMap::Provenance fingerprint_provenance;
    fingerprint_provenance.nominal_electric_field_vcm = 400.0;
    fingerprint_provenance.space_charge_scale_side0 = 1.0;
    fingerprint_provenance.space_charge_scale_side1 = 1.0;
    fingerprint_provenance.magnetic_field_scale = 1.0;
    fingerprint_provenance.zero_magnetic_field = true;
    fingerprint_provenance.gas_file = "synthetic.gas";
    fingerprint_provenance.gas_file_md5 =
        "00000000000000000000000000000000";
    fingerprint_provenance.maximum_integration_steps = 1000;
    fingerprint_provenance.reach_probe_scheme =
        "corners+center+faces+cage-biased-v1";
    const double nominal_reference_time_ns =
        (ReadoutAbsZCm - CentralMembraneHalfWidthCm) /
        fingerprint_config.nominal_phg4_drift_velocity_cm_per_ns;
    fingerprint_provenance.scale_zero_reference_time_side0_ns =
        nominal_reference_time_ns;
    fingerprint_provenance.scale_zero_reference_time_side1_ns =
        nominal_reference_time_ns;
    fingerprint_provenance.scale_zero_maximum_time_mismatch = 0.0;
    if (extended)
    {
      fingerprint_provenance.field_model = "per_side_3d";
      fingerprint_provenance.electric_field_map_3d = {{"south.root", "north.root"}};
      fingerprint_provenance.electric_field_map_3d_md5 = {{std::string(32, 'a'), std::string(32, 'b')}};
      fingerprint_provenance.frame_field_map_3d = {{"frames_south.root", "frames_north.root"}};
      fingerprint_provenance.frame_field_map_3d_md5 = {{std::string(32, 'c'), std::string(32, 'd')}};
      fingerprint_provenance.frame_charge_scale = {{-180.0, -180.0}};
    }
    fingerprint_provenance.reach_probe_scheme = adaptive ? PHGarfieldAdaptiveMap::ProbeScheme : fingerprint_provenance.reach_probe_scheme;
    fingerprint_config.interpolation_tolerance_um = adaptive ? adaptive->policy.guard_um : numerical ? 3. : 0.;
    const int schema_version = adaptive ? PHGarfieldAdaptiveMap::SchemaVersion : numerical ? PHGarfieldForwardMap::NumericalSchemaVersion : extended ? PHGarfieldForwardMap::ExtendedSchemaVersion : PHGarfieldForwardMap::SchemaVersion;
    const std::string canonical_configuration =
        PHGarfieldForwardMap::BuildCanonicalConfiguration(
            fingerprint_config,
            fingerprint_provenance, schema_version);
    write_named("schema_name", PHGarfieldForwardMap::SchemaName);
    write_parameter(
        "schema_version", schema_version);
    write_named(
        "transport_contract", PHGarfieldForwardMap::TransportContract);
    write_named("coordinate_frame", "TPC_ENVELOPE_LOCAL");
    write_named("axis_order_units", "phi_rad,r_cm,z_cm");
    write_named(
        "transverse_quantity", "delta_r_cm,r0_delta_phi_cm");
    write_named("time_quantity", "absolute_drift_duration_ns");
    write_named("integrator", "explicit_euler");
    write_named(
        "magnetic_field_domain_policy",
        PHGarfieldForwardMap::MagneticFieldDomainPolicy);
    write_named(
        "nominal_time_reference_policy",
        PHGarfieldForwardMap::NominalTimeReferencePolicy);
    write_named("reach_encoding", (numerical || adaptive) ? "0=AllLost,1=Mixed,2=AllReach,3=NumericallyUnresolved" : "0=AllLost,1=Mixed,2=AllReach");
    if (numerical || adaptive)
    {
      write_parameter("interpolation_tolerance_um", fingerprint_config.interpolation_tolerance_um);
      write_named("interpolation_error_probe_scheme", adaptive ? PHGarfieldAdaptiveMap::ProbeScheme : "center_and_face_centers_transverse_time_equivalent_v1");
    }
    write_named(
        "reach_probe_scheme", fingerprint_provenance.reach_probe_scheme);
    write_named("configuration_label", "synthetic reader fixture");
    if (extended)
    {
      write_named("field_model", fingerprint_provenance.field_model);
      write_named("frame_field_map", "");
      write_named("frame_field_map_md5", "");
      for (std::size_t side = 0; side < 2; ++side)
      {
        const std::string side_tag = "_side" + std::to_string(side);
        write_named(("electric_field_map_3d" + side_tag).c_str(), fingerprint_provenance.electric_field_map_3d[side]);
        write_named(("electric_field_map_3d_md5" + side_tag).c_str(), fingerprint_provenance.electric_field_map_3d_md5[side]);
        write_named(("frame_field_map_3d" + side_tag).c_str(), fingerprint_provenance.frame_field_map_3d[side]);
        write_named(("frame_field_map_3d_md5" + side_tag).c_str(), fingerprint_provenance.frame_field_map_3d_md5[side]);
        write_parameter(("frame_charge_scale" + side_tag).c_str(), -180.0);
        write_parameter(("cm_field_factor" + side_tag).c_str(), 1.0);
      }
      write_parameter("use_ifc_voltage_distortion", 0);
      write_parameter("use_ofc_voltage_distortion", 0);
      for (std::size_t index = 0; index < 4; ++index)
      {
        write_parameter(("field_cage_offset_v_" + std::to_string(index)).c_str(), 0.0);
      }
    }
    write_named("electric_field_map", "");
    write_named("magnetic_field_map", "");
    write_named("gas_file", "synthetic.gas");
    write_named("electric_field_map_md5", "");
    write_named("magnetic_field_map_md5", "");
    write_named("gas_file_md5", "00000000000000000000000000000000");
    write_named("cdb_global_tag", "SYNTHETIC");
    write_named("canonical_configuration", canonical_configuration);
    write_named(
        "configuration_fingerprint_fnv1a64",
        hexadecimal(fnv1a64(canonical_configuration)));

    write_parameter("phi_nodes", PhiNodes);
    write_parameter("radial_nodes", RadialNodes);
    write_parameter("longitudinal_nodes", LongitudinalNodes);
    write_parameter("gas_inner_radius_cm", InnerRadiusCm);
    write_parameter("gas_outer_radius_cm", OuterRadiusCm);
    write_parameter(
        "central_membrane_half_width_cm",
        CentralMembraneHalfWidthCm);
    write_parameter("readout_abs_z_cm", ReadoutAbsZCm);
    write_parameter("integration_step_ns", 12.5);
    write_parameter(
        "nominal_phg4_drift_velocity_cm_per_ns", 0.00755);
    write_parameter("maximum_nominal_time_mismatch", 0.02);
    write_parameter(
        "scale_zero_reference_time_side0_ns",
        nominal_reference_time_ns);
    write_parameter(
        "scale_zero_reference_time_side1_ns",
        nominal_reference_time_ns);
    write_parameter("scale_zero_maximum_time_mismatch", 0.0);
    write_parameter("maximum_drift_time_ns", adaptive ? adaptive->maximum_drift_time_ns : maximum_drift_time_ns);
    if (adaptive)
    {
      write_parameter("generation_root_first", static_cast<Long64_t>(adaptive->first_root));
      write_parameter("generation_root_end", static_cast<Long64_t>(adaptive->end_root));
    }
    write_parameter("tpc_place_x_cm", 0.0);
    write_parameter("tpc_place_y_cm", 0.0);
    write_parameter("tpc_place_z_cm", 0.0);
    write_parameter("tpc_rotation_x_rad", 0.0);
    write_parameter("tpc_rotation_y_rad", 0.0);
    write_parameter("tpc_rotation_z_rad", 0.0);
    write_parameter("nominal_electric_field_vcm", 400.0);
    write_parameter("space_charge_scale_side0", 1.0);
    write_parameter("space_charge_scale_side1", 1.0);
    write_parameter("magnetic_field_scale", 1.0);
    write_parameter("zero_magnetic_field", 1);
    write_parameter("maximum_integration_steps", 1000);
    write_parameter("cell_center_probes", 1);
    write_named("reach_probe_trajectories", "0");
    write_named("cdb_timestamp", "1");

    output->Write();
    output->Flush();
    const bool write_ok =
        !output->TestBit(TFile::kWriteError) &&
        output->GetBytesWritten() > 0;
    output->Close();
    return write_ok;
  }

  bool valid_sample(
      const PHG4TpcDistortion::GarfieldTransportSample& sample)
  {
    if (sample.reach ==
        PHG4TpcDistortion::ReachCellClass::AllReach)
    {
      return std::isfinite(sample.delta_r_cm) &&
             std::isfinite(sample.delta_rphi_cm) &&
             std::isfinite(sample.drift_time_ns) &&
             sample.drift_time_ns >= 0.0;
    }
    return (sample.reach ==
                PHG4TpcDistortion::ReachCellClass::AllLost ||
            sample.reach ==
                PHG4TpcDistortion::ReachCellClass::Mixed) &&
           sample.delta_r_cm == 0.0 &&
           sample.delta_rphi_cm == 0.0 &&
           sample.drift_time_ns == 0.0;
  }

  bool close_sample(
      const PHG4TpcDistortion::GarfieldTransportSample& first,
      const PHG4TpcDistortion::GarfieldTransportSample& second)
  {
    constexpr double tolerance = 1.0e-11;
    return first.reach == second.reach &&
           std::abs(first.delta_r_cm - second.delta_r_cm) < tolerance &&
           std::abs(first.delta_rphi_cm - second.delta_rphi_cm) <
               tolerance &&
           std::abs(first.drift_time_ns - second.drift_time_ns) <
               tolerance;
  }

  PHGarfieldForwardMap::Config adaptive_geometry_configuration()
  {
    PHGarfieldForwardMap::Config config;
    config.phi_nodes = PhiNodes;
    config.radial_nodes = RadialNodes;
    config.longitudinal_nodes = LongitudinalNodes;
    config.inner_radius_cm = InnerRadiusCm;
    config.outer_radius_cm = OuterRadiusCm;
    config.central_membrane_half_width_cm = CentralMembraneHalfWidthCm;
    config.readout_abs_z_cm = ReadoutAbsZCm;
    return config;
  }

  // An explicit analytic numerical fixture, not a gas/detector response model.
  // Its quadratic radial term requires a radial split at the10um guard; its
  // drift duration is linear in signed z within each side and zero at readout.
  bool build_adaptive_fixture(PHGarfieldAdaptiveMap::Data& data,
                              const unsigned int query_budget = 2048,
                              const bool partial = false)
  {
    const auto geometry = PHGarfieldAdaptiveMap::FromConfig(adaptive_geometry_configuration());
    PHGarfieldAdaptiveMap::Policy policy;
    policy.guard_um = 10.;
    policy.max_queries_per_root = query_budget;
    const auto sampler = [](const std::array<double, 3>& coordinate,
                            PHGarfieldAdaptiveMap::Point& point, std::ostream&)
    {
      point.coordinate = coordinate;
      point.reaches = 1;
      point.forward_status = static_cast<int>(PHGarfield::ForwardDriftStatus::ReachedReadout);
      const double drift = ReadoutAbsZCm - std::abs(coordinate[2]);
      const double fraction = drift / (ReadoutAbsZCm - CentralMembraneHalfWidthCm);
      point.value = {{1.e-4 * (coordinate[1] - InnerRadiusCm) *
                         (OuterRadiusCm - coordinate[1]) * fraction,
                      1.e-5 * std::sin(coordinate[0]) * fraction, 10. * drift}};
      return true;
    };
    PHGarfieldAdaptiveMap::BuildSummary summary;
    std::ostringstream log;
    const auto end = PHGarfieldAdaptiveMap::RootCount(geometry) - (partial ? 1U : 0U);
    return PHGarfieldAdaptiveMap::Build(geometry, policy, 0, end, sampler, data, summary, log, std::cerr);
  }

  bool rejects_initialization(const std::string& filename)
  {
    const auto child = fork();
    if (child < 0) { return false; }
    if (child == 0)
    {
      PHG4TpcDistortion map;
      map.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
      map.set_do_static_distortions(true);
      map.set_static_distortion_filename(filename);
      map.Init();
      _exit(0);
    }
    int status = 0;
    return waitpid(child, &status, 0) == child && WIFEXITED(status) && WEXITSTATUS(status) != 0;
  }

  bool adaptive_reader_fixtures()
  {
    using Map = PHGarfieldAdaptiveMap;
    using Status = PHG4TpcDistortion::GarfieldTransportQueryStatus;
    Map::Data data;
    TemporaryFile file;
    if (!build_adaptive_fixture(data) || !reserve_temporary_filename(file) ||
        !write_synthetic_map(file.filename, true, false, &data))
    { std::cerr << "Could not build analytic adaptive ROOT fixture\n"; return false; }
    unsigned int splits = 0;
    for (const auto& cell : data.cells)
    {
      if (cell.split_axis >= 0)
      {
        if (cell.split_axis != 1)
        { std::cerr << "Analytic radial curvature refined another axis\n"; return false; }
        ++splits;
      }
    }
    if (!splits) { std::cerr << "Adaptive fixture never refined radial curvature\n"; return false; }
    PHG4TpcDistortion map;
    map.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
    map.set_do_static_distortions(true);
    map.set_static_distortion_filename(file.filename);
    map.Init();
    const auto& metadata = map.garfield_lut_metadata();
    if (metadata.schema_version != 4 || metadata.interpolation_tolerance_um != 10. ||
        metadata.adaptive_cell_count != data.cells.size() ||
        metadata.adaptive_root_count != data.roots.size() ||
        metadata.adaptive_all_lost_policy != Map::AllLostPolicy ||
        gROOT->GetListOfFiles()->FindObject(file.filename.c_str()))
    { std::cerr << "Adaptive metadata or owning-reader lifetime mismatch\n"; return false; }
    // The reader must remain usable without an open ROOT file or live TTree.
    if (std::remove(file.filename.c_str()) != 0) { return false; }
    unsigned long long leaf_queries = 0, readout_queries = 0;
    for (const double z : {-99., -50.5, -1., 1., 50.5, 99.})
    {
      for (const double radius : {InnerRadiusCm, 22.5, 25., 27.5, OuterRadiusCm})
      {
        for (const double phi : {.2, 2.4, 5.1})
        {
          PHG4TpcDistortion::GarfieldTransportSample sample;
          if (map.query_garfield_transport(radius, phi, z, sample) != Status::Success ||
              !sample.used_adaptive_lookup || sample.used_direct_fallback ||
              sample.adaptive_leaf_class != static_cast<int>(Map::LeafClass::Certified))
          { std::cerr << "Certified adaptive query failed\n"; return false; }
          ++leaf_queries;
          const auto found = std::find_if(data.cells.begin(), data.cells.end(),
              [&sample](const auto& cell) { return cell.id == sample.adaptive_cell_id; });
          if (found == data.cells.end()) { return false; }
          // Independently reconstruct the tensor weights from persisted corner
          // values. Do not call the implementation's Interpolate for expected values.
          std::array<double, 3> expected{{0., 0., 0.}};
          const std::array<double, 3> coordinate{{phi, radius, z}};
          for (unsigned int corner = 0; corner < 8; ++corner)
          {
            double weight = 1.;
            for (unsigned int axis = 0; axis < 3; ++axis)
            {
              const double fraction = (coordinate[axis] - found->lo[axis]) / (found->hi[axis] - found->lo[axis]);
              weight *= (corner & (1U << (2U - axis))) ? fraction : 1. - fraction;
            }
            const auto point = std::find_if(data.points.begin(), data.points.end(),
                [id = found->corner_id[corner]](const auto& value) { return value.id == static_cast<std::uint64_t>(id); });
            if (point == data.points.end()) { return false; }
            for (unsigned int component = 0; component < 3; ++component)
            { expected[component] += weight * point->value[component]; }
          }
          if (std::abs(sample.delta_r_cm - expected[0]) > 1.e-12 ||
              std::abs(sample.delta_rphi_cm - expected[1]) > 1.e-12 ||
              std::abs(sample.drift_time_ns - expected[2]) > 1.e-10)
          { std::cerr << "Adaptive tensor interpolation disagrees with independent weights\n"; return false; }
          const double truth = 1.e-4 * (radius - InnerRadiusCm) * (OuterRadiusCm - radius) *
              (ReadoutAbsZCm - std::abs(z)) / (ReadoutAbsZCm - CentralMembraneHalfWidthCm);
          if (1.e4 * std::abs(sample.delta_r_cm - truth) > 10. + 1.e-8)
          { std::cerr << "Analytic adaptive interpolation exceeds fixture guard\n"; return false; }
          // Midpoint equality must select the upper radial child, even though
          // the shared face has continuous values from either child.
          if (radius == 25.)
          {
            const auto root = std::find_if(data.roots.begin(), data.roots.end(),
                [&found](const auto& value) { return value.root_id == found->root_id; });
            if (root == data.roots.end()) { return false; }
            const auto parent = std::find_if(data.cells.begin(), data.cells.end(),
                [&root](const auto& value) { return value.id == root->cell_id; });
            if (parent == data.cells.end() || parent->child_id[1] != static_cast<std::int64_t>(sample.adaptive_cell_id))
            { std::cerr << "Adaptive split equality did not choose upper child\n"; return false; }
          }
        }
      }
    }
    for (const double z : {-1., 1.})
    {
      for (const auto& pair : {std::array<double, 2>{{-1.e-8, TwoPi - 1.e-8}},
                              std::array<double, 2>{{1.e-8, TwoPi + 1.e-8}},
                              std::array<double, 2>{{0., TwoPi}}})
      {
        PHG4TpcDistortion::GarfieldTransportSample first, second;
        if (map.query_garfield_transport(25., pair[0], z, first) != Status::Success ||
            map.query_garfield_transport(25., pair[1], z, second) != Status::Success ||
            !close_sample(first, second) || first.adaptive_cell_id != second.adaptive_cell_id)
        { std::cerr << "Adaptive periodic seam mismatch\n"; return false; }
        leaf_queries += 2;
      }
    }
    for (const double z : {-ReadoutAbsZCm, ReadoutAbsZCm})
    {
      for (const double radius : {InnerRadiusCm, 25., OuterRadiusCm})
      {
        PHG4TpcDistortion::GarfieldTransportSample sample;
        if (map.query_garfield_transport(radius, TwoPi, z, sample) != Status::Success ||
            sample.reach != PHG4TpcDistortion::ReachCellClass::AllReach ||
            sample.used_adaptive_lookup || sample.used_direct_fallback ||
            sample.delta_r_cm != 0. || sample.delta_rphi_cm != 0. || sample.drift_time_ns != 0.)
        { std::cerr << "Adaptive readout endpoint is not exact zero\n"; return false; }
        ++readout_queries;
      }
    }
    PHG4TpcDistortion::GarfieldTransportSample invalid;
    if (map.query_garfield_transport(25., 0., 0., invalid) != Status::OutsideDriftVolume ||
        map.query_garfield_transport(25., std::numeric_limits<double>::quiet_NaN(), 1., invalid) != Status::Invalid)
    { return false; }
    const auto& stats = map.boundary_fallback_stats();
    if (stats.adaptive_leaf_queries != leaf_queries || stats.adaptive_interpolated_queries != leaf_queries ||
        stats.adaptive_readout_queries != readout_queries || stats.adaptive_queries != leaf_queries + readout_queries ||
        stats.total_queries != leaf_queries + readout_queries + 2 || stats.outside_queries != 1 || stats.invalid_queries != 1 ||
        stats.lut_reaching_queries != leaf_queries + readout_queries || stats.calls != 0 || stats.adaptive_max_depth_seen != 1)
    { std::cerr << "Adaptive query/source counter partition mismatch\n"; return false; }
    // A legal budget leaf is never treated as interpolated or sampled loss
    // without the explicitly enabled pinned direct model.
    Map::Data budget;
    TemporaryFile budget_file;
    if (!build_adaptive_fixture(budget, 27) || !reserve_temporary_filename(budget_file) ||
        !write_synthetic_map(budget_file.filename, true, false, &budget)) { return false; }
    PHG4TpcDistortion budget_map;
    budget_map.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
    budget_map.set_do_static_distortions(true); budget_map.set_static_distortion_filename(budget_file.filename);
    budget_map.Init();
    PHG4TpcDistortion::GarfieldTransportSample unresolved;
    if (budget_map.query_garfield_transport(25., .2, 1., unresolved) != Status::Invalid ||
        !unresolved.used_adaptive_lookup || unresolved.used_direct_fallback ||
        unresolved.reach != PHG4TpcDistortion::ReachCellClass::Invalid ||
        unresolved.adaptive_leaf_class == static_cast<int>(Map::LeafClass::Certified) ||
        budget_map.boundary_fallback_stats().adaptive_interpolated_queries != 0 ||
        budget_map.boundary_fallback_stats().invalid_queries != 1)
    { std::cerr << "Uncertified adaptive leaf did not fail closed\n"; return false; }
    TemporaryFile missing;
    if (!reserve_temporary_filename(missing) || !write_synthetic_map(missing.filename, true, false, &data)) { return false; }
    {
      TFile changed(missing.filename.c_str(), "UPDATE");
      auto* transport = changed.GetDirectory("GarfieldFullTransport");
      transport->Delete("Adaptive;*"); changed.Close();
    }
    if (!rejects_initialization(missing.filename))
    { std::cerr << "Schema4 accepted missing adaptive payload despite regular histograms\n"; return false; }
    // The duplicate top range must be a signed ROOT scalar and must describe
    // the same complete partition as the independently validated adaptive data.
    // Test both endpoints, signedness/type, and absence rather than accepting a
    // correct payload with misleading or silently defaulted top metadata.
    for (unsigned int mutation = 0; mutation < 5; ++mutation)
    {
      TemporaryFile changed_range;
      if (!reserve_temporary_filename(changed_range) ||
          !write_synthetic_map(changed_range.filename, true, false, &data)) { return false; }
      {
        TFile changed(changed_range.filename.c_str(), "UPDATE");
        auto* changed_metadata = changed.GetDirectory("GarfieldFullTransport/Metadata");
        if (changed.IsZombie() || !changed_metadata) { return false; }
        changed_metadata->cd();
        switch (mutation)
        {
        case 0:
          TParameter<Long64_t>("generation_root_first", 1).Write("", TObject::kOverwrite); break;
        case 1:
          TParameter<Long64_t>("generation_root_end", static_cast<Long64_t>(data.end_root - 1)).Write("", TObject::kOverwrite); break;
        case 2:
          TParameter<Long64_t>("generation_root_first", -1).Write("", TObject::kOverwrite); break;
        case 3:
          TParameter<int>("generation_root_first", 0).Write("", TObject::kOverwrite); break;
        case 4:
          changed_metadata->Delete("generation_root_end;*"); break;
        }
        changed.Close();
      }
      if (!rejects_initialization(changed_range.filename))
      { std::cerr << "Schema4 accepted invalid top root range mutation " << mutation << '\n'; return false; }
    }
    Map::Data partial;
    TemporaryFile partial_file;
    if (!build_adaptive_fixture(partial, 2048, true) || !reserve_temporary_filename(partial_file) ||
        !write_synthetic_map(partial_file.filename, true, false, &partial) || !rejects_initialization(partial_file.filename))
    { std::cerr << "Schema4 did not reject incomplete root coverage\n"; return false; }
    return true;
  }

  // Exercise the actual detector reader with the same valid physical envelope
  // used by schema4, while the pure compact tests cover all source classes.
  bool compact_reader_fixtures()
  {
    using Map = PHGarfieldCompactMap;
    using Stats = PHG4TpcDistortion::BoundaryFallbackStats;
    using Caps = std::array<std::uint64_t, 5>;
    const Caps large_caps{{1500000000ULL, 250000000ULL, 150000000ULL,
                           1500000000ULL, 128ULL*1024*1024*1024}};
    const auto set_caps = [](PHG4TpcDistortion& map, const Caps& caps)
    { map.set_garfield_compact_read_limits(caps[0], caps[1], caps[2], caps[3], caps[4]); };
    const auto observed_caps = [](const PHG4TpcDistortion::GarfieldCompactReadLimits& caps)
    { return Caps{{caps.max_source_points, caps.max_nodes, caps.max_corners, caps.max_values, caps.max_bytes}}; };
    const Map::Limits defaults;
    const Caps default_caps{{defaults.max_source_points, defaults.max_nodes, defaults.max_corners,
                             defaults.max_values, defaults.max_bytes}};
    PHG4TpcDistortion admission;
    if (observed_caps(admission.garfield_compact_read_limits()) != default_caps ||
        admission.garfield_compact_read_limits_explicit()) { return false; }
    // Setter validation is scalar-only. Bad input leaves the accepted request
    // unchanged and must not attempt an allocation at these upper bounds.
    set_caps(admission, large_caps);
    for (unsigned int axis = 0; axis < large_caps.size(); ++axis)
    for (const bool excessive : {false, true})
    {
      auto bad = large_caps; bad[axis] = excessive ? bad[axis]+1 : 0;
      bool rejected = false;
      try { set_caps(admission, bad); }
      catch (const std::invalid_argument&) { rejected = true; }
      if (!rejected || observed_caps(admission.garfield_compact_read_limits()) != large_caps)
      { std::cerr << "Invalid compact admission limits accepted or changed state\n"; return false; }
    }
    static_assert(std::has_unique_object_representations_v<Stats>);
    for (const unsigned int budget : {2048U, 27U})
    {
      PHGarfieldAdaptiveMap::Data training;
      TemporaryFile source_file, runtime_file;
      if (!build_adaptive_fixture(training, budget) || !reserve_temporary_filename(source_file) ||
          !reserve_temporary_filename(runtime_file) ||
          !write_synthetic_map(source_file.filename, true, false, &training)) { return false; }
      Map::SourceBinding binding;
      binding.file_sha256 = Map::FileSha256(source_file.filename);
      binding.content_fingerprint_fnv1a64 = PHGarfieldAdaptiveMap::ContentFingerprint(training);
      binding.canonical_policy = PHGarfieldAdaptiveMap::BuildCanonicalPolicy(training);
      binding.converter_sha256 = std::string(64, 'a');  // Explicit synthetic converter identity.
      Caps required_caps{};
      {
        TFile input(source_file.filename.c_str(), "READ");
        auto* metadata = input.GetDirectory("GarfieldFullTransport/Metadata");
        auto* fingerprint = metadata ? dynamic_cast<TNamed*>(metadata->Get("configuration_fingerprint_fnv1a64")) : nullptr;
        if (!fingerprint) { return false; }
        binding.configuration_fingerprint_fnv1a64 = fingerprint->GetTitle();
        Map::Data packed;
        if (!Map::Convert(training, binding, packed, Map::Limits{}, std::cerr, true)) { return false; }
        required_caps = {{training.points.size(), packed.nodes.size(), packed.corners.size(), packed.values.size(), 1025}};
        TFile output(runtime_file.filename.c_str(), "CREATE");
        if (output.IsZombie()) { return false; }
        auto* top = output.mkdir("GarfieldFullTransport");
        auto* identity = top->mkdir("Metadata");
        identity->cd();
        TParameter<int>("schema_version", Map::SchemaVersion).Write();
        TNamed("schema_name", Map::SchemaName).Write();
        TParameter<int>("source_schema_version", 4).Write();
        TNamed("compact_layout", Map::Layout).Write();
        TNamed("source_file_sha256", binding.file_sha256.c_str()).Write();
        TNamed("source_configuration_fingerprint_fnv1a64", binding.configuration_fingerprint_fnv1a64.c_str()).Write();
        TNamed("source_content_fingerprint_fnv1a64", binding.content_fingerprint_fnv1a64.c_str()).Write();
        auto* copied = top->mkdir("SourceMetadata");
        copied->cd();
        TIter next(metadata->GetListOfKeys());
        while (auto* key = static_cast<TKey*>(next()))
        {
          std::unique_ptr<TObject> object(key->ReadObj());
          if (!object || object->Write(key->GetName()) <= 0) { return false; }
        }
        if (!Map::Write(top->mkdir("Runtime"), packed, Map::Limits{}, std::cerr)) { return false; }
        output.Write();
        if (output.TestBit(TFile::kWriteError)) { return false; }
        output.Close();
      }
      const std::string runtime_hash = Map::FileSha256(runtime_file.filename);
      if (budget == 2048)
      {
        // Each of the five caps must reach the actual native reader. The sixth
        // mutation rejects caps on schema4; the seventh rejects legacy mode.
        for (unsigned int axis = 0; axis < 7; ++axis)
        {
          if (axis < 5 && required_caps[axis] <= 1) { return false; }
          const auto child = fork();
          if (child < 0) { return false; }
          if (child == 0)
          {
            PHG4TpcDistortion candidate;
            if (axis != 6)
              candidate.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
            candidate.set_do_static_distortions(true);
            candidate.set_static_distortion_filename(axis == 5 ? source_file.filename : runtime_file.filename);
            candidate.set_garfield_compact_source_binding(binding.file_sha256,
                binding.configuration_fingerprint_fnv1a64, binding.content_fingerprint_fnv1a64, runtime_hash);
            auto caps = large_caps;
            if (axis < 5) { caps[axis] = required_caps[axis]-1; }
            set_caps(candidate, caps);
            candidate.Init();
            _exit(0);
          }
          int status = 0;
          if (waitpid(child, &status, 0) != child || !WIFEXITED(status) || WEXITSTATUS(status) == 0)
          { std::cerr << "Compact admission cap did not reject input: " << axis << '\n'; return false; }
        }
      }
      // Missing or wrong external attestations must fail even for an otherwise
      // valid file. Do not accept self-asserted training metadata as approval.
      for (int mutation = -1; mutation < 4; ++mutation)
      {
        const auto child = fork();
        if (child < 0) { return false; }
        if (child == 0)
        {
          PHG4TpcDistortion candidate;
          candidate.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
          candidate.set_do_static_distortions(true);
          candidate.set_static_distortion_filename(runtime_file.filename);
          if (mutation >= 0)
          {
            candidate.set_garfield_compact_source_binding(
                mutation == 0 ? std::string(64, '0') : binding.file_sha256,
                mutation == 1 ? std::string(16, '0') : binding.configuration_fingerprint_fnv1a64,
                mutation == 2 ? std::string(16, '0') : binding.content_fingerprint_fnv1a64,
                mutation == 3 ? std::string(64, '0') : runtime_hash);
          }
          candidate.Init();
          _exit(0);
        }
        int status = 0;
        if (waitpid(child, &status, 0) != child || !WIFEXITED(status) || WEXITSTATUS(status) == 0)
        { std::cerr << "Compact reader accepted missing/wrong external binding " << mutation << '\n'; return false; }
      }
      PHG4TpcDistortion full, compact;
      for (auto* map : {&full, &compact})
      {
        map->set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
        map->set_do_static_distortions(true);
      }
      full.set_static_distortion_filename(source_file.filename);
      compact.set_static_distortion_filename(runtime_file.filename);
      compact.set_garfield_compact_source_binding(binding.file_sha256,
          binding.configuration_fingerprint_fnv1a64, binding.content_fingerprint_fnv1a64, runtime_hash);
      if (budget == 2048) { set_caps(compact, large_caps); }
      full.Init(); compact.Init();
      const auto& metadata = compact.garfield_lut_metadata();
      if (metadata.schema_version != 5 || metadata.compact_source_sha256 != binding.file_sha256 ||
          metadata.compact_file_sha256 != runtime_hash || metadata.adaptive_point_count != training.points.size() ||
          metadata.adaptive_cell_count != training.cells.size() || metadata.adaptive_max_queries_per_root != budget ||
          metadata.compact_read_limits_explicit != (budget == 2048) ||
          observed_caps(metadata.compact_read_limits) != (budget == 2048 ? large_caps : default_caps) ||
          gROOT->GetListOfFiles()->FindObject(runtime_file.filename.c_str())) { return false; }
      bool rejected_after_init = false;
      try { set_caps(compact, large_caps); }
      catch (const std::logic_error&) { rejected_after_init = true; }
      if (!rejected_after_init) { return false; }
      if (std::remove(source_file.filename.c_str()) || std::remove(runtime_file.filename.c_str())) { return false; }
      for (double z : {-101., -100., -50.5, -1., 0., 1., 50.5, 100., 101.})
      for (double radius : {19., 20., 22.5, 25., 27.5, 30., 31.})
      for (double phi : {-TwoPi, -1.e-8, 0., 1.e-8, .2, 2.4, 5.1, TwoPi, TwoPi + 1.e-8})
      {
        PHG4TpcDistortion::GarfieldTransportSample a, b;
        const auto sa = full.query_garfield_transport(radius, phi, z, a);
        const auto sb = compact.query_garfield_transport(radius, phi, z, b);
        if (sa != sb || a.reach != b.reach || a.delta_r_cm != b.delta_r_cm ||
            a.delta_rphi_cm != b.delta_rphi_cm || a.drift_time_ns != b.drift_time_ns ||
            a.used_adaptive_lookup != b.used_adaptive_lookup || a.used_direct_fallback != b.used_direct_fallback ||
            a.direct_fallback_reason != b.direct_fallback_reason || a.adaptive_cell_id != b.adaptive_cell_id ||
            a.adaptive_depth != b.adaptive_depth || a.adaptive_leaf_class != b.adaptive_leaf_class ||
            a.adaptive_fallback_reason != b.adaptive_fallback_reason)
        { std::cerr << "Compact detector query differs from full source\n"; return false; }
      }
      if (std::memcmp(&full.boundary_fallback_stats(), &compact.boundary_fallback_stats(), sizeof(Stats)))
      { std::cerr << "Compact detector counters differ\n"; return false; }
    }
    return true;
  }

}  // namespace

int main()
{
  using QueryStatus =
      PHG4TpcDistortion::GarfieldTransportQueryStatus;

  if (!adaptive_reader_fixtures()) { return 1; }
  if (!compact_reader_fixtures()) { return 1; }

  PHG4TpcDistortion uninitialized_map;
  PHG4TpcDistortion::GarfieldTransportSample uninitialized_sample;
  if (uninitialized_map.query_garfield_transport(
          25.0, 0.0, 50.0, uninitialized_sample) !=
      QueryStatus::Invalid)
  {
    std::cerr << "Uninitialized LUT query was not invalid\n";
    return 1;
  }

  // A second fixture exercises the extended schema independently of the
  // historical schema fixture and checks that metadata edits are rejected.
  TemporaryFile extended_file;
  if (!reserve_temporary_filename(extended_file) || !write_synthetic_map(extended_file.filename, true))
  {
    std::cerr << "Failed to create extended fixture\n";
    return 1;
  }
  {
    PHG4TpcDistortion extended_map;
    extended_map.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
    extended_map.set_do_static_distortions(true);
    extended_map.set_static_distortion_filename(extended_file.filename);
    extended_map.Init();
    const auto& metadata = extended_map.garfield_lut_metadata();
    if (metadata.schema_version != 2 || metadata.electric_field_map_3d[0] != "south.root" ||
        metadata.frame_field_map_3d_md5[1] != std::string(32, 'd') || metadata.frame_charge_scale[1] != -180.0)
    {
      std::cerr << "Extended provenance did not round-trip\n";
      return 1;
    }
  }
  {
    TFile changed(extended_file.filename.c_str(), "UPDATE");
    changed.cd("GarfieldFullTransport/Metadata");
    TParameter<double>("frame_charge_scale_side0", -179.0).Write("", TObject::kOverwrite);
    changed.Close();
  }
  const pid_t child = fork();
  if (child < 0) { return 1; }
  if (child == 0)
  {
    PHG4TpcDistortion changed_map;
    changed_map.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
    changed_map.set_do_static_distortions(true);
    changed_map.set_static_distortion_filename(extended_file.filename);
    changed_map.Init();
    _exit(0);
  }
  int child_status = 0;
  if (waitpid(child, &child_status, 0) != child || !WIFEXITED(child_status) || WEXITSTATUS(child_status) == 0)
  {
    std::cerr << "Tampered extended field provenance was accepted\n";
    return 1;
  }

  TemporaryFile numerical_file;
  if (!reserve_temporary_filename(numerical_file) || !write_synthetic_map(numerical_file.filename, true, true)) { return 1; }
  {
    PHG4TpcDistortion numerical_map;
    numerical_map.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
    numerical_map.set_do_static_distortions(true);
    numerical_map.set_static_distortion_filename(numerical_file.filename);
    numerical_map.Init();
    PHG4TpcDistortion::GarfieldTransportSample sample;
    const double radius = .5*(InnerRadiusCm+OuterRadiusCm), height=.5*(CentralMembraneHalfWidthCm+ReadoutAbsZCm);
    if (numerical_map.garfield_lut_metadata().schema_version != 3 ||
        numerical_map.query_garfield_transport(radius, .5*TwoPi/PhiNodes, height, sample) != QueryStatus::Invalid ||
        numerical_map.query_garfield_transport(radius, 1.5*TwoPi/PhiNodes, height, sample) != QueryStatus::Success ||
        numerical_map.boundary_fallback_stats().total_queries != 2 ||
        numerical_map.boundary_fallback_stats().invalid_queries != 1 ||
        numerical_map.boundary_fallback_stats().lut_reaching_queries != 1)
    { std::cerr << "Numerical mask did not fail closed or query counters disagree\n"; return 1; }
  }
  {
    TFile changed(numerical_file.filename.c_str(), "UPDATE");
    auto* errors = dynamic_cast<TH3*>(changed.Get("GarfieldFullTransport/hInterpolationTransverseError_posz"));
    errors->SetBinContent(1,1,1,0.);
    changed.cd("GarfieldFullTransport");errors->Write("",TObject::kOverwrite);changed.Close();
  }
  const pid_t numerical_child=fork();
  if(numerical_child<0) { return 1; }
  if(numerical_child==0)
  {
    PHG4TpcDistortion changed_map;changed_map.set_application_mode(PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
    changed_map.set_do_static_distortions(true);changed_map.set_static_distortion_filename(numerical_file.filename);changed_map.Init();_exit(0);
  }
  int numerical_status=0;
  if(waitpid(numerical_child,&numerical_status,0)!=numerical_child || !WIFEXITED(numerical_status) || WEXITSTATUS(numerical_status)==0)
  { std::cerr << "Tampered numerical cell mask was accepted\n";return 1; }

  TemporaryFile temporary;
  const char* environment_filename =
      std::getenv("TPC_GARFIELD_LUT_TEST_FILE");
  const bool using_synthetic_map =
      !environment_filename || environment_filename[0] == '\0';
  std::string filename;
  if (using_synthetic_map)
  {
    if (!reserve_temporary_filename(temporary) ||
        !write_synthetic_map(temporary.filename))
    {
      std::cerr << "Failed to create synthetic Garfield LUT fixture\n";
      return 1;
    }
    filename = temporary.filename;
  }
  else
  {
    filename = environment_filename;
  }

  PHG4TpcDistortion map;
  map.set_application_mode(
      PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport);
  map.set_do_static_distortions(true);
  map.set_static_distortion_filename(filename);
  map.Init();
  if (!map.is_initialized())
  {
    std::cerr << "LUT reader did not initialize\n";
    return 1;
  }

  const auto& metadata = map.garfield_lut_metadata();
  // Init validates every stored node, support cell and numerical mask.
  // Independent query/precision certification is performed by the study
  // harness; a masked file must not silently query without explicit fallback.
  if (std::getenv("TPC_GARFIELD_LUT_TEST_INITIALIZATION_ONLY"))
  {
    std::cout << "PASS: complete LUT structure/provenance validation; queries require independent study validation\n";
    return 0;
  }
  bool found_all_reach = false;
  bool checked_readout_zero = false;
  for (unsigned int side = 0; side < 2; ++side)
  {
    const double z_min =
        side == 0
            ? -metadata.readout_abs_z_cm
            : metadata.central_membrane_half_width_cm;
    const double z_max =
        side == 0
            ? -metadata.central_membrane_half_width_cm
            : metadata.readout_abs_z_cm;
    for (int iphi = 0; iphi < metadata.phi_nodes; ++iphi)
    {
      const double phi =
          (static_cast<double>(iphi) + 0.5) *
          TwoPi /
          static_cast<double>(metadata.phi_nodes);
      for (int ir = 0; ir + 1 < metadata.radial_nodes; ++ir)
      {
        const double r =
            metadata.gas_inner_radius_cm +
            (static_cast<double>(ir) + 0.5) *
                (metadata.gas_outer_radius_cm -
                 metadata.gas_inner_radius_cm) /
                static_cast<double>(metadata.radial_nodes - 1);
        for (int iz = 0; iz + 1 < metadata.longitudinal_nodes; ++iz)
        {
          const double z =
              z_min +
              (static_cast<double>(iz) + 0.5) *
                  (z_max - z_min) /
                  static_cast<double>(
                      metadata.longitudinal_nodes - 1);
          PHG4TpcDistortion::GarfieldTransportSample sample;
          if (!map.get_garfield_transport(r, phi, z, sample) ||
              !valid_sample(sample))
          {
            std::cerr << "Valid in-volume LUT query failed\n";
            return 1;
          }
          if (sample.reach ==
              PHG4TpcDistortion::ReachCellClass::AllReach)
          {
            found_all_reach = true;
          }
          if (using_synthetic_map &&
              sample.reach !=
                  PHG4TpcDistortion::ReachCellClass::AllReach)
          {
            std::cerr << "Synthetic fixture cell is not AllReach\n";
            return 1;
          }
        }
      }
    }
  }

  // The physical r/z boundaries are queryable on both drift sides.
  for (const double z :
       {-metadata.central_membrane_half_width_cm,
        metadata.central_membrane_half_width_cm})
  {
    for (const double r :
         {metadata.gas_inner_radius_cm,
          metadata.gas_outer_radius_cm})
    {
      PHG4TpcDistortion::GarfieldTransportSample sample;
      if (!map.get_garfield_transport(r, 0.37, z, sample) ||
          !valid_sample(sample))
      {
        std::cerr << "Physical r/z boundary query failed\n";
        return 1;
      }
      if (using_synthetic_map &&
          sample.reach !=
              PHG4TpcDistortion::ReachCellClass::AllReach)
      {
        std::cerr << "Synthetic boundary cell is not AllReach\n";
        return 1;
      }
    }
  }

  // Queries immediately on either representation of the phi seam agree.
  constexpr double seam_epsilon = 1.0e-8;
  const double test_r =
      0.5 *
      (metadata.gas_inner_radius_cm + metadata.gas_outer_radius_cm);
  for (const double z :
       {-metadata.central_membrane_half_width_cm,
        metadata.central_membrane_half_width_cm})
  {
    PHG4TpcDistortion::GarfieldTransportSample below_zero;
    PHG4TpcDistortion::GarfieldTransportSample below_two_pi;
    PHG4TpcDistortion::GarfieldTransportSample above_zero;
    PHG4TpcDistortion::GarfieldTransportSample above_two_pi;
    if (!map.get_garfield_transport(
            test_r, -seam_epsilon, z, below_zero) ||
        !map.get_garfield_transport(
            test_r, TwoPi - seam_epsilon, z, below_two_pi) ||
        !map.get_garfield_transport(
            test_r, seam_epsilon, z, above_zero) ||
        !map.get_garfield_transport(
            test_r, TwoPi + seam_epsilon, z, above_two_pi) ||
        !close_sample(below_zero, below_two_pi) ||
        !close_sample(above_zero, above_two_pi))
    {
      std::cerr << "Periodic phi-seam query failed\n";
      return 1;
    }
  }

  // Every in-range readout-plane query is exactly zero transport, including
  // radial boundaries and phi wrapping.
  for (const double z :
       {-metadata.readout_abs_z_cm, metadata.readout_abs_z_cm})
  {
    for (const double r :
         {metadata.gas_inner_radius_cm,
          test_r,
          metadata.gas_outer_radius_cm})
    {
      for (const double phi : {-seam_epsilon, 0.0, TwoPi})
      {
        PHG4TpcDistortion::GarfieldTransportSample sample;
        if (!map.get_garfield_transport(r, phi, z, sample) ||
            sample.reach !=
                PHG4TpcDistortion::ReachCellClass::AllReach ||
            sample.delta_r_cm != 0.0 ||
            sample.delta_rphi_cm != 0.0 ||
            sample.drift_time_ns != 0.0)
        {
          std::cerr << "Readout endpoint is not exact zero transport\n";
          return 1;
        }
        checked_readout_zero = true;
      }
    }
  }

  PHG4TpcDistortion::GarfieldTransportSample outside;
  const std::array<std::array<double, 3>, 7> outside_queries = {{
      {{metadata.gas_inner_radius_cm - 1.0, 0.0,
        metadata.central_membrane_half_width_cm}},
      {{metadata.gas_outer_radius_cm + 1.0, 0.0,
        metadata.central_membrane_half_width_cm}},
      {{test_r, 0.0, 0.0}},
      {{test_r, 0.0,
        metadata.readout_abs_z_cm + 1.0}},
      {{test_r, 0.0,
        -metadata.readout_abs_z_cm - 1.0}},
      {{test_r, 0.0,
        0.5 * metadata.central_membrane_half_width_cm}},
      {{test_r, 0.0,
        -0.5 * metadata.central_membrane_half_width_cm}},
  }};
  for (const auto& query : outside_queries)
  {
    if (map.query_garfield_transport(
            query[0], query[1], query[2], outside) !=
            QueryStatus::OutsideDriftVolume ||
        map.get_garfield_transport(
            query[0], query[1], query[2], outside))
    {
      std::cerr
          << "Out-of-volume query was not classified and rejected\n";
      return 1;
    }
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (map.query_garfield_transport(
          test_r,
          nan,
          metadata.central_membrane_half_width_cm,
          outside) != QueryStatus::Invalid)
  {
    std::cerr << "Nonfinite LUT query was not invalid\n";
    return 1;
  }

  if (!found_all_reach || !checked_readout_zero)
  {
    std::cerr << "Map did not provide the expected reaching test cell\n";
    return 1;
  }
  return 0;
}
