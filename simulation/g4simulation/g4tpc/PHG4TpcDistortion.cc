// $Id: $

/*!
 * \file PHG4TpcDistortion.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>, Henry Klest <henry.klest@stonybrook.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4TpcDistortion.h"

#include <phgarfield/PHGarfieldForwardMap.h>
#include <phgarfield/PHGarfieldAdaptiveMap.h>
#include <phgarfield/PHGarfieldCompactMap.h>

#include <TFile.h>
#include <TH3.h>
#include <TNamed.h>
#include <TMD5.h>
#include <TParameter.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>    // for sqrt, fabs, NAN
#include <cstdlib>  // for exit
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
  PHGarfieldCompactMap::Limits compact_limits(const PHG4TpcDistortion::GarfieldCompactReadLimits& limits)
  {
    PHGarfieldCompactMap::Limits result;
    result.max_source_points = limits.max_source_points;
    result.max_nodes = limits.max_nodes;
    result.max_corners = limits.max_corners;
    result.max_values = limits.max_values;
    result.max_bytes = limits.max_bytes;
    return result;
  }

  template <class T>
  constexpr T square(const T& x)
  {
    return x * x;
  }

  // check boundaries in axis
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries(const TAxis* axis, double value)
  {
    const auto bin = axis->FindBin(value);
    return (bin >= 2 && bin < axis->GetNbins());
  }

  // check boundaries in histogram, before interpolation
  /* for the interpolation to work, the value must be within the range of the provided axis, and not into the first and last bin */
  inline bool check_boundaries(const TH3* h, double r, double phi, double z)
  {
    return check_boundaries(h->GetXaxis(), r) && check_boundaries(h->GetYaxis(), phi) && check_boundaries(h->GetZaxis(), z);
  }

  bool same_axis(const TAxis* first, const TAxis* second)
  {
    if (!first || !second || first->GetNbins() != second->GetNbins())
    {
      return false;
    }
    constexpr double tolerance = 1.0e-12;
    for (int bin = 1; bin <= first->GetNbins() + 1; ++bin)
    {
      const double first_edge = first->GetBinLowEdge(bin);
      const double second_edge = second->GetBinLowEdge(bin);
      const double scale =
          std::max({1.0, std::abs(first_edge), std::abs(second_edge)});
      if (std::abs(first_edge - second_edge) > tolerance * scale)
      {
        return false;
      }
    }
    return true;
  }

  bool same_node_axes(const TH3* first, const TH3* second)
  {
    return first && second &&
           same_axis(first->GetXaxis(), second->GetXaxis()) &&
           same_axis(first->GetYaxis(), second->GetYaxis()) &&
           same_axis(first->GetZaxis(), second->GetZaxis());
  }

  bool close_value(
      const double first,
      const double second,
      const double tolerance = 1.0e-12)
  {
    const double scale =
        std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <= tolerance * scale;
  }

  bool cell_edges_match_node_centers(
      const TAxis* cell_axis,
      const TAxis* node_axis)
  {
    if (!cell_axis || !node_axis ||
        node_axis->GetNbins() < cell_axis->GetNbins() + 2)
    {
      return false;
    }
    for (int edge = 1; edge <= cell_axis->GetNbins() + 1; ++edge)
    {
      if (!close_value(
              cell_axis->GetBinLowEdge(edge),
              node_axis->GetBinCenter(edge + 1)))
      {
        return false;
      }
    }
    return true;
  }

  bool guard_nodes_match(const TH3* histogram)
  {
    if (!histogram)
    {
      return false;
    }
    const int nx = histogram->GetNbinsX();
    const int ny = histogram->GetNbinsY();
    const int nz = histogram->GetNbinsZ();
    if (nx < 3 || ny < 3 || nz < 3)
    {
      return false;
    }

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
          if (histogram->GetBinContent(ix, iy, iz) !=
              histogram->GetBinContent(source_x, source_y, source_z))
          {
            return false;
          }
        }
      }
    }
    return true;
  }

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

  bool is_lowercase_hex_digest(
      const std::string& value,
      const std::size_t length)
  {
    return value.size() == length &&
           std::all_of(
               value.begin(),
               value.end(),
               [](const unsigned char character)
               {
                 return std::isdigit(character) ||
                        (character >= 'a' && character <= 'f');
               });
  }

  template <class T>
  bool read_parameter(TFile* file, const char* path, T& value)
  {
    const auto* parameter =
        dynamic_cast<TParameter<T>*>(file ? file->Get(path) : nullptr);
    if (!parameter)
    {
      return false;
    }
    value = parameter->GetVal();
    return true;
  }

  bool read_named(TFile* file, const char* path, std::string& value)
  {
    const auto* named =
        dynamic_cast<TNamed*>(file ? file->Get(path) : nullptr);
    if (!named)
    {
      return false;
    }
    value = named->GetTitle();
    return true;
  }

  bool finite_metadata(
      const PHG4TpcDistortion::GarfieldLutMetadata& metadata)
  {
    return metadata.schema_version > 0 &&
           metadata.phi_nodes >= 3 &&
           metadata.radial_nodes >= 2 &&
           metadata.longitudinal_nodes >= 2 &&
           std::isfinite(metadata.gas_inner_radius_cm) &&
           std::isfinite(metadata.gas_outer_radius_cm) &&
           std::isfinite(metadata.central_membrane_half_width_cm) &&
           std::isfinite(metadata.readout_abs_z_cm) &&
           std::isfinite(metadata.integration_step_ns) &&
           std::isfinite(
               metadata.nominal_phg4_drift_velocity_cm_per_ns) &&
           std::isfinite(metadata.maximum_nominal_time_mismatch) &&
           std::isfinite(
               metadata.scale_zero_reference_time_side0_ns) &&
           std::isfinite(
               metadata.scale_zero_reference_time_side1_ns) &&
           std::isfinite(
               metadata.scale_zero_maximum_time_mismatch) &&
           std::isfinite(metadata.maximum_drift_time_ns) &&
           std::isfinite(metadata.tpc_place_x_cm) &&
           std::isfinite(metadata.tpc_place_y_cm) &&
           std::isfinite(metadata.tpc_place_z_cm) &&
           std::isfinite(metadata.tpc_rotation_x_rad) &&
           std::isfinite(metadata.tpc_rotation_y_rad) &&
           std::isfinite(metadata.tpc_rotation_z_rad) &&
           std::isfinite(metadata.nominal_electric_field_vcm) &&
           std::isfinite(metadata.space_charge_scale_side0) &&
           std::isfinite(metadata.space_charge_scale_side1) &&
           std::isfinite(metadata.magnetic_field_scale) &&
           metadata.gas_inner_radius_cm > 0.0 &&
           metadata.gas_outer_radius_cm >
               metadata.gas_inner_radius_cm &&
           metadata.central_membrane_half_width_cm > 0.0 &&
           metadata.readout_abs_z_cm >
               metadata.central_membrane_half_width_cm &&
           metadata.integration_step_ns > 0.0 &&
           metadata.nominal_phg4_drift_velocity_cm_per_ns > 0.0 &&
           metadata.maximum_nominal_time_mismatch >= 0.0 &&
           metadata.scale_zero_reference_time_side0_ns > 0.0 &&
           metadata.scale_zero_reference_time_side1_ns > 0.0 &&
           metadata.scale_zero_maximum_time_mismatch >= 0.0 &&
           metadata.scale_zero_maximum_time_mismatch <=
               metadata.maximum_nominal_time_mismatch &&
           metadata.maximum_drift_time_ns >= 0.0 &&
           metadata.nominal_electric_field_vcm > 0.0 &&
           metadata.maximum_integration_steps > 0;
  }

  // print histogram
  [[maybe_unused]] void print_histogram(TH3* h)
  {
    std::cout << "PHG4TpcDistortion::print_histogram - name: " << h->GetName() << std::endl;
    for (const auto& axis : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()})
    {
      std::cout
          << "  " << axis->GetName()
          << " bins: " << axis->GetNbins()
          << " min: " << axis->GetXmin()
          << " max: " << axis->GetXmax()
          << std::endl;
    }
    std::cout << std::endl;
  }

}  // namespace

struct PHG4TpcDistortion::AdaptiveTransportData
{
  PHGarfieldAdaptiveMap::Data data;
};

struct PHG4TpcDistortion::CompactTransportData
{
  PHGarfieldCompactMap::Data data;
};

PHG4TpcDistortion::PHG4TpcDistortion()
{
  const PHGarfieldCompactMap::Limits defaults;
  m_compactReadLimits = {defaults.max_source_points, defaults.max_nodes,
      defaults.max_corners, defaults.max_values, defaults.max_bytes};
}
PHG4TpcDistortion::~PHG4TpcDistortion() = default;

void PHG4TpcDistortion::set_garfield_compact_read_limits(
    const std::uint64_t max_source_points, const std::uint64_t max_nodes,
    const std::uint64_t max_corners, const std::uint64_t max_values,
    const std::uint64_t max_bytes)
{
  if (m_initialized) { throw std::logic_error("compact read limits must be set before Init"); }
  const GarfieldCompactReadLimits requested{max_source_points, max_nodes, max_corners, max_values, max_bytes};
  std::ostringstream error;
  if (!PHGarfieldCompactMap::ValidateLimits(compact_limits(requested), error))
  { throw std::invalid_argument("invalid compact read limits: " + error.str()); }
  m_compactReadLimits = requested;
  m_compactReadLimitsExplicit = true;
}

//__________________________________________________________________________________________________________
void PHG4TpcDistortion::Init()
{
  m_initialized = false;
  m_garfield_metadata = GarfieldLutMetadata{};
  m_garfieldAdaptive.reset();
  m_garfieldCompact.reset();
  m_garfield_interpolation_transverse_error = {{nullptr, nullptr}};
  m_garfield_interpolation_time_error = {{nullptr, nullptr}};
  m_boundaryGarfield.reset();
  m_boundaryFallbackStats = BoundaryFallbackStats{};
  if (m_compactReadLimitsExplicit && m_application_mode != ApplicationMode::GarfieldFullTransport)
  {
    std::cerr << "Explicit compact read limits require GarfieldFullTransport application mode" << std::endl;
    exit(1);
  }
  std::cout << "PHG4TpcDistortion::Init - m_phi_hist_in_radians: " << m_phi_hist_in_radians << std::endl;

  if (m_application_mode == ApplicationMode::GarfieldFullTransport)
  {
    if (!m_do_static_distortions || m_do_time_ordered_distortions)
    {
      std::cerr
          << "PHG4TpcDistortion::Init - Garfield full transport requires "
             "exactly one static map and forbids time-ordered map addition"
          << std::endl;
      exit(1);
    }

    std::cout
        << "PHG4TpcDistortion::Init - Garfield full transport map: "
        << m_static_distortion_filename << std::endl;
    m_static_tfile = std::make_unique<TFile>(
        m_static_distortion_filename.c_str(), "READ");
    if (!m_static_tfile->IsOpen() || m_static_tfile->IsZombie() ||
        !initialize_garfield_full_transport())
    {
      std::cerr
          << "PHG4TpcDistortion::Init - invalid Garfield full transport map"
          << std::endl;
      exit(1);
    }
    if (m_garfieldBoundaryFallback && !initialize_garfield_boundary_fallback())
    {
      std::cerr << "PHG4TpcDistortion::Init - invalid direct boundary fallback configuration" << std::endl;
      exit(1);
    }
    m_initialized = true;
    return;
  }

  if (m_do_static_distortions)
  {
    std::cout << "PHG4TpcDistortion::Init - m_static_distortion_filename: " << m_static_distortion_filename << std::endl;
    m_static_tfile.reset(new TFile(m_static_distortion_filename.c_str()));
    if (!m_static_tfile->IsOpen())
    {
      std::cout << "PHG4TpcDistortion::Init - Static distortion file could not be opened!" << std::endl;
      exit(1);
    }

    // Open Static Space Charge Maps
    hDRint[0] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionR_negz"));
    hDRint[1] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionR_posz"));
    hDPint[0] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionP_negz"));
    hDPint[1] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionP_posz"));
    hDZint[0] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionZ_negz"));
    hDZint[1] = dynamic_cast<TH3*>(m_static_tfile->Get("hIntDistortionZ_posz"));
    if (m_do_ReachesReadout)
    {
      hReach[0] = dynamic_cast<TH3*>(m_static_tfile->Get("hReachesReadout_negz"));
      hReach[1] = dynamic_cast<TH3*>(m_static_tfile->Get("hReachesReadout_posz"));
    }
  }

  if (m_do_time_ordered_distortions)
  {
    std::cout << "PHG4TpcDistortion::Init - m_time_ordered_distortion_filename: " << m_time_ordered_distortion_filename << std::endl;
    m_time_ordered_tfile.reset(new TFile(m_time_ordered_distortion_filename.c_str()));
    if (!m_time_ordered_tfile->IsOpen())
    {
      std::cout << "PHG4TpcDistortion::Init - TimeOrdered distortion file could not be opened!" << std::endl;
      exit(1);
    }

    TimeTree = static_cast<TTree*>(m_time_ordered_tfile->Get("TimeDists"));
    if (!TimeTree)
    {
      std::cout << "PHG4TpcDistortion::Init - TimeOrdered distortion tree could not be found!" << std::endl;
      exit(1);
    }

    // create histograms
    TimehDR[0] = new TH3F();
    TimehDR[1] = new TH3F();
    TimehDP[0] = new TH3F();
    TimehDP[1] = new TH3F();
    TimehDZ[0] = new TH3F();
    TimehDZ[1] = new TH3F();
    TimehRR[0] = new TH3F();  // RR stands for ReachesReadout
    TimehRR[1] = new TH3F();  // RR stands for ReachesReadout

    // assign to tree branches
    TimeTree->SetBranchAddress("hIntDistortionR_negz", &(TimehDR[0]));
    TimeTree->SetBranchAddress("hIntDistortionR_posz", &(TimehDR[1]));
    TimeTree->SetBranchAddress("hIntDistortionP_negz", &(TimehDP[0]));
    TimeTree->SetBranchAddress("hIntDistortionP_posz", &(TimehDP[1]));
    TimeTree->SetBranchAddress("hIntDistortionZ_negz", &(TimehDZ[0]));
    TimeTree->SetBranchAddress("hIntDistortionZ_posz", &(TimehDZ[1]));
    if (m_do_ReachesReadout)
    {
      TimeTree->SetBranchAddress("hReachesReadout_negz", &(TimehRR[0]));
      TimeTree->SetBranchAddress("hReachesReadout_posz", &(TimehRR[1]));
    }
  }
  m_initialized = true;
}

//__________________________________________________________________________________________________________
bool PHG4TpcDistortion::initialize_garfield_boundary_fallback()
{
  const auto& m = m_garfield_metadata;
  if (m.schema_version < PHGarfieldForwardMap::ExtendedSchemaVersion || m.field_model != "per_side_3d")
  {
    std::cerr << "Direct boundary fallback requires explicit schema-2 per-side 3D fields" << std::endl;
    return false;
  }
  const auto verify = [](const std::string& path, const std::string& expected)
  {
    if (path.empty()) { return expected.empty(); }
    std::unique_ptr<TMD5> actual(TMD5::FileChecksum(path.c_str()));
    return actual && expected == actual->AsString();
  };
  bool payloads_valid = verify(m.electric_field_map, m.electric_field_map_md5) &&
      verify(m.frame_field_map, m.frame_field_map_md5) && verify(m.gas_file, m.gas_file_md5) &&
      verify(m.magnetic_field_map, m.magnetic_field_map_md5);
  for (std::size_t side = 0; side < 2; ++side)
  {
    payloads_valid = payloads_valid && verify(m.electric_field_map_3d[side], m.electric_field_map_3d_md5[side]) &&
        verify(m.frame_field_map_3d[side], m.frame_field_map_3d_md5[side]);
  }
  if (!payloads_valid)
  {
    std::cerr << "Direct boundary fallback payload checksum mismatch" << std::endl;
    return false;
  }
  auto g = std::make_unique<PHGarfield>("GarfieldLutBoundaryFallback", m.electric_field_map,
      m.space_charge_scale_side0, m.space_charge_scale_side1, m.electric_field_map_3d[0], m.electric_field_map_3d[1]);
  g->SetFrameElectricFieldMap(m.frame_field_map);
  g->SetFrameElectricFieldMap3D(m.frame_field_map_3d[0], m.frame_field_map_3d[1]);
  g->SetFrameChargeScaleSide0(m.frame_charge_scale[0]);
  g->SetFrameChargeScaleSide1(m.frame_charge_scale[1]);
  g->SetCMVoltageDefault(m.nominal_electric_field_vcm);
  g->SetCMFieldRunFactor(m.cm_field_factor[0], m.cm_field_factor[1]);
  g->SetUsePTCorrection(false);
  g->SetUseBCOkEffs(false);
  g->SetUseSurveyGeometry(false);
  g->SetUseCDBDriftCorrections(false);
  g->SetDeltaT0Tbins(0., 0.);
  g->SetUseIFCVoltageDistortion(m.use_ifc_voltage_distortion);
  g->SetUseOFCVoltageDistortion(m.use_ofc_voltage_distortion);
  g->SetFieldCageVoltageOffsets(m.field_cage_offsets_v[0], m.field_cage_offsets_v[1], m.field_cage_offsets_v[2], m.field_cage_offsets_v[3]);
  g->SetLoadTpcFeeMap(false);
  g->SetPublishToNodeTree(false);
  g->SetGasFile(m.gas_file);
  g->SetMagneticFieldMap(m.magnetic_field_map, m.magnetic_field_scale);
  g->SetZeroMagneticField(m.zero_magnetic_field);
  g->SetCentralMembraneHalfWidth(m.central_membrane_half_width_cm);
  g->SetDriftVolumeRadii(m.gas_inner_radius_cm, m.gas_outer_radius_cm);
  g->MoveTpc(m.tpc_place_x_cm, m.tpc_place_y_cm, m.tpc_place_z_cm);
  g->SetTpcRotation(m.tpc_rotation_x_rad, m.tpc_rotation_y_rad, m.tpc_rotation_z_rad);
  if (g->InitRun(nullptr) != 0 || !g->IsReady() ||
      !g->HasLoadedElectricFieldMap3D(0) || !g->HasLoadedElectricFieldMap3D(1))
  { return false; }
  m_boundaryGarfield = std::move(g);
  std::cout << "PHG4TpcDistortion: direct Mixed-cell fallback enabled, Euler step " << m.integration_step_ns << " ns" << std::endl;
  return true;
}

PHG4TpcDistortion::GarfieldTransportQueryStatus PHG4TpcDistortion::query_garfield_boundary_fallback(
    const double r, const double phi, const double z, GarfieldTransportSample& sample, const bool numerical) const
{
  sample = GarfieldTransportSample{};
  sample.used_direct_fallback = true;
  sample.direct_fallback_reason = numerical ? 2U : 1U;
  ++m_boundaryFallbackStats.calls;
  if (numerical) { ++m_boundaryFallbackStats.numerical_calls; }
  else { ++m_boundaryFallbackStats.mixed_calls; }
  if (!m_boundaryGarfield)
  {
    ++m_boundaryFallbackStats.failures;
    return GarfieldTransportQueryStatus::Invalid;
  }
  PHGarfield::ForwardDriftResult result;
  const bool ok = m_boundaryGarfield->ForwardDrift(r * std::cos(phi), r * std::sin(phi), z,
      m_garfield_metadata.readout_abs_z_cm, m_garfield_metadata.integration_step_ns, result, true, true);
  if (ok && result.reached_readout)
  {
    PHGarfieldForwardMap::NodeValue encoded;
    if (PHGarfieldForwardMap::EncodeNode(r, phi, result, encoded))
    {
      sample.reach = ReachCellClass::AllReach;
      sample.delta_r_cm = encoded.delta_r_cm;
      sample.delta_rphi_cm = encoded.delta_rphi_cm;
      sample.drift_time_ns = encoded.drift_time_ns;
      ++m_boundaryFallbackStats.reaching;
      return GarfieldTransportQueryStatus::Success;
    }
  }
  if (result.status == PHGarfield::ForwardDriftStatus::LeftDriftVolume)
  {
    sample.reach = ReachCellClass::AllLost;
    ++m_boundaryFallbackStats.physical_losses;
    return GarfieldTransportQueryStatus::Success;
  }
  ++m_boundaryFallbackStats.failures;
  return GarfieldTransportQueryStatus::Invalid;
}

bool PHG4TpcDistortion::initialize_garfield_full_transport()
{
  constexpr const char* directory = "GarfieldFullTransport/";
  const char* metadata_directory =
      "GarfieldFullTransport/Metadata/";
  const std::array<std::string, 2> suffix = {"negz", "posz"};

  int input_schema = 0;
  if (!read_parameter(m_static_tfile.get(),
                      "GarfieldFullTransport/Metadata/schema_version", input_schema))
  { return false; }

  const bool collection = input_schema == PHGarfieldCompactMap::CollectionSchemaVersion;
  const bool compact = input_schema == PHGarfieldCompactMap::SchemaVersion || collection;
  if (!compact && m_compactReadLimitsExplicit)
  { std::cerr << "Explicit compact read limits require a schema5 transport LUT" << std::endl; return false; }
  if (compact)
  {
    std::string compact_name, layout, training_hash, configuration_hash, content_hash;
    int source_schema = 0, declared_source_schema = 0;
    if (m_static_tfile->TestBit(TFile::kRecovered) ||
        !read_named(m_static_tfile.get(), "GarfieldFullTransport/Metadata/schema_name", compact_name) ||
        compact_name != (collection ? PHGarfieldCompactMap::CollectionSchemaName : PHGarfieldCompactMap::SchemaName) || m_compactFileSha256.empty() ||
        PHGarfieldCompactMap::FileSha256(m_static_distortion_filename) != m_compactFileSha256 ||
        !read_parameter(m_static_tfile.get(), "GarfieldFullTransport/Metadata/source_schema_version", declared_source_schema) ||
        declared_source_schema != PHGarfieldAdaptiveMap::SchemaVersion ||
        !read_named(m_static_tfile.get(), "GarfieldFullTransport/Metadata/compact_layout", layout) ||
        layout != PHGarfieldCompactMap::Layout ||
        !read_named(m_static_tfile.get(), "GarfieldFullTransport/Metadata/source_file_sha256", training_hash) ||
        training_hash != m_compactTrainingSha256 ||
        !read_named(m_static_tfile.get(), "GarfieldFullTransport/Metadata/source_configuration_fingerprint_fnv1a64", configuration_hash) ||
        configuration_hash != m_compactConfigurationFingerprint ||
        !read_named(m_static_tfile.get(), "GarfieldFullTransport/Metadata/source_content_fingerprint_fnv1a64", content_hash) ||
        content_hash != m_compactContentFingerprint ||
        !read_parameter(m_static_tfile.get(), "GarfieldFullTransport/SourceMetadata/schema_version", source_schema) ||
        source_schema != PHGarfieldAdaptiveMap::SchemaVersion)
    { std::cerr << "Invalid or externally unapproved compact transport payload" << std::endl; return false; }
    if (collection)
    {
      std::string source_kind, manifest;
      if (!read_named(m_static_tfile.get(), "GarfieldFullTransport/Metadata/source_kind", source_kind) ||
          source_kind != PHGarfieldCompactMap::CollectionSourceKind ||
          !read_named(m_static_tfile.get(), "GarfieldFullTransport/Metadata/source_manifest_tsv", manifest) ||
          PHGarfieldCompactMap::Sha256(manifest) != training_hash ||
          hexadecimal(fnv1a64(manifest)) != content_hash)
      { std::cerr << "Invalid compact shard-collection binding" << std::endl; return false; }
    }
    // All physical/payload/canonical-configuration checks below apply to the
    // immutable source schema4 metadata. Runtime identity is set only after Read.
    metadata_directory = "GarfieldFullTransport/SourceMetadata/";
  }

  if (input_schema != PHGarfieldAdaptiveMap::SchemaVersion && !compact)
  {
    for (unsigned int side = 0; side < 2; ++side)
    {
      m_garfield_delta_r[side] = dynamic_cast<TH3*>(
          m_static_tfile->Get(
              (std::string(directory) + "hDeltaR_" + suffix[side]).c_str()));
      m_garfield_delta_rphi[side] = dynamic_cast<TH3*>(
          m_static_tfile->Get(
              (std::string(directory) + "hDeltaRPhi_" + suffix[side]).c_str()));
      m_garfield_drift_time[side] = dynamic_cast<TH3*>(
          m_static_tfile->Get(
              (std::string(directory) + "hDriftTime_" + suffix[side]).c_str()));
      m_garfield_reach_node[side] = dynamic_cast<TH3*>(
          m_static_tfile->Get(
              (std::string(directory) +
               "hReachesReadoutNode_" + suffix[side])
                  .c_str()));
      m_garfield_reach_cell_class[side] = dynamic_cast<TH3*>(
          m_static_tfile->Get(
              (std::string(directory) +
               "hReachCellClass_" + suffix[side])
                  .c_str()));

      if (!m_garfield_delta_r[side] ||
          !m_garfield_delta_rphi[side] ||
          !m_garfield_drift_time[side] ||
          !m_garfield_reach_node[side] ||
          !m_garfield_reach_cell_class[side])
      {
        std::cerr
            << "PHG4TpcDistortion::initialize_garfield_full_transport - "
               "missing histogram for "
            << suffix[side] << std::endl;
        return false;
      }
    }

  }

  std::string schema_name;
  std::string transport_contract;
  std::string coordinate_frame;
  std::string axis_order_units;
  std::string transverse_quantity;
  std::string time_quantity;
  std::string reach_encoding;
  std::string canonical_configuration;
  std::string cdb_timestamp_text;
  if (!read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "schema_name").c_str(),
          schema_name) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "transport_contract").c_str(),
          transport_contract) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "coordinate_frame").c_str(),
          coordinate_frame) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "axis_order_units").c_str(),
          axis_order_units) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "transverse_quantity").c_str(),
          transverse_quantity) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "time_quantity").c_str(),
          time_quantity) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "reach_encoding").c_str(),
          reach_encoding) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "integrator").c_str(),
          m_garfield_metadata.integrator) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) +
           "magnetic_field_domain_policy")
              .c_str(),
          m_garfield_metadata.magnetic_field_domain_policy) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) +
           "nominal_time_reference_policy")
              .c_str(),
          m_garfield_metadata.nominal_time_reference_policy) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "reach_probe_scheme").c_str(),
          m_garfield_metadata.reach_probe_scheme) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "electric_field_map").c_str(),
          m_garfield_metadata.electric_field_map) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "magnetic_field_map").c_str(),
          m_garfield_metadata.magnetic_field_map) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "gas_file").c_str(),
          m_garfield_metadata.gas_file) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "electric_field_map_md5").c_str(),
          m_garfield_metadata.electric_field_map_md5) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "magnetic_field_map_md5").c_str(),
          m_garfield_metadata.magnetic_field_map_md5) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "gas_file_md5").c_str(),
          m_garfield_metadata.gas_file_md5) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "cdb_global_tag").c_str(),
          m_garfield_metadata.cdb_global_tag) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "configuration_label").c_str(),
          m_garfield_metadata.configuration_label) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "cdb_timestamp").c_str(),
          cdb_timestamp_text) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) + "canonical_configuration").c_str(),
          canonical_configuration) ||
      !read_named(
          m_static_tfile.get(),
          (std::string(metadata_directory) +
           "configuration_fingerprint_fnv1a64")
              .c_str(),
          m_garfield_metadata.configuration_fingerprint))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "missing required metadata"
        << std::endl;
    return false;
  }

  if (schema_name != (collection ? PHGarfieldCompactMap::CollectionSourceName : PHGarfieldForwardMap::SchemaName) ||
      transport_contract != PHGarfieldForwardMap::TransportContract ||
      coordinate_frame != "TPC_ENVELOPE_LOCAL" ||
      axis_order_units != "phi_rad,r_cm,z_cm" ||
      transverse_quantity != "delta_r_cm,r0_delta_phi_cm" ||
      time_quantity != "absolute_drift_duration_ns" ||
      (reach_encoding != "0=AllLost,1=Mixed,2=AllReach" && reach_encoding != "0=AllLost,1=Mixed,2=AllReach,3=NumericallyUnresolved") ||
      m_garfield_metadata.integrator != "explicit_euler" ||
      m_garfield_metadata.magnetic_field_domain_policy !=
          PHGarfieldForwardMap::MagneticFieldDomainPolicy ||
      m_garfield_metadata.nominal_time_reference_policy !=
          PHGarfieldForwardMap::NominalTimeReferencePolicy ||
      canonical_configuration.empty())
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "unsupported schema semantics"
        << std::endl;
    return false;
  }

  const auto metadata_path =
      [metadata_directory](const char* name)
  {
    return std::string(metadata_directory) + name;
  };
  int zero_magnetic_field = -1;
  int cell_center_probes = -1;
  if (!read_parameter(
          m_static_tfile.get(),
          metadata_path("schema_version").c_str(),
          m_garfield_metadata.schema_version) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("phi_nodes").c_str(),
          m_garfield_metadata.phi_nodes) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("radial_nodes").c_str(),
          m_garfield_metadata.radial_nodes) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("longitudinal_nodes").c_str(),
          m_garfield_metadata.longitudinal_nodes) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("gas_inner_radius_cm").c_str(),
          m_garfield_metadata.gas_inner_radius_cm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("gas_outer_radius_cm").c_str(),
          m_garfield_metadata.gas_outer_radius_cm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("central_membrane_half_width_cm").c_str(),
          m_garfield_metadata.central_membrane_half_width_cm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("readout_abs_z_cm").c_str(),
          m_garfield_metadata.readout_abs_z_cm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("integration_step_ns").c_str(),
          m_garfield_metadata.integration_step_ns) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("nominal_phg4_drift_velocity_cm_per_ns").c_str(),
          m_garfield_metadata.nominal_phg4_drift_velocity_cm_per_ns) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("maximum_nominal_time_mismatch").c_str(),
          m_garfield_metadata.maximum_nominal_time_mismatch) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("scale_zero_reference_time_side0_ns").c_str(),
          m_garfield_metadata.scale_zero_reference_time_side0_ns) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("scale_zero_reference_time_side1_ns").c_str(),
          m_garfield_metadata.scale_zero_reference_time_side1_ns) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("scale_zero_maximum_time_mismatch").c_str(),
          m_garfield_metadata.scale_zero_maximum_time_mismatch) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("maximum_drift_time_ns").c_str(),
          m_garfield_metadata.maximum_drift_time_ns) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("tpc_place_x_cm").c_str(),
          m_garfield_metadata.tpc_place_x_cm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("tpc_place_y_cm").c_str(),
          m_garfield_metadata.tpc_place_y_cm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("tpc_place_z_cm").c_str(),
          m_garfield_metadata.tpc_place_z_cm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("tpc_rotation_x_rad").c_str(),
          m_garfield_metadata.tpc_rotation_x_rad) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("tpc_rotation_y_rad").c_str(),
          m_garfield_metadata.tpc_rotation_y_rad) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("tpc_rotation_z_rad").c_str(),
          m_garfield_metadata.tpc_rotation_z_rad) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("nominal_electric_field_vcm").c_str(),
          m_garfield_metadata.nominal_electric_field_vcm) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("space_charge_scale_side0").c_str(),
          m_garfield_metadata.space_charge_scale_side0) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("space_charge_scale_side1").c_str(),
          m_garfield_metadata.space_charge_scale_side1) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("magnetic_field_scale").c_str(),
          m_garfield_metadata.magnetic_field_scale) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("maximum_integration_steps").c_str(),
          m_garfield_metadata.maximum_integration_steps) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("zero_magnetic_field").c_str(),
          zero_magnetic_field) ||
      !read_parameter(
          m_static_tfile.get(),
          metadata_path("cell_center_probes").c_str(),
          cell_center_probes))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "missing numeric metadata"
        << std::endl;
    return false;
  }
  try
  {
    std::size_t parsed = 0;
    m_garfield_metadata.cdb_timestamp =
        std::stoull(cdb_timestamp_text, &parsed);
    if (parsed != cdb_timestamp_text.size())
    {
      throw std::invalid_argument("trailing timestamp characters");
    }
  }
  catch (const std::exception&)
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "invalid CDB timestamp metadata"
        << std::endl;
    return false;
  }

  if ((m_garfield_metadata.schema_version != PHGarfieldForwardMap::SchemaVersion &&
       m_garfield_metadata.schema_version != PHGarfieldForwardMap::ExtendedSchemaVersion &&
       m_garfield_metadata.schema_version != PHGarfieldForwardMap::NumericalSchemaVersion &&
       m_garfield_metadata.schema_version != PHGarfieldAdaptiveMap::SchemaVersion) ||
      (zero_magnetic_field != 0 && zero_magnetic_field != 1) ||
      (cell_center_probes != 0 && cell_center_probes != 1) ||
      !finite_metadata(m_garfield_metadata))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "invalid or unsupported metadata values"
        << std::endl;
    return false;
  }
  if (reach_encoding != (m_garfield_metadata.schema_version >= PHGarfieldForwardMap::NumericalSchemaVersion ?
      "0=AllLost,1=Mixed,2=AllReach,3=NumericallyUnresolved" : "0=AllLost,1=Mixed,2=AllReach"))
  { std::cerr << "Reach encoding does not match schema" << std::endl; return false; }
  m_garfield_metadata.zero_magnetic_field =
      zero_magnetic_field != 0;
  m_garfield_metadata.cell_center_probes =
      cell_center_probes != 0;
  const std::string expected_probe_scheme =
      m_garfield_metadata.schema_version == PHGarfieldAdaptiveMap::SchemaVersion
          ? PHGarfieldAdaptiveMap::ProbeScheme
          : m_garfield_metadata.cell_center_probes
          ? "corners+center+faces+cage-biased-v1"
          : "corners-only-debug";
  if (m_garfield_metadata.reach_probe_scheme != expected_probe_scheme)
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "reach-probe metadata is inconsistent"
        << std::endl;
    return false;
  }
  const double nominal_reference_time_ns =
      (m_garfield_metadata.readout_abs_z_cm -
       m_garfield_metadata.central_membrane_half_width_cm) /
      m_garfield_metadata.nominal_phg4_drift_velocity_cm_per_ns;
  const double observed_reference_mismatch =
      std::max(
          std::abs(
              m_garfield_metadata.scale_zero_reference_time_side0_ns -
              nominal_reference_time_ns),
          std::abs(
              m_garfield_metadata.scale_zero_reference_time_side1_ns -
              nominal_reference_time_ns)) /
      nominal_reference_time_ns;
  if (!close_value(
          observed_reference_mismatch,
          m_garfield_metadata.scale_zero_maximum_time_mismatch))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "scale-zero timing metadata is internally inconsistent"
        << std::endl;
    return false;
  }
  if (m_garfield_metadata.gas_file.empty() ||
      !is_lowercase_hex_digest(
          m_garfield_metadata.gas_file_md5, 32) ||
      (m_garfield_metadata.electric_field_map.empty() !=
       m_garfield_metadata.electric_field_map_md5.empty()) ||
      (!m_garfield_metadata.electric_field_map_md5.empty() &&
       !is_lowercase_hex_digest(
           m_garfield_metadata.electric_field_map_md5, 32)) ||
      (!m_garfield_metadata.magnetic_field_map_md5.empty() &&
       !is_lowercase_hex_digest(
           m_garfield_metadata.magnetic_field_map_md5, 32)) ||
      (!m_garfield_metadata.zero_magnetic_field &&
       (m_garfield_metadata.magnetic_field_map.empty() ||
        !is_lowercase_hex_digest(
            m_garfield_metadata.magnetic_field_map_md5, 32))))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "invalid input-payload provenance"
        << std::endl;
    return false;
  }

  if (m_garfield_metadata.schema_version >= PHGarfieldForwardMap::ExtendedSchemaVersion)
  {
    auto& metadata = m_garfield_metadata;
    int use_ifc = -1;
    int use_ofc = -1;
    bool valid = read_named(m_static_tfile.get(), metadata_path("field_model").c_str(), metadata.field_model) &&
        read_named(m_static_tfile.get(), metadata_path("frame_field_map").c_str(), metadata.frame_field_map) &&
        read_named(m_static_tfile.get(), metadata_path("frame_field_map_md5").c_str(), metadata.frame_field_map_md5) &&
        read_parameter(m_static_tfile.get(), metadata_path("use_ifc_voltage_distortion").c_str(), use_ifc) &&
        read_parameter(m_static_tfile.get(), metadata_path("use_ofc_voltage_distortion").c_str(), use_ofc);
    const auto payload_valid = [](const std::string& path, const std::string& digest)
    {
      return path.empty() ? digest.empty() : is_lowercase_hex_digest(digest, 32);
    };
    valid = valid && payload_valid(metadata.frame_field_map, metadata.frame_field_map_md5) &&
        (use_ifc == 0 || use_ifc == 1) && (use_ofc == 0 || use_ofc == 1) &&
        (metadata.field_model == "per_side_3d" || metadata.field_model == "axisymmetric_2d" || metadata.field_model == "nominal");
    for (std::size_t side = 0; side < 2; ++side)
    {
      const std::string side_tag = "_side" + std::to_string(side);
      const auto path = [&metadata_path, &side_tag](const char* key) { return metadata_path((std::string(key) + side_tag).c_str()); };
      valid = read_named(m_static_tfile.get(), path("electric_field_map_3d").c_str(), metadata.electric_field_map_3d[side]) && valid;
      valid = read_named(m_static_tfile.get(), path("electric_field_map_3d_md5").c_str(), metadata.electric_field_map_3d_md5[side]) && valid;
      valid = read_named(m_static_tfile.get(), path("frame_field_map_3d").c_str(), metadata.frame_field_map_3d[side]) && valid;
      valid = read_named(m_static_tfile.get(), path("frame_field_map_3d_md5").c_str(), metadata.frame_field_map_3d_md5[side]) && valid;
      valid = read_parameter(m_static_tfile.get(), path("frame_charge_scale").c_str(), metadata.frame_charge_scale[side]) && valid;
      valid = read_parameter(m_static_tfile.get(), path("cm_field_factor").c_str(), metadata.cm_field_factor[side]) && valid;
      valid = valid && payload_valid(metadata.electric_field_map_3d[side], metadata.electric_field_map_3d_md5[side]) &&
          payload_valid(metadata.frame_field_map_3d[side], metadata.frame_field_map_3d_md5[side]) &&
          std::isfinite(metadata.frame_charge_scale[side]) && std::isfinite(metadata.cm_field_factor[side]) && metadata.cm_field_factor[side] > 0.0;
      valid = valid && ((metadata.field_model == "per_side_3d") != metadata.electric_field_map_3d[side].empty());
    }
    valid = valid && (metadata.frame_field_map_3d[0].empty() == metadata.frame_field_map_3d[1].empty());
    for (std::size_t index = 0; index < 4; ++index)
    {
      valid = read_parameter(m_static_tfile.get(), metadata_path(("field_cage_offset_v_" + std::to_string(index)).c_str()).c_str(),
                             metadata.field_cage_offsets_v[index]) && valid;
      valid = valid && std::isfinite(metadata.field_cage_offsets_v[index]);
    }
    if (!valid)
    {
      std::cerr << "PHG4TpcDistortion: missing or invalid schema-2 field provenance" << std::endl;
      return false;
    }
    metadata.use_ifc_voltage_distortion = use_ifc != 0;
    metadata.use_ofc_voltage_distortion = use_ofc != 0;
  }

  if (m_garfield_metadata.schema_version >= PHGarfieldForwardMap::NumericalSchemaVersion)
  {
    auto& metadata = m_garfield_metadata;
    if (!read_parameter(m_static_tfile.get(), metadata_path("interpolation_tolerance_um").c_str(), metadata.interpolation_tolerance_um) ||
        !read_named(m_static_tfile.get(), metadata_path("interpolation_error_probe_scheme").c_str(), metadata.interpolation_error_probe_scheme) ||
        !std::isfinite(metadata.interpolation_tolerance_um) || metadata.interpolation_tolerance_um <= 0. || !metadata.cell_center_probes ||
        metadata.interpolation_error_probe_scheme !=
            (metadata.schema_version == PHGarfieldAdaptiveMap::SchemaVersion
                 ? PHGarfieldAdaptiveMap::ProbeScheme
                 : "center_and_face_centers_transverse_time_equivalent_v1"))
    { std::cerr << "Invalid numerical-mask provenance" << std::endl; return false; }
    if (metadata.schema_version == PHGarfieldForwardMap::NumericalSchemaVersion)
    {
      for (unsigned int side = 0; side < 2; ++side)
      {
        m_garfield_interpolation_transverse_error[side] = dynamic_cast<TH3*>(m_static_tfile->Get((std::string(directory)+"hInterpolationTransverseError_"+suffix[side]).c_str()));
        m_garfield_interpolation_time_error[side] = dynamic_cast<TH3*>(m_static_tfile->Get((std::string(directory)+"hInterpolationTimeEquivalentError_"+suffix[side]).c_str()));
        if (!m_garfield_interpolation_transverse_error[side] || !m_garfield_interpolation_time_error[side])
        { std::cerr << "Missing numerical-mask histograms" << std::endl; return false; }
      }
    }
  }

  PHGarfieldForwardMap::Config fingerprint_config;
  fingerprint_config.interpolation_tolerance_um = m_garfield_metadata.interpolation_tolerance_um;
  fingerprint_config.phi_nodes =
      static_cast<unsigned int>(m_garfield_metadata.phi_nodes);
  fingerprint_config.radial_nodes =
      static_cast<unsigned int>(m_garfield_metadata.radial_nodes);
  fingerprint_config.longitudinal_nodes =
      static_cast<unsigned int>(m_garfield_metadata.longitudinal_nodes);
  fingerprint_config.inner_radius_cm =
      m_garfield_metadata.gas_inner_radius_cm;
  fingerprint_config.outer_radius_cm =
      m_garfield_metadata.gas_outer_radius_cm;
  fingerprint_config.central_membrane_half_width_cm =
      m_garfield_metadata.central_membrane_half_width_cm;
  fingerprint_config.readout_abs_z_cm =
      m_garfield_metadata.readout_abs_z_cm;
  fingerprint_config.integration_step_ns =
      m_garfield_metadata.integration_step_ns;
  fingerprint_config.nominal_phg4_drift_velocity_cm_per_ns =
      m_garfield_metadata.nominal_phg4_drift_velocity_cm_per_ns;
  fingerprint_config.maximum_nominal_time_mismatch =
      m_garfield_metadata.maximum_nominal_time_mismatch;
  fingerprint_config.tpc_place_x_cm =
      m_garfield_metadata.tpc_place_x_cm;
  fingerprint_config.tpc_place_y_cm =
      m_garfield_metadata.tpc_place_y_cm;
  fingerprint_config.tpc_place_z_cm =
      m_garfield_metadata.tpc_place_z_cm;
  fingerprint_config.tpc_rotation_x_rad =
      m_garfield_metadata.tpc_rotation_x_rad;
  fingerprint_config.tpc_rotation_y_rad =
      m_garfield_metadata.tpc_rotation_y_rad;
  fingerprint_config.tpc_rotation_z_rad =
      m_garfield_metadata.tpc_rotation_z_rad;
  fingerprint_config.cdb_global_tag =
      m_garfield_metadata.cdb_global_tag;
  fingerprint_config.cdb_timestamp =
      m_garfield_metadata.cdb_timestamp;
  fingerprint_config.probe_cell_centers =
      m_garfield_metadata.cell_center_probes;
  fingerprint_config.configuration_label =
      m_garfield_metadata.configuration_label;

  PHGarfieldForwardMap::Provenance fingerprint_provenance;
  fingerprint_provenance.nominal_electric_field_vcm =
      m_garfield_metadata.nominal_electric_field_vcm;
  fingerprint_provenance.space_charge_scale_side0 =
      m_garfield_metadata.space_charge_scale_side0;
  fingerprint_provenance.space_charge_scale_side1 =
      m_garfield_metadata.space_charge_scale_side1;
  fingerprint_provenance.magnetic_field_scale =
      m_garfield_metadata.magnetic_field_scale;
  fingerprint_provenance.zero_magnetic_field =
      m_garfield_metadata.zero_magnetic_field;
  fingerprint_provenance.electric_field_map =
      m_garfield_metadata.electric_field_map;
  fingerprint_provenance.magnetic_field_map =
      m_garfield_metadata.magnetic_field_map;
  fingerprint_provenance.gas_file =
      m_garfield_metadata.gas_file;
  fingerprint_provenance.electric_field_map_md5 =
      m_garfield_metadata.electric_field_map_md5;
  fingerprint_provenance.magnetic_field_map_md5 =
      m_garfield_metadata.magnetic_field_map_md5;
  fingerprint_provenance.gas_file_md5 =
      m_garfield_metadata.gas_file_md5;
  fingerprint_provenance.maximum_integration_steps =
      static_cast<unsigned int>(
          m_garfield_metadata.maximum_integration_steps);
  fingerprint_provenance.reach_probe_scheme =
      m_garfield_metadata.reach_probe_scheme;
  fingerprint_provenance.scale_zero_reference_time_side0_ns =
      m_garfield_metadata.scale_zero_reference_time_side0_ns;
  fingerprint_provenance.scale_zero_reference_time_side1_ns =
      m_garfield_metadata.scale_zero_reference_time_side1_ns;
  fingerprint_provenance.scale_zero_maximum_time_mismatch =
      m_garfield_metadata.scale_zero_maximum_time_mismatch;

  fingerprint_provenance.field_model = m_garfield_metadata.field_model;
  fingerprint_provenance.electric_field_map_3d = m_garfield_metadata.electric_field_map_3d;
  fingerprint_provenance.electric_field_map_3d_md5 = m_garfield_metadata.electric_field_map_3d_md5;
  fingerprint_provenance.frame_field_map = m_garfield_metadata.frame_field_map;
  fingerprint_provenance.frame_field_map_md5 = m_garfield_metadata.frame_field_map_md5;
  fingerprint_provenance.frame_field_map_3d = m_garfield_metadata.frame_field_map_3d;
  fingerprint_provenance.frame_field_map_3d_md5 = m_garfield_metadata.frame_field_map_3d_md5;
  fingerprint_provenance.frame_charge_scale = m_garfield_metadata.frame_charge_scale;
  fingerprint_provenance.cm_field_factor = m_garfield_metadata.cm_field_factor;
  fingerprint_provenance.field_cage_offsets_v = m_garfield_metadata.field_cage_offsets_v;
  fingerprint_provenance.use_ifc_voltage_distortion = m_garfield_metadata.use_ifc_voltage_distortion;
  fingerprint_provenance.use_ofc_voltage_distortion = m_garfield_metadata.use_ofc_voltage_distortion;

  const std::string reconstructed_configuration =
      PHGarfieldForwardMap::BuildCanonicalConfiguration(
          fingerprint_config,
          fingerprint_provenance, m_garfield_metadata.schema_version);
  if (canonical_configuration != reconstructed_configuration ||
      m_garfield_metadata.configuration_fingerprint !=
          hexadecimal(fnv1a64(reconstructed_configuration)))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "configuration fingerprint is not bound to the stored metadata"
        << std::endl;
    return false;
  }

  if (!m_garfield_metadata.cell_center_probes &&
      !m_allow_unprobed_garfield_lut_for_debug)
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "corner-only reach maps are disabled for production"
        << std::endl;
    return false;
  }

  if (compact)
  {
    const bool ready = initialize_garfield_compact_transport();
    if (ready) { m_garfield_metadata.schema_version = input_schema; }
    return ready;
  }
  if (m_garfield_metadata.schema_version == PHGarfieldAdaptiveMap::SchemaVersion)
  { return initialize_garfield_adaptive_transport(); }

  constexpr double axis_tolerance = 1.0e-9;
  const auto close =
      [](const double first, const double second)
  {
    const double scale =
        std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <= axis_tolerance * scale;
  };

  double observed_maximum_drift_time_ns = 0.0;
  for (unsigned int side = 0; side < 2; ++side)
  {
    const TH3* reference = m_garfield_delta_r[side];
    if (!same_node_axes(reference, m_garfield_delta_rphi[side]) ||
        !same_node_axes(reference, m_garfield_drift_time[side]) ||
        !same_node_axes(reference, m_garfield_reach_node[side]) ||
        reference->GetNbinsX() != m_garfield_metadata.phi_nodes + 2 ||
        reference->GetNbinsY() != m_garfield_metadata.radial_nodes + 2 ||
        reference->GetNbinsZ() !=
            m_garfield_metadata.longitudinal_nodes + 2)
    {
      std::cerr
          << "PHG4TpcDistortion::initialize_garfield_full_transport - "
             "inconsistent node histogram axes for "
          << suffix[side] << std::endl;
      return false;
    }

    const TH3* cell_class = m_garfield_reach_cell_class[side];
    if (cell_class->GetNbinsX() != m_garfield_metadata.phi_nodes ||
        cell_class->GetNbinsY() !=
            m_garfield_metadata.radial_nodes - 1 ||
        cell_class->GetNbinsZ() !=
            m_garfield_metadata.longitudinal_nodes - 1 ||
        !close(cell_class->GetXaxis()->GetXmin(), 0.0) ||
        !close(
            cell_class->GetXaxis()->GetXmax(),
            2.0 * std::numbers::pi) ||
        !close(
            cell_class->GetYaxis()->GetXmin(),
            m_garfield_metadata.gas_inner_radius_cm) ||
        !close(
            cell_class->GetYaxis()->GetXmax(),
            m_garfield_metadata.gas_outer_radius_cm) ||
        !cell_edges_match_node_centers(
            cell_class->GetXaxis(), reference->GetXaxis()) ||
        !cell_edges_match_node_centers(
            cell_class->GetYaxis(), reference->GetYaxis()) ||
        !cell_edges_match_node_centers(
            cell_class->GetZaxis(), reference->GetZaxis()))
    {
      std::cerr
          << "PHG4TpcDistortion::initialize_garfield_full_transport - "
             "invalid reach-cell axes for "
          << suffix[side] << std::endl;
      return false;
    }

    if (m_garfield_metadata.schema_version >= PHGarfieldForwardMap::NumericalSchemaVersion &&
        (!same_node_axes(cell_class, m_garfield_interpolation_transverse_error[side]) ||
         !same_node_axes(cell_class, m_garfield_interpolation_time_error[side])))
    { std::cerr << "Numerical-mask cell axes mismatch" << std::endl; return false; }

    const double expected_z_min =
        side == 0
            ? -m_garfield_metadata.readout_abs_z_cm
            : m_garfield_metadata.central_membrane_half_width_cm;
    const double expected_z_max =
        side == 0
            ? -m_garfield_metadata.central_membrane_half_width_cm
            : m_garfield_metadata.readout_abs_z_cm;
    if (!close(cell_class->GetZaxis()->GetXmin(), expected_z_min) ||
        !close(cell_class->GetZaxis()->GetXmax(), expected_z_max) ||
        !close(
            reference->GetYaxis()->GetBinCenter(2),
            m_garfield_metadata.gas_inner_radius_cm) ||
        !close(
            reference->GetYaxis()->GetBinCenter(
                reference->GetNbinsY() - 1),
            m_garfield_metadata.gas_outer_radius_cm) ||
        !close(reference->GetZaxis()->GetBinCenter(2), expected_z_min) ||
        !close(
            reference->GetZaxis()->GetBinCenter(
                reference->GetNbinsZ() - 1),
            expected_z_max))
    {
      std::cerr
          << "PHG4TpcDistortion::initialize_garfield_full_transport - "
             "node centers do not match physical metadata for "
          << suffix[side] << std::endl;
      return false;
    }

    if (!guard_nodes_match(m_garfield_delta_r[side]) ||
        !guard_nodes_match(m_garfield_delta_rphi[side]) ||
        !guard_nodes_match(m_garfield_drift_time[side]) ||
        !guard_nodes_match(m_garfield_reach_node[side]))
    {
      std::cerr
          << "PHG4TpcDistortion::initialize_garfield_full_transport - "
             "invalid periodic or copied guard nodes for "
          << suffix[side] << std::endl;
      return false;
    }

    for (int ix = 1; ix <= reference->GetNbinsX(); ++ix)
    {
      for (int iy = 1; iy <= reference->GetNbinsY(); ++iy)
      {
        for (int iz = 1; iz <= reference->GetNbinsZ(); ++iz)
        {
          const double reach =
              m_garfield_reach_node[side]->GetBinContent(ix, iy, iz);
          const double delta_r =
              m_garfield_delta_r[side]->GetBinContent(ix, iy, iz);
          const double delta_rphi =
              m_garfield_delta_rphi[side]->GetBinContent(ix, iy, iz);
          const double time =
              m_garfield_drift_time[side]->GetBinContent(ix, iy, iz);
          if (reach != 0.0 && reach != 1.0)
          {
            std::cerr
                << "PHG4TpcDistortion::initialize_garfield_full_transport - "
                   "node reach value is not binary"
                << std::endl;
            return false;
          }
          if (!std::isfinite(delta_r) ||
              !std::isfinite(delta_rphi) ||
              !std::isfinite(time) ||
              time < 0.0)
          {
            std::cerr
                << "PHG4TpcDistortion::initialize_garfield_full_transport - "
                   "invalid continuous node value"
                << std::endl;
            return false;
          }
          if (reach == 0.0 &&
              (delta_r != 0.0 ||
               delta_rphi != 0.0 ||
               time != 0.0))
          {
            std::cerr
                << "PHG4TpcDistortion::initialize_garfield_full_transport - "
                   "lost node has nonzero transport values"
                << std::endl;
            return false;
          }

          const bool physical_node =
              ix >= 2 && ix < reference->GetNbinsX() &&
              iy >= 2 && iy < reference->GetNbinsY() &&
              iz >= 2 && iz < reference->GetNbinsZ();
          if (reach == 1.0 && physical_node)
          {
            const double source_radius_cm =
                reference->GetYaxis()->GetBinCenter(iy);
            const double endpoint_radius_cm =
                source_radius_cm + delta_r;
            if (endpoint_radius_cm <
                    m_garfield_metadata.gas_inner_radius_cm -
                        axis_tolerance ||
                endpoint_radius_cm >
                    m_garfield_metadata.gas_outer_radius_cm +
                        axis_tolerance ||
                source_radius_cm <= 0.0 ||
                std::abs(delta_rphi / source_radius_cm) >
                    std::numbers::pi + axis_tolerance ||
                time >
                    m_garfield_metadata.maximum_drift_time_ns +
                        axis_tolerance)
            {
              std::cerr
                  << "PHG4TpcDistortion::initialize_garfield_full_transport - "
                     "physical node violates endpoint or timing semantics"
                  << std::endl;
              return false;
            }
            observed_maximum_drift_time_ns =
                std::max(observed_maximum_drift_time_ns, time);
          }
        }
      }
    }

    for (int ix = 1; ix <= cell_class->GetNbinsX(); ++ix)
    {
      const int node_x_low = ix + 1;
      const int node_x_high = node_x_low + 1;
      for (int iy = 1; iy <= cell_class->GetNbinsY(); ++iy)
      {
        const int node_y_low = iy + 1;
        const int node_y_high = node_y_low + 1;
        for (int iz = 1; iz <= cell_class->GetNbinsZ(); ++iz)
        {
          const int node_z_low = iz + 1;
          const int node_z_high = node_z_low + 1;
          const double encoded =
              cell_class->GetBinContent(ix, iy, iz);
          if (encoded != 0.0 && encoded != 1.0 && encoded != 2.0 &&
              !(encoded == 3.0 && m_garfield_metadata.schema_version >= PHGarfieldForwardMap::NumericalSchemaVersion))
          {
            std::cerr
                << "PHG4TpcDistortion::initialize_garfield_full_transport - "
                   "invalid reach cell class"
                << std::endl;
            return false;
          }

          if (m_garfield_metadata.schema_version >= PHGarfieldForwardMap::NumericalSchemaVersion)
          {
            const double transverse_error = m_garfield_interpolation_transverse_error[side]->GetBinContent(ix, iy, iz);
            const double time_error = m_garfield_interpolation_time_error[side]->GetBinContent(ix, iy, iz);
            if (!std::isfinite(transverse_error) || !std::isfinite(time_error) || transverse_error < 0. || time_error < 0. ||
                ((encoded == 2.0 || encoded == 3.0) &&
                 ((std::max(transverse_error, time_error) > m_garfield_metadata.interpolation_tolerance_um) != (encoded == 3.0))))
            { std::cerr << "Numerical-mask class disagrees with stored probe error" << std::endl; return false; }
          }
          unsigned int reaching_corners = 0;
          double minimum_corner_delta_phi =
              std::numeric_limits<double>::infinity();
          double maximum_corner_delta_phi =
              -std::numeric_limits<double>::infinity();
          for (const int xbin : {node_x_low, node_x_high})
          {
            for (const int ybin : {node_y_low, node_y_high})
            {
              for (const int zbin : {node_z_low, node_z_high})
              {
                if (m_garfield_reach_node[side]->GetBinContent(
                        xbin, ybin, zbin) == 1.0)
                {
                  ++reaching_corners;
                  const double source_radius_cm =
                      m_garfield_reach_node[side]
                          ->GetYaxis()
                          ->GetBinCenter(ybin);
                  const double delta_phi =
                      m_garfield_delta_rphi[side]->GetBinContent(
                          xbin, ybin, zbin) /
                      source_radius_cm;
                  minimum_corner_delta_phi =
                      std::min(minimum_corner_delta_phi, delta_phi);
                  maximum_corner_delta_phi =
                      std::max(maximum_corner_delta_phi, delta_phi);
                }
              }
            }
          }

          if ((encoded == 0.0 && reaching_corners != 0) ||
              ((encoded == 2.0 || encoded == 3.0) &&
               (reaching_corners != 8 ||
                maximum_corner_delta_phi -
                        minimum_corner_delta_phi >
                    std::numbers::pi)))
          {
            std::cerr
                << "PHG4TpcDistortion::initialize_garfield_full_transport - "
                   "reach cell class disagrees with its support nodes"
                << std::endl;
            return false;
          }
        }
      }
    }

    const int readout_node_bin =
        side == 0 ? 2 : reference->GetNbinsZ() - 1;
    for (int ix = 2; ix < reference->GetNbinsX(); ++ix)
    {
      for (int iy = 2; iy < reference->GetNbinsY(); ++iy)
      {
        if (m_garfield_reach_node[side]->GetBinContent(
                ix, iy, readout_node_bin) != 1.0 ||
            std::abs(m_garfield_delta_r[side]->GetBinContent(
                ix, iy, readout_node_bin)) > 1.0e-12 ||
            std::abs(m_garfield_delta_rphi[side]->GetBinContent(
                ix, iy, readout_node_bin)) > 1.0e-12 ||
            std::abs(m_garfield_drift_time[side]->GetBinContent(
                ix, iy, readout_node_bin)) > 1.0e-12)
        {
          std::cerr
              << "PHG4TpcDistortion::initialize_garfield_full_transport - "
                 "physical readout node is not an exact zero-displacement "
                 "endpoint"
              << std::endl;
          return false;
        }
      }
    }
  }

  if (!same_axis(
          m_garfield_delta_r[0]->GetXaxis(),
          m_garfield_delta_r[1]->GetXaxis()) ||
      !same_axis(
          m_garfield_delta_r[0]->GetYaxis(),
          m_garfield_delta_r[1]->GetYaxis()))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "positive/negative transverse axes differ"
        << std::endl;
    return false;
  }
  if (!close(
          observed_maximum_drift_time_ns,
          m_garfield_metadata.maximum_drift_time_ns))
  {
    std::cerr
        << "PHG4TpcDistortion::initialize_garfield_full_transport - "
           "maximum drift-time metadata does not match the node payload"
        << std::endl;
    return false;
  }

  std::cout
      << "PHG4TpcDistortion::Init - loaded Garfield full transport schema v"
      << m_garfield_metadata.schema_version
      << ", Euler step " << m_garfield_metadata.integration_step_ns
      << " ns, configuration "
      << m_garfield_metadata.configuration_fingerprint << std::endl;
  return true;
}

// Schema4 uses owning primitive arrays. The shared loader rejects partial,
// overlapping, cyclic or otherwise ambiguous partitions before any query.
bool PHG4TpcDistortion::initialize_garfield_adaptive_transport()
{
  if (!m_static_tfile || m_static_tfile->TestBit(TFile::kRecovered))
  { return false; }
  Long64_t declared_first_root = -1;
  Long64_t declared_end_root = -1;
  if (!read_parameter(m_static_tfile.get(),
          "GarfieldFullTransport/Metadata/generation_root_first", declared_first_root) ||
      !read_parameter(m_static_tfile.get(),
          "GarfieldFullTransport/Metadata/generation_root_end", declared_end_root) ||
      declared_first_root < 0 || declared_end_root < 0)
  {
    std::cerr << "PHG4TpcDistortion: missing, incorrectly typed or negative adaptive top root range" << std::endl;
    return false;
  }
  PHGarfieldForwardMap::Config configuration;
  configuration.phi_nodes = static_cast<unsigned int>(m_garfield_metadata.phi_nodes);
  configuration.radial_nodes = static_cast<unsigned int>(m_garfield_metadata.radial_nodes);
  configuration.longitudinal_nodes = static_cast<unsigned int>(m_garfield_metadata.longitudinal_nodes);
  configuration.inner_radius_cm = m_garfield_metadata.gas_inner_radius_cm;
  configuration.outer_radius_cm = m_garfield_metadata.gas_outer_radius_cm;
  configuration.central_membrane_half_width_cm = m_garfield_metadata.central_membrane_half_width_cm;
  configuration.readout_abs_z_cm = m_garfield_metadata.readout_abs_z_cm;
  configuration.nominal_phg4_drift_velocity_cm_per_ns = m_garfield_metadata.nominal_phg4_drift_velocity_cm_per_ns;
  auto adaptive = std::make_unique<AdaptiveTransportData>();
  if (!PHGarfieldAdaptiveMap::Read(
          m_static_tfile->GetDirectory("GarfieldFullTransport/Adaptive"),
          PHGarfieldAdaptiveMap::FromConfig(configuration), adaptive->data, std::cerr, true) ||
      static_cast<std::uint64_t>(declared_first_root) != adaptive->data.first_root ||
      static_cast<std::uint64_t>(declared_end_root) != adaptive->data.end_root ||
      adaptive->data.policy.guard_um != m_garfield_metadata.interpolation_tolerance_um ||
      adaptive->data.maximum_drift_time_ns != m_garfield_metadata.maximum_drift_time_ns)
  {
    std::cerr << "PHG4TpcDistortion: adaptive topology, root range, guard or duration metadata mismatch" << std::endl;
    return false;
  }
  m_garfield_metadata.adaptive_layout = PHGarfieldAdaptiveMap::SplitPolicy;
  m_garfield_metadata.adaptive_policy = PHGarfieldAdaptiveMap::BuildCanonicalPolicy(adaptive->data);
  m_garfield_metadata.adaptive_all_lost_policy = PHGarfieldAdaptiveMap::AllLostPolicy;
  m_garfield_metadata.adaptive_max_depth = static_cast<int>(adaptive->data.policy.max_depth);
  m_garfield_metadata.adaptive_max_nodes_per_root = adaptive->data.policy.max_nodes_per_root;
  m_garfield_metadata.adaptive_max_queries_per_root = adaptive->data.policy.max_queries_per_root;
  m_garfield_metadata.adaptive_point_count = adaptive->data.points.size();
  m_garfield_metadata.adaptive_cell_count = adaptive->data.cells.size();
  m_garfield_metadata.adaptive_root_count = adaptive->data.roots.size();
  m_garfieldAdaptive = std::move(adaptive);
  // No histogram, TTree, branch-address or shard-file lifetime is required by
  // subsequent adaptive queries, including direct fallback initialization.
  m_static_tfile.reset();
  std::cout << "PHG4TpcDistortion: loaded adaptive schema4, "
            << m_garfield_metadata.adaptive_cell_count << " cells, max-probe guard "
            << m_garfield_metadata.interpolation_tolerance_um << " um; configuration "
            << m_garfield_metadata.configuration_fingerprint << std::endl;
  return true;
}

bool PHG4TpcDistortion::initialize_garfield_compact_transport()
{
  if (!m_static_tfile || m_static_tfile->TestBit(TFile::kRecovered) ||
      m_compactConfigurationFingerprint != m_garfield_metadata.configuration_fingerprint)
  { return false; }
  Long64_t first = -1, end = -1;
  if (!read_parameter(m_static_tfile.get(),
          "GarfieldFullTransport/SourceMetadata/generation_root_first", first) ||
      !read_parameter(m_static_tfile.get(),
          "GarfieldFullTransport/SourceMetadata/generation_root_end", end) || first < 0 || end < first)
  { return false; }
  PHGarfieldForwardMap::Config configuration;
  configuration.phi_nodes = static_cast<unsigned int>(m_garfield_metadata.phi_nodes);
  configuration.radial_nodes = static_cast<unsigned int>(m_garfield_metadata.radial_nodes);
  configuration.longitudinal_nodes = static_cast<unsigned int>(m_garfield_metadata.longitudinal_nodes);
  configuration.inner_radius_cm = m_garfield_metadata.gas_inner_radius_cm;
  configuration.outer_radius_cm = m_garfield_metadata.gas_outer_radius_cm;
  configuration.central_membrane_half_width_cm = m_garfield_metadata.central_membrane_half_width_cm;
  configuration.readout_abs_z_cm = m_garfield_metadata.readout_abs_z_cm;
  configuration.nominal_phg4_drift_velocity_cm_per_ns = m_garfield_metadata.nominal_phg4_drift_velocity_cm_per_ns;
  const PHGarfieldCompactMap::ExpectedSource expected{
      m_compactTrainingSha256, m_compactConfigurationFingerprint, m_compactContentFingerprint};
  auto compact = std::make_unique<CompactTransportData>();
  if (!PHGarfieldCompactMap::Read(m_static_tfile->GetDirectory("GarfieldFullTransport/Runtime"),
          PHGarfieldAdaptiveMap::FromConfig(configuration), expected, compact->data,
          compact_limits(m_compactReadLimits), std::cerr, true) ||
      compact->data.first_root != static_cast<std::uint64_t>(first) ||
      compact->data.end_root != static_cast<std::uint64_t>(end) ||
      compact->data.guard_um != m_garfield_metadata.interpolation_tolerance_um ||
      compact->data.maximum_drift_time_ns != m_garfield_metadata.maximum_drift_time_ns)
  { std::cerr << "Compact source, topology, range, guard or duration mismatch" << std::endl; return false; }
  const auto& data = compact->data;
  m_garfield_metadata.adaptive_layout = PHGarfieldCompactMap::Layout;
  m_garfield_metadata.adaptive_policy = data.source.canonical_policy;
  m_garfield_metadata.adaptive_all_lost_policy = PHGarfieldAdaptiveMap::AllLostPolicy;
  m_garfield_metadata.adaptive_max_depth = static_cast<int>(data.max_depth);
  m_garfield_metadata.adaptive_max_nodes_per_root = data.max_nodes_per_root;
  m_garfield_metadata.adaptive_max_queries_per_root = data.max_queries_per_root;
  m_garfield_metadata.adaptive_point_count = data.proof.source_points;
  m_garfield_metadata.adaptive_cell_count = data.nodes.size();
  m_garfield_metadata.adaptive_root_count = data.end_root - data.first_root;
  m_garfield_metadata.compact_source_sha256 = data.source.file_sha256;
  m_garfield_metadata.compact_file_sha256 = m_compactFileSha256;
  m_garfield_metadata.compact_converter_sha256 = data.source.converter_sha256;
  m_garfield_metadata.compact_proof_sha256 = data.proof.sha256;
  m_garfield_metadata.compact_value_count = data.values.size();
  m_garfield_metadata.compact_corner_count = data.corners.size();
  m_garfield_metadata.compact_read_limits = m_compactReadLimits;
  m_garfield_metadata.compact_read_limits_explicit = m_compactReadLimitsExplicit;
  m_garfieldCompact = std::move(compact);
  m_static_tfile.reset();
  std::cout << "PHG4TpcDistortion: loaded compact schema5, " << data.nodes.size()
            << " nodes and " << data.values.size() << " corner values; source "
            << data.source.file_sha256 << std::endl;
  return true;
}

PHG4TpcDistortion::GarfieldTransportQueryStatus
PHG4TpcDistortion::query_garfield_adaptive_transport(
    const unsigned int side, const double r, const double phi, const double z,
    GarfieldTransportSample& sample) const
{
  if (!m_garfieldAdaptive && !m_garfieldCompact) { return GarfieldTransportQueryStatus::Invalid; }
  PHGarfieldCompactMap::Lookup compact_lookup;
  PHGarfieldAdaptiveMap::Cell compact_source;
  const PHGarfieldAdaptiveMap::Cell* cell = nullptr;
  if (m_garfieldCompact)
  {
    if (!PHGarfieldCompactMap::FindLeaf(m_garfieldCompact->data, side, phi, r, z, compact_lookup))
    { return GarfieldTransportQueryStatus::Invalid; }
    compact_source.id = compact_lookup.cell_id;
    compact_source.depth = compact_lookup.depth;
    compact_source.leaf_class = compact_lookup.leaf_class;
    compact_source.fallback_reason = compact_lookup.fallback_reason;
    cell = &compact_source;
  }
  else { cell = PHGarfieldAdaptiveMap::FindLeaf(m_garfieldAdaptive->data, side, phi, r, z); }
  if (!cell) { return GarfieldTransportQueryStatus::Invalid; }
  ++m_boundaryFallbackStats.adaptive_leaf_queries;
  m_boundaryFallbackStats.adaptive_max_depth_seen = std::max(
      m_boundaryFallbackStats.adaptive_max_depth_seen,
      static_cast<unsigned long long>(cell->depth));
  const auto preserve_source = [&sample, cell]()
  {
    sample.used_adaptive_lookup = true;
    sample.adaptive_cell_id = cell->id;
    sample.adaptive_depth = cell->depth;
    sample.adaptive_leaf_class = cell->leaf_class;
    sample.adaptive_fallback_reason = cell->fallback_reason;
  };
  preserve_source();
  const auto source = static_cast<PHGarfieldAdaptiveMap::LeafClass>(cell->leaf_class);
  if (source != PHGarfieldAdaptiveMap::LeafClass::Certified)
  {
    switch (source)
    {
    case PHGarfieldAdaptiveMap::LeafClass::Mixed:
      ++m_boundaryFallbackStats.adaptive_mixed_queries; break;
    case PHGarfieldAdaptiveMap::LeafClass::AllLost:
      ++m_boundaryFallbackStats.adaptive_lost_queries; break;
    case PHGarfieldAdaptiveMap::LeafClass::NumericallyUnresolved:
      ++m_boundaryFallbackStats.adaptive_numerical_queries; break;
    case PHGarfieldAdaptiveMap::LeafClass::Uncertified:
      ++m_boundaryFallbackStats.adaptive_uncertified_queries; break;
    default:
      return GarfieldTransportQueryStatus::Invalid;
    }
    // Finite probes cannot certify an all-lost continuous volume. In schema4
    // every noncertified leaf, including sampled-AllLost, requires the pinned
    // physical direct model. Without it, fail closed rather than report a loss.
    if (!m_garfieldBoundaryFallback)
    { return GarfieldTransportQueryStatus::Invalid; }
    const auto reason = static_cast<PHGarfieldAdaptiveMap::FallbackReason>(cell->fallback_reason);
    const bool numerical = reason != PHGarfieldAdaptiveMap::FallbackReason::PhysicalMixed &&
                           reason != PHGarfieldAdaptiveMap::FallbackReason::SampledAllLost;
    if (numerical) { ++m_boundaryFallbackStats.adaptive_numerical_direct_queries; }
    else { ++m_boundaryFallbackStats.adaptive_physical_direct_queries; }
    if (cell->fallback_reason <= 0 ||
        static_cast<std::size_t>(cell->fallback_reason) >= m_boundaryFallbackStats.adaptive_direct_reason_queries.size())
    { return GarfieldTransportQueryStatus::Invalid; }
    ++m_boundaryFallbackStats.adaptive_direct_reason_queries[static_cast<std::size_t>(cell->fallback_reason)];
    const auto status = query_garfield_boundary_fallback(r, phi, z, sample, numerical);
    if (status != GarfieldTransportQueryStatus::Success) { ++m_boundaryFallbackStats.adaptive_direct_failures; }
    else if (sample.reach == ReachCellClass::AllReach) { ++m_boundaryFallbackStats.adaptive_direct_reaching_queries; }
    else if (sample.reach == ReachCellClass::AllLost) { ++m_boundaryFallbackStats.adaptive_direct_physical_losses; }
    else { ++m_boundaryFallbackStats.adaptive_direct_failures; }
    // The legacy fallback resets the sample; retain the actual adaptive source
    // independently of the subsequently observed reaching/loss classification.
    preserve_source();
    return status;
  }
  PHGarfieldForwardMap::NodeValue value;
  constexpr double tolerance = 1.0e-9;
  const bool interpolated = m_garfieldCompact
      ? PHGarfieldCompactMap::Interpolate(m_garfieldCompact->data, compact_lookup, phi, r, z, value)
      : PHGarfieldAdaptiveMap::Interpolate(m_garfieldAdaptive->data, *cell, phi, r, z, value);
  if (!interpolated ||
      !value.reaches_readout || !std::isfinite(value.delta_r_cm) ||
      !std::isfinite(value.delta_rphi_cm) || !std::isfinite(value.drift_time_ns) ||
      value.drift_time_ns < 0.0 ||
      value.drift_time_ns > m_garfield_metadata.maximum_drift_time_ns + tolerance ||
      r + value.delta_r_cm < m_garfield_metadata.gas_inner_radius_cm - tolerance ||
      r + value.delta_r_cm > m_garfield_metadata.gas_outer_radius_cm + tolerance ||
      r <= 0.0 || std::abs(value.delta_rphi_cm / r) > std::numbers::pi + tolerance)
  { return GarfieldTransportQueryStatus::Invalid; }
  sample.reach = ReachCellClass::AllReach;
  sample.delta_r_cm = value.delta_r_cm;
  sample.delta_rphi_cm = value.delta_rphi_cm;
  sample.drift_time_ns = value.drift_time_ns;
  ++m_boundaryFallbackStats.adaptive_interpolated_queries;
  return GarfieldTransportQueryStatus::Success;
}

//__________________________________________________________________________________________________________
void PHG4TpcDistortion::load_event(int event_num)
{
  if (TimeTree)
  {
    int nentries = TimeTree->GetEntries();
    if (event_num > nentries)
    {
      event_num = event_num % nentries;
    }
    if (event_num % nentries == 0 && event_num != 0)
    {
      std::cout << "Distortion map sequence repeating as of event number " << event_num << std::endl;
    }
    TimeTree->GetEntry(event_num);
  }

  return;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_x_distortion_cartesian(double x, double y, double z) const
{
  double r = sqrt(x * x + y * y);
  double phi = std::atan2(y, x);

  // get components
  double dr = get_distortion('r', r, phi, z);
  double dphi = get_distortion('p', r, phi, z);

  // rotate into cartesian based on local r phi:
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double dx = dr * cosphi - dphi * sinphi;
  return dx;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_y_distortion_cartesian(double x, double y, double z) const
{
  double r = sqrt(x * x + y * y);
  double phi = std::atan2(y, x);

  // get components
  double dr = get_distortion('r', r, phi, z);
  double dphi = get_distortion('p', r, phi, z);

  // rotate into cartesian based on local r phi:
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double dy = dphi * cosphi + dr * sinphi;
  return dy;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_z_distortion_cartesian(double x, double y, double z) const
{
  double r = sqrt(x * x + y * y);
  double phi = std::atan2(y, x);

  // get components
  double dz = get_distortion('z', r, phi, z);

  return dz;
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_r_distortion(double r, double phi, double z) const
{
  return get_distortion('r', r, phi, z);
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_rphi_distortion(double r, double phi, double z) const
{
  if (m_phi_hist_in_radians)
  {  // if the hist is in radians, multiply by r to get the rphi distortion
    return r * get_distortion('p', r, phi, z);
  }

  return get_distortion('p', r, phi, z);
}

//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_z_distortion(double r, double phi, double z) const
{
  return get_distortion('z', r, phi, z);
}
//__________________________________________________________________________________________________________
double PHG4TpcDistortion::get_reaches_readout(double r, double phi, double z) const
{
  if (m_do_ReachesReadout)
  {
    return get_distortion('R', r, phi, z);
  }

  return 1;
}

//__________________________________________________________________________________________________________
PHG4TpcDistortion::GarfieldTransportQueryStatus
PHG4TpcDistortion::query_garfield_transport(
    const double r, const double phi, const double z, GarfieldTransportSample& sample) const
{
  ++m_boundaryFallbackStats.total_queries;
  const auto status = query_garfield_transport_impl(r, phi, z, sample);
  if (status == GarfieldTransportQueryStatus::OutsideDriftVolume) { ++m_boundaryFallbackStats.outside_queries; }
  else if (status == GarfieldTransportQueryStatus::Invalid) { ++m_boundaryFallbackStats.invalid_queries; }
  else if (sample.used_direct_fallback && sample.direct_fallback_reason == 1) { ++m_boundaryFallbackStats.mixed_direct_queries; }
  else if (sample.used_direct_fallback && sample.direct_fallback_reason == 2) { ++m_boundaryFallbackStats.numerical_direct_queries; }
  else if (sample.reach == ReachCellClass::AllReach) { ++m_boundaryFallbackStats.lut_reaching_queries; }
  else if (sample.reach == ReachCellClass::AllLost) { ++m_boundaryFallbackStats.lut_lost_queries; }
  else if (sample.reach == ReachCellClass::Mixed) { ++m_boundaryFallbackStats.unresolved_mixed_queries; }
  else { ++m_boundaryFallbackStats.invalid_queries; return GarfieldTransportQueryStatus::Invalid; }
  return status;
}

PHG4TpcDistortion::GarfieldTransportQueryStatus
PHG4TpcDistortion::query_garfield_transport_impl(
    const double r,
    double phi,
    double z,
    GarfieldTransportSample& sample) const
{
  sample = GarfieldTransportSample{};
  if (!m_initialized ||
      m_application_mode != ApplicationMode::GarfieldFullTransport ||
      !std::isfinite(r) ||
      !std::isfinite(phi) ||
      !std::isfinite(z))
  {
    return GarfieldTransportQueryStatus::Invalid;
  }

  constexpr double tolerance = 1.0e-9;
  const auto in_closed_range =
      [](const double value, const double minimum, const double maximum)
  {
    const double scale =
        std::max({1.0, std::abs(value), std::abs(minimum), std::abs(maximum)});
    return value >= minimum - tolerance * scale &&
           value <= maximum + tolerance * scale;
  };

  if (!in_closed_range(
          r,
          m_garfield_metadata.gas_inner_radius_cm,
          m_garfield_metadata.gas_outer_radius_cm) ||
      std::abs(z) <
          m_garfield_metadata.central_membrane_half_width_cm -
              tolerance ||
      std::abs(z) >
          m_garfield_metadata.readout_abs_z_cm + tolerance)
  {
    return GarfieldTransportQueryStatus::OutsideDriftVolume;
  }

  const double two_pi = 2.0 * std::numbers::pi;
  phi = std::fmod(phi, two_pi);
  if (phi < 0.0)
  {
    phi += two_pi;
  }
  if (phi >= two_pi)
  {
    phi = 0.0;
  }

  const unsigned int side = z > 0.0 ? 1U : 0U;
  const double r_query = std::clamp(
      r,
      m_garfield_metadata.gas_inner_radius_cm,
      m_garfield_metadata.gas_outer_radius_cm);
  const double z_min =
      side == 0
          ? -m_garfield_metadata.readout_abs_z_cm
          : m_garfield_metadata.central_membrane_half_width_cm;
  const double z_max =
      side == 0
          ? -m_garfield_metadata.central_membrane_half_width_cm
          : m_garfield_metadata.readout_abs_z_cm;
  z = std::clamp(z, z_min, z_max);

  const bool adaptive = m_garfield_metadata.schema_version == PHGarfieldAdaptiveMap::SchemaVersion ||
      m_garfield_metadata.schema_version == PHGarfieldCompactMap::SchemaVersion ||
      m_garfield_metadata.schema_version == PHGarfieldCompactMap::CollectionSchemaVersion;
  if (adaptive) { ++m_boundaryFallbackStats.adaptive_queries; }

  // Transport from the physical readout plane is exactly zero for every
  // in-range r/phi. Do not let an adjacent conservative Mixed cell turn this
  // known endpoint into a false loss.
  if (std::abs(
          std::abs(z) -
          m_garfield_metadata.readout_abs_z_cm) <= tolerance)
  {
    sample.reach = ReachCellClass::AllReach;
    if (adaptive) { ++m_boundaryFallbackStats.adaptive_readout_queries; }
    return GarfieldTransportQueryStatus::Success;
  }

  if (adaptive)
  { return query_garfield_adaptive_transport(side, r_query, phi, z, sample); }

  const TH3* cell_class = m_garfield_reach_cell_class[side];
  int phi_bin = cell_class->GetXaxis()->FindFixBin(phi);
  int radial_bin = cell_class->GetYaxis()->FindFixBin(r_query);
  int longitudinal_bin = cell_class->GetZaxis()->FindFixBin(z);
  phi_bin = std::clamp(phi_bin, 1, cell_class->GetNbinsX());
  radial_bin = std::clamp(radial_bin, 1, cell_class->GetNbinsY());
  longitudinal_bin =
      std::clamp(longitudinal_bin, 1, cell_class->GetNbinsZ());

  const double encoded_class =
      cell_class->GetBinContent(phi_bin, radial_bin, longitudinal_bin);
  if (encoded_class == 0.0)
  {
    sample.reach = ReachCellClass::AllLost;
    return GarfieldTransportQueryStatus::Success;
  }
  if (encoded_class == 1.0)
  {
    if (m_garfieldBoundaryFallback)
    {
      return query_garfield_boundary_fallback(r_query, phi, z, sample);
    }
    sample.reach = ReachCellClass::Mixed;
    return GarfieldTransportQueryStatus::Success;
  }
  if (encoded_class == 3.0 && m_garfield_metadata.schema_version >= PHGarfieldForwardMap::NumericalSchemaVersion)
  {
    if (!m_garfieldBoundaryFallback) { return GarfieldTransportQueryStatus::Invalid; }
    return query_garfield_boundary_fallback(r_query, phi, z, sample, true);
  }
  if (encoded_class != 2.0)
  {
    return GarfieldTransportQueryStatus::Invalid;
  }

  sample.reach = ReachCellClass::AllReach;
  sample.delta_r_cm =
      m_garfield_delta_r[side]->Interpolate(phi, r_query, z);
  sample.delta_rphi_cm =
      m_garfield_delta_rphi[side]->Interpolate(phi, r_query, z);
  sample.drift_time_ns =
      m_garfield_drift_time[side]->Interpolate(phi, r_query, z);

  const double endpoint_radius_cm = r_query + sample.delta_r_cm;
  if (!std::isfinite(sample.delta_r_cm) ||
      !std::isfinite(sample.delta_rphi_cm) ||
      !std::isfinite(sample.drift_time_ns) ||
      sample.drift_time_ns < 0.0 ||
      sample.drift_time_ns >
          m_garfield_metadata.maximum_drift_time_ns + tolerance ||
      endpoint_radius_cm <
          m_garfield_metadata.gas_inner_radius_cm - tolerance ||
      endpoint_radius_cm >
          m_garfield_metadata.gas_outer_radius_cm + tolerance ||
      r_query <= 0.0 ||
      std::abs(sample.delta_rphi_cm / r_query) >
          std::numbers::pi + tolerance)
  {
    sample = GarfieldTransportSample{};
    return GarfieldTransportQueryStatus::Invalid;
  }
  return GarfieldTransportQueryStatus::Success;
}

bool PHG4TpcDistortion::get_garfield_transport(
    const double r,
    const double phi,
    const double z,
    GarfieldTransportSample& sample) const
{
  return query_garfield_transport(r, phi, z, sample) ==
         GarfieldTransportQueryStatus::Success;
}

double PHG4TpcDistortion::get_distortion(char axis, double r, double phi, double z) const
{
  if (phi < 0)
  {
    phi += 2 * M_PI;
  }
  const int zpart = (z > 0 ? 1 : 0);  // z<0 corresponds to the negative side, which is element 0.

  TH3* hdistortion = nullptr;

  if (axis != 'r' && axis != 'p' && axis != 'z' && axis != 'R')
  {
    std::cout << "Distortion Requested along axis " << axis << " which is invalid.  Exiting.\n"
              << std::endl;
    exit(1);
  }

  double _distortion = 0.;

  // select the appropriate histogram:
  if (m_do_static_distortions)
  {
    if (axis == 'r')
    {
      hdistortion = hDRint[zpart];
    }
    else if (axis == 'p')
    {
      hdistortion = hDPint[zpart];
    }
    else if (axis == 'z')
    {
      hdistortion = hDZint[zpart];
    }
    else if (axis == 'R')
    {
      hdistortion = hReach[zpart];
    }
    if (hdistortion)
    {
      if (check_boundaries(hdistortion, phi, r, z))
      {
        _distortion += hdistortion->Interpolate(phi, r, z);
      }
    }
    else
    {
      std::cout << "Static Distortion Requested along axis " << axis << ", but distortion map does not exist.  Exiting.\n"
                << std::endl;
      exit(1);
    }
  }

  if (m_do_time_ordered_distortions)
  {
    if (axis == 'r')
    {
      hdistortion = TimehDR[zpart];
    }
    else if (axis == 'p')
    {
      hdistortion = TimehDP[zpart];
    }
    else if (axis == 'z')
    {
      hdistortion = TimehDZ[zpart];
    }
    else if (axis == 'R')
    {
      hdistortion = TimehRR[zpart];
    }
    if (hdistortion)
    {
      if (check_boundaries(hdistortion, phi, r, z))
      {
        _distortion += hdistortion->Interpolate(phi, r, z);
      }
    }
    else
    {
      std::cout << "Time Series Distortion Requested along axis " << axis << ", but distortion map does not exist.  Exiting.\n"
                << std::endl;
      exit(1);
    }
  }

  return _distortion;
}
