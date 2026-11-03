// $Id: $

/*!
 * \file PHG4TPCDistortion.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>, Henry Klest <henry.klest@stonybrook.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>, Ross Corliss <ross.corliss@stonybrook.edu>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4TPC_PHG4TPCDISTORTION_H
#define G4TPC_PHG4TPCDISTORTION_H
#define PHG4TPCDISTORTION_GARFIELD_BOUNDARY_FALLBACK 1
#define PHG4TPCDISTORTION_GARFIELD_ADAPTIVE_LUT 1
#define PHG4TPCDISTORTION_GARFIELD_COMPACT_LUT 1
#define PHG4TPCDISTORTION_GARFIELD_COMPACT_LIMITS 1

#include <array>
#include <cstdint>
#include <memory>
#include <string>

class TFile;
class TH3;
class PHGarfield;
class TTree;

//! handle distortions (static and time-ordered)
class PHG4TpcDistortion
{
 public:
  enum class ApplicationMode
  {
    LegacyDisplacement,
    GarfieldFullTransport
  };

  enum class ReachCellClass : std::uint8_t
  {
    AllLost = 0,
    Mixed = 1,
    AllReach = 2,
    NumericallyUnresolved = 3,
    Invalid = 255
  };

  enum class GarfieldTransportQueryStatus : std::uint8_t
  {
    Success = 0,
    OutsideDriftVolume = 1,
    Invalid = 255
  };

  struct GarfieldTransportSample
  {
    ReachCellClass reach{ReachCellClass::Invalid};
    bool used_direct_fallback{false};
    unsigned int direct_fallback_reason{0};  // 0 none,1 physical Mixed,2 numerical mask
    bool used_adaptive_lookup{false};
    std::uint64_t adaptive_cell_id{0};  // meaningful only when used_adaptive_lookup
    unsigned int adaptive_depth{0};
    // Schema4 source classifications are PHGarfieldAdaptiveMap enums, not
    // ReachCellClass. They survive a direct query that changes reach/loss.
    int adaptive_leaf_class{-1};
    int adaptive_fallback_reason{-1};
    double delta_r_cm{0.0};
    double delta_rphi_cm{0.0};
    double drift_time_ns{0.0};
  };

  struct GarfieldCompactReadLimits
  {
    std::uint64_t max_source_points{0};
    std::uint64_t max_nodes{0};
    std::uint64_t max_corners{0};
    std::uint64_t max_values{0};
    std::uint64_t max_bytes{0};
  };

  struct GarfieldLutMetadata
  {
    double interpolation_tolerance_um{0.0};
    std::string interpolation_error_probe_scheme;
    std::string adaptive_layout;
    std::string adaptive_policy;
    std::string adaptive_all_lost_policy;
    int adaptive_max_depth{0};
    unsigned int adaptive_max_nodes_per_root{0};
    unsigned int adaptive_max_queries_per_root{0};
    unsigned long long adaptive_point_count{0};
    unsigned long long adaptive_cell_count{0};
    unsigned long long adaptive_root_count{0};
    std::string compact_source_sha256;
    std::string compact_file_sha256;
    std::string compact_converter_sha256;
    std::string compact_proof_sha256;
    unsigned long long compact_value_count{0};
    unsigned long long compact_corner_count{0};
    GarfieldCompactReadLimits compact_read_limits;
    bool compact_read_limits_explicit{false};
    std::string field_model;
    std::array<std::string, 2> electric_field_map_3d{};
    std::array<std::string, 2> electric_field_map_3d_md5{};
    std::string frame_field_map;
    std::string frame_field_map_md5;
    std::array<std::string, 2> frame_field_map_3d{};
    std::array<std::string, 2> frame_field_map_3d_md5{};
    std::array<double, 2> frame_charge_scale{{0.0, 0.0}};
    std::array<double, 2> cm_field_factor{{1.0, 1.0}};
    std::array<double, 4> field_cage_offsets_v{{0.0, 0.0, 0.0, 0.0}};
    bool use_ifc_voltage_distortion{false};
    bool use_ofc_voltage_distortion{false};

    int schema_version{0};
    int phi_nodes{0};
    int radial_nodes{0};
    int longitudinal_nodes{0};
    double gas_inner_radius_cm{0.0};
    double gas_outer_radius_cm{0.0};
    double central_membrane_half_width_cm{0.0};
    double readout_abs_z_cm{0.0};
    double integration_step_ns{0.0};
    double nominal_phg4_drift_velocity_cm_per_ns{0.0};
    double maximum_nominal_time_mismatch{0.0};
    double scale_zero_reference_time_side0_ns{0.0};
    double scale_zero_reference_time_side1_ns{0.0};
    double scale_zero_maximum_time_mismatch{0.0};
    double maximum_drift_time_ns{0.0};
    double tpc_place_x_cm{0.0};
    double tpc_place_y_cm{0.0};
    double tpc_place_z_cm{0.0};
    double tpc_rotation_x_rad{0.0};
    double tpc_rotation_y_rad{0.0};
    double tpc_rotation_z_rad{0.0};
    double nominal_electric_field_vcm{0.0};
    double space_charge_scale_side0{0.0};
    double space_charge_scale_side1{0.0};
    double magnetic_field_scale{0.0};
    int maximum_integration_steps{0};
    bool zero_magnetic_field{false};
    bool cell_center_probes{false};
    unsigned long long cdb_timestamp{0};
    std::string integrator;
    std::string electric_field_map;
    std::string magnetic_field_map;
    std::string gas_file;
    std::string electric_field_map_md5;
    std::string magnetic_field_map_md5;
    std::string gas_file_md5;
    std::string cdb_global_tag;
    std::string configuration_label;
    std::string reach_probe_scheme;
    std::string magnetic_field_domain_policy;
    std::string nominal_time_reference_policy;
    std::string configuration_fingerprint;
  };

  //! constructor
  PHG4TpcDistortion();
  ~PHG4TpcDistortion();

  //!@name accessors
  //@{

  //! x distortion for a given truth location of the primary ionization
  double get_x_distortion_cartesian(double x, double y, double z) const;

  //! y distortion for a given truth location of the primary ionization
  double get_y_distortion_cartesian(double x, double y, double z) const;

  //! z distortion for a given truth location of the primary ionization
  double get_z_distortion_cartesian(double x, double y, double z) const;

  //! radial distortion for a given cylindrical truth location of the primary ionization
  double get_r_distortion(double r, double phi, double z) const;

  //! R*phi (hence the unitful phi-hat) distortion for a given cylindrical truth location of the primary ionization
  double get_rphi_distortion(double r, double phi, double z) const;

  //! z distortion for a given cylindrical truth location of the primary ionization
  double get_z_distortion(double r, double phi, double z) const;

  // The ReachesReadout serves as a fourth axis in the distortion histogram
  double get_reaches_readout(double r, double phi, double z) const;

  //! Atomically sample a full deterministic Garfield transport LUT while
  //! distinguishing physical source coordinates outside the drift volume from
  //! invalid configuration or LUT contents.
  GarfieldTransportQueryStatus query_garfield_transport(
      double r,
      double phi,
      double z,
      GarfieldTransportSample &sample) const;

  //! Backward-compatible boolean query. False combines an outside-volume
  //! source with an invalid query; callers that need to distinguish them must
  //! use query_garfield_transport.
  bool get_garfield_transport(double r,
                              double phi,
                              double z,
                              GarfieldTransportSample &sample) const;

  ApplicationMode application_mode() const
  {
    return m_application_mode;
  }

  const GarfieldLutMetadata &garfield_lut_metadata() const
  {
    return m_garfield_metadata;
  }

  bool is_initialized() const
  {
    return m_initialized;
  }

  //! Gets the verbosity of this module.
  int Verbosity() const
  {
    return verbosity;
  }

  //@}

  //!@name modifiers
  //@{

  //! enable distortions from a single fixed source
  void set_do_static_distortions(bool value)
  {
    m_do_static_distortions = value;
  }

  //! set the filename for the single fixed distortion
  void set_static_distortion_filename(const std::string &value)
  {
    m_static_distortion_filename = value;
  }

  //! enable time ordered distortions
  void set_do_time_ordered_distortions(bool value)
  {
    m_do_time_ordered_distortions = value;
  }

  //! time ordered distortion filename
  void set_time_ordered_distortion_filename(const std::string &value)
  {
    m_time_ordered_distortion_filename = value;
  }

  //! enable reaches readout
  void set_do_ReachesReadout(bool value)
  {
    m_do_ReachesReadout = value;
  }

  void set_read_phi_as_radians(bool flag)
  {
    m_phi_hist_in_radians = flag;
  }

  void set_application_mode(ApplicationMode mode)
  {
    m_application_mode = mode;
  }

  struct BoundaryFallbackStats
  {
    unsigned long long total_queries{0};
    // Adaptive in-volume queries partition into exact-readout and located
    // leaves. A Certified interpolation failure is additionally Invalid; no
    // direct/error result is ever relabeled as interpolated.
    unsigned long long adaptive_queries{0};
    unsigned long long adaptive_readout_queries{0};
    unsigned long long adaptive_leaf_queries{0};
    unsigned long long adaptive_interpolated_queries{0};
    unsigned long long adaptive_lost_queries{0};  // sampled-AllLost source leaf, still queried directly
    unsigned long long adaptive_mixed_queries{0};
    unsigned long long adaptive_uncertified_queries{0};
    unsigned long long adaptive_numerical_queries{0};
    unsigned long long adaptive_max_depth_seen{0};
    unsigned long long adaptive_physical_direct_queries{0};
    unsigned long long adaptive_numerical_direct_queries{0};
    // Indices are PHGarfieldAdaptiveMap::FallbackReason. None(index0) stays0.
    std::array<unsigned long long, 8> adaptive_direct_reason_queries{};
    unsigned long long adaptive_direct_reaching_queries{0};
    unsigned long long adaptive_direct_physical_losses{0};
    unsigned long long adaptive_direct_failures{0};
    unsigned long long lut_reaching_queries{0};
    unsigned long long lut_lost_queries{0};
    unsigned long long mixed_direct_queries{0};
    unsigned long long numerical_direct_queries{0};
    unsigned long long unresolved_mixed_queries{0};
    unsigned long long outside_queries{0};
    unsigned long long invalid_queries{0};
    unsigned long long calls{0};
    unsigned long long mixed_calls{0};
    unsigned long long numerical_calls{0};
    unsigned long long reaching{0};
    unsigned long long physical_losses{0};
    unsigned long long failures{0};
  };
  //! Opt in before Init: schema-2 Mixed cells, schema-3 numerical masks and
  //! schema-4 noncertified leaves use the same pinned direct model. Schema-4
  //! finite-probe AllLost leaves are also queried directly. Defaults stay disabled.
  void set_garfield_boundary_fallback(bool value) { m_garfieldBoundaryFallback = value; }
  bool garfield_boundary_fallback_enabled() const { return m_garfieldBoundaryFallback; }
  const BoundaryFallbackStats& boundary_fallback_stats() const { return m_boundaryFallbackStats; }

  //! Schema5 requires the externally approved converter manifest, including
  //! the runtime file hash. Self-reported training hashes are insufficient.
  void set_garfield_compact_source_binding(const std::string& training_sha256,
      const std::string& configuration_fingerprint, const std::string& content_fingerprint,
      const std::string& compact_file_sha256)
  {
    m_compactTrainingSha256 = training_sha256;
    m_compactConfigurationFingerprint = configuration_fingerprint;
    m_compactContentFingerprint = content_fingerprint;
    m_compactFileSha256 = compact_file_sha256;
  }

  //! Optional, explicit admission limits before Init. These bound allocations;
  //! they neither preallocate storage nor change the transport/lookup policy.
  //! Omitting this setter retains PHGarfieldCompactMap's default limits.
  void set_garfield_compact_read_limits(std::uint64_t max_source_points,
      std::uint64_t max_nodes, std::uint64_t max_corners,
      std::uint64_t max_values, std::uint64_t max_bytes);
  const GarfieldCompactReadLimits& garfield_compact_read_limits() const { return m_compactReadLimits; }
  bool garfield_compact_read_limits_explicit() const { return m_compactReadLimitsExplicit; }

  //! Permit corner-only reach maps for dedicated debugging. Production keeps
  //! this false and requires the generator's interior cell probes.
  void set_allow_unprobed_garfield_lut_for_debug(bool value)
  {
    m_allow_unprobed_garfield_lut_for_debug = value;
  }

  //! initialize
  void Init();

  //! get relevant histogram from time ordered TTrees
  void load_event(int event_num);

  //! Sets the verbosity of this module (0 by default=quiet).
  void Verbosity(const int ival)
  {
    verbosity = ival;
  }

  //@}

 private:
  //! get distortion for a set of histogram and an input momentum distribution
  double get_distortion(char axis, double r, double phi, double z) const;
  bool initialize_garfield_full_transport();
  bool initialize_garfield_adaptive_transport();
  bool initialize_garfield_compact_transport();
  GarfieldTransportQueryStatus query_garfield_adaptive_transport(unsigned int side, double r, double phi, double z, GarfieldTransportSample& sample) const;
  GarfieldTransportQueryStatus query_garfield_transport_impl(double r, double phi, double z, GarfieldTransportSample& sample) const;
  bool initialize_garfield_boundary_fallback();
  GarfieldTransportQueryStatus query_garfield_boundary_fallback(double r, double phi, double z, GarfieldTransportSample& sample, bool numerical = false) const;

  //! The verbosity level. 0 means not verbose at all.
  int verbosity = 0;
  bool m_initialized = false;
  ApplicationMode m_application_mode{ApplicationMode::LegacyDisplacement};
  bool m_allow_unprobed_garfield_lut_for_debug{false};

  //! Flag controls whether to assume the phi hist units are radians or cm.
  bool m_phi_hist_in_radians = true;

  bool m_do_ReachesReadout = false;

  //!@name static histograms
  //@{
  bool m_do_static_distortions = false;
  std::string m_static_distortion_filename;
  std::unique_ptr<TFile> m_static_tfile;
  TH3 *hDRint[2] = {nullptr, nullptr};
  TH3 *hDPint[2] = {nullptr, nullptr};
  TH3 *hDZint[2] = {nullptr, nullptr};
  TH3 *hReach[2] = {nullptr, nullptr};
  //@}

  //!@name full deterministic Garfield transport
  //@{
  std::array<TH3 *, 2> m_garfield_delta_r{{nullptr, nullptr}};
  std::array<TH3 *, 2> m_garfield_delta_rphi{{nullptr, nullptr}};
  std::array<TH3 *, 2> m_garfield_drift_time{{nullptr, nullptr}};
  std::array<TH3 *, 2> m_garfield_reach_node{{nullptr, nullptr}};
  std::array<TH3 *, 2> m_garfield_reach_cell_class{{nullptr, nullptr}};
  std::array<TH3 *, 2> m_garfield_interpolation_transverse_error{{nullptr, nullptr}};
  std::array<TH3 *, 2> m_garfield_interpolation_time_error{{nullptr, nullptr}};
  struct AdaptiveTransportData;
  std::unique_ptr<AdaptiveTransportData> m_garfieldAdaptive;
  struct CompactTransportData;
  std::unique_ptr<CompactTransportData> m_garfieldCompact;
  std::string m_compactTrainingSha256;
  std::string m_compactConfigurationFingerprint;
  std::string m_compactContentFingerprint;
  std::string m_compactFileSha256;
  GarfieldCompactReadLimits m_compactReadLimits;
  bool m_compactReadLimitsExplicit{false};
  GarfieldLutMetadata m_garfield_metadata;
  bool m_garfieldBoundaryFallback{false};
  std::unique_ptr<PHGarfield> m_boundaryGarfield;
  mutable BoundaryFallbackStats m_boundaryFallbackStats;
  //@}

  //!@name time ordered histograms
  //@{
  bool m_do_time_ordered_distortions = false;
  std::string m_time_ordered_distortion_filename;
  std::unique_ptr<TFile> m_time_ordered_tfile;
  TTree *TimeTree = nullptr;
  TH3 *TimehDR[2] = {nullptr, nullptr};
  TH3 *TimehDP[2] = {nullptr, nullptr};
  TH3 *TimehDZ[2] = {nullptr, nullptr};
  TH3 *TimehRR[2] = {nullptr, nullptr};
  //@}
};

#endif /* G4TPC_PHG4TPCDISTORTION_H */
