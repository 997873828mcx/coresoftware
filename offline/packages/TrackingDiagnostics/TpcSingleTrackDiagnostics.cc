#include "TpcSingleTrackDiagnostics.h"

#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>
#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4main/PHG4EventHeader.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/getClass.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <tuple>

namespace
{
  constexpr double quiet_nan = std::numeric_limits<double>::quiet_NaN();

  template <std::size_t N> void invalidate(std::array<double, N> &values) { values.fill(quiet_nan); }

  template <typename Row> void branches_identity(TTree *tree, Row &row)
  {
    tree->Branch("sample_id", &row.sample_id, "sample_id/l");
    tree->Branch("entry", &row.entry, "entry/l");
    tree->Branch("run", &row.run, "run/I");
    tree->Branch("event", &row.event, "event/I");
  }

  std::array<double, 3> transform_point(const std::array<double, 3> &point,
                                       const std::array<double, 9> &basis,
                                       const std::array<double, 3> &origin,
                                       const bool to_world)
  {
    std::array<double, 3> output{};
    for (int row = 0; row < 3; ++row)
    {
      for (int column = 0; column < 3; ++column)
      {
        output[row] += to_world ? basis[column * 3 + row] * point[column]
                                : basis[row * 3 + column] * (point[column] - origin[column]);
      }
      if (to_world) { output[row] += origin[row]; }
    }
    return output;
  }
}  // namespace

TpcSingleTrackDiagnostics::TpcSingleTrackDiagnostics(const std::string &name,
                                                     const std::string &filename)
  : SubsysReco(name), m_filename(filename)
{
}

TpcSingleTrackDiagnostics::~TpcSingleTrackDiagnostics()
{
  if (m_file) { m_file->Close(); delete m_file; }
}

int TpcSingleTrackDiagnostics::Init(PHCompositeNode *)
{
  if ((m_require_field_map && !m_config.magnetic_field) || !(m_readout_z_cm > 0.0))
  {
    std::cerr << Name() << ": explicit mapped field and positive readout z are required\n";
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_file = TFile::Open(m_filename.c_str(), "CREATE");
  if (!m_file || m_file->IsZombie()) { return Fun4AllReturnCodes::ABORTRUN; }
  create_trees();
  TNamed("schema_name", "TpcSingleTrackDiagnostics").Write();
  TParameter<int>("schema_version", 1).Write();
  TNamed("provenance", m_provenance.c_str()).Write();
  TNamed("coordinate_contract",
      "Fit state/world positions: cm; momentum GeV/c; covariance row-major (x,y,z,phi,signedInversePt,tanLambda). "
      "Physical-state signedInversePt=-charge/pT in this fitter's curvature convention; pT=1/abs(signedInversePt). "
      "Cluster x/y/z and layer cylinders: TPC envelope. envelope_basis_world rows are envelope axes in world. "
      "Residual = measurement minus prediction, local (radial,tangential,longitudinal); basis uses predicted "
      "envelope phi. Detector side is raw TpcDefs side 0/1. Path is transverse arc length cm. "
      "Poly state is stored unchanged at its reconstructed beamline PCA; Kalman PCA is an unconstrained "
      "mapped-field transverse minimum; inner state is at the declared envelope cylinder.").Write();
  TNamed("residual_contract",
      "type1: included-smoothed state propagated to geometry layer cylinder, also retained for nominal "
      "rejections; never an unbiased residual; pull_valid=0. type2: held point and all raw-hit-overlap "
      "points removed before ordering/seed/fit/robust decisions; pattern associations fixed; prediction "
      "uses geometry radius and remaining fit only. Physical LOO charge/orientation is recomputed from "
      "retained measurements with the production poly charge rule; full-poly charge is ignored. "
      "independent_fit_charge records that hypothesis; independent_charge_ambiguous flags radial-order "
      "reversals or near-collinear orientation in the remaining points and is not a hidden selection. "
      "pull covariance=R+P assumes independent measurement "
      "and prediction under linearized Gaussian model; frozen robust/global-seed covariance is approximate. "
      "Multiple cylinder intersections are counted; branch nearest geometry radial-distance hint is stored. "
      "Neighboring response correlations and track-finding/selection correlations can remain. "
      "Truth diagnostics are computed only after all event fits and never enter reconstruction.").Write();
  std::ostringstream configuration;
  configuration.precision(17);
#define CONFIG(field) configuration << #field << '=' << m_config.field << '\n'
  CONFIG(bfield_t); CONFIG(analytic_uniform_propagation); CONFIG(rkn_max_step_cm);
  CONFIG(rkn_step_tolerance); CONFIG(rkn_max_step_trials); CONFIG(rkn_max_total_steps);
  CONFIG(rkn_variable_field_within_step); CONFIG(rkn_fast_field_jacobian); CONFIG(rkn_fast_field_pca);
  CONFIG(rkn_field_pca_refine_iterations); CONFIG(meas_sigma_rphi_cm); CONFIG(meas_sigma_r_cm);
  CONFIG(meas_sigma_z_cm); CONFIG(min_measurement_sigma_cm); CONFIG(initial_sigma_pos_cm);
  CONFIG(initial_sigma_phi); CONFIG(initial_sigma_qop_t); CONFIG(use_relative_initial_sigma_qop_t);
  CONFIG(initial_sigma_qop_t_relative); CONFIG(initial_sigma_qop_t_floor); CONFIG(initial_sigma_tanl);
  CONFIG(process_sigma_pos_cm); CONFIG(process_sigma_phi); CONFIG(process_sigma_qop_t);
  CONFIG(process_sigma_tanl); CONFIG(collect_innovation_components); CONFIG(robust_measurement_updates);
  CONFIG(robust_chi2_soft); CONFIG(robust_chi2_hard); CONFIG(fit_iterations); CONFIG(material_x0_per_cm);
  CONFIG(multiple_scattering_scale); CONFIG(energy_loss_gev_per_cm); CONFIG(energy_loss_sigma_fraction);
  CONFIG(min_pt_gev);
#undef CONFIG
  configuration << "point_order=radius-descending\nregion_model=Yuri_R1_.040_.200_R2_.025_.050_R3_.025_.200_z_.100_cm\n"
                   "transition_layers=22,23,38,39\ntransition_sigma_rphi_r_z=.200,.500,.200_cm\n"
                   "exclude_transition_layers=true\nvertex_constraint=false\nV0_pairing=false\n"
                << "write_leave_one_out=" << m_write_leave_one_out << "\nwrite_truth=" << m_write_truth
                << "\nreadout_z_cm=" << m_readout_z_cm << '\n';
  TNamed("resolved_configuration", configuration.str().c_str()).Write();
  return m_file->TestBit(TFile::kWriteError) ? Fun4AllReturnCodes::ABORTRUN : Fun4AllReturnCodes::EVENT_OK;
}

void TpcSingleTrackDiagnostics::create_trees()
{
  m_file->cd();
  m_event_tree = new TTree("eventTree", "All processed event identities and diagnostic counters");
  m_track_tree = new TTree("trackTree", "Original polytracks and independent nominal mapped Kalman fits, including failures");
  m_cluster_tree = new TTree("clusterTree", "Complete original cluster catalog, including unassociated clusters and raw-hit keys");
  m_residual_tree = new TTree("residualTree", "All associated clusters, included and independent held-out residuals, including failures");
  m_truth_reference_tree = new TTree("truthReferenceTree", "Optional simulation truth surface diagnostics, evaluated after fitting");
  branches_identity(m_event_tree, m_event); branches_identity(m_track_tree, m_track);
  branches_identity(m_cluster_tree, m_cluster); branches_identity(m_residual_tree, m_residual);
  branches_identity(m_truth_reference_tree, m_truth_reference);
#define B(tree, row, field, type) tree->Branch(#field, &row.field, #field "/" type)
#define V(tree, row, field) tree->Branch(#field, &row.field)
#define A(tree, row, field, size) tree->Branch(#field, row.field.data(), #field "[" #size "]/D")
#define E(field) B(m_event_tree, m_event, field, "i")
  E(n_input_tracks); E(n_input_clusters); E(n_valid_polytracks); E(n_kalman_success); E(n_kalman_failure);
  E(n_pca_failure); E(n_missing_cluster_keys); E(n_duplicate_track_ids); E(n_duplicate_cluster_ids);
  E(n_unassociated_clusters); E(n_residual_rows); E(n_loo_success); E(n_loo_failure);
  B(m_event_tree, m_event, has_geometry, "I"); B(m_event_tree, m_event, has_crossing_decisions, "I");
  B(m_event_tree, m_event, has_truth, "I"); B(m_event_tree, m_event, inner_reference_radius_cm, "D");
#undef E
#define T(field, type) B(m_track_tree, m_track, field, type)
  T(track_id,"i"); T(source_assembled_track_id,"i"); T(nclusters,"i"); T(charge,"I"); T(poly_valid,"I");
  T(poly_fit_status,"I"); T(kalman_success,"I"); T(crossing_valid,"I"); T(crossing,"I");
  T(crossing_status,"I"); T(crossing_tier,"I"); T(crossing_score,"D"); T(poly_chi2,"D"); T(poly_ndof,"D");
  T(pca_success,"I"); T(inner_success,"I"); T(pca_path_cm,"D"); T(inner_path_cm,"D");
  T(inner_reference_radius_cm,"D"); T(inner_intersections,"i"); T(naccepted,"i"); T(nrejected,"i"); T(chi2,"D"); T(ndof,"I");
#undef T
  A(m_track_tree,m_track,poly_position,3); A(m_track_tree,m_track,poly_momentum,3); A(m_track_tree,m_track,poly_covariance,36);
  A(m_track_tree,m_track,kalman_endpoint_state,6); A(m_track_tree,m_track,kalman_pca_state,6); A(m_track_tree,m_track,kalman_inner_state,6);
  A(m_track_tree,m_track,kalman_endpoint_covariance,36); A(m_track_tree,m_track,kalman_pca_covariance,36); A(m_track_tree,m_track,kalman_inner_covariance,36);
  V(m_track_tree,m_track,fit_message); V(m_track_tree,m_track,pca_message); V(m_track_tree,m_track,inner_message);
  V(m_track_tree,m_track,cluster_ids); V(m_track_tree,m_track,measurement_used); V(m_track_tree,m_track,measurement_rejection_reason);
  V(m_track_tree,m_track,measurement_chi2_raw); V(m_track_tree,m_track,measurement_weight_scale);
#define C(field, type) B(m_cluster_tree,m_cluster,field,type)
  C(cluster_id,"l"); C(trkr_cluster_key,"l"); C(source_assembled_track_id,"i"); C(n_associated_tracks,"i");
  C(valid,"I"); C(layer,"I"); C(sector,"I"); C(side,"I"); C(geometry_valid,"I");
  C(x,"D"); C(y,"D"); C(z,"D"); C(layer_radius_cm,"D"); C(drift_length_cm,"D");
  C(adc,"D"); C(rms_x,"D"); C(rms_y,"D"); C(rms_z,"D"); C(nhits,"i"); C(phi_width,"i"); C(time_width,"i");
#undef C
  A(m_cluster_tree,m_cluster,world_position,3); A(m_cluster_tree,m_cluster,envelope_basis_world,9); A(m_cluster_tree,m_cluster,envelope_origin_world,3);
  V(m_cluster_tree,m_cluster,hitset_keys); V(m_cluster_tree,m_cluster,hit_keys);
#define R(field,type) B(m_residual_tree,m_residual,field,type)
  R(track_id,"i"); R(cluster_id,"l"); R(trkr_cluster_key,"l"); R(residual_type,"I"); R(layer,"I"); R(sector,"I");
  R(side,"I"); R(charge,"I"); R(success,"I"); R(nominal_used,"I"); R(nominal_rejection_reason,"I"); R(pull_valid,"I");
  R(independent_fit_charge,"I"); R(independent_charge_ambiguous,"I");
  R(excluded_measurements,"i"); R(excluded_shared_raw_hit_measurements,"i"); R(intersections,"i");
  R(layer_radius_cm,"D"); R(path_cm,"D"); R(tangent_angle,"D"); R(longitudinal_angle,"D");
#undef R
  A(m_residual_tree,m_residual,measurement,3); A(m_residual_tree,m_residual,residual,3); A(m_residual_tree,m_residual,pull,3);
  A(m_residual_tree,m_residual,prediction_state,6); A(m_residual_tree,m_residual,prediction_covariance,36);
  A(m_residual_tree,m_residual,measurement_covariance_local,9); A(m_residual_tree,m_residual,prediction_covariance_local,9);
  A(m_residual_tree,m_residual,residual_covariance_local,9); V(m_residual_tree,m_residual,retained_cluster_ids); V(m_residual_tree,m_residual,message);
#define Q(field,type) B(m_truth_reference_tree,m_truth_reference,field,type)
  Q(cluster_id,"l"); Q(reference_type,"I"); Q(truth_track_id,"I"); Q(pid,"I"); Q(is_primary,"I"); Q(layer,"I");
  Q(success,"I"); Q(shared_raw_hits,"i"); Q(ambiguous_raw_hits,"i"); Q(intersections,"i"); Q(interpolation_method,"I");
  Q(g4hit_key,"l"); Q(reference_radius_cm,"D"); Q(interpolation_span_cm,"D"); Q(time_ns,"D");
#undef Q
  A(m_truth_reference_tree,m_truth_reference,position,3); A(m_truth_reference_tree,m_truth_reference,momentum,3);
  A(m_truth_reference_tree,m_truth_reference,residual,3); V(m_truth_reference_tree,m_truth_reference,message);
#undef B
#undef V
#undef A
}

TpcSingleTrackDiagnostics::Identity TpcSingleTrackDiagnostics::identity() const
{
  return {m_sample_id, m_entry, m_event.run, m_event.event};
}

int TpcSingleTrackDiagnostics::process_event(PHCompositeNode *topNode)
{
  auto *clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_cluster_node);
  auto *tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_track_node);
  auto *geometry = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  auto *acts = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  auto *crossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, "TPC_CROSSING_DECISIONS");
  if (!clusters || !tracks || (m_require_field_map && !acts)) { return Fun4AllReturnCodes::ABORTRUN; }
  m_event = {};
  m_event.sample_id = m_sample_id; m_event.entry = m_entry; m_event.event = static_cast<int>(m_entry);
  if (auto *header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  { m_event.run = header->get_RunNumber(); m_event.event = header->get_EvtSequence(); }
  else if (auto *g4header = findNode::getClass<PHG4EventHeader>(topNode, "EventHeader"))
  { m_event.run = 1; m_event.event = g4header->get_EvtSequence(); }
  m_event.n_input_tracks = tracks->size(); m_event.n_input_clusters = clusters->size();
  m_event.has_geometry = geometry && acts; m_event.has_crossing_decisions = crossings != nullptr;
  double inner_radius = m_inner_radius_cm;
  if (!(inner_radius > 0.0) && geometry)
  {
    const auto range = geometry->get_begin_end();
    if (range.first != range.second && range.first->second) { inner_radius = range.first->second->get_radius(); }
  }
  m_event.inner_reference_radius_cm = inner_radius;
  std::array<double, 3> origin{};
  std::array<double, 9> basis{{1,0,0,0,1,0,0,0,1}};
  if (acts)
  {
    const auto world_origin = acts->transformTpcEnvelopeToWorld(Acts::Vector3::Zero());
    for (int coordinate = 0; coordinate < 3; ++coordinate) { origin[coordinate] = world_origin[coordinate]; }
    for (int axis = 0; axis < 3; ++axis)
    {
      Acts::Vector3 unit = Acts::Vector3::Zero(); unit[axis] = 1.0;
      const auto world_axis = acts->transformTpcEnvelopeToWorld(unit) - world_origin;
      for (int coordinate = 0; coordinate < 3; ++coordinate) { basis[axis * 3 + coordinate] = world_axis[coordinate]; }
    }
  }
  const TpcTrackVec3 cylinder_origin{origin[0],origin[1],origin[2]};
  const TpcTrackVec3 cylinder_axis{basis[6],basis[7],basis[8]};
  std::vector<ClusterRow> catalog;
  std::map<std::uint64_t, std::size_t> by_key;
  std::set<std::uint64_t> cluster_ids;
  for (unsigned int index = 0; index < clusters->size(); ++index)
  {
    const auto *cluster = clusters->get_cluster(index);
    if (!cluster) { continue; }
    ClusterRow row; static_cast<Identity &>(row) = identity();
    row.cluster_id = cluster->get_cluster_id(); row.trkr_cluster_key = cluster->get_trkr_cluster_key();
    row.source_assembled_track_id = cluster->get_source_assembled_track_id(); row.valid = cluster->isValid();
    row.x = cluster->get_centroid_x(); row.y = cluster->get_centroid_y(); row.z = cluster->get_centroid_z();
    row.world_position = transform_point({row.x,row.y,row.z}, basis, origin, true);
    row.envelope_basis_world = basis; row.envelope_origin_world = origin;
    row.side = cluster->get_side(); row.nhits = cluster->size_hits(); row.adc = cluster->get_adc();
    row.rms_x = cluster->get_rms_x(); row.rms_y = cluster->get_rms_y(); row.rms_z = cluster->get_rms_z();
    row.phi_width = cluster->get_phi_width(); row.time_width = cluster->get_time_width();
    row.drift_length_cm = m_readout_z_cm - std::abs(row.z);
    for (const auto &hit : cluster->get_hit_indices())
    { row.hitset_keys.push_back(hit.first); row.hit_keys.push_back(hit.second); }
    if (!row.hitset_keys.empty())
    { row.layer = TrkrDefs::getLayer(row.hitset_keys.front()); row.sector = TpcDefs::getSectorId(row.hitset_keys.front()); }
    if (geometry && row.layer >= 0)
    {
      const auto *layer = geometry->GetLayerCellGeom(row.layer);
      if (layer) { row.layer_radius_cm = layer->get_radius(); row.geometry_valid = std::isfinite(row.layer_radius_cm) && row.layer_radius_cm>0.0; }
    }
    if (!cluster_ids.insert(row.cluster_id).second) { ++m_event.n_duplicate_cluster_ids; }
    if (!by_key.emplace(row.trkr_cluster_key, catalog.size()).second) { ++m_event.n_duplicate_cluster_ids; }
    catalog.push_back(std::move(row));
  }
  std::set<unsigned int> track_ids;
  for (unsigned int index = 0; index < tracks->size(); ++index)
  {
    const auto *track = tracks->get_track(index);
    if (!track) { continue; }
    m_track = {}; static_cast<Identity &>(m_track) = identity();
    m_track.track_id = track->get_track_id(); m_track.source_assembled_track_id = track->get_source_assembled_track_id();
    if (!track_ids.insert(m_track.track_id).second) { ++m_event.n_duplicate_track_ids; continue; }
    m_track.poly_valid = track->isValid(); m_track.poly_fit_status = track->get_fit_status();
    m_track.charge = track->get_charge() > 0.0 ? 1 : track->get_charge() < 0.0 ? -1 : 0;
    m_track.poly_position = {track->get_x(),track->get_y(),track->get_z()};
    m_track.poly_momentum = {track->get_px(),track->get_py(),track->get_pz()};
    m_track.poly_chi2 = track->get_chi2(); m_track.poly_ndof = track->get_ndf();
    for (unsigned int row = 0; row < 6; ++row)
    { for (unsigned int column = 0; column < 6; ++column) { m_track.poly_covariance[row*6+column] = track->get_cov(row,column); } }
    invalidate(m_track.kalman_endpoint_state); invalidate(m_track.kalman_pca_state); invalidate(m_track.kalman_inner_state);
    invalidate(m_track.kalman_endpoint_covariance); invalidate(m_track.kalman_pca_covariance); invalidate(m_track.kalman_inner_covariance);
    m_track.inner_reference_radius_cm = inner_radius;
    if (const auto *decision = crossings ? crossings->get_decision(m_track.source_assembled_track_id) : nullptr)
    {
      m_track.crossing_valid = decision->isValid(); m_track.crossing = decision->get_selected_crossing();
      m_track.crossing_status = static_cast<int>(decision->get_status());
      m_track.crossing_tier = decision->get_selected_tier(); m_track.crossing_score = decision->get_selected_score();
    }
    std::vector<TpcTrackDiagnosticMeasurement> measurements;
    std::vector<std::size_t> catalog_indices;
    std::set<std::uint64_t> unique_keys;
    for (const auto key : track->get_cluster_keys())
    {
      if (!unique_keys.insert(key).second) { ++m_event.n_missing_cluster_keys; continue; }
      const auto found = by_key.find(key);
      if (found == by_key.end()) { ++m_event.n_missing_cluster_keys; continue; }
      auto &cluster = catalog[found->second]; ++cluster.n_associated_tracks;
      TpcTrackDiagnosticMeasurement measurement;
      measurement.cluster_id = cluster.cluster_id; measurement.layer_radius_cm = cluster.layer_radius_cm;
      measurement.cylinder_origin = cylinder_origin; measurement.cylinder_axis = cylinder_axis;
      measurement.point.original_index = measurements.size(); measurement.point.layer = cluster.layer;
      measurement.point.position = {cluster.world_position[0],cluster.world_position[1],cluster.world_position[2]};
      TpcTrackDiagnosticFitter::assign_yuri_measurement(measurement.point);
      for (std::size_t hit = 0; hit < cluster.hitset_keys.size(); ++hit)
      { measurement.raw_hits.emplace_back(cluster.hitset_keys[hit],cluster.hit_keys[hit]); }
      measurements.push_back(std::move(measurement)); catalog_indices.push_back(found->second);
      m_track.cluster_ids.push_back(cluster.cluster_id);
    }
    m_track.nclusters = measurements.size();
    TpcKalmanResult fit;
    if (m_track.poly_valid && m_track.poly_fit_status != 0 && m_track.charge != 0)
    {
      ++m_event.n_valid_polytracks;
      m_track.kalman_success = TpcTrackDiagnosticFitter::fit(measurements,m_track.charge,m_config,fit);
      m_track.fit_message = fit.message;
    }
    else { m_track.fit_message = "input polytrack invalid, failed, or has zero charge"; }
    if (m_track.kalman_success)
    {
      ++m_event.n_kalman_success;
      m_track.naccepted = fit.naccepted; m_track.nrejected = fit.nrejected; m_track.chi2 = fit.chi2; m_track.ndof = fit.ndof;
      TpcTrackKalmanFitter::physical_endpoint(fit,m_track.kalman_endpoint_state,m_track.kalman_endpoint_covariance);
      TpcKalmanSurfaceResult pca;
      m_track.pca_success = TpcTrackKalmanFitter::predict_to_beamline_pca(fit,{},m_config,pca);
      m_track.pca_message = pca.message;
      if (pca.success) { m_track.kalman_pca_state = pca.state; m_track.kalman_pca_covariance = pca.covariance; m_track.pca_path_cm = pca.path_cm; }
      else { ++m_event.n_pca_failure; }
      TpcKalmanSurfaceResult inner;
      m_track.inner_success = TpcTrackKalmanFitter::predict_to_cylinder(m_track.kalman_endpoint_state,
          m_track.kalman_endpoint_covariance,inner_radius,0.0,-80.0,120.0,m_config,inner,fit.mass_gev,cylinder_origin,cylinder_axis);
      m_track.inner_message = inner.message; m_track.inner_intersections = inner.intersections;
      if (inner.success) { m_track.kalman_inner_state=inner.state; m_track.kalman_inner_covariance=inner.covariance; m_track.inner_path_cm=inner.path_cm; }
    }
    else { ++m_event.n_kalman_failure; }
    // Reindex the fitter's sorted sequence using its explicit original index.
    std::map<std::size_t,std::size_t> state_for_input;
    for (std::size_t state = 0; state < fit.measurement_original_index.size(); ++state)
    { state_for_input[fit.measurement_original_index[state]] = state; }
    for (std::size_t measurement_index=0; measurement_index<measurements.size(); ++measurement_index)
    {
      const auto &measurement=measurements[measurement_index]; const auto &cluster=catalog[catalog_indices[measurement_index]];
      const auto state_found=state_for_input.find(measurement_index);
      const bool has_state=state_found!=state_for_input.end(); const std::size_t state_index=has_state?state_found->second:0;
      const bool used=has_state && fit.measurement_used[state_index];
      const int reason=has_state?fit.measurement_rejection_reason[state_index]:-1;
      m_track.measurement_used.push_back(used); m_track.measurement_rejection_reason.push_back(reason<0?255:reason);
      m_track.measurement_chi2_raw.push_back(has_state?fit.measurement_chi2_raw[state_index]:quiet_nan);
      m_track.measurement_weight_scale.push_back(has_state?fit.measurement_weight_scale[state_index]:quiet_nan);
      TpcKalmanSurfaceResult included; included.message="nominal fit unavailable";
      if (has_state)
      {
        TpcKalmanResult local; local.success=fit.success; local.charge=fit.charge;
        local.bfield_t=fit.bfield_t; local.magnetic_field=fit.magnetic_field;
        local.analytic_uniform_propagation=fit.analytic_uniform_propagation;
        local.states_smoothed={fit.states_smoothed[state_index]}; local.covs_smoothed={fit.covs_smoothed[state_index]};
        std::array<double,6> local_state{}; std::array<double,36> local_cov{};
        if (TpcTrackKalmanFitter::physical_endpoint(local,local_state,local_cov))
        { TpcTrackKalmanFitter::predict_to_cylinder(local_state,local_cov,measurement.layer_radius_cm,0.0,-5.0,5.0,m_config,included,fit.mass_gev,cylinder_origin,cylinder_axis); }
      }
      auto included_measurement=measurement;
      if (has_state)
      {
        included_measurement.point.sigma_r_cm=fit.measurement_sigma_r_used[state_index];
        included_measurement.point.sigma_rphi_cm=fit.measurement_sigma_rphi_used[state_index];
        included_measurement.point.sigma_z_cm=fit.measurement_sigma_z_used[state_index];
      }
      fill_residual(included_measurement,cluster,m_track.track_id,m_track.charge,1,used,reason,included,nullptr);
      if (m_write_leave_one_out)
      {
        TpcTrackLeaveOneOutResult loo;
        TpcTrackDiagnosticFitter::leave_one_out(measurements,measurement_index,m_track.charge,m_config,loo,
            m_track.poly_valid && m_track.poly_fit_status && m_track.charge);
        if (loo.success) { ++m_event.n_loo_success; } else { ++m_event.n_loo_failure; }
        fill_residual(measurement,cluster,m_track.track_id,m_track.charge,2,used,reason,loo.prediction,&loo);
      }
    }
    m_track_tree->Fill();
  }
  for (const auto &cluster : catalog)
  { m_cluster=cluster; m_cluster_tree->Fill(); if (!cluster.n_associated_tracks) { ++m_event.n_unassociated_clusters; } }
  // Deliberately read no truth node until all event reconstruction is finished.
  const int truth_status = m_write_truth ? fill_truth(topNode,catalog) : Fun4AllReturnCodes::EVENT_OK;
  m_event_tree->Fill(); ++m_entry;
  if (m_file->TestBit(TFile::kWriteError)) { return Fun4AllReturnCodes::ABORTRUN; }
  return truth_status;
}

void TpcSingleTrackDiagnostics::fill_residual(
    const TpcTrackDiagnosticMeasurement &measurement, const ClusterRow &cluster,
    const unsigned int track_id, const int charge, const int type, const bool used,
    const int rejection_reason, const TpcKalmanSurfaceResult &prediction,
    const TpcTrackLeaveOneOutResult *loo)
{
  m_residual={}; static_cast<Identity &>(m_residual)=identity();
  m_residual.track_id=track_id; m_residual.cluster_id=cluster.cluster_id; m_residual.trkr_cluster_key=cluster.trkr_cluster_key;
  m_residual.residual_type=type; m_residual.layer=cluster.layer; m_residual.sector=cluster.sector;
  m_residual.side=cluster.side; m_residual.charge=charge; m_residual.nominal_used=used;
  m_residual.nominal_rejection_reason=rejection_reason; m_residual.layer_radius_cm=cluster.layer_radius_cm;
  m_residual.measurement=cluster.world_position; m_residual.message=loo?loo->message:prediction.message;
  invalidate(m_residual.residual); invalidate(m_residual.pull); invalidate(m_residual.prediction_state);
  invalidate(m_residual.prediction_covariance); invalidate(m_residual.prediction_covariance_local);
  invalidate(m_residual.residual_covariance_local); m_residual.tangent_angle=quiet_nan; m_residual.longitudinal_angle=quiet_nan;
  if (loo)
  {
    m_residual.excluded_measurements=loo->excluded_measurements;
    m_residual.excluded_shared_raw_hit_measurements=loo->excluded_shared_raw_hit_measurements;
    m_residual.retained_cluster_ids=loo->retained_cluster_ids;
    m_residual.independent_fit_charge=loo->independent_charge;
    m_residual.independent_charge_ambiguous=loo->independent_charge_ambiguous;
  }
  m_residual.success=prediction.success;
  m_residual.intersections=prediction.intersections;
  if (prediction.success)
  {
    m_residual.path_cm=prediction.path_cm; m_residual.prediction_state=prediction.state;
    m_residual.prediction_covariance=prediction.covariance;
    const std::array<double,3> world_prediction{{prediction.state[0],prediction.state[1],prediction.state[2]}};
    const auto envelope_prediction=transform_point(world_prediction,cluster.envelope_basis_world,cluster.envelope_origin_world,false);
    const double phi=std::atan2(envelope_prediction[1],envelope_prediction[0]);
    const double cosine=std::cos(phi), sine=std::sin(phi);
    std::array<double,9> rotation{};
    for (int coordinate=0;coordinate<3;++coordinate)
    {
      rotation[coordinate]=cosine*cluster.envelope_basis_world[coordinate]+sine*cluster.envelope_basis_world[3+coordinate];
      rotation[3+coordinate]=-sine*cluster.envelope_basis_world[coordinate]+cosine*cluster.envelope_basis_world[3+coordinate];
      rotation[6+coordinate]=cluster.envelope_basis_world[6+coordinate];
    }
    const auto momentum=TpcTrackKalmanFitter::state_momentum(prediction.state);
    const std::array<double,3> world_momentum{{momentum.x,momentum.y,momentum.z}};
    std::array<double,3> local_momentum{};
    for (int row=0;row<3;++row)
    {
      m_residual.residual[row]=0.0;
      for (int column=0;column<3;++column)
      {
        m_residual.residual[row]+=rotation[row*3+column]*(cluster.world_position[column]-world_prediction[column]);
        local_momentum[row]+=rotation[row*3+column]*world_momentum[column];
      }
    }
    m_residual.tangent_angle=std::atan2(local_momentum[1],local_momentum[0]);
    m_residual.longitudinal_angle=std::atan2(local_momentum[2],std::hypot(local_momentum[0],local_momentum[1]));
    const double delta_phi=std::atan2(cluster.y,cluster.x)-phi;
    const double cd=std::cos(delta_phi),sd=std::sin(delta_phi);
    const double vr=measurement.point.sigma_r_cm*measurement.point.sigma_r_cm;
    const double vt=measurement.point.sigma_rphi_cm*measurement.point.sigma_rphi_cm;
    const double vz=measurement.point.sigma_z_cm*measurement.point.sigma_z_cm;
    m_residual.measurement_covariance_local={{cd*cd*vr+sd*sd*vt,cd*sd*(vr-vt),0,
                                             cd*sd*(vr-vt),sd*sd*vr+cd*cd*vt,0,0,0,vz}};
    for (int row=0;row<3;++row)
    {
      for (int column=0;column<3;++column)
      {
        double value=0.0;
        for (int i=0;i<3;++i)
        { for (int j=0;j<3;++j) { value+=rotation[row*3+i]*prediction.covariance[i*6+j]*rotation[column*3+j]; } }
        m_residual.prediction_covariance_local[row*3+column]=value;
        if (type==2)
        { m_residual.residual_covariance_local[row*3+column]=value+m_residual.measurement_covariance_local[row*3+column]; }
      }
    }
    if (type==2)
    {
      m_residual.pull_valid=1;
      for (int coordinate=0;coordinate<3;++coordinate)
      {
        const double variance=m_residual.residual_covariance_local[coordinate*3+coordinate];
        if (std::isfinite(variance) && variance>0.0)
        { m_residual.pull[coordinate]=m_residual.residual[coordinate]/std::sqrt(variance); }
        else { m_residual.pull_valid=0; }
      }
    }
  }
  else
  {
    // Measurement errors remain available even for failed predictions.
    m_residual.measurement_covariance_local={{measurement.point.sigma_r_cm*measurement.point.sigma_r_cm,0,0,
        0,measurement.point.sigma_rphi_cm*measurement.point.sigma_rphi_cm,0,0,0,measurement.point.sigma_z_cm*measurement.point.sigma_z_cm}};
  }
  m_residual_tree->Fill(); ++m_event.n_residual_rows;
}

int TpcSingleTrackDiagnostics::fill_truth(PHCompositeNode *topNode,
                                         const std::vector<ClusterRow> &catalog)
{
  auto *truth=findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  auto *g4hits=findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_TPC");
  auto *associations=findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  if (!truth || !g4hits || !associations) { return Fun4AllReturnCodes::ABORTRUN; }
  m_event.has_truth=1;
  std::map<int,std::vector<const PHG4Hit *>> by_truth;
  const auto hit_range=g4hits->getHits();
  for (auto iter=hit_range.first;iter!=hit_range.second;++iter)
  { if (iter->second) { by_truth[iter->second->get_trkid()].push_back(iter->second); } }
  struct Crossing
  {
    std::array<double,3> position{},momentum{};
    double time{0.0},span{0.0}; std::uint64_t key{0}; int method{-1};
  };
  auto crossings=[&](const int truth_id,const double radius,const ClusterRow &frame,const bool outward)
  {
    std::vector<Crossing> result;
    for (const auto *hit:by_truth[truth_id])
    {
      const std::array<double,3> world0{{hit->get_x(0),hit->get_y(0),hit->get_z(0)}};
      const std::array<double,3> world1{{hit->get_x(1),hit->get_y(1),hit->get_z(1)}};
      const auto local0=transform_point(world0,frame.envelope_basis_world,frame.envelope_origin_world,false);
      const auto local1=transform_point(world1,frame.envelope_basis_world,frame.envelope_origin_world,false);
      const double dx=local1[0]-local0[0],dy=local1[1]-local0[1],dz=local1[2]-local0[2];
      const double aa=dx*dx+dy*dy,bb=2*(local0[0]*dx+local0[1]*dy),cc=local0[0]*local0[0]+local0[1]*local0[1]-radius*radius;
      const double discriminant=bb*bb-4*aa*cc;
      if (!(aa>0.0) || discriminant<0.0 || !(radius>0.0)) { continue; }
      for (const int sign:{-1,1})
      {
        const double fraction=(-bb+sign*std::sqrt(discriminant))/(2*aa);
        if (fraction < -1e-7 || fraction > 1+1e-7) { continue; }
        const double f=std::clamp(fraction,0.0,1.0);
        if (outward && (local0[0]+f*dx)*dx+(local0[1]+f*dy)*dy<=0.0) { continue; }
        Crossing crossing;
        for (int coordinate=0;coordinate<3;++coordinate)
        { crossing.position[coordinate]=world0[coordinate]+f*(world1[coordinate]-world0[coordinate]); }
        crossing.momentum={{hit->get_px(0)+f*(hit->get_px(1)-hit->get_px(0)),
                            hit->get_py(0)+f*(hit->get_py(1)-hit->get_py(0)),
                            hit->get_pz(0)+f*(hit->get_pz(1)-hit->get_pz(0))}};
        crossing.time=hit->get_t(0)+f*(hit->get_t(1)-hit->get_t(0));
        crossing.span=std::sqrt(dx*dx+dy*dy+dz*dz); crossing.key=hit->get_hit_id();
        crossing.method=(f<1e-7 || f>1-1e-7)?0:1;
        result.push_back(crossing);
      }
    }
    std::sort(result.begin(),result.end(),[](const Crossing &a,const Crossing &b) { return std::tie(a.time,a.key)<std::tie(b.time,b.key); });
    result.erase(std::unique(result.begin(),result.end(),[](const Crossing &a,const Crossing &b) { return std::abs(a.time-b.time)<1e-7; }),result.end());
    return result;
  };
  ClusterRow default_frame;
  default_frame.envelope_basis_world={{1,0,0,0,1,0,0,0,1}};
  if (!catalog.empty()) { default_frame=catalog.front(); }
  else if (const auto *acts=findNode::getClass<ActsGeometry>(topNode,"ActsGeometry"))
  {
    const auto origin=acts->transformTpcEnvelopeToWorld(Acts::Vector3::Zero());
    for (int coordinate=0;coordinate<3;++coordinate) { default_frame.envelope_origin_world[coordinate]=origin[coordinate]; }
    for (int axis=0;axis<3;++axis)
    {
      Acts::Vector3 unit=Acts::Vector3::Zero();unit[axis]=1;
      const auto direction=acts->transformTpcEnvelopeToWorld(unit)-origin;
      for (int coordinate=0;coordinate<3;++coordinate) { default_frame.envelope_basis_world[axis*3+coordinate]=direction[coordinate]; }
    }
  }
  auto begin_row=[&]()
  {
    m_truth_reference={}; static_cast<Identity &>(m_truth_reference)=identity();
    invalidate(m_truth_reference.position); invalidate(m_truth_reference.momentum); invalidate(m_truth_reference.residual);
    m_truth_reference.interpolation_span_cm=quiet_nan; m_truth_reference.time_ns=quiet_nan;
  };
  auto fill_crossing=[&](const std::vector<Crossing> &candidates,const double time_hint)
  {
    m_truth_reference.intersections=candidates.size();
    if (candidates.empty()) { m_truth_reference.message="no recorded G4-segment crossing"; return; }
    const auto selected=std::min_element(candidates.begin(),candidates.end(),[time_hint](const Crossing &a,const Crossing &b)
        { return std::make_pair(std::abs(a.time-time_hint),a.key)<std::make_pair(std::abs(b.time-time_hint),b.key); });
    m_truth_reference.success=1; m_truth_reference.position=selected->position; m_truth_reference.momentum=selected->momentum;
    m_truth_reference.time_ns=selected->time; m_truth_reference.g4hit_key=selected->key;
    m_truth_reference.interpolation_method=selected->method; m_truth_reference.interpolation_span_cm=selected->span;
    m_truth_reference.message=selected->method?"linear interpolation within saved Geant step":"saved Geant step endpoint";
  };
  const auto primary_range=truth->GetPrimaryParticleRange();
  for (auto iter=primary_range.first;iter!=primary_range.second;++iter)
  {
    const auto *particle=iter->second;
    if (!particle || std::abs(particle->get_pid())!=211) { continue; }
    begin_row(); m_truth_reference.reference_type=1; m_truth_reference.truth_track_id=particle->get_track_id();
    m_truth_reference.pid=particle->get_pid(); m_truth_reference.is_primary=1;
    m_truth_reference.reference_radius_cm=m_event.inner_reference_radius_cm;
    const auto primary_crossings=crossings(particle->get_track_id(),m_event.inner_reference_radius_cm,default_frame,true);
    fill_crossing(primary_crossings,primary_crossings.empty()?0.0:primary_crossings.front().time);
    m_truth_reference_tree->Fill();
  }
  for (const auto &cluster:catalog)
  {
    begin_row(); m_truth_reference.reference_type=2; m_truth_reference.cluster_id=cluster.cluster_id;
    m_truth_reference.layer=cluster.layer; m_truth_reference.reference_radius_cm=cluster.layer_radius_cm;
    std::map<int,unsigned int> contributors;
    std::map<int,std::vector<double>> associated_times;
    for (std::size_t raw=0;raw<cluster.hitset_keys.size();++raw)
    {
      TrkrHitTruthAssoc::MMap links;
      associations->getG4Hits(cluster.hitset_keys[raw],cluster.hit_keys[raw],links);
      std::set<int> distinct_truth;
      std::set<PHG4HitDefs::keytype> distinct_hits;
      for (const auto &link:links) { distinct_hits.insert(link.second.second); }
      for (const auto key:distinct_hits)
      {
        const auto *hit=g4hits->findHit(key);
        if (hit) { distinct_truth.insert(hit->get_trkid()); associated_times[hit->get_trkid()].push_back(0.5*(hit->get_t(0)+hit->get_t(1))); }
      }
      if (distinct_truth.size()>1) { ++m_truth_reference.ambiguous_raw_hits; }
      for (const int id:distinct_truth) { ++contributors[id]; }
    }
    if (contributors.empty()) { m_truth_reference.message="no raw-hit truth contributor"; m_truth_reference_tree->Fill(); continue; }
    const auto best=std::max_element(contributors.begin(),contributors.end(),[](const auto &a,const auto &b)
        { return a.second<b.second || (a.second==b.second && a.first>b.first); });
    const int id=best->first; m_truth_reference.truth_track_id=id; m_truth_reference.shared_raw_hits=best->second;
    if (const auto *particle=truth->GetParticle(id)) { m_truth_reference.pid=particle->get_pid(); m_truth_reference.is_primary=truth->is_primary(particle); }
    const auto &times=associated_times[id];
    double time_hint=0.0; for (const double time:times) { time_hint+=time; } time_hint/=times.size();
    fill_crossing(crossings(id,cluster.layer_radius_cm,cluster,false),time_hint);
    if (m_truth_reference.success)
    {
      const auto predicted=transform_point(m_truth_reference.position,cluster.envelope_basis_world,cluster.envelope_origin_world,false);
      const double phi=std::atan2(predicted[1],predicted[0]),dx=cluster.x-predicted[0],dy=cluster.y-predicted[1];
      m_truth_reference.residual={{std::cos(phi)*dx+std::sin(phi)*dy,-std::sin(phi)*dx+std::cos(phi)*dy,cluster.z-predicted[2]}};
    }
    m_truth_reference_tree->Fill();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSingleTrackDiagnostics::End(PHCompositeNode *)
{
  if (!m_file) { return Fun4AllReturnCodes::EVENT_OK; }
  m_file->cd(); bool failed=false;
  for (auto *tree:{m_event_tree,m_track_tree,m_cluster_tree,m_residual_tree,m_truth_reference_tree})
  { failed=tree->Write("",TObject::kOverwrite)<=0 || failed; }
  m_file->Flush(); failed=failed || m_file->TestBit(TFile::kWriteError);
  m_file->Close(); delete m_file; m_file=nullptr;
  return failed?Fun4AllReturnCodes::ABORTRUN:Fun4AllReturnCodes::EVENT_OK;
}
