// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef TRACKINGDIAGNOSTICS_TPCPOLYTRACKTRUTHEVALUATOR_H
#define TRACKINGDIAGNOSTICS_TPCPOLYTRACKTRUTHEVALUATOR_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>

class PHCompositeNode;
class TFile;
class TTree;

/**
 * Write hit-level truth associations for TPC pattern-recognition tracks.
 *
 * The evaluator deliberately runs after Tpc_PolyTrackReco, while the original
 * TRKR_HITSET, TRKR_HITTRUTHASSOC, G4HIT_TPC, and G4TruthInfo nodes are still
 * present.  It joins a polytrack to its polyclusters by the cluster keys stored
 * on the track, then joins each unique raw hit to its contributing G4 hits and
 * truth particles.
 *
 * A raw hit shared by several truth particles is counted once for each distinct
 * truth track, never once per G4-hit association.  Such hits are also marked as
 * ambiguous and receive equal fractional bookkeeping weight among their
 * distinct contributors.  The fractional weight is a diagnostic only; the
 * association object does not contain per-contributor ADC fractions.
 */
class TpcPolyTrackTruthEvaluator : public SubsysReco
{
 public:
  TpcPolyTrackTruthEvaluator(
      const std::string &name = "TpcPolyTrackTruthEvaluator",
      const std::string &filename = "TpcPolyTrackTruth.root");
  ~TpcPolyTrackTruthEvaluator() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_output_file(const std::string &filename) { m_filename = filename; }
  void set_sample_id(std::uint64_t value) { m_sample_id = value; }
  // Select generated primary charged pions, including those producing no
  // TPC G4/raw hits. Raw-hit truth ranking and overlap metrics are unchanged.
  void set_primary_pion_mode(bool value) { m_primary_pion_mode = value; }
  void set_decay_chain(int parent_pdg,
                       int positive_daughter_pdg,
                       int negative_daughter_pdg)
  {
    m_parent_pdg = parent_pdg;
    m_positive_daughter_pdg = positive_daughter_pdg;
    m_negative_daughter_pdg = negative_daughter_pdg;
  }

  void set_event_header_node(const std::string &name) { m_event_header_node = name; }
  void set_truth_info_node(const std::string &name) { m_truth_info_node = name; }
  void set_g4hit_tpc_node(const std::string &name) { m_g4hit_tpc_node = name; }
  void set_hitset_node(const std::string &name) { m_hitset_node = name; }
  void set_hit_truth_assoc_node(const std::string &name) { m_hit_truth_assoc_node = name; }
  void set_poly_cluster_node(const std::string &name) { m_poly_cluster_node = name; }
  void set_poly_track_node(const std::string &name) { m_poly_track_node = name; }

 private:
  struct EventRow
  {
    unsigned int schema_version{1};
    std::uint64_t sample_id{0};
    std::uint64_t entry{0};
    int run{0};
    int event{0};
    int has_event_header{0};
    unsigned int n_truth_particles{0};
    unsigned int n_ionizing_g4hits{0};
    unsigned int n_raw_hits{0};
    unsigned int n_truth_associated_raw_hits{0};
    unsigned int n_ambiguous_raw_hits{0};
    unsigned int n_polytracks{0};
    unsigned int n_valid_polytracks{0};
    unsigned int n_generated_primary_particles{0};
    unsigned int n_generated_primary_pions{0};
    int primary_pion_truth_track_id{0};
    unsigned int n_primary_parent{0};
    unsigned int n_direct_positive_daughter{0};
    unsigned int n_direct_negative_daughter{0};
    int parent_truth_track_id{0};
    int positive_daughter_truth_track_id{0};
    int negative_daughter_truth_track_id{0};
    int positive_daughter_overlap_poly_track_id{-1};
    int negative_daughter_overlap_poly_track_id{-1};
    int positive_daughter_has_overlap_assignment{0};
    int negative_daughter_has_overlap_assignment{0};
    // Legacy Kshort aliases retained for existing analysis files.
    unsigned int n_primary_kshort{0};
    unsigned int n_direct_piplus{0};
    unsigned int n_direct_piminus{0};
    unsigned int n_missing_cluster_keys{0};
    unsigned int n_cluster_source_mismatches{0};
    unsigned int n_missing_raw_hit_references{0};
    unsigned int n_missing_g4hits{0};
    unsigned int n_missing_truth_particles{0};
    int kshort_truth_track_id{0};
    int piplus_truth_track_id{0};
    int piminus_truth_track_id{0};
    int piplus_overlap_poly_track_id{-1};
    int piminus_overlap_poly_track_id{-1};
    int piplus_has_overlap_assignment{0};
    int piminus_has_overlap_assignment{0};
    int both_daughters_present{0};
    int both_daughters_have_g4hits{0};
    int both_daughters_have_digitized_hits{0};
    int both_daughters_overlap_assigned_distinct{0};
  };

  struct TruthRow
  {
    std::uint64_t sample_id{0};
    std::uint64_t entry{0};
    int run{0};
    int event{0};
    int truth_track_id{0};
    int pid{0};
    int parent_id{0};
    int parent_pid{0};
    int primary_id{0};
    int vertex_id{0};
    int barcode{0};
    int embed_id{0};
    int is_primary{0};
    int is_selected_primary_pion{0};
    int is_selected_parent{0};
    int is_selected_kshort{0};
    int direct_daughter_role{0};
    float px{0.F};
    float py{0.F};
    float pz{0.F};
    float e{0.F};
    float pt{0.F};
    float eta{0.F};
    float phi{0.F};
    float rapidity{0.F};
    float vx{0.F};
    float vy{0.F};
    float vz{0.F};
    float vt{0.F};
    float vertex_radius{0.F};
    unsigned int n_ionizing_g4hits{0};
    unsigned int n_ionizing_g4hit_layers{0};
    int first_ionizing_g4hit_layer{-1};
    int last_ionizing_g4hit_layer{-1};
    unsigned int n_available_raw_hits{0};
    unsigned int n_available_layers{0};
    int first_layer{-1};
    int last_layer{-1};
    int overlap_assigned_poly_track_id{-1};
  };

  struct TrackRow
  {
    std::uint64_t sample_id{0};
    std::uint64_t entry{0};
    int run{0};
    int event{0};
    unsigned int poly_track_id{0};
    unsigned int source_assembled_track_id{0};
    int fit_status{0};
    int charge{0};
    unsigned int nclusters{0};
    float px{0.F};
    float py{0.F};
    float pz{0.F};
    float pt{0.F};
    float eta{0.F};
    float phi{0.F};
    unsigned int n_raw_hits{0};
    unsigned int n_layers{0};
    unsigned int n_associated_raw_hits{0};
    unsigned int n_unassociated_raw_hits{0};
    unsigned int n_ambiguous_raw_hits{0};
    unsigned int n_truth_contributors{0};
    unsigned int n_missing_cluster_keys{0};
    unsigned int n_cluster_source_mismatches{0};
    unsigned int n_missing_raw_hit_references{0};
    int best_truth_track_id{0};
    int best_pid{0};
    int best_parent_id{0};
    int best_parent_pid{0};
    unsigned int best_shared_raw_hits{0};
    unsigned int best_shared_layers{0};
    unsigned int best_exclusive_raw_hits{0};
    float best_fractional_shared_weight{0.F};
    float best_hit_overlap_fraction{0.F};
    float best_associated_hit_overlap_fraction{0.F};
    float best_fractional_purity{0.F};
    float best_completeness{0.F};
    float best_layer_completeness{0.F};
    int best_direct_daughter_role{0};
    int overlap_assigned_daughter_role{0};
  };

  struct MatchRow
  {
    std::uint64_t sample_id{0};
    std::uint64_t entry{0};
    int run{0};
    int event{0};
    unsigned int poly_track_id{0};
    unsigned int source_assembled_track_id{0};
    int truth_track_id{0};
    int pid{0};
    int parent_id{0};
    int parent_pid{0};
    unsigned int reco_raw_hits{0};
    unsigned int reco_associated_raw_hits{0};
    unsigned int truth_available_raw_hits{0};
    unsigned int truth_available_layers{0};
    unsigned int shared_raw_hits{0};
    unsigned int shared_layers{0};
    unsigned int exclusive_raw_hits{0};
    float fractional_shared_weight{0.F};
    float hit_overlap_fraction{0.F};
    float associated_hit_overlap_fraction{0.F};
    float fractional_purity{0.F};
    float completeness{0.F};
    float layer_completeness{0.F};
    int is_best_match{0};
    int direct_daughter_role{0};
    int overlap_assigned_daughter_role{0};
  };

  void create_trees();
  void reset_rows();

  std::string m_filename;
  std::string m_event_header_node{"EventHeader"};
  std::string m_truth_info_node{"G4TruthInfo"};
  std::string m_g4hit_tpc_node{"G4HIT_TPC"};
  std::string m_hitset_node{"TRKR_HITSET"};
  std::string m_hit_truth_assoc_node{"TRKR_HITTRUTHASSOC"};
  std::string m_poly_cluster_node{"TPC_POLYCLUSTERS"};
  std::string m_poly_track_node{"TPC_POLYTRACKS"};

  std::uint64_t m_sample_id{0};
  std::uint64_t m_entry{0};
  int m_parent_pdg{310};
  int m_positive_daughter_pdg{211};
  int m_negative_daughter_pdg{-211};
  bool m_primary_pion_mode{false};
  TFile *m_file{nullptr};
  TTree *m_event_tree{nullptr};
  TTree *m_truth_tree{nullptr};
  TTree *m_track_tree{nullptr};
  TTree *m_match_tree{nullptr};

  EventRow m_event;
  TruthRow m_truth;
  TrackRow m_track;
  MatchRow m_match;
};

#endif
