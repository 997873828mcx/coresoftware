// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_PHG4TPCTRUTHPOINTBUILDER_H
#define G4TPC_PHG4TPCTRUTHPOINTBUILDER_H

#include <fun4all/SubsysReco.h>

#include <cstddef>
#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHG4TpcGeomContainer;

class PHG4TpcTruthPointBuilder : public SubsysReco
{
 public:
  explicit PHG4TpcTruthPointBuilder(const std::string &name = "PHG4TpcTruthPointBuilder");
  ~PHG4TpcTruthPointBuilder() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_input_node(const std::string &name) { m_input_node_name = name; }
  void set_output_node(const std::string &name) { m_output_node_name = name; }
  void set_geometry_node(const std::string &name) { m_geometry_node_name = name; }
  void set_max_intersections_per_track_layer(const std::size_t n) { m_max_intersections_per_track_layer = n; }
  void set_radial_tolerance_cm(const double value) { m_radial_tolerance_cm = value; }

 private:
  struct TpcLayerRadius
  {
    unsigned int layer{0};
    double radius{0.0};
  };

  struct TruthPoint
  {
    int track_id{0};
    int shower_id{0};
    unsigned int layer{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double t{0.0};
    double px{0.0};
    double py{0.0};
    double pz{0.0};
    double path{0.0};
  };

  bool create_output_node(PHCompositeNode *topNode);
  bool load_layer_radii(PHCompositeNode *topNode);
  void add_intersections(const PHG4Hit *hit, double track_segment_start,
                         std::vector<TruthPoint> &points) const;

  std::string m_input_node_name{"G4HIT_TPC"};
  std::string m_output_node_name{"G4HIT_TPC_TRUECLUSTER"};
  std::string m_geometry_node_name{"TPCGEOMCONTAINER"};

  PHG4HitContainer *m_output_hits{nullptr};
  PHG4TpcGeomContainer *m_geom_container{nullptr};
  std::vector<TpcLayerRadius> m_tpc_layer_radii;

  std::size_t m_max_intersections_per_track_layer{0};
  double m_radial_tolerance_cm{5e-3};
};

#endif
