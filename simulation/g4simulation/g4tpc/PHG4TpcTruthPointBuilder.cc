#include "PHG4TpcTruthPointBuilder.h"

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <utility>

namespace
{
  template <class T>
  constexpr T square(const T &x)
  {
    return x * x;
  }
}  // namespace

PHG4TpcTruthPointBuilder::PHG4TpcTruthPointBuilder(const std::string &name)
  : SubsysReco(name)
{
}

int PHG4TpcTruthPointBuilder::InitRun(PHCompositeNode *topNode)
{
  if (!load_layer_radii(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!create_output_node(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": loaded " << m_tpc_layer_radii.size()
              << " TPC layer radii from " << m_geometry_node_name
              << ", writing truth layer points to " << m_output_node_name
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcTruthPointBuilder::process_event(PHCompositeNode *topNode)
{
  if (!m_output_hits && !create_output_node(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_output_hits->Reset();

  auto *input_hits = findNode::getClass<PHG4HitContainer>(topNode, m_input_node_name);
  if (!input_hits)
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << Name() << ": missing input node " << m_input_node_name << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::map<int, double> track_path_offset;
  std::map<std::pair<int, unsigned int>, std::vector<TruthPoint>> points_by_track_layer;

  PHG4HitContainer::ConstRange range = input_hits->getHits();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    const PHG4Hit *hit = iter->second;
    if (!hit)
    {
      continue;
    }

    const double dx = hit->get_x(1) - hit->get_x(0);
    const double dy = hit->get_y(1) - hit->get_y(0);
    const double dz = hit->get_z(1) - hit->get_z(0);
    const double step_length = std::sqrt(square(dx) + square(dy) + square(dz));
    if (step_length <= 0.0 || !std::isfinite(step_length))
    {
      continue;
    }

    const int track_id = hit->get_trkid();
    const double track_segment_start = track_path_offset[track_id];

    std::vector<TruthPoint> hit_points;
    add_intersections(hit, track_segment_start, hit_points);
    for (const auto &point : hit_points)
    {
      points_by_track_layer[{point.track_id, point.layer}].push_back(point);
    }

    track_path_offset[track_id] = track_segment_start + step_length;
  }

  std::size_t n_written = 0;
  for (auto &entry : points_by_track_layer)
  {
    auto &points = entry.second;
    std::sort(points.begin(), points.end(),
              [](const TruthPoint &lhs, const TruthPoint &rhs)
              { return lhs.path < rhs.path; });

    const std::size_t n_to_write =
        (m_max_intersections_per_track_layer == 0) ? points.size() : std::min(points.size(), m_max_intersections_per_track_layer);

    for (std::size_t i = 0; i < n_to_write; ++i)
    {
      const auto &point = points[i];

      auto *truth_hit = new PHG4Hitv1();
      truth_hit->set_trkid(point.track_id);
      truth_hit->set_shower_id(point.shower_id);
      truth_hit->set_layer(point.layer);
      truth_hit->set_x(0, static_cast<float>(point.x));
      truth_hit->set_y(0, static_cast<float>(point.y));
      truth_hit->set_z(0, static_cast<float>(point.z));
      truth_hit->set_t(0, static_cast<float>(point.t));
      truth_hit->set_x(1, static_cast<float>(point.x));
      truth_hit->set_y(1, static_cast<float>(point.y));
      truth_hit->set_z(1, static_cast<float>(point.z));
      truth_hit->set_t(1, static_cast<float>(point.t));
      truth_hit->set_px(0, static_cast<float>(point.px));
      truth_hit->set_py(0, static_cast<float>(point.py));
      truth_hit->set_pz(0, static_cast<float>(point.pz));
      truth_hit->set_px(1, static_cast<float>(point.px));
      truth_hit->set_py(1, static_cast<float>(point.py));
      truth_hit->set_pz(1, static_cast<float>(point.pz));
      truth_hit->set_path_length(static_cast<float>(point.path));

      m_output_hits->AddHit(point.layer, truth_hit);
      ++n_written;
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << Name() << ": wrote " << n_written
              << " truth layer points to " << m_output_node_name
              << " from " << input_hits->size() << " input hits" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHG4TpcTruthPointBuilder::create_output_node(PHCompositeNode *topNode)
{
  m_output_hits = findNode::getClass<PHG4HitContainer>(topNode, m_output_node_name);
  if (m_output_hits)
  {
    return true;
  }

  PHNodeIterator iter(topNode);
  auto *dst_node = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cout << PHWHERE << Name() << ": DST node missing" << std::endl;
    return false;
  }

  m_output_hits = new PHG4HitContainer(m_output_node_name);
  dst_node->addNode(new PHIODataNode<PHObject>(m_output_hits, m_output_node_name, "PHObject"));
  return true;
}

bool PHG4TpcTruthPointBuilder::load_layer_radii(PHCompositeNode *topNode)
{
  m_geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, m_geometry_node_name);
  if (!m_geom_container)
  {
    std::cout << PHWHERE << Name() << ": missing geometry node " << m_geometry_node_name << std::endl;
    return false;
  }

  m_tpc_layer_radii.clear();
  const auto layer_range = m_geom_container->get_begin_end();
  for (auto layer_iter = layer_range.first; layer_iter != layer_range.second; ++layer_iter)
  {
    const auto *layer_geom = layer_iter->second;
    if (!layer_geom)
    {
      continue;
    }

    const double radius = layer_geom->get_radius();
    if (!std::isfinite(radius) || radius <= 0.0)
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": skipping TPC geometry layer " << layer_geom->get_layer()
                  << " with invalid radius " << radius << std::endl;
      }
      continue;
    }

    m_tpc_layer_radii.push_back({static_cast<unsigned int>(layer_geom->get_layer()), radius});
  }

  std::sort(m_tpc_layer_radii.begin(), m_tpc_layer_radii.end(),
            [](const TpcLayerRadius &lhs, const TpcLayerRadius &rhs)
            { return lhs.layer < rhs.layer; });

  return !m_tpc_layer_radii.empty();
}

void PHG4TpcTruthPointBuilder::add_intersections(const PHG4Hit *hit,
                                                 const double track_segment_start,
                                                 std::vector<TruthPoint> &points) const
{
  if (!hit || m_tpc_layer_radii.empty())
  {
    return;
  }

  const double x0 = hit->get_x(0);
  const double y0 = hit->get_y(0);
  const double z0 = hit->get_z(0);
  const double t0 = hit->get_t(0);
  const double dx = hit->get_x(1) - x0;
  const double dy = hit->get_y(1) - y0;
  const double dz = hit->get_z(1) - z0;
  const double dt = hit->get_t(1) - t0;
  const double step_length = std::sqrt(square(dx) + square(dy) + square(dz));
  if (step_length <= 0.0 || !std::isfinite(step_length))
  {
    return;
  }

  const double a = square(dx) + square(dy);
  const double b = 2.0 * (dx * x0 + dy * y0);
  const double c0 = square(x0) + square(y0);
  static constexpr double t_tolerance = 1e-6;

  if (a <= 1e-12)
  {
    return;
  }

  const double r0 = std::sqrt(c0);
  const double r1 = std::sqrt(square(hit->get_x(1)) + square(hit->get_y(1)));
  const double closest_fraction = std::clamp(-b / (2.0 * a), 0.0, 1.0);
  const double closest_r2 = c0 + b * closest_fraction + a * square(closest_fraction);
  const double closest_r = std::sqrt(std::max(0.0, closest_r2));
  const double rmin = std::min({r0, r1, closest_r});
  const double rmax = std::max(r0, r1);

  for (const auto &layer_info : m_tpc_layer_radii)
  {
    const double radius = layer_info.radius;
    if (!std::isfinite(radius) || radius <= 0.0)
    {
      continue;
    }

    if (radius < rmin - m_radial_tolerance_cm || radius > rmax + m_radial_tolerance_cm)
    {
      continue;
    }

    const double c = c0 - square(radius);
    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < -1e-10)
    {
      continue;
    }

    discriminant = std::max(0.0, discriminant);
    const double sqrt_discriminant = std::sqrt(discriminant);
    const std::array<double, 2> candidates = {
        (-b - sqrt_discriminant) / (2.0 * a),
        (-b + sqrt_discriminant) / (2.0 * a)};

    std::vector<double> hit_fractions;

    for (const double candidate : candidates)
    {
      if (candidate < -t_tolerance || candidate > 1.0 + t_tolerance)
      {
        continue;
      }

      const double hit_fraction = std::clamp(candidate, 0.0, 1.0);
      const double xi = x0 + hit_fraction * dx;
      const double yi = y0 + hit_fraction * dy;
      const double residual = std::abs(std::sqrt(square(xi) + square(yi)) - radius);
      if (residual > m_radial_tolerance_cm)
      {
        continue;
      }

      const bool duplicate = std::any_of(hit_fractions.begin(), hit_fractions.end(),
                                         [hit_fraction](const double previous)
                                         { return std::abs(previous - hit_fraction) < 1e-5; });
      if (!duplicate)
      {
        hit_fractions.push_back(hit_fraction);
      }
    }

    for (const double hit_fraction : hit_fractions)
    {
      TruthPoint point;
      point.track_id = hit->get_trkid();
      point.shower_id = hit->get_shower_id();
      point.layer = layer_info.layer;
      point.x = x0 + hit_fraction * dx;
      point.y = y0 + hit_fraction * dy;
      point.z = z0 + hit_fraction * dz;
      point.t = t0 + hit_fraction * dt;
      point.px = hit->get_px(0) + hit_fraction * (hit->get_px(1) - hit->get_px(0));
      point.py = hit->get_py(0) + hit_fraction * (hit->get_py(1) - hit->get_py(0));
      point.pz = hit->get_pz(0) + hit_fraction * (hit->get_pz(1) - hit->get_pz(0));
      point.path = track_segment_start + hit_fraction * step_length;

      if (std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z))
      {
        points.push_back(point);
      }
    }
  }
}
