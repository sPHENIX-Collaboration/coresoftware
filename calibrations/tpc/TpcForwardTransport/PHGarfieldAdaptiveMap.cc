#include "PHGarfieldAdaptiveMap.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace
{
  using Map = PHGarfieldAdaptiveMap;
  constexpr double TwoPi = 6.283185307179586476925286766559;
  struct CoordinateHash
  {
    std::size_t operator()(const std::array<double, 3>& coordinate) const
    {
      std::size_t value = 0;
      for (double x : coordinate)
      {
        value ^= std::hash<double>{}(x) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
      }
      return value;
    }
  };

  bool require(bool condition, std::ostream& error, const std::string& message)
  {
    if (!condition)
    {
      error << "Adaptive LUT: " << message << '\n';
    }
    return condition;
  }

  bool geometry_valid(const Map::Geometry& g)
  {
    return g.phi_nodes >= 3 && g.radial_nodes >= 2 && g.longitudinal_nodes >= 2 &&
           g.phi_nodes <= 100000 && g.radial_nodes <= 100000 && g.longitudinal_nodes <= 100000 &&
           std::isfinite(g.inner_radius_cm) && std::isfinite(g.outer_radius_cm) &&
           g.inner_radius_cm > 0. && g.outer_radius_cm > g.inner_radius_cm &&
           std::isfinite(g.central_membrane_half_width_cm) && std::isfinite(g.readout_abs_z_cm) &&
           g.central_membrane_half_width_cm > 0. && g.readout_abs_z_cm > g.central_membrane_half_width_cm &&
           std::isfinite(g.nominal_velocity_cm_ns) && g.nominal_velocity_cm_ns > 0.;
  }

  double node(unsigned int i, unsigned int count, double lo, double hi)
  {
    if (!i)
    {
      return lo;
    }
    if (i + 1 == count)
    {
      return hi;
    }
    return lo + static_cast<double>(i) * (hi - lo) / static_cast<double>(count - 1);
  }

  bool policy_valid(const Map::Policy& p)
  {
    return std::isfinite(p.guard_um) && p.guard_um > 0. && p.max_depth <= 63 &&
           p.max_nodes_per_root >= 1 && p.max_queries_per_root >= 1 &&
           p.max_total_points >= 1 && p.max_total_points <= static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) &&
           p.max_total_cells >= 1 && p.max_total_cells <= static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) &&
           p.max_estimated_bytes >= 1024 && p.max_estimated_bytes <= Map::HardMaximumEstimatedBytes;
  }

  bool validate_impl(const Map::Data& data, std::ostream& error, bool complete)
  {
    if (!require(geometry_valid(data.geometry) && policy_valid(data.policy), error, "invalid geometry/policy"))
    {
      return false;
    }
    const auto count = Map::RootCount(data.geometry);
    std::uint64_t bytes = 0;
    if (!require(count && data.first_root <= data.end_root && data.end_root <= count &&
                     (!complete || (data.first_root == 0 && data.end_root == count)) &&
                     data.roots.size() == data.end_root - data.first_root &&
                     data.points.size() <= data.policy.max_total_points && data.cells.size() <= data.policy.max_total_cells &&
                     Map::EstimatedBytes(data.points.size(), data.cells.size(), count, bytes) &&
                     bytes <= data.policy.max_estimated_bytes &&
                     std::isfinite(data.maximum_drift_time_ns) && data.maximum_drift_time_ns >= 0.,
                 error,
                 "incomplete root range or allocation/count limit"))
    {
      return false;
    }
    std::unordered_map<std::array<double, 3>, std::uint64_t, CoordinateHash> coordinates;
    coordinates.reserve(data.points.size());
    for (std::size_t i = 0; i < data.points.size(); ++i)
    {
      const auto& p = data.points[i];
      if (!require(p.id == i, error, "invalid point identity") ||
          !Map::ValidatePoint(p, data.geometry, data.maximum_drift_time_ns, error))
      {
        return false;
      }
      if (!require(coordinates.emplace(p.coordinate, p.id).second, error, "duplicate sampled coordinate"))
      {
        return false;
      }
    }
    std::vector<unsigned char> visited(data.cells.size(), 0);
    std::uint64_t visited_count = 0;
    std::function<bool(std::uint64_t, std::uint64_t, unsigned int, const Map::Bounds&)> visit;
    visit = [&](std::uint64_t id, std::uint64_t root, unsigned int depth, const Map::Bounds& bounds)
    {
      if (!require(id < data.cells.size() && !visited[id], error, "missing/overlapping/cyclic child"))
      {
        return false;
      }
      visited[id] = 1;
      ++visited_count;
      const auto& c = data.cells[id];
      if (!require(c.id == id && c.root_id == root && c.depth == depth && depth <= data.policy.max_depth &&
                       c.lo == bounds.lo && c.hi == bounds.hi && c.probe_count <= 19 && c.probe_reaching <= c.probe_count &&
                       (c.angular_valid == 0 || c.angular_valid == 1) &&
                       std::isfinite(c.max_transverse_error_um) && c.max_transverse_error_um >= -1. &&
                       std::isfinite(c.max_time_error_um) && c.max_time_error_um >= -1. &&
                       std::all_of(c.directional_error_um.begin(), c.directional_error_um.end(), [](double x)
                                   { return std::isfinite(x) && x >= -1.; }) &&
                       std::all_of(c.directional_reach_changes.begin(), c.directional_reach_changes.end(), [](unsigned int x)
                                   { return x <= 18; }),
                   error, "invalid cell bounds/depth/probe diagnostics"))
      {
        return false;
      }
      for (unsigned int axis = 0; axis < 3; ++axis)
      {
        if (!require(c.hi[axis] > c.lo[axis], error, "empty/reversed cell extent"))
        {
          return false;
        }
      }
      if (c.split_axis >= 0)
      {
        if (!require(c.split_axis < 3 && c.child_id[0] >= 0 && c.child_id[1] >= 0 && depth < data.policy.max_depth,
                     error, "invalid internal split"))
        {
          return false;
        }
        const double midpoint = .5 * (c.lo[c.split_axis] + c.hi[c.split_axis]);
        if (!require(midpoint > c.lo[c.split_axis] && midpoint < c.hi[c.split_axis], error, "dyadic midpoint collapses at numerical precision"))
        {
          return false;
        }
        return visit(c.child_id[0], root, depth + 1, Map::ChildBounds(bounds, c.split_axis, 0)) &&
               visit(c.child_id[1], root, depth + 1, Map::ChildBounds(bounds, c.split_axis, 1));
      }
      if (!require(c.split_axis == -1 && c.child_id[0] == -1 && c.child_id[1] == -1 &&
                       c.leaf_class >= 0 && c.leaf_class <= 4 && c.fallback_reason >= 0 && c.fallback_reason <= 7,
                   error, "invalid leaf classification/children"))
      {
        return false;
      }
      unsigned int corners = 0;
      unsigned int reaching = 0;
      double minimum_angle = std::numeric_limits<double>::infinity();
      double maximum_angle = -minimum_angle;
      for (unsigned int corner = 0; corner < 8; ++corner)
      {
        const auto point_id = c.corner_id[corner];
        if (point_id == -1)
        {
          continue;
        }
        if (!require(point_id >= 0 && static_cast<std::uint64_t>(point_id) < data.points.size(), error, "invalid corner index"))
        {
          return false;
        }
        const auto& point = data.points[point_id];
        if (!require(point.coordinate == Map::CornerCoordinate(bounds, corner), error, "corner coordinate does not match dyadic cell"))
        {
          return false;
        }
        ++corners;
        reaching += point.reaches;
        if (point.reaches)
        {
          const double angle = point.value[1] / point.coordinate[1];
          minimum_angle = std::min(minimum_angle, angle);
          maximum_angle = std::max(maximum_angle, angle);
        }
      }
      const auto type = static_cast<Map::LeafClass>(c.leaf_class);
      if (type == Map::LeafClass::Certified)
      {
        if (!require(corners == 8 && reaching == 8 && c.probe_count == 19 && c.probe_reaching == 19 && c.angular_valid == 1 &&
                         maximum_angle - minimum_angle <= .5 * TwoPi && c.fallback_reason == 0 &&
                         c.max_transverse_error_um >= 0. && c.max_time_error_um >= 0. &&
                         std::max(c.max_transverse_error_um, c.max_time_error_um) <= data.policy.guard_um,
                     error, "uncertified interpolation leaf"))
        {
          return false;
        }
        // The complete stored3^3 stencil is independently checked. Rewriting a
        // guard diagnostic and its checksum cannot make an inaccurate leaf pass.
        double maximum_transverse = 0.;
        double maximum_time = 0.;
        double minimum_probe_angle = minimum_angle;
        double maximum_probe_angle = maximum_angle;
        for (unsigned int iphi = 0; iphi < 3; ++iphi)
        {
          for (unsigned int ir = 0; ir < 3; ++ir)
          {
            for (unsigned int iz = 0; iz < 3; ++iz)
            {
              const std::array<unsigned int, 3> index{{iphi, ir, iz}};
              std::array<double, 3> unwrapped{};
              std::array<double, 3> coordinate{};
              std::array<double, 3> fraction{};
              std::array<double, 3> prediction{};
              for (unsigned int axis = 0; axis < 3; ++axis)
              {
                if (index[axis] == 0)
                {
                  unwrapped[axis] = c.lo[axis];
                }
                else if (index[axis] == 2)
                {
                  unwrapped[axis] = c.hi[axis];
                }
                else
                {
                  unwrapped[axis] = .5 * (c.lo[axis] + c.hi[axis]);
                }
                coordinate[axis] = unwrapped[axis];
                fraction[axis] = .5 * index[axis];
              }
              coordinate[0] = Map::CanonicalPhi(coordinate[0]);
              const auto found = coordinates.find(coordinate);
              if (!require(found != coordinates.end() && data.points[found->second].reaches, error, "certified stencil has missing/lost sampled point"))
              {
                return false;
              }
              const auto& actual = data.points[found->second];
              const double angle = actual.value[1] / actual.coordinate[1];
              minimum_probe_angle = std::min(minimum_probe_angle, angle);
              maximum_probe_angle = std::max(maximum_probe_angle, angle);
              for (unsigned int corner = 0; corner < 8; ++corner)
              {
                const double weight = ((corner & 4U) ? fraction[0] : 1. - fraction[0]) * ((corner & 2U) ? fraction[1] : 1. - fraction[1]) * ((corner & 1U) ? fraction[2] : 1. - fraction[2]);
                for (unsigned int component = 0; component < 3; ++component)
                {
                  prediction[component] += weight * data.points[c.corner_id[corner]].value[component];
                }
              }
              const double r = actual.coordinate[1];
              const double ra = r + prediction[0];
              const double aa = prediction[1] / r;
              const double rb = r + actual.value[0];
              const double ab = actual.value[1] / r;
              maximum_transverse = std::max(maximum_transverse, 1.e4 * std::hypot(ra * std::cos(aa) - rb * std::cos(ab), ra * std::sin(aa) - rb * std::sin(ab)));
              maximum_time = std::max(maximum_time, 1.e4 * data.geometry.nominal_velocity_cm_ns * std::abs(prediction[2] - actual.value[2]));
            }
          }
        }
        const double tolerance = 1.e-7;
        if (!require(maximum_probe_angle - minimum_probe_angle <= .5 * TwoPi &&
                         std::abs(maximum_transverse - c.max_transverse_error_um) <= tolerance &&
                         std::abs(maximum_time - c.max_time_error_um) <= tolerance &&
                         std::max(maximum_transverse, maximum_time) <= data.policy.guard_um,
                     error, "certified error diagnostics disagree with stored stencil"))
        {
          return false;
        }
      }
      else
      {
        if (!require(c.fallback_reason != 0, error, "direct leaf without explicit reason"))
        {
          return false;
        }
        if (type == Map::LeafClass::Uncertified && !require(c.fallback_reason == static_cast<int>(Map::FallbackReason::QueryBudget), error, "invalid uncertified reason"))
        {
          return false;
        }
        if (type == Map::LeafClass::AllLost && !require(c.fallback_reason == static_cast<int>(Map::FallbackReason::SampledAllLost), error, "invalid sampled-AllLost reason"))
        {
          return false;
        }
        if (type == Map::LeafClass::Mixed && !require(c.fallback_reason == 1 || c.fallback_reason == 3 || c.fallback_reason == 4, error, "invalid Mixed reason"))
        {
          return false;
        }
        if (type == Map::LeafClass::NumericallyUnresolved && !require(c.fallback_reason == 2 || c.fallback_reason == 3 || c.fallback_reason == 4 || c.fallback_reason == 6, error, "invalid numerical reason"))
        {
          return false;
        }
        if (type == Map::LeafClass::AllLost && !require(corners == 8 && !reaching && c.probe_count == 19 && !c.probe_reaching,
                                                        error, "invalid sampled-AllLost leaf"))
        {
          return false;
        }
        if (type == Map::LeafClass::Mixed && !require(corners == 8 && c.probe_count == 19 &&
                                                          reaching + c.probe_reaching > 0 && reaching + c.probe_reaching < 27,
                                                      error, "invalid Mixed leaf"))
        {
          return false;
        }
        if (type == Map::LeafClass::NumericallyUnresolved && !require(corners == 8 && reaching == 8 && c.probe_count == 19 && c.probe_reaching == 19 &&
                                                                          (!c.angular_valid || std::max(c.max_transverse_error_um, c.max_time_error_um) > data.policy.guard_um),
                                                                      error, "invalid numerical fallback leaf"))
        {
          return false;
        }
      }
      return true;
    };
    for (std::size_t i = 0; i < data.roots.size(); ++i)
    {
      const auto& root = data.roots[i];
      if (!require(root.root_id == data.first_root + i && root.cell_id < data.cells.size(), error, "duplicate/missing root ownership"))
      {
        return false;
      }
      const auto before = visited_count;
      if (!visit(root.cell_id, root.root_id, 0, Map::RootBounds(data.geometry, root.root_id)))
      {
        return false;
      }
      const auto after = visited_count;
      if (!require(after - before <= data.policy.max_nodes_per_root, error, "per-root node budget exceeded"))
      {
        return false;
      }
    }
    return require(std::all_of(visited.begin(), visited.end(), [](unsigned char v)
                               { return v == 1; }),
                   error, "orphan/disconnected cells");
  }

}  // namespace

bool PHGarfieldAdaptiveMap::ValidatePoint(const Point& p, const Geometry& g,
                                          double maximum_drift_time_ns, std::ostream& error)
{
  const auto fail = [&](const char* reason)
  {
    const auto flags = error.flags();
    const auto precision = error.precision();
    error << "Adaptive LUT: " << reason << std::setprecision(17)
          << " point_id=" << p.id << " coordinate=(" << p.coordinate[0] << ',' << p.coordinate[1] << ',' << p.coordinate[2]
          << ") reaches=" << p.reaches << " status=" << p.forward_status
          << " value=(" << p.value[0] << ',' << p.value[1] << ',' << p.value[2] << ')'
          << " endpoint_radius_cm=" << p.coordinate[1] + p.value[0]
          << " allowed_radii_cm=(" << g.inner_radius_cm << ',' << g.outer_radius_cm << ")\n";
    error.flags(flags);
    error.precision(precision);
    return false;
  };
  if (!(std::all_of(p.coordinate.begin(), p.coordinate.end(), [](double x)
                    { return std::isfinite(x); }) &&
        p.coordinate[0] >= 0. && p.coordinate[0] < TwoPi &&
        p.coordinate[1] >= g.inner_radius_cm && p.coordinate[1] <= g.outer_radius_cm &&
        std::abs(p.coordinate[2]) >= g.central_membrane_half_width_cm &&
        std::abs(p.coordinate[2]) <= g.readout_abs_z_cm &&
        (p.reaches == 0 || p.reaches == 1) &&
        p.forward_status == static_cast<int>(p.reaches ? PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout : PHGarfieldForwardTransport::ForwardDriftStatus::LeftDriftVolume) &&
        std::all_of(p.value.begin(), p.value.end(), [](double x)
                    { return std::isfinite(x); }) &&
        p.value[2] >= 0. && (!p.reaches || p.value[2] <= maximum_drift_time_ns)))
  {
    return fail("invalid point encoding/status/coordinate");
  }
  const bool readout = std::abs(p.coordinate[2]) == g.readout_abs_z_cm;
  if (!((p.reaches || p.value == std::array<double, 3>{{0., 0., 0.}}) &&
        (!p.reaches || (p.coordinate[1] + p.value[0] >= g.inner_radius_cm - EndpointRadialToleranceCm &&
                        p.coordinate[1] + p.value[0] <= g.outer_radius_cm + EndpointRadialToleranceCm &&
                        std::abs(p.value[1] / p.coordinate[1]) <= .5 * TwoPi)) &&
        (!readout || (p.reaches && p.value == std::array<double, 3>{{0., 0., 0.}}))))
  {
    return fail("invalid lost/endpoint/angular/readout point values");
  }
  return true;
}

PHGarfieldAdaptiveMap::Geometry PHGarfieldAdaptiveMap::FromConfig(const PHGarfieldForwardMap::Config& c)
{
  return {c.phi_nodes, c.radial_nodes, c.longitudinal_nodes, c.inner_radius_cm, c.outer_radius_cm,
          c.central_membrane_half_width_cm, c.readout_abs_z_cm, c.nominal_phg4_drift_velocity_cm_per_ns};
}

double PHGarfieldAdaptiveMap::CanonicalPhi(double phi)
{
  phi = std::fmod(phi, TwoPi);
  if (phi < 0.)
  {
    phi += TwoPi;
  }
  if (phi >= TwoPi)
  {
    phi = 0.;
  }
  return phi == 0. ? 0. : phi;
}

std::uint64_t PHGarfieldAdaptiveMap::RootCount(const Geometry& g)
{
  if (!geometry_valid(g))
  {
    return 0;
  }
  std::uint64_t value = 2;
  for (const auto factor : {g.phi_nodes, g.radial_nodes - 1, g.longitudinal_nodes - 1})
  {
    if (value > std::numeric_limits<std::uint64_t>::max() / factor)
    {
      return 0;
    }
    value *= factor;
  }
  return value;
}

std::uint64_t PHGarfieldAdaptiveMap::RootIndex(const Geometry& g, unsigned int side, double phi, double r, double z)
{
  // Compare against the same exact node expressions as RootBounds. Division
  // followed by integer truncation can assign a representable face to its
  // lower neighbour through roundoff, violating the declared half-open rule.
  auto interval = [](double value, unsigned int intervals, const auto& coordinate)
  {
    unsigned int lo = 0;
    unsigned int hi = intervals;
    while (lo + 1 < hi)
    {
      const auto mid = lo + (hi - lo) / 2;
      if (value >= coordinate(mid))
      {
        lo = mid;
      }
      else
      {
        hi = mid;
      }
    }
    return lo;
  };
  const auto iphi = interval(phi, g.phi_nodes, [&](unsigned int i)
                             { return static_cast<double>(i) * TwoPi / g.phi_nodes; });
  const auto ir = interval(r, g.radial_nodes - 1, [&](unsigned int i)
                           { return node(i, g.radial_nodes, g.inner_radius_cm, g.outer_radius_cm); });
  const double zlo = side ? g.central_membrane_half_width_cm : -g.readout_abs_z_cm;
  const double zhi = side ? g.readout_abs_z_cm : -g.central_membrane_half_width_cm;
  const auto iz = interval(z, g.longitudinal_nodes - 1, [&](unsigned int i)
                           { return node(i, g.longitudinal_nodes, zlo, zhi); });
  return ((static_cast<std::uint64_t>(side) * g.phi_nodes + iphi) * (g.radial_nodes - 1) + ir) * (g.longitudinal_nodes - 1) + iz;
}

PHGarfieldAdaptiveMap::Bounds PHGarfieldAdaptiveMap::RootBounds(const Geometry& g, std::uint64_t id)
{
  if (id >= RootCount(g))
  {
    throw std::out_of_range("adaptive root identity");
  }
  const auto iz = id % (g.longitudinal_nodes - 1);
  id /= g.longitudinal_nodes - 1;
  const auto ir = id % (g.radial_nodes - 1);
  id /= g.radial_nodes - 1;
  const auto iphi = id % g.phi_nodes;
  const auto side = id / g.phi_nodes;
  const double zlo = side ? g.central_membrane_half_width_cm : -g.readout_abs_z_cm;
  const double zhi = side ? g.readout_abs_z_cm : -g.central_membrane_half_width_cm;
  return {{{static_cast<double>(iphi) * TwoPi / g.phi_nodes, node(ir, g.radial_nodes, g.inner_radius_cm, g.outer_radius_cm), node(iz, g.longitudinal_nodes, zlo, zhi)}},
          {{static_cast<double>(iphi + 1) * TwoPi / g.phi_nodes, node(ir + 1, g.radial_nodes, g.inner_radius_cm, g.outer_radius_cm), node(iz + 1, g.longitudinal_nodes, zlo, zhi)}}};
}

PHGarfieldAdaptiveMap::Bounds PHGarfieldAdaptiveMap::ChildBounds(const Bounds& parent, int axis, unsigned int child)
{
  if (axis < 0 || axis > 2 || child > 1)
  {
    throw std::out_of_range("adaptive child identity");
  }
  auto result = parent;
  const double midpoint = .5 * (parent.lo[axis] + parent.hi[axis]);
  if (child)
  {
    result.lo[axis] = midpoint;
  }
  else
  {
    result.hi[axis] = midpoint;
  }
  return result;
}

std::array<double, 3> PHGarfieldAdaptiveMap::CornerCoordinate(const Bounds& b, unsigned int corner)
{
  if (corner >= 8)
  {
    throw std::out_of_range("adaptive corner identity");
  }
  return {{CanonicalPhi((corner & 4U) ? b.hi[0] : b.lo[0]), (corner & 2U) ? b.hi[1] : b.lo[1], (corner & 1U) ? b.hi[2] : b.lo[2]}};
}

bool PHGarfieldAdaptiveMap::EstimatedBytes(std::uint64_t points, std::uint64_t cells, std::uint64_t roots, std::uint64_t& bytes)
{
  bytes = 1024;
  // Conservative peak allowance includes vector growth, builder/cache and
  // independent validation indexes; it is intentionally above sizeof records.
  for (const auto& item : {std::make_pair(points, 384ULL), std::make_pair(cells, 768ULL), std::make_pair(roots, 32ULL)})
  {
    if (item.first > (std::numeric_limits<std::uint64_t>::max() - bytes) / item.second)
    {
      return false;
    }
    bytes += item.first * item.second;
  }
  return true;
}

bool PHGarfieldAdaptiveMap::Validate(const Data& data, std::ostream& error, bool complete)
{
  return validate_impl(data, error, complete);
}

bool PHGarfieldAdaptiveMap::HasUsableCorners(const Data& data, const Cell& cell)
{
  if (cell.split_axis != -1 ||
      (cell.leaf_class != static_cast<int>(LeafClass::Certified) &&
       cell.leaf_class != static_cast<int>(LeafClass::Uncertified) &&
       cell.leaf_class != static_cast<int>(LeafClass::NumericallyUnresolved)) ||
      cell.fallback_reason == static_cast<int>(FallbackReason::AngularBranch))
  {
    return false;
  }
  if (cell.leaf_class == static_cast<int>(LeafClass::NumericallyUnresolved) && cell.angular_valid != 1)
  {
    return false;
  }
  const Bounds bounds{cell.lo, cell.hi};
  return UsableCornerValues(data.geometry, bounds, data.maximum_drift_time_ns,
                            [&](unsigned int corner) -> const std::array<double, 3>*
                            {
                              const auto id = cell.corner_id[corner];
                              if (id < 0 || static_cast<std::uint64_t>(id) >= data.points.size())
                              {
                                return nullptr;
                              }
                              const auto& point = data.points[id];
                              return point.reaches == 1 && point.coordinate == CornerCoordinate(bounds, corner) ? &point.value : nullptr;
                            });
}

std::string PHGarfieldAdaptiveMap::BuildCanonicalPolicy(const Data& d)
{
  std::ostringstream s;
  s.imbue(std::locale::classic());
  s << std::setprecision(17);
  s << "adaptive_policy=1;probe=" << ProbeScheme << ";split=" << SplitPolicy << ";all_lost=" << AllLostPolicy
    << ";leaf_encoding=" << LeafEncoding << ";fallback_encoding=" << FallbackEncoding
    << ";phi=" << d.geometry.phi_nodes << ";r=" << d.geometry.radial_nodes << ";z=" << d.geometry.longitudinal_nodes
    << ";rmin=" << d.geometry.inner_radius_cm << ";rmax=" << d.geometry.outer_radius_cm
    << ";cm=" << d.geometry.central_membrane_half_width_cm << ";readout=" << d.geometry.readout_abs_z_cm
    << ";velocity=" << d.geometry.nominal_velocity_cm_ns << ";guard=" << d.policy.guard_um
    << ";depth=" << d.policy.max_depth << ";nodes_per_root=" << d.policy.max_nodes_per_root
    << ";queries_per_root=" << d.policy.max_queries_per_root << ";max_points=" << d.policy.max_total_points
    << ";max_cells=" << d.policy.max_total_cells << ";max_bytes=" << d.policy.max_estimated_bytes;
  return s.str();
}
