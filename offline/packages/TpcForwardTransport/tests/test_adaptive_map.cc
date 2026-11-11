#include "../PHGarfieldAdaptiveMap.h"
#include "../PHGarfieldCompactMap.h"
#include "PHGarfieldTestUtils.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

namespace
{
  using Map = PHGarfieldAdaptiveMap;
  using PHGarfieldTest::check;
  using Stored = PHGarfieldCompactMap;
  Stored::Data pack(const Map::Data& data)
  {
    Stored::Data result;
    std::ostringstream errors;
    const bool valid = Stored::Pack(data, "0123456789abcdef", result, Stored::Limits{}, errors, false);
    check(valid, errors.str());
    return result;
  }

  const Map::Cell* find_leaf(const Map::Data& source, unsigned int side, double phi, double r, double z)
  {
    const auto stored = pack(source);
    Stored::Lookup lookup;
    return Stored::FindLeaf(stored, side, phi, r, z, lookup) ? &source.cells.at(lookup.cell_id) : nullptr;
  }

  bool interpolate(const Map::Data& source, const Map::Cell& cell, double phi, double r, double z,
                   PHGarfieldForwardMap::NodeValue& prediction)
  {
    const auto stored = pack(source);
    Stored::Lookup lookup;
    return Stored::FindLeaf(stored, z > 0., phi, r, z, lookup) && lookup.cell_id == cell.id &&
           Stored::Interpolate(stored, lookup, phi, r, z, prediction);
  }

  Map::Geometry geometry()
  {
    Map::Geometry g;
    g.phi_nodes = 3;
    g.radial_nodes = 2;
    g.longitudinal_nodes = 2;
    g.inner_radius_cm = 10.;
    g.outer_radius_cm = 12.;
    g.central_membrane_half_width_cm = 1.;
    g.readout_abs_z_cm = 3.;
    g.nominal_velocity_cm_ns = .1;
    return g;
  }

  Map::Sampler sampler(unsigned int axis, bool mixed = false)
  {
    return [axis, mixed](const std::array<double, 3>& x, Map::Point& p, std::ostream&)
    {
      p.reaches = !mixed || x[1] < 11.1 || std::abs(x[2]) == 3.;
      p.forward_status = static_cast<int>(p.reaches ? PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout : PHGarfieldForwardTransport::ForwardDriftStatus::LeftDriftVolume);
      p.value = {{0., 0., 0.}};
      if (p.reaches)
      {
        p.value[2] = (3. - std::abs(x[2])) / .1;
        if (!mixed)
        {
          const double fraction = (3. - std::abs(x[2])) / 2.;
          if (axis == 0)
          {
            p.value[1] = .01 * x[1] * std::cos(x[0]) * fraction;
          }
          if (axis == 1)
          {
            p.value[0] = .02 * (x[1] - 10.) * (12. - x[1]) * fraction;
          }
          if (axis == 2)
          {
            p.value[2] += .1 * (3. - std::abs(x[2])) * (3. - std::abs(x[2]));
          }
        }
      }
      return true;
    };
  }

  double volume(const Map::Cell& c)
  {
    return .5 * (c.hi[1] * c.hi[1] - c.lo[1] * c.lo[1]) * (c.hi[0] - c.lo[0]) * (c.hi[2] - c.lo[2]);
  }

  Map::Sampler localized_corner_sampler()
  {
    return [](const std::array<double, 3>& x, Map::Point& p, std::ostream&)
    {
      p.reaches = 1;
      p.forward_status = static_cast<int>(PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout);
      const double remaining = 3. - std::abs(x[2]);
      const double r = x[1] - 10.;
      // A narrow difficult inner-radius corner and smooth bulk that needs only
      // two radial subdivisions. DFS used the last split in the inner corner,
      // leaving the larger outer half unrefined under the five-node budget.
      p.value = {{0., 0., remaining / .1 + .5 * remaining * (.08 * r * (2. - r) + .5 * std::exp(-200. * r * r))}};
      return true;
    };
  }

  Map::Data build(unsigned int axis, const Map::Policy& policy = Map::Policy{}, bool mixed = false)
  {
    Map::Data d;
    Map::BuildSummary s;
    std::ostringstream log;
    std::ostringstream error;
    const auto g = geometry();
    check(Map::Build(g, policy, 0, Map::RootCount(g), sampler(axis, mixed), d, s, log, error), error.str());
    check(Map::Validate(d, error, true), error.str());
    return d;
  }

}  // namespace

int main()
{
  try
  {
    const auto g = geometry();
    {
      Map::Point p;
      p.coordinate = {{2.5951546084829742, 21.620715725806456, -98.63753761574074}};
      p.reaches = 1;
      p.forward_status = static_cast<int>(PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout);
      p.value = {{-0.020715726070122287, -0.10039665469166942, 536.14044537822394}};
      std::ostringstream error;
      check(Map::ValidatePoint(p, Map::Geometry{}, 537., error), "actual failed-shard endpoint violates the direct solver's accepted tolerance");
      p.value[0] -= 2. * Map::EndpointRadialToleranceCm;
      check(!Map::ValidatePoint(p, Map::Geometry{}, 537., error), "endpoint beyond the direct solver tolerance accepted");
      p.coordinate = {{0., 76.4, 50.}};
      p.value = {{.5 * Map::EndpointRadialToleranceCm, 0., 1.}};
      check(Map::ValidatePoint(p, Map::Geometry{}, 2., error), "outer endpoint within native tolerance rejected");
      p.value[0] = 2. * Map::EndpointRadialToleranceCm;
      check(!Map::ValidatePoint(p, Map::Geometry{}, 2., error), "outer endpoint beyond native tolerance accepted");
    }
    {
      unsigned int calls = 0;
      Map::Data invalid;
      Map::BuildSummary summary;
      std::ostringstream log;
      std::ostringstream error;
      const auto bad = [&](const std::array<double, 3>&, Map::Point& p, std::ostream&)
      {++calls;p.reaches=1;p.forward_status=static_cast<int>(PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout);p.value={{-0.001,0.,0.}};return true; };
      check(!Map::Build(g, Map::Policy{}, 0, 1, bad, invalid, summary, log, error), "invalid endpoint accepted during sampling");
      check(calls == 1 && summary.actual_sampler_queries == 1 && invalid.points.empty(), "invalid point did not abort immediately without partial publication");
      check(error.str().find("point_id=0 coordinate=(") != std::string::npos &&
                error.str().find("endpoint_radius_cm=") != std::string::npos &&
                error.str().find("root=0 after queries=1") != std::string::npos,
            "failure omits reproducible point/root context");
    }
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
      const auto data = build(axis);
      for (const auto& root : data.roots)
      {
        check(data.cells[root.cell_id].split_axis == static_cast<int>(axis), "split did not follow actual directional curvature");
      }
      for (const auto& cell : data.cells)
      {
        check(cell.split_axis < 0 || cell.split_axis == static_cast<int>(axis), "unnecessary isotropic/other-axis split");
      }
      for (unsigned int side = 0; side < 2; ++side)
      {
        for (unsigned int i = 0; i < 50; ++i)
        {
          const double phi = (i + .37) * 6.283185307179586 / 50.;
          const double r = 10. + 2. * (i + .13) / 50.;
          const double z = (side ? 1. : -1.) * (1. + 2. * (i + .71) / 50.);
          const auto* cell = find_leaf(data, side, phi, r, z);
          check(cell && cell->leaf_class == static_cast<int>(Map::LeafClass::Certified), "analytic smooth sample not certified");
          PHGarfieldForwardMap::NodeValue prediction;
          check(interpolate(data, *cell, phi, r, z, prediction), "analytic interpolation failed");
          Map::Point actual;
          std::ostringstream error;
          sampler(axis)({{phi, r, z}}, actual, error);
          check(1.e4 * std::abs(prediction.delta_r_cm - actual.value[0]) <= 10.01, "off-stencil transverse error exceeds guard");
          check(1.e4 * r * std::abs(prediction.delta_rphi_cm / r - actual.value[1] / r) <= 10.01, "off-stencil angular error exceeds guard");
          check(1.e4 * .1 * std::abs(prediction.drift_time_ns - actual.value[2]) <= 10.01, "off-stencil time error exceeds guard");
        }
      }
      check(Stored::ContentFingerprint(pack(data)) == Stored::ContentFingerprint(pack(build(axis))), "nondeterministic generation");
    }
    auto radial = build(1);
    const auto& root = radial.cells[radial.roots[0].cell_id];
    const double face = .5 * (root.lo[1] + root.hi[1]);
    const auto* chosen = find_leaf(radial, 0, .3, face, -2.);
    check(chosen && chosen->lo[1] == face, "split face equality did not select upper child");
    const auto* first = find_leaf(radial, 0, 0., 10.5, -2.);
    const auto* periodic = find_leaf(radial, 0, 6.283185307179586, 10.5, -2.);
    check(first && periodic && first->id == periodic->id, "phi2pi was not canonicalphi0");
    check(Map::CanonicalPhi(-1.e-30) < 6.283185307179586, "negative tinyphi rounded outside canonical domain");
    const auto* seam = find_leaf(radial, 0, std::nextafter(6.283185307179586, 0.), 10.5, -2.);
    check(seam && seam->root_id != first->root_id, "near2pi did not select last phi cell");

    Map::Policy boundary_policy;
    boundary_policy.max_depth = 10;
    boundary_policy.max_nodes_per_root = 511;
    boundary_policy.max_queries_per_root = 8192;
    auto boundary = build(0, boundary_policy, true);
    for (const auto& r : boundary.roots)
    {
      check(boundary.cells[r.cell_id].split_axis == 1, "reach transitions did not select radial refinement");
    }
    const auto* interior = find_leaf(boundary, 0, .3, 10.5, -2.);
    const auto* edge = find_leaf(boundary, 0, .3, 11.1, -2.);
    check(interior && interior->leaf_class == static_cast<int>(Map::LeafClass::Certified), "mixed root failed to recover reachable interior");
    check(edge && edge->leaf_class == static_cast<int>(Map::LeafClass::Mixed) && edge->hi[1] - edge->lo[1] <= 2. / 64., "physical boundary did not shrink with anisotropic refinement");

    Map::Policy corner_policy;
    corner_policy.max_nodes_per_root = 5;
    Map::Data corner;
    Map::BuildSummary corner_summary;
    std::ostringstream corner_log;
    std::ostringstream corner_errors;
    check(Map::Build(g, corner_policy, 0, Map::RootCount(g), localized_corner_sampler(), corner, corner_summary, corner_log, corner_errors), corner_errors.str());
    for (const auto& r : corner.roots)
    {
      const auto& parent = corner.cells[r.cell_id];
      check(parent.split_axis == 1, "localized fixture did not measure radial curvature");
      const auto& lower = corner.cells[parent.child_id[0]];
      const auto& upper = corner.cells[parent.child_id[1]];
      check(lower.split_axis == -1 && lower.fallback_reason == static_cast<int>(Map::FallbackReason::NodeBudget), "small difficult corner consumed the final split");
      check(upper.split_axis == 1, "larger unresolved cylinder volume was not refined first");
      double certified_volume = 0.;
      for (const auto& c : corner.cells)
      {
        if (c.root_id == r.root_id && c.split_axis == -1 && c.leaf_class == static_cast<int>(Map::LeafClass::Certified))
        {
          certified_volume += volume(c);
        }
      }
      check(certified_volume / volume(parent) > .52, "volume-first refinement did not recover the smooth outer half under a five-node budget");
      Map::Data shard;
      Map::BuildSummary shard_summary;
      check(Map::Build(g, corner_policy, r.root_id, r.root_id + 1, localized_corner_sampler(), shard, shard_summary, corner_log, corner_errors), corner_errors.str());
      std::map<std::array<double, 3>, std::array<double, 3>> full_samples;
      for (const auto& p : corner.points)
      {
        full_samples.emplace(p.coordinate, p.value);
      }
      for (const auto& p : shard.points)
      {
        check(full_samples.contains(p.coordinate) && full_samples.at(p.coordinate) == p.value, "sharding changed the sampled physical response");
      }
      for (unsigned int i = 0; i < 20; ++i)
      {
        const double phi = parent.lo[0] + (parent.hi[0] - parent.lo[0]) * (i + .31) / 20.;
        const double radius = parent.lo[1] + (parent.hi[1] - parent.lo[1]) * (i + .47) / 20.;
        const double z = parent.lo[2] + (parent.hi[2] - parent.lo[2]) * (i + .59) / 20.;
        const unsigned int side = z > 0.;
        const auto* a = find_leaf(corner, side, phi, radius, z);
        const auto* b = find_leaf(shard, side, phi, radius, z);
        check(a && b && a->lo == b->lo && a->hi == b->hi && a->leaf_class == b->leaf_class && a->fallback_reason == b->fallback_reason, "volume-priority routing depends on shard boundaries");
      }
    }
    check(Stored::ContentFingerprint(pack(corner)) == [&]()
              {Map::Data repeated;Map::BuildSummary repeated_summary;check(Map::Build(g,corner_policy,0,Map::RootCount(g),localized_corner_sampler(),repeated,repeated_summary,corner_log,corner_errors),corner_errors.str());return Stored::ContentFingerprint(pack(repeated)); }(),
          "equal-volume tie ordering is not deterministic");

    Map::Policy query_budget;
    query_budget.max_queries_per_root = 4;
    const auto limited = build(1, query_budget);
    for (const auto& c : limited.cells)
    {
      check(c.leaf_class == static_cast<int>(Map::LeafClass::Uncertified) && c.fallback_reason == static_cast<int>(Map::FallbackReason::QueryBudget), "query-budget leaf was falsely certified");
    }
    Map::Policy node_budget;
    node_budget.max_nodes_per_root = 1;
    const auto unsplit = build(1, node_budget);
    for (const auto& c : unsplit.cells)
    {
      check(c.fallback_reason == static_cast<int>(Map::FallbackReason::NodeBudget), "node budget not preserved explicitly");
    }
    std::ostringstream errors;
    std::ostringstream log;
    Map::Data failed;
    Map::BuildSummary summary;
    Map::Policy allocation;
    allocation.max_total_points = 1;
    check(!Map::Build(g, allocation, 0, Map::RootCount(g), sampler(1), failed, summary, log, errors), "global allocation limit was ignored");
    std::uint64_t bytes = 0;
    check(!Map::EstimatedBytes(std::numeric_limits<std::uint64_t>::max(), 1, 1, bytes), "allocation arithmetic overflow was not caught");
    constexpr std::uint64_t gib = 1024ULL * 1024 * 1024;
    check(Map::Policy{}.max_estimated_bytes == 2 * gib && Map::Policy{}.max_total_points == 5000000 && Map::Policy{}.max_total_cells == 2000000,
          "raising the engineering ceiling changed nominal default budgets");
    // Historical capacities and the engineering ceiling must all admit a tiny
    // shard without reserving their declared maximum allocation in advance.
    const std::array<std::uint64_t, 4> capacities[] = {
        {{600000000, 100000000, 320 * gib, 128 * gib}},
        {{1500000000, 250000000, 768 * gib, 512 * gib}},
        {{200000000, 40000000, 128 * gib, 0}},
        {{100000000, 20000000, 64 * gib, 0}},
        {{600000000, 100000000, Map::HardMaximumEstimatedBytes, 0}},
    };
    Map::Policy capacity;
    for (const auto& limits : capacities)
    {
      capacity.max_total_points = limits[0];
      capacity.max_total_cells = limits[1];
      capacity.max_estimated_bytes = limits[2];
      check(Map::EstimatedBytes(limits[0], limits[1], Map::RootCount(Map::Geometry{}), bytes) && bytes > limits[3] && bytes < limits[2],
            "declared footprint does not fit its explicit capacity");
      Map::Data tiny_capacity;
      Map::BuildSummary capacity_summary;
      check(Map::Build(g, capacity, 0, 1, sampler(1), tiny_capacity, capacity_summary, log, errors), errors.str());
      check(tiny_capacity.points.size() < 1000 && tiny_capacity.cells.size() < 1000 &&
                tiny_capacity.points.capacity() < 1000 && tiny_capacity.cells.capacity() < 1000,
            "declared global capacity caused eager allocation/growth");
    }
    auto excess = capacity;
    excess.max_estimated_bytes = Map::HardMaximumEstimatedBytes + 1;
    std::uint64_t sampler_calls = 0;
    const auto counting = [&](const std::array<double, 3>& x, Map::Point& p, std::ostream& e)
    {++sampler_calls;return sampler(1)(x,p,e); };
    check(!Map::Build(g, excess, 0, 1, counting, failed, summary, log, errors) && sampler_calls == 0, "oversized hard-cap policy reached the physical sampler");
    auto corrupt = radial;
    corrupt.cells[corrupt.roots[0].cell_id].child_id[0] = corrupt.roots[0].cell_id;
    check(!Map::Validate(corrupt, errors, true), "cycle accepted");
    corrupt = radial;
    corrupt.roots[1].root_id = corrupt.roots[0].root_id;
    check(!Map::Validate(corrupt, errors, true), "duplicate root accepted");
    corrupt = radial;
    corrupt.cells.back().lo[1] += .01;
    check(!Map::Validate(corrupt, errors, true), "overlap/gap accepted");
    corrupt = radial;
    corrupt.points.push_back(corrupt.points.front());
    corrupt.points.back().id = corrupt.points.size() - 1;
    check(!Map::Validate(corrupt, errors, true), "duplicate physical point accepted");
    corrupt = radial;
    for (auto& point : corrupt.points)
    {
      if (std::abs(point.coordinate[2]) < g.readout_abs_z_cm)
      {
        point.value[2] += 1.;
        break;
      }
    }
    corrupt.maximum_drift_time_ns += 1.;
    check(!Map::Validate(corrupt, errors, true), "changed actual stencil passed stored guard diagnostics");
    corrupt = radial;
    corrupt.points.front().value[0] = g.outer_radius_cm;
    check(!Map::Validate(corrupt, errors, true), "endpoint outside gas passed");
    corrupt = radial;
    corrupt.roots.pop_back();
    --corrupt.end_root;
    check(!Map::Validate(corrupt, errors, true), "incomplete runtime tree accepted");

    std::cout << "PASS adaptive axes, off-stencil errors, boundary ties, reach refinement, volume-priority budgets, shard invariance and capacity limits\n";
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "FAIL adaptive map: " << e.what() << '\n';
    return 1;
  }
}
