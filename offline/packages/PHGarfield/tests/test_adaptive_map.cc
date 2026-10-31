#include "../PHGarfieldAdaptiveMap.h"
#include "../PHGarfieldRootOutput.h"

#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TParameter.h>
#include <TTree.h>

#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
  using Map = PHGarfieldAdaptiveMap;
  void check(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
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
      p.forward_status = static_cast<int>(p.reaches ? PHGarfield::ForwardDriftStatus::ReachedReadout : PHGarfield::ForwardDriftStatus::LeftDriftVolume);
      p.value = {{0., 0., 0.}};
      if (p.reaches)
      {
        p.value[2] = (3. - std::abs(x[2])) / .1;
        if (!mixed)
        {
          const double fraction = (3. - std::abs(x[2])) / 2.;
          if (axis == 0) p.value[1] = .01 * x[1] * std::cos(x[0]) * fraction;
          if (axis == 1) p.value[0] = .02 * (x[1] - 10.) * (12. - x[1]) * fraction;
          if (axis == 2) p.value[2] += .1 * (3. - std::abs(x[2])) * (3. - std::abs(x[2]));
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
      p.forward_status = static_cast<int>(PHGarfield::ForwardDriftStatus::ReachedReadout);
      const double remaining = 3. - std::abs(x[2]), r = x[1] - 10.;
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
    std::ostringstream log, error;
    const auto g = geometry();
    check(Map::Build(g, policy, 0, Map::RootCount(g), sampler(axis, mixed), d, s, log, error), error.str());
    check(Map::Validate(d, error, true), error.str());
    return d;
  }

  void check_persisted_multibasket_tree(TFile& file, TDirectory* directory, const char* name)
  {
    auto* tree = dynamic_cast<TTree*>(directory->Get(name));
    check(tree && tree->GetEntries() > 0, "missing nonempty streamed primitive tree");
    auto* branches = tree->GetListOfBranches();
    for (int branch_index = 0; branch_index < branches->GetEntriesFast(); ++branch_index)
    {
      auto* branch = dynamic_cast<TBranch*>(branches->At(branch_index));
      check(branch && branch->GetEntries() == tree->GetEntries(), "streamed primitive branch count mismatch");
      const auto* starts = branch->GetBasketEntry();
      const auto* sizes = branch->GetBasketBytes();
      const int last = branch->GetWriteBasket();
      check(starts && sizes && last >= 0 && last < branch->GetMaxBaskets(), "missing streamed basket metadata");
      Long64_t covered = 0;
      unsigned int persisted = 0;
      // Check persisted coverage before GetEntry can materialize any basket.
      // Inline-only baskets would let the old detached writer falsely pass a
      // value-only roundtrip at small size but cannot scale to a full map.
      for (int basket_index = 0; basket_index <= last; ++basket_index)
      {
        const auto begin = starts[basket_index];
        const auto end = basket_index < last ? starts[basket_index + 1] : tree->GetEntries();
        if (end == begin) continue;
        const auto seek = branch->GetBasketSeek(basket_index);
        check(begin == covered && end > begin && end <= tree->GetEntries(), "persisted baskets do not tile every primitive row");
        check(seek > 0 && sizes[basket_index] > 0 && seek <= file.GetSize() && sizes[basket_index] <= file.GetSize() - seek,
              "primitive payload was retained inline instead of streamed into the declared output");
        covered = end;
        ++persisted;
      }
      check(covered == tree->GetEntries() && persisted >= 2, "fixture did not exercise multiple persisted baskets per branch");
      for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
        check(branch->GetEntry(entry) > 0, "streamed primitive branch row is unreadable");
    }
  }

  void test_streamed_multibasket_roundtrip(const std::string& filename)
  {
    auto g = geometry();
    g.phi_nodes = 36;
    g.radial_nodes = 13;
    g.longitudinal_nodes = 13;
    Map::Policy policy;
    policy.max_nodes_per_root = 1;
    Map::Data original;
    Map::BuildSummary summary;
    std::ostringstream errors, log;
    check(Map::Build(g, policy, 0, Map::RootCount(g), sampler(3, true), original, summary, log, errors), errors.str());
    check(original.points.size() > 50000 && original.cells.size() > 8192, "streamed fixture did not exceed basket thresholds");
    std::array<std::uint64_t, 5> classes{};
    for (const auto& cell : original.cells) ++classes.at(cell.leaf_class);
    check(classes[static_cast<int>(Map::LeafClass::Certified)] && classes[static_cast<int>(Map::LeafClass::Mixed)] &&
              classes[static_cast<int>(Map::LeafClass::AllLost)],
          "streamed fixture lacks distinct interpolation/direct source classes");
    const std::string foreign_name = filename + ".foreign";
    {
      TFile output(filename.c_str(), "CREATE");
      check(!output.IsZombie(), "cannot create streamed ROOT fixture");
      auto* directory = output.mkdir("Adaptive");
      TFile foreign(foreign_name.c_str(), "CREATE");
      check(!foreign.IsZombie(), "cannot create foreign current-directory fixture");
      auto* previous = static_cast<TDirectory*>(gDirectory);
      check(Map::Write(directory, original, errors), errors.str());
      check(static_cast<TDirectory*>(gDirectory) == previous, "adaptive writer did not restore caller's current directory");
      check(foreign.Get("Adaptive/Points") == nullptr, "streamed baskets were written into a foreign current file");
      check(PHGarfieldFlushAndCloseRootOutput(output, errors), errors.str());
      check(PHGarfieldFlushAndCloseRootOutput(foreign, errors), errors.str());
    }
    Map::Data readback;
    {
      TFile file(filename.c_str(), "READ");
      check(!file.IsZombie() && !file.TestBit(TFile::kRecovered), "streamed fixture reopened recovered/unreadable");
      auto* directory = file.GetDirectory("Adaptive");
      check(directory != nullptr, "missing streamed output directory");
      for (const char* name : {"Points", "Cells", "Roots"}) check_persisted_multibasket_tree(file, directory, name);
      check(Map::Read(directory, g, readback, errors, true), errors.str());
      check(!Map::Write(directory, original, errors), "read-only adaptive output accepted writes");
      file.Close();
    }
    std::remove(filename.c_str());
    std::remove(foreign_name.c_str());
    check(Map::ContentFingerprint(original) == Map::ContentFingerprint(readback), "streamed primitive/metadata bytes changed after full readback");
    for (const auto& root : original.roots)
    {
      const auto bounds = Map::RootBounds(g, root.root_id);
      const double phi = .5 * (bounds.lo[0] + bounds.hi[0]), radius = .5 * (bounds.lo[1] + bounds.hi[1]), z = .5 * (bounds.lo[2] + bounds.hi[2]);
      const unsigned int side = z > 0.;
      const auto* before = Map::FindLeaf(original, side, phi, radius, z);
      const auto* after = Map::FindLeaf(readback, side, phi, radius, z);
      check(before && after && before->id == after->id && before->leaf_class == after->leaf_class &&
                before->fallback_reason == after->fallback_reason,
            "streaming changed query routing/source classification");
      PHGarfieldForwardMap::NodeValue a, b;
      const bool accepted_before = Map::Interpolate(original, *before, phi, radius, z, a);
      const bool accepted_after = Map::Interpolate(readback, *after, phi, radius, z, b);
      check(accepted_before == accepted_after, "streaming changed interpolation/direct decision");
      if (accepted_before)
        check(a.delta_r_cm == b.delta_r_cm && a.delta_rphi_cm == b.delta_rphi_cm && a.drift_time_ns == b.drift_time_ns,
              "streaming changed exact primitive interpolation values");
    }
    std::cout << "PASS persisted multibasket roundtrip points=" << readback.points.size() << " cells=" << readback.cells.size()
              << " roots=" << readback.roots.size() << " exact primitives and Certified/Mixed/AllLost routing\n";
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
      p.forward_status = static_cast<int>(PHGarfield::ForwardDriftStatus::ReachedReadout);
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
      std::ostringstream log, error;
      const auto bad = [&](const std::array<double, 3>&, Map::Point& p, std::ostream&)
      {++calls;p.reaches=1;p.forward_status=static_cast<int>(PHGarfield::ForwardDriftStatus::ReachedReadout);p.value={{-0.001,0.,0.}};return true; };
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
        check(data.cells[root.cell_id].split_axis == static_cast<int>(axis), "split did not follow actual directional curvature");
      for (const auto& cell : data.cells)
        check(cell.split_axis < 0 || cell.split_axis == static_cast<int>(axis), "unnecessary isotropic/other-axis split");
      for (unsigned int side = 0; side < 2; ++side)
        for (unsigned int i = 0; i < 50; ++i)
        {
          const double phi = (i + .37) * 6.283185307179586 / 50.;
          const double r = 10. + 2. * (i + .13) / 50.;
          const double z = (side ? 1. : -1.) * (1. + 2. * (i + .71) / 50.);
          const auto* cell = Map::FindLeaf(data, side, phi, r, z);
          check(cell && cell->leaf_class == static_cast<int>(Map::LeafClass::Certified), "analytic smooth sample not certified");
          PHGarfieldForwardMap::NodeValue prediction;
          check(Map::Interpolate(data, *cell, phi, r, z, prediction), "analytic interpolation failed");
          Map::Point actual;
          std::ostringstream error;
          sampler(axis)({{phi, r, z}}, actual, error);
          check(1.e4 * std::abs(prediction.delta_r_cm - actual.value[0]) <= 10.01, "off-stencil transverse error exceeds guard");
          check(1.e4 * r * std::abs(prediction.delta_rphi_cm / r - actual.value[1] / r) <= 10.01, "off-stencil angular error exceeds guard");
          check(1.e4 * .1 * std::abs(prediction.drift_time_ns - actual.value[2]) <= 10.01, "off-stencil time error exceeds guard");
        }
      check(Map::ContentFingerprint(data) == Map::ContentFingerprint(build(axis)), "nondeterministic generation");
    }
    auto radial = build(1);
    const auto& root = radial.cells[radial.roots[0].cell_id];
    const double face = .5 * (root.lo[1] + root.hi[1]);
    const auto* chosen = Map::FindLeaf(radial, 0, .3, face, -2.);
    check(chosen && chosen->lo[1] == face, "split face equality did not select upper child");
    const auto* first = Map::FindLeaf(radial, 0, 0., 10.5, -2.);
    const auto* periodic = Map::FindLeaf(radial, 0, 6.283185307179586, 10.5, -2.);
    check(first && periodic && first->id == periodic->id, "phi2pi was not canonicalphi0");
    check(Map::CanonicalPhi(-1.e-30) < 6.283185307179586, "negative tinyphi rounded outside canonical domain");
    const auto* seam = Map::FindLeaf(radial, 0, std::nextafter(6.283185307179586, 0.), 10.5, -2.);
    check(seam && seam->root_id != first->root_id, "near2pi did not select last phi cell");

    Map::Policy boundary_policy;
    boundary_policy.max_depth = 10;
    boundary_policy.max_nodes_per_root = 511;
    boundary_policy.max_queries_per_root = 8192;
    auto boundary = build(0, boundary_policy, true);
    for (const auto& r : boundary.roots)
      check(boundary.cells[r.cell_id].split_axis == 1, "reach transitions did not select radial refinement");
    const auto* interior = Map::FindLeaf(boundary, 0, .3, 10.5, -2.);
    const auto* edge = Map::FindLeaf(boundary, 0, .3, 11.1, -2.);
    check(interior && interior->leaf_class == static_cast<int>(Map::LeafClass::Certified), "mixed root failed to recover reachable interior");
    check(edge && edge->leaf_class == static_cast<int>(Map::LeafClass::Mixed) && edge->hi[1] - edge->lo[1] <= 2. / 64., "physical boundary did not shrink with anisotropic refinement");

    Map::Policy corner_policy;
    corner_policy.max_nodes_per_root = 5;
    Map::Data corner;
    Map::BuildSummary corner_summary;
    std::ostringstream corner_log, corner_errors;
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
        if (c.root_id == r.root_id && c.split_axis == -1 && c.leaf_class == static_cast<int>(Map::LeafClass::Certified)) certified_volume += volume(c);
      check(certified_volume / volume(parent) > .52, "volume-first refinement did not recover the smooth outer half under a five-node budget");
      Map::Data shard;
      Map::BuildSummary shard_summary;
      check(Map::Build(g, corner_policy, r.root_id, r.root_id + 1, localized_corner_sampler(), shard, shard_summary, corner_log, corner_errors), corner_errors.str());
      std::map<std::array<double, 3>, std::array<double, 3>> full_samples;
      for (const auto& p : corner.points) full_samples.emplace(p.coordinate, p.value);
      for (const auto& p : shard.points)
        check(full_samples.count(p.coordinate) && full_samples.at(p.coordinate) == p.value, "sharding changed the sampled physical response");
      for (unsigned int i = 0; i < 20; ++i)
      {
        const double phi = parent.lo[0] + (parent.hi[0] - parent.lo[0]) * (i + .31) / 20.;
        const double radius = parent.lo[1] + (parent.hi[1] - parent.lo[1]) * (i + .47) / 20.;
        const double z = parent.lo[2] + (parent.hi[2] - parent.lo[2]) * (i + .59) / 20.;
        const unsigned int side = z > 0.;
        const auto* a = Map::FindLeaf(corner, side, phi, radius, z);
        const auto* b = Map::FindLeaf(shard, side, phi, radius, z);
        check(a && b && a->lo == b->lo && a->hi == b->hi && a->leaf_class == b->leaf_class && a->fallback_reason == b->fallback_reason, "volume-priority routing depends on shard boundaries");
      }
    }
    check(Map::ContentFingerprint(corner) == [&]()
              {Map::Data repeated;Map::BuildSummary repeated_summary;check(Map::Build(g,corner_policy,0,Map::RootCount(g),localized_corner_sampler(),repeated,repeated_summary,corner_log,corner_errors),corner_errors.str());return Map::ContentFingerprint(repeated); }(),
          "equal-volume tie ordering is not deterministic");

    Map::Policy query_budget;
    query_budget.max_queries_per_root = 4;
    const auto limited = build(1, query_budget);
    for (const auto& c : limited.cells)
      check(c.leaf_class == static_cast<int>(Map::LeafClass::Uncertified) && c.fallback_reason == static_cast<int>(Map::FallbackReason::QueryBudget), "query-budget leaf was falsely certified");
    Map::Policy node_budget;
    node_budget.max_nodes_per_root = 1;
    const auto unsplit = build(1, node_budget);
    for (const auto& c : unsplit.cells)
      check(c.fallback_reason == static_cast<int>(Map::FallbackReason::NodeBudget), "node budget not preserved explicitly");
    std::ostringstream errors, log;
    Map::Data failed;
    Map::BuildSummary summary;
    Map::Policy allocation;
    allocation.max_total_points = 1;
    check(!Map::Build(g, allocation, 0, Map::RootCount(g), sampler(1), failed, summary, log, errors), "global allocation limit was ignored");
    std::uint64_t bytes = 0;
    check(!Map::EstimatedBytes(std::numeric_limits<std::uint64_t>::max(), 1, 1, bytes), "allocation arithmetic overflow was not caught");
    Map::Policy capacity;
    capacity.max_total_points = 600000000;
    capacity.max_total_cells = 100000000;
    capacity.max_estimated_bytes = 320ULL * 1024 * 1024 * 1024;
    check(Map::Policy{}.max_estimated_bytes == 2ULL * 1024 * 1024 * 1024 && Map::Policy{}.max_total_points == 5000000 && Map::Policy{}.max_total_cells == 2000000, "raising the engineering ceiling changed nominal default budgets");
    check(Map::EstimatedBytes(capacity.max_total_points, capacity.max_total_cells, Map::RootCount(Map::Geometry{}), bytes) &&
              bytes > 128ULL * 1024 * 1024 * 1024 && bytes < capacity.max_estimated_bytes,
          "declared large-capacity footprint does not fit its explicit ceiling");
    Map::Data tiny_capacity;
    Map::BuildSummary capacity_summary;
    check(Map::Build(g, capacity, 0, 1, sampler(1), tiny_capacity, capacity_summary, log, errors), errors.str());
    check(tiny_capacity.points.size() < 1000 && tiny_capacity.cells.size() < 1000 &&
              tiny_capacity.points.capacity() < 1000 && tiny_capacity.cells.capacity() < 1000,
          "declared global capacity caused eager allocation/growth");
    auto aggregate_capacity = capacity;
    aggregate_capacity.max_total_points = 1500000000ULL;
    aggregate_capacity.max_total_cells = 250000000ULL;
    aggregate_capacity.max_estimated_bytes = 768ULL * 1024 * 1024 * 1024;
    check(Map::EstimatedBytes(aggregate_capacity.max_total_points, aggregate_capacity.max_total_cells, Map::RootCount(Map::Geometry{}), bytes) &&
              bytes > 512ULL * 1024 * 1024 * 1024 && bytes < aggregate_capacity.max_estimated_bytes,
          "explicit aggregate footprint does not fit768GiB");
    check(Map::Build(g, aggregate_capacity, 0, 1, sampler(1), tiny_capacity, capacity_summary, log, errors), "explicit1.5B/250M/768GiB policy rejected");
    check(tiny_capacity.points.capacity() < 1000 && tiny_capacity.cells.capacity() < 1000, "aggregate policy caused eager allocation");
    auto previous_capacity = capacity;
    previous_capacity.max_total_points = 200000000;
    previous_capacity.max_total_cells = 40000000;
    previous_capacity.max_estimated_bytes = 128ULL * 1024 * 1024 * 1024;
    check(Map::Build(g, previous_capacity, 0, 1, sampler(1), tiny_capacity, capacity_summary, log, errors), "previous128GiB policy stopped working");
    previous_capacity.max_total_points = 100000000;
    previous_capacity.max_total_cells = 20000000;
    previous_capacity.max_estimated_bytes = 64ULL * 1024 * 1024 * 1024;
    check(Map::Build(g, previous_capacity, 0, 1, sampler(1), tiny_capacity, capacity_summary, log, errors), "previous64GiB policy stopped working");
    auto hard_capacity = capacity;
    hard_capacity.max_estimated_bytes = Map::HardMaximumEstimatedBytes;
    check(Map::Build(g, hard_capacity, 0, 1, sampler(1), tiny_capacity, capacity_summary, log, errors), "engineering ceiling boundary was rejected");
    check(tiny_capacity.points.capacity() < 1000 && tiny_capacity.cells.capacity() < 1000, "engineering ceiling boundary eagerly reserved capacity");
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

    const std::string filename = "adaptive-map-test-" + std::to_string(getpid()) + ".root";
    {
      TFile file(filename.c_str(), "CREATE");
      check(!file.IsZombie(), "cannot create analytic ROOT fixture");
      check(Map::Write(file.mkdir("Adaptive"), radial, errors), errors.str());
      // Map::Write streams and destroys its scoped trees before returning.
      // A zero final file-level Write is valid and must not discard the file.
      check(file.Write() == 0, "scoped-tree fixture did not exercise zero final Write");
      check(PHGarfieldFlushAndCloseRootOutput(file, errors), "zero final Write was treated as an I/O failure");
    }
    Map::Data readback;
    {
      TFile file(filename.c_str(), "READ");
      check(Map::Read(file.GetDirectory("Adaptive"), g, readback, errors, true), errors.str());
      file.Close();
    }
    std::remove(filename.c_str());
    check(Map::ContentFingerprint(readback) == Map::ContentFingerprint(radial), "primitive ROOT roundtrip changed data");
    check(Map::FindLeaf(readback, 0, .3, 10.6, -2.) != nullptr, "owning data depends on closed/unlinked ROOT input");
    test_streamed_multibasket_roundtrip(filename + ".streamed");
    {
      TFile file(filename.c_str(), "CREATE");
      auto* directory = file.mkdir("Adaptive");
      check(Map::Write(directory, radial, errors), errors.str());
      directory->cd();
      TParameter<Long64_t>("point_count", std::numeric_limits<Long64_t>::max()).Write("point_count", TObject::kOverwrite);
      file.Close();
    }
    {
      TFile file(filename.c_str(), "READ");
      Map::Data rejected;
      check(!Map::Read(file.GetDirectory("Adaptive"), g, rejected, errors, true), "absurd declared allocation accepted");
      file.Close();
    }
    std::remove(filename.c_str());
    // Reproduce the historical QA ownership failure: large branches flushed
    // while attached to one file, then only the TTree header written elsewhere.
    // The reader must reject missing persisted baskets before stale/default
    // values or a positive aggregate GetEntry can be mistaken for full rows.
    const std::string owner_name = filename + ".owner";
    {
      TFile owner(owner_name.c_str(), "CREATE");
      TTree points("Points", "foreign-owned primitive payload");
      ULong64_t id = 0;
      double coordinate[3] = {0., 10.5, -2.}, value[3] = {0., 0., 10.};
      int reaches = 1, status = 1;
      points.Branch("id", &id, "id/l");
      points.Branch("coordinate", coordinate, "coordinate[3]/D");
      points.Branch("reaches", &reaches, "reaches/I");
      points.Branch("forward_status", &status, "forward_status/I");
      points.Branch("value", value, "value[3]/D");
      points.SetBasketSize("*", 256);
      for (id = 0; id < 10000; ++id) points.Fill();
      TFile broken(filename.c_str(), "CREATE");
      auto* directory = broken.mkdir("Adaptive");
      check(Map::Write(directory, radial, errors), errors.str());
      directory->cd();
      points.Write("Points", TObject::kOverwrite);
      TParameter<Long64_t>("point_count", 10000).Write("point_count", TObject::kOverwrite);
      broken.Close();
      points.SetDirectory(nullptr);
      owner.Close();
    }
    {
      TFile file(filename.c_str(), "READ");
      Map::Data rejected;
      std::ostringstream structural;
      check(!Map::Read(file.GetDirectory("Adaptive"), g, rejected, structural, true), "foreign-owned basket omission accepted");
      check(structural.str().find("basket") != std::string::npos, "missing-basket fixture was not rejected structurally");
      file.Close();
    }
    std::remove(filename.c_str());
    std::remove(owner_name.c_str());
    std::cout << "PASS adaptive analytic axes, off-stencil errors, periodicity, exact ties, Mixed refinement, volume-priority budgets, shard invariance, checked capacity limits, malformed trees and primitive ROOT ownership\n";
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "FAIL adaptive map: " << e.what() << '\n';
    return 1;
  }
}
