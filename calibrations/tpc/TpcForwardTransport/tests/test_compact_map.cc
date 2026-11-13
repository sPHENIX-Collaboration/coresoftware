#include <TBasket.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <unistd.h>
#include <algorithm>
#include <bit>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include "../PHGarfieldCompactMap.h"
#include "../PHGarfieldRootOutput.h"
#include "PHGarfieldTestUtils.h"

namespace
{
  using Map = PHGarfieldCompactMap;
  using Source = PHGarfieldAdaptiveMap;
  using PHGarfieldTest::check;
  using PHGarfieldTest::check_persisted_multibasket_tree;
  bool bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
  Source::Geometry geometry()
  {
    Source::Geometry g;
    g.phi_nodes = 4;
    g.radial_nodes = 4;
    g.longitudinal_nodes = 4;
    g.inner_radius_cm = 10.;
    g.outer_radius_cm = 12.;
    g.central_membrane_half_width_cm = 1.;
    g.readout_abs_z_cm = 3.;
    g.nominal_velocity_cm_ns = .1;
    return g;
  }
  Source::Sampler sampler(bool curved = false, bool lost = false)
  {
    return [curved, lost](const std::array<double, 3>& x, Source::Point& p, std::ostream&)
    {
      p.reaches = !lost || x[1] < 11.1 || std::abs(x[2]) == 3.;
      p.forward_status = static_cast<int>(p.reaches ? PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout : PHGarfieldForwardTransport::ForwardDriftStatus::LeftDriftVolume);
      p.value = {{0., 0., 0.}};
      if (p.reaches)
      {
        const double remaining = 3. - std::abs(x[2]);
        p.value[2] = remaining / .1;
        p.value[1] = .003 * x[1] * remaining;
        if (curved)
        {
          p.value[0] = .02 * (x[1] - 10.) * (12. - x[1]) * remaining / 2.;
        }
      }
      return true;
    };
  }
  Source::Data training(const Source::Policy& policy = Source::Policy{}, bool curved = false, bool lost = false, const Source::Geometry& g = geometry())
  {
    Source::Data source;
    Source::BuildSummary summary;
    std::ostringstream log;
    std::ostringstream error;
    check(Source::Build(g, policy, 0, Source::RootCount(g), sampler(curved, lost), source, summary, log, error), error.str());
    return source;
  }
  constexpr const char* Configuration = "0123456789abcdef";
  Map::Data pack(const Source::Data& source, bool complete = true)
  {
    Map::Data data;
    std::ostringstream error;
    check(Map::Pack(source, Configuration, data, Map::Limits{}, error, complete), error.str());
    return data;
  }
  void invalid(const Map::Data& data, const std::string& label)
  {
    std::ostringstream error;
    check(!Map::Validate(data, Map::Limits{}, error), "corruption accepted: " + label);
  }
  // An independent linear scan of leaf bounds checks runtime tree routing,
  // including the half-open split-plane convention and closed gas boundary.
  const Source::Cell* oracle(const Source::Data& source, unsigned int side, double phi, double r, double z)
  {
    const auto& g = source.geometry;
    if (side > 1 || !std::isfinite(phi) || !std::isfinite(r) || !std::isfinite(z) || r < g.inner_radius_cm || r > g.outer_radius_cm ||
        (side ? z < g.central_membrane_half_width_cm || z > g.readout_abs_z_cm : z < -g.readout_abs_z_cm || z > -g.central_membrane_half_width_cm))
    {
      return nullptr;
    }
    const std::array<double, 3> x{{Source::CanonicalPhi(phi), r, z}};
    const std::array<double, 3> upper{{2. * std::numbers::pi, g.outer_radius_cm, side ? g.readout_abs_z_cm : -g.central_membrane_half_width_cm}};
    const Source::Cell* found = nullptr;
    for (const auto& c : source.cells)
    {
      if (c.split_axis != -1)
      {
        continue;
      }
      bool inside = true;
      for (unsigned int axis = 0; axis < 3; ++axis)
      {
        inside &= x[axis] >= c.lo[axis] && (x[axis] < c.hi[axis] || (x[axis] == upper[axis] && c.hi[axis] == upper[axis]));
      }
      if (inside)
      {
        check(!found, "overlapping oracle leaf bounds");
        found = &c;
      }
    }
    return found;
  }
  void match(const Source::Data& source, const Map::Data& data, unsigned int side, double phi, double r, double z)
  {
    const auto* a = oracle(source, side, phi, r, z);
    Map::Lookup b;
    check(Map::FindLeaf(data, side, phi, r, z, b) == (a != nullptr), "lookup validity changed");
    if (!a)
    {
      return;
    }
    check(a->id == b.cell_id && a->depth == b.depth && a->leaf_class == b.leaf_class && a->fallback_reason == b.fallback_reason, "routing/class/depth changed");
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
      check(bits(a->lo[axis], b.bounds.lo[axis]) && bits(a->hi[axis], b.bounds.hi[axis]), "bounds changed bits");
    }
    PHGarfieldForwardMap::NodeValue expected;
    PHGarfieldForwardMap::NodeValue strict;
    PHGarfieldForwardMap::NodeValue relaxed;
    const bool certified = a->leaf_class == static_cast<int>(Source::LeafClass::Certified);
    const bool usable = Source::HasUsableCorners(source, *a);
    bool flag = false;
    check(Map::Interpolate(data, b, phi, r, z, strict) == certified, "strict policy changed");
    check(Map::HasUsableCorners(data, b) == usable && Map::InterpolateRelaxed(data, b, phi, r, z, relaxed, flag) == usable && flag == certified,
          "relaxed policy promoted certification or invented corners");
    if (usable)
    {
      check(Source::InterpolateCorners({a->lo, a->hi}, {{Source::CanonicalPhi(phi), r, z}}, [&](unsigned int i)
                                       { return &source.points[a->corner_id[i]].value; }, expected),
            "source corners rejected");
      check(bits(expected.delta_r_cm, relaxed.delta_r_cm) && bits(expected.delta_rphi_cm, relaxed.delta_rphi_cm) && bits(expected.drift_time_ns, relaxed.drift_time_ns), "interpolation bits changed");
    }
  }
  void population(const Source::Data& source, const Map::Data& data)
  {
    std::uint64_t direct = 0;
    std::vector<std::uint32_t> mapping(source.points.size(), Map::Invalid);
    std::vector<std::uint32_t> reverse(data.values.size(), Map::Invalid);
    for (const auto& c : source.cells)
    {
      const auto& n = data.nodes.at(c.id);
      if (c.split_axis >= 0)
      {
        check(n.tag == static_cast<unsigned int>(c.split_axis) && n.first == c.child_id[0] && n.second == c.child_id[1], "split routing changed");
      }
      else
      {
        const bool certified = c.leaf_class == static_cast<int>(Source::LeafClass::Certified);
        const bool retain = Source::HasUsableCorners(source, c);
        direct += !certified;
        check(n.tag == (3U | (static_cast<unsigned int>(c.leaf_class) << 2U) | (static_cast<unsigned int>(c.fallback_reason) << 5U)) && n.second == 0,
              "leaf class/reason changed");
        check((n.first != Map::Invalid) == retain, "wrong leaf corner retention");
        if (retain)
        {
          for (unsigned int i = 0; i < 8; ++i)
          {
            const auto old_id = c.corner_id[i];
            const auto new_id = data.corners.at(n.first)[i];
            if (mapping.at(old_id) == Map::Invalid)
            {
              mapping[old_id] = new_id;
            }
            if (reverse.at(new_id) == Map::Invalid)
            {
              reverse[new_id] = old_id;
            }
            check(mapping[old_id] == new_id && reverse[new_id] == old_id, "point identity deduplication changed");
            for (unsigned int component = 0; component < 3; ++component)
            {
              check(bits(source.points[old_id].value[component], data.values[new_id][component]), "retained corner changed double bits");
            }
          }
        }
      }
      const double phi = .5 * (c.lo[0] + c.hi[0]);
      const double r = .5 * (c.lo[1] + c.hi[1]);
      const double z = .5 * (c.lo[2] + c.hi[2]);
      match(source, data, z > 0., phi, r, z);
      if (c.split_axis >= 0)
      {
        std::array<double, 3> x{{phi, r, z}};
        const double mid = x[c.split_axis];
        for (double v : {std::nextafter(mid, -INFINITY), mid, std::nextafter(mid, INFINITY)})
        {
          x[c.split_axis] = v;
          match(source, data, z > 0., x[0], x[1], x[2]);
        }
      }
    }
    check(Map::CountNonCertifiedLeaves(data) == direct && std::all_of(reverse.begin(), reverse.end(), [](auto v)
                                                                      { return v != Map::Invalid; }),
          "population/deduplication count differs");
    std::mt19937_64 rng(2026101809);
    std::uniform_real_distribution<double> u(0., 1.);
    const auto& g = source.geometry;
    for (unsigned int i = 0; i < 300; ++i)
    {
      match(source, data, i % 2, (8. * u(rng) - 4.) * std::numbers::pi, g.inner_radius_cm + (g.outer_radius_cm - g.inner_radius_cm) * u(rng),
            (i % 2 ? 1. : -1.) * (g.central_membrane_half_width_cm + (g.readout_abs_z_cm - g.central_membrane_half_width_cm) * u(rng)));
    }
    for (unsigned int side = 0; side < 2; ++side)
    {
      for (double phi : {-2. * std::numbers::pi, -0., 0., std::nextafter(2. * std::numbers::pi, 0.), 2. * std::numbers::pi, 4. * std::numbers::pi})
      {
        for (double r : {g.inner_radius_cm, g.outer_radius_cm})
        {
          for (double z : {g.central_membrane_half_width_cm, g.readout_abs_z_cm})
          {
            match(source, data, side, phi, r, (side ? 1. : -1.) * z);
          }
        }
      }
      match(source, data, side, 0., g.inner_radius_cm - .01, (side ? 1. : -1.) * 2.);
      match(source, data, side, 0., 11., 0.);
      match(source, data, side, NAN, 11., (side ? 1. : -1.) * 2.);
    }
  }
  void write(const std::filesystem::path& path, const Map::Data& data)
  {
    std::ostringstream error;
    TFile file(path.c_str(), "CREATE");
    check(!file.IsZombie() && Map::Write(file.mkdir("Runtime"), data, Map::Limits{}, error), error.str());
    // Every primitive tree is scoped and explicitly written: final Write may
    // validly write no objects. Reopening below verifies the actual payload.
    check(file.Write() == 0, "fixture no longer exercises zero final Write");
    check(PHGarfieldFlushAndCloseRootOutput(file, error), error.str());
  }
  Map::Data read(const std::filesystem::path& path, const Map::Data& original, bool complete = true)
  {
    TFile file(path.c_str(), "READ");
    check(!file.IsZombie() && !file.TestBit(TFile::kRecovered), "fixture unreadable/recovered");
    Map::Data d;
    std::ostringstream error;
    check(Map::Read(file.GetDirectory("Runtime"), original.geometry, d, Map::Limits{}, error, complete), error.str());
    return d;
  }
  void retention()
  {
    Source::Policy policy;
    policy.max_depth = 0;
    policy.max_nodes_per_root = 1;
    policy.guard_um = .001;
    auto source = training(policy, true);
    population(source, pack(source));
    const auto target = std::find_if(source.cells.begin(), source.cells.end(), [](const auto& c)
                                     { return c.leaf_class == static_cast<int>(Source::LeafClass::NumericallyUnresolved); });
    check(target != source.cells.end(), "no numerical retention fixture");
    const auto id = target->id;
    auto unknown = source;
    unknown.cells[id].leaf_class = static_cast<int>(Source::LeafClass::Uncertified);
    unknown.cells[id].fallback_reason = static_cast<int>(Source::FallbackReason::QueryBudget);
    unknown.cells[id].angular_valid = 0;
    check(pack(unknown).nodes[id].first != Map::Invalid, "complete unverified corners dropped");
    unknown.cells[id].corner_id[7] = -1;
    check(pack(unknown).nodes[id].first == Map::Invalid, "partial corners retained");
    auto angular = source;
    angular.cells[id].angular_valid = 0;
    check(pack(angular).nodes[id].first == Map::Invalid, "known angular-invalid cell retained");
    angular.cells[id].fallback_reason = static_cast<int>(Source::FallbackReason::AngularBranch);
    check(pack(angular).nodes[id].first == Map::Invalid, "angular-branch cell retained");
  }
  void capacity(const Source::Data& source, const Map::Data& data)
  {
    std::ostringstream error;
    std::uint64_t bytes = 0;
    check(!Map::EstimatedBytes(UINT64_MAX, 0, 0, 0, bytes), "size overflow accepted");
    using Bound = std::pair<std::uint64_t Map::Limits::*, std::uint64_t>;
    for (const auto& [member, maximum] : std::array<Bound, 5>{{{&Map::Limits::max_nodes, Map::HardMaximumNodes}, {&Map::Limits::max_values, Map::HardMaximumValues}, {&Map::Limits::max_corners, Map::HardMaximumCorners}, {&Map::Limits::max_roots, Map::HardMaximumRoots}, {&Map::Limits::max_bytes, Map::HardMaximumBytes}}})
    {
      for (auto v : std::array<std::uint64_t, 2>{{0, maximum + 1}})
      {
        Map::Limits limits;
        limits.*member = v;
        Map::Data rejected;
        check(!Map::ValidateLimits(limits, error) && !Map::Pack(source, Configuration, rejected, limits, error) && rejected.nodes.capacity() == 0 &&
                  !Map::Read(nullptr, source.geometry, rejected, limits, error) && rejected.nodes.capacity() == 0,
              "invalid limits accepted or allocated");
      }
    }
    Map::Limits tiny;
    tiny.max_nodes = 1;
    Map::Data rejected;
    check(!Map::Pack(source, Configuration, rejected, tiny, error) && !Map::Validate(data, tiny, error), "node budget ignored");
    Map::Limits wide{Map::HardMaximumNodes, Map::HardMaximumValues, Map::HardMaximumCorners, Map::HardMaximumRoots, Map::HardMaximumBytes};
    Map::Data compact;
    check(Map::Pack(source, Configuration, compact, wide, error) && Map::ContentFingerprint(compact) == Map::ContentFingerprint(data), "explicit wide limits change output");
    check(compact.nodes.capacity() < 10000 && compact.values.capacity() < 10000, "limits eagerly reserved arrays");
  }
  void shards()
  {
    std::array<Map::Data, 2> pieces;
    const auto g = geometry();
    const auto middle = Source::RootCount(g) / 2;
    std::ostringstream log;
    std::ostringstream error;
    for (unsigned int i = 0; i < 2; ++i)
    {
      Source::Data source;
      Source::BuildSummary summary;
      check(Source::Build(g, Source::Policy{}, i ? middle : 0, i ? Source::RootCount(g) : middle, sampler(true), source, summary, log, error), error.str());
      pieces[i] = pack(source, false);
    }
    Map::Data merged;
    for (const auto& piece : pieces)
    {
      check(Map::Append(merged, piece, Map::Limits{}, error), error.str());
    }
    check(Map::Validate(merged, Map::Limits{}, error), error.str());
    const auto node_offset = pieces[0].nodes.size();
    const auto corner_offset = pieces[0].corners.size();
    const auto value_offset = pieces[0].values.size();
    for (unsigned int side = 0; side < 2; ++side)
    {
      const auto& original = pieces[side];
      for (double phi : {.01, 2., 6.1})
      {
        Map::Lookup before;
        Map::Lookup after;
        PHGarfieldForwardMap::NodeValue a;
        PHGarfieldForwardMap::NodeValue b;
        bool ca;
        bool cb;
        check(Map::FindLeaf(original, side, phi, 11., side ? 2. : -2., before) && Map::FindLeaf(merged, side, phi, 11., side ? 2. : -2., after), "merged route missing");
        check(after.cell_id == before.cell_id + (side ? node_offset : 0) && after.corner_record == before.corner_record + (side ? corner_offset : 0), "merged route offsets differ");
        check(Map::InterpolateRelaxed(original, before, phi, 11., side ? 2. : -2., a, ca) && Map::InterpolateRelaxed(merged, after, phi, 11., side ? 2. : -2., b, cb) && ca == cb && bits(a.delta_r_cm, b.delta_r_cm) && bits(a.drift_time_ns, b.drift_time_ns), "merged query changed");
      }
    }
    check(std::equal(pieces[1].values.begin(), pieces[1].values.end(), merged.values.begin() + value_offset), "merged corner values differ");
    for (auto member : {&Map::Data::canonical_policy, &Map::Data::configuration_fingerprint_fnv1a64})
    {
      auto changed = pieces[1];
      changed.*member += "0";
      auto rejected = pieces[0];
      check(!Map::Append(rejected, changed, Map::Limits{}, error), "incompatible shard accepted");
    }
    check(!Map::Append(merged, pieces[1], Map::Limits{}, error) && !Map::Append(merged, merged, Map::Limits{}, error), "overlap/self append accepted");
    Map::Limits tiny;
    tiny.max_nodes = merged.nodes.size() - 1;
    auto rejected = pieces[0];
    check(!Map::Append(rejected, pieces[1], tiny, error) && rejected.nodes.size() == pieces[0].nodes.size(), "merge allocation budget ignored");
  }
}  // namespace

int main()
{
  char pattern[] = "/tmp/phgarfield-compact-test-XXXXXX";
  const auto* tmp = mkdtemp(pattern);
  if (!tmp)
  {
    return 1;
  }
  const std::filesystem::path directory = tmp;
  try
  {
    retention();
    shards();
    auto source = training(Source::Policy{}, true);
    auto data = pack(source);
    population(source, data);
    capacity(source, data);
    check(data.values.size() < source.points.size(), "construction-only samples retained");
    std::array<bool, 5> classes{};
    for (auto settings : std::array<std::array<unsigned int, 3>, 5>{{{{127, 12, 2048}}, {{1, 0, 2048}}, {{127, 0, 2048}}, {{127, 12, 1}}, {{3, 12, 2048}}}})
    {
      Source::Policy policy;
      policy.max_nodes_per_root = settings[0];
      policy.max_depth = settings[1];
      policy.max_queries_per_root = settings[2];
      auto fixture = training(policy, settings[2] != 1, true);
      population(fixture, pack(fixture));
      for (const auto& c : fixture.cells)
      {
        if (c.split_axis == -1)
        {
          classes[c.leaf_class] = true;
        }
      }
    }
    Source::Policy budget;
    budget.max_nodes_per_root = 1;
    budget.max_depth = 0;
    auto numeric = training(budget, true);
    population(numeric, pack(numeric));
    for (const auto& c : numeric.cells)
    {
      if (c.split_axis == -1)
      {
        classes[c.leaf_class] = true;
      }
    }
    check(std::all_of(classes.begin(), classes.end(), [](bool x)
                      { return x; }),
          "not every leaf class exercised");
    for (bool outer : {false, true})
    {
      const auto g = geometry();
      Source::Data boundary_source;
      Source::BuildSummary summary;
      std::ostringstream log;
      std::ostringstream error;
      check(Source::Build(g, Source::Policy{}, 0, Source::RootCount(g), [&](const auto& x, Source::Point& point, std::ostream&)
                          {
        const double remaining = g.readout_abs_z_cm - std::abs(x[2]);
        const double endpoint = outer ? g.outer_radius_cm + .5 * Source::EndpointRadialToleranceCm : g.inner_radius_cm - .5 * Source::EndpointRadialToleranceCm;
        point.reaches = 1; point.forward_status = static_cast<int>(PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout);
        point.value = {{(endpoint - x[1]) * remaining / (g.readout_abs_z_cm - g.central_membrane_half_width_cm), 0., remaining / g.nominal_velocity_cm_ns}};
        return true; }, boundary_source, summary, log, error),
            error.str());
      population(boundary_source, pack(boundary_source));
    }
    const auto internal = std::find_if(data.nodes.begin(), data.nodes.end(), [](const auto& n)
                                       { return (n.tag & 3U) < 3; });
    check(internal != data.nodes.end(), "no internal node fixture");
    const auto index = internal - data.nodes.begin();
    using Mutator = std::function<void(Map::Data&)>;
    const std::vector<std::pair<std::string, Mutator>> mutations{
        {"missing root", [](auto& d)
         { d.root_lookup[0] = Map::Invalid; }},
        {"cycle", [index](auto& d)
         { d.nodes[index].first = index; }},
        {"overlap", [index](auto& d)
         { d.nodes[index].second = d.nodes[index].first; }},
        {"reserved bits", [index](auto& d)
         { d.nodes[index].tag |= 256U; }},
        {"dangling corner", [](auto& d)
         { d.corners[0][0] = Map::Invalid; }},
        {"nonfinite", [](auto& d)
         { d.values[0][2] = NAN; }},
        {"orphan node", [](auto& d)
         { d.nodes.push_back(d.nodes.back()); }},
        {"orphan value", [](auto& d)
         { d.values.push_back(d.values.back()); }},
        {"wrong physical endpoint", [](auto& d)
         { d.values[0][0] = 100.; }}};
    for (const auto& [label, mutate] : mutations)
    {
      auto bad = data;
      mutate(bad);
      invalid(bad, label);
    }
    std::ostringstream errors;
    const auto root = directory / "runtime.root";
    write(root, data);
    auto owned = read(root, data);
    check(Map::ContentFingerprint(owned) == Map::ContentFingerprint(data), "ROOT roundtrip differs");
    std::filesystem::remove(root);
    population(source, owned);  // No ROOT object lifetime dependency.
    for (const auto& [key, value] : std::array<std::pair<const char*, Long64_t>, 2>{{{"node_count", INT64_MAX}, {"schema_version", 8}}})
    {
      const auto path = directory / (std::string(key) + ".root");
      write(path, data);
      {
        TFile output(path.c_str(), "UPDATE");
        output.cd("Runtime");
        if (std::string(key) == "schema_version")
        {
          TParameter<int>(key, value).Write(key, TObject::kOverwrite);
        }
        else
        {
          TParameter<Long64_t>(key, value).Write(key, TObject::kOverwrite);
        }
      }
      TFile input(path.c_str(), "READ");
      Map::Data rejected;
      check(!Map::Read(input.GetDirectory("Runtime"), data.geometry, rejected, Map::Limits{}, errors) && rejected.nodes.empty(), "invalid schema/count allocated or accepted");
    }
    const auto changed_path = directory / "changed.root";
    write(changed_path, data);
    {
      TFile output(changed_path.c_str(), "UPDATE");
      output.cd("Runtime");
      TNamed("content_fingerprint_fnv1a64", "0000000000000000").Write("content_fingerprint_fnv1a64", TObject::kOverwrite);
    }
    {
      TFile input(changed_path.c_str(), "READ");
      Map::Data rejected;
      check(!Map::Read(input.GetDirectory("Runtime"), data.geometry, rejected, Map::Limits{}, errors), "changed fingerprint accepted");
    }
    auto large_geometry = geometry();
    large_geometry.phi_nodes = 24;
    large_geometry.radial_nodes = 20;
    large_geometry.longitudinal_nodes = 20;
    auto large_source = training(Source::Policy{}, false, false, large_geometry);
    auto large = pack(large_source);
    const auto streamed = directory / "streamed.root";
    write(streamed, large);
    {
      TFile input(streamed.c_str(), "READ");
      for (const char* name : {"Nodes", "Corners", "Values", "Roots"})
      {
        check_persisted_multibasket_tree(input, input.GetDirectory("Runtime"), name);
      }
    }
    check(Map::ContentFingerprint(read(streamed, large)) == Map::ContentFingerprint(large), "multibasket roundtrip changed");
    const auto damaged = directory / "missing-baskets.root";
    const auto foreign_path = directory / "foreign.root";
    write(damaged, large);
    {
      TFile foreign(foreign_path.c_str(), "CREATE");
      TNamed("marker", "foreign fixture").Write();
    }
    const auto destination_path = directory / "foreign-current-writer.root";
    {
      TFile output(destination_path.c_str(), "CREATE");
      auto* destination = output.mkdir("Runtime");
      TFile foreign(foreign_path.c_str(), "READ");
      foreign.cd();
      check(Map::Write(destination, data, Map::Limits{}, errors) && gDirectory == &foreign, "writer used or changed foreign current directory");
      check(PHGarfieldFlushAndCloseRootOutput(output, errors), errors.str());
    }
    check(Map::ContentFingerprint(read(destination_path, data)) == Map::ContentFingerprint(data), "foreign-current destination payload differs");
    {
      TFile foreign(foreign_path.c_str(), "READ");
      foreign.cd();
      TTree tree("Nodes", "foreign readonly lost baskets");
      Map::Node node;
      tree.Branch("first", &node.first, "first/i", 256);
      tree.Branch("second", &node.second, "second/i", 10000000);
      tree.Branch("tag", &node.tag, "tag/i", 10000000);
      const auto saved = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kFatal;
      for (const auto& n : large.nodes)
      {
        node = n;
        tree.Fill();
      }
      TFile output(damaged.c_str(), "UPDATE");
      output.cd("Runtime");
      check(tree.Write("Nodes", TObject::kOverwrite) > 0, "cannot write damaged fixture");
      tree.SetDirectory(nullptr);
      gErrorIgnoreLevel = saved;
    }
    {
      TFile input(damaged.c_str(), "READ");
      auto* tree = dynamic_cast<TTree*>(input.Get("Runtime/Nodes"));
      check(tree && tree->GetBranch("first")->GetWriteBasket() > 1, "damaged fixture lacked lost baskets");
      Map::Data rejected;
      std::ostringstream error;
      check(!Map::Read(input.GetDirectory("Runtime"), large.geometry, rejected, Map::Limits{}, error) && error.str().find("basket") != std::string::npos, "missing persisted baskets accepted");
    }
    std::cout << "PASS direct compact packing/routing/allclasses/retention/shards/limits/corruption/streaming/ownership\n";
    std::filesystem::remove_all(directory);
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "FAIL compact fixture preserved at " << directory << ": " << e.what() << '\n';
    return 1;
  }
}
