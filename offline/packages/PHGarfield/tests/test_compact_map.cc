#include <TBasket.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <sys/wait.h>
#include <unistd.h>
#include <algorithm>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include "../PHGarfieldCompactMap.h"
#include "../PHGarfieldRootOutput.h"

namespace
{
  using Map = PHGarfieldCompactMap;
  using Source = PHGarfieldAdaptiveMap;
  void check(bool value, const std::string& text)
  {
    if (!value) throw std::runtime_error(text);
  }
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
      p.forward_status = static_cast<int>(p.reaches ? PHGarfield::ForwardDriftStatus::ReachedReadout : PHGarfield::ForwardDriftStatus::LeftDriftVolume);
      p.value = {{0., 0., 0.}};
      if (p.reaches)
      {
        const double remaining = 3. - std::abs(x[2]);
        p.value[2] = remaining / .1;
        if (curved) p.value[0] = .02 * (x[1] - 10.) * (12. - x[1]) * remaining / 2.;
      }
      return true;
    };
  }
  Source::Data training(const Source::Policy& policy = Source::Policy{}, bool curved = false, bool lost = false, Source::Geometry g = geometry())
  {
    Source::Data source;
    Source::BuildSummary summary;
    std::ostringstream log, error;
    check(Source::Build(g, policy, 0, Source::RootCount(g), sampler(curved, lost), source, summary, log, error), error.str());
    return source;
  }
  Map::SourceBinding binding(const Source::Data& s)
  {
    return {Map::Sha256("software-fixture-not-physical-training"), "0123456789abcdef", Source::ContentFingerprint(s), Source::BuildCanonicalPolicy(s), Map::Sha256("software-fixture-converter")};
  }
  Map::ExpectedSource expected(const Map::Data& d)
  {
    return {d.source.file_sha256, d.source.configuration_fingerprint_fnv1a64, d.source.content_fingerprint_fnv1a64};
  }
  Map::Data convert(Source::Data& source)
  {
    Map::Data data;
    std::ostringstream errors;
    check(Map::Convert(source, binding(source), data, Map::Limits{}, errors, true), errors.str());
    return data;
  }
  void reseal(Map::Data& data) { data.proof.sha256 = Map::Sha256(Map::ProofCanonical(data)); }
  void invalid(Map::Data data, const std::string& label)
  {
    reseal(data);
    std::ostringstream error;
    check(!Map::Validate(data, Map::Limits{}, error, true), "corruption accepted: " + label);
  }
  void match(const Source::Data& source, const Map::Data& data, unsigned int side, double phi, double r, double z)
  {
    const auto* a = Source::FindLeaf(source, side, phi, r, z);
    Map::Lookup b;
    check(Map::FindLeaf(data, side, phi, r, z, b) == (a != nullptr), "lookup validity changed");
    if (!a) return;
    check(a->id == b.cell_id && a->depth == b.depth && a->leaf_class == b.leaf_class && a->fallback_reason == b.fallback_reason, "routing/source/depth changed");
    for (unsigned int axis = 0; axis < 3; ++axis) check(bits(a->lo[axis], b.bounds.lo[axis]) && bits(a->hi[axis], b.bounds.hi[axis]), "lookup bounds changed bits");
    PHGarfieldForwardMap::NodeValue x, y;
    const bool before = Source::Interpolate(source, *a, phi, r, z, x), after = Map::Interpolate(data, b, phi, r, z, y);
    check(before == after, "interpolation/direct routing changed");
    if (before) check(x.reaches_readout == y.reaches_readout && bits(x.delta_r_cm, y.delta_r_cm) && bits(x.delta_rphi_cm, y.delta_rphi_cm) && bits(x.drift_time_ns, y.drift_time_ns), "interpolation double bits changed");
  }
  void population(const Source::Data& source, const Map::Data& data)
  {
    std::uint64_t expected_direct = 0;
    for (const auto& cell : source.cells)
    {
      if (cell.split_axis == -1 && cell.leaf_class != static_cast<int>(Source::LeafClass::Certified))
      {
        ++expected_direct;
      }
    }
    check(Map::CountNonCertifiedLeaves(data) == expected_direct,
          "compact fallback-leaf count disagrees with the adaptive source");
    const double pi = std::numbers::pi;
    std::mt19937_64 rng(2026101809);
    std::uniform_real_distribution<double> u(0., 1.);
    for (unsigned int i = 0; i < 2000; ++i)
    {
      const auto side = i % 2;
      match(source, data, side, (8. * u(rng) - 4.) * pi, 10. + 2. * u(rng), (side ? 1. : -1.) * (1. + 2. * u(rng)));
    }
    for (const auto& cell : source.cells)
    {
      const double phi = .5 * (cell.lo[0] + cell.hi[0]), r = .5 * (cell.lo[1] + cell.hi[1]), z = .5 * (cell.lo[2] + cell.hi[2]);
      const unsigned int side = z > 0.;
      match(source, data, side, phi, r, z);
      if (cell.split_axis >= 0)
      {
        std::array<double, 3> x{{phi, r, z}};
        const unsigned int axis = cell.split_axis;
        const double mid = .5 * (cell.lo[axis] + cell.hi[axis]);
        for (const double value : {std::nextafter(mid, -INFINITY), mid, std::nextafter(mid, INFINITY)})
        {
          x[axis] = value;
          match(source, data, side, x[0], x[1], x[2]);
        }
      }
    }
    for (unsigned int side = 0; side < 2; ++side)
    {
      for (const double phi : {-2. * pi, -0., 0., std::nextafter(2. * pi, 0.), 2. * pi, 4. * pi})
        for (const double r : {10., 11., 12.})
          for (const double z : {1., 2., 3.}) match(source, data, side, phi, r, (side ? 1. : -1.) * z);
      match(source, data, side, 0., 9.99, (side ? 1. : -1.) * 2.);
      match(source, data, side, 0., 11., 0.);
      match(source, data, side, NAN, 11., (side ? 1. : -1.) * 2.);
    }
  }
  void write(const std::filesystem::path& path, const Map::Data& data)
  {
    std::ostringstream error;
    TFile file(path.c_str(), "CREATE");
    check(!file.IsZombie(), "cannot create fixture");
    check(Map::Write(file.mkdir("Runtime"), data, Map::Limits{}, error), error.str());
    check(PHGarfieldFlushAndCloseRootOutput(file, error), error.str());
  }
  Map::Data read(const std::filesystem::path& path, const Map::Data& original)
  {
    TFile file(path.c_str(), "READ");
    check(!file.IsZombie() && !file.TestBit(TFile::kRecovered), "fixture unreadable/recovered");
    Map::Data d;
    std::ostringstream errors;
    check(Map::Read(file.GetDirectory("Runtime"), original.geometry, expected(original), d, Map::Limits{}, errors, true), errors.str());
    return d;
  }
  void persisted(TDirectory* directory, const char* name)
  {
    auto* tree = dynamic_cast<TTree*>(directory->Get(name));
    check(tree && tree->GetEntries() > 0, "missing nonempty streamed tree");
    for (auto* object : *tree->GetListOfBranches())
    {
      auto* branch = dynamic_cast<TBranch*>(object);
      const auto last = branch->GetWriteBasket();
      Long64_t covered = 0;
      unsigned int count = 0;
      for (int i = 0; i <= last; ++i)
      {
        const auto begin = branch->GetBasketEntry()[i], end = i < last ? branch->GetBasketEntry()[i + 1] : tree->GetEntries();
        if (begin == end) continue;
        check(begin == covered && end > begin && end <= tree->GetEntries(), "streamed ranges incomplete");
        const auto seek = branch->GetBasketSeek(i);
        const auto bytes = branch->GetBasketBytes()[i];
        check(seek > 0 && bytes > 0 && seek <= directory->GetFile()->GetSize() && bytes <= directory->GetFile()->GetSize() - seek, "basket not persisted in destination");
        covered = end;
        ++count;
      }
      check(covered == tree->GetEntries() && count >= 2, "fixture did not span multiple persisted baskets");
    }
  }
  void write_training(const std::filesystem::path& path, const Source::Data& source)
  {
    std::ostringstream errors;
    TFile file(path.c_str(), "CREATE");
    auto* top = file.mkdir("GarfieldFullTransport");
    auto* m = top->mkdir("Metadata");
    {
      TDirectory::TContext context(m);
      TNamed("schema_name", PHGarfieldForwardMap::SchemaName).Write();
      TParameter<int>("schema_version", 4).Write();
      TNamed("configuration_fingerprint_fnv1a64", "0123456789abcdef").Write();
      TParameter<int>("phi_nodes", source.geometry.phi_nodes).Write();
      TParameter<int>("radial_nodes", source.geometry.radial_nodes).Write();
      TParameter<int>("longitudinal_nodes", source.geometry.longitudinal_nodes).Write();
      TParameter<double>("gas_inner_radius_cm", source.geometry.inner_radius_cm).Write();
      TParameter<double>("gas_outer_radius_cm", source.geometry.outer_radius_cm).Write();
      TParameter<double>("central_membrane_half_width_cm", source.geometry.central_membrane_half_width_cm).Write();
      TParameter<double>("readout_abs_z_cm", source.geometry.readout_abs_z_cm).Write();
      TParameter<double>("nominal_phg4_drift_velocity_cm_per_ns", source.geometry.nominal_velocity_cm_ns).Write();
    }
    check(Source::Write(top->mkdir("Adaptive"), source, errors), errors.str());
    check(PHGarfieldFlushAndCloseRootOutput(file, errors), errors.str());
  }
  int execute(const std::filesystem::path& converter, const std::vector<std::string>& arguments)
  {
    const auto pid = fork();
    check(pid >= 0, "fork failed");
    if (!pid)
    {
      std::vector<char*> args;
      std::string binary = converter;
      args.push_back(binary.data());
      for (const auto& value : arguments) args.push_back(const_cast<char*>(value.c_str()));
      args.push_back(nullptr);
      execv(binary.c_str(), args.data());
      _exit(127);
    }
    int status = 0;
    check(waitpid(pid, &status, 0) == pid, "waitpid failed");
    return WIFEXITED(status) ? WEXITSTATUS(status) : 128;
  }
  Map::Limits wide_limits()
  {
    Map::Limits limits;
    limits.max_source_points = 1500000000;
    limits.max_nodes = 250000000;
    limits.max_values = 1500000000;
    limits.max_corners = 150000000;
    limits.max_bytes = 128ULL * 1024 * 1024 * 1024;
    return limits;
  }
  void capacity(const std::filesystem::path& directory)
  {
    const Map::Limits defaults;
    check(defaults.max_source_points == 600000000 && defaults.max_nodes == 100000000 && defaults.max_values == 600000000 &&
              defaults.max_corners == 50000000 && defaults.max_roots == 1000000 && defaults.max_bytes == 32ULL * 1024 * 1024 * 1024,
          "compact defaults changed");
    const auto limits = wide_limits();
    std::ostringstream errors;
    check(Map::ValidateLimits(defaults, errors) && Map::ValidateLimits(limits, errors), "bounded compact limits rejected");
    std::uint64_t estimated = 0;
    check(Map::EstimatedBytes(limits.max_nodes, limits.max_corners, limits.max_values, defaults.max_roots, estimated) && estimated <= limits.max_bytes,
          "explicit compact ceilings do not fit128GiB conservative byte allowance");
    auto g = geometry();
    g.radial_nodes = 2;
    g.longitudinal_nodes = 2;
    Source::Policy policy;
    policy.max_total_points = 1500000000;
    policy.max_total_cells = 250000000;
    policy.max_estimated_bytes = 768ULL * 1024 * 1024 * 1024;
    auto source = training(policy, false, false, g);
    Map::Data compact;
    check(Map::Convert(source, binding(source), compact, limits, errors, true), "explicit large limits failed tiny conversion");
    check(source.points.size() < 1000 && source.cells.size() < 1000 && source.points.capacity() < 1000 && source.cells.capacity() < 1000,
          "declared source global caps eagerly reserved vectors");
    check(compact.nodes.size() < 1000 && compact.corners.size() < 1000 && compact.values.size() < 1000 && compact.root_lookup.size() < 1000 &&
              compact.nodes.capacity() < 1000 && compact.corners.capacity() < 1000 && compact.values.capacity() < 1000 && compact.root_lookup.capacity() < 1000,
          "declared compact caps eagerly reserved vectors");
    using Bound = std::pair<std::uint64_t Map::Limits::*, std::uint64_t>;
    for (const auto& [member, maximum] : std::array<Bound, 6>{{{&Map::Limits::max_source_points, Map::HardMaximumSourcePoints}, {&Map::Limits::max_nodes, Map::HardMaximumNodes}, {&Map::Limits::max_values, Map::HardMaximumValues}, {&Map::Limits::max_corners, Map::HardMaximumCorners}, {&Map::Limits::max_roots, Map::HardMaximumRoots}, {&Map::Limits::max_bytes, Map::HardMaximumBytes}}})
      for (const auto value : std::array<std::uint64_t, 2>{{0, maximum + 1}})
      {
        auto invalid_limits = limits;
        invalid_limits.*member = value;
        std::ostringstream diagnostic;
        check(!Map::ValidateLimits(invalid_limits, diagnostic), "zero/excess compact limit accepted");
        Map::Data rejected;
        check(!Map::Convert(source, binding(source), rejected, invalid_limits, diagnostic, true) && rejected.nodes.capacity() == 0,
              "invalid compact limits allocated conversion output");
        check(!Map::Read(nullptr, g, expected(compact), rejected, invalid_limits, diagnostic, true) && rejected.nodes.capacity() == 0,
              "invalid compact limits reached ROOT allocation");
      }
    auto source_limit = limits;
    source_limit.max_source_points = source.points.size() - 1;
    Map::Data rejected;
    check(!Map::Convert(source, binding(source), rejected, source_limit, errors, true) && rejected.nodes.empty(), "source point admission limit ignored");
    // Synthetic declaration fixture only: extra source probes need not survive
    // compact conversion. This checks admission/serialization without allocating
    // hundreds of millions of rows. It is explicitly not a conversion proof for
    // this tiny source; VerifyConversion must reject that false claim below.
    auto declared = compact;
    declared.proof.source_points = 700000001;
    reseal(declared);
    check(!Map::Validate(declared, defaults, errors, true) && Map::Validate(declared, limits, errors, true), "source proof cap did not follow explicit limits");
    check(!Map::VerifyConversion(source, declared, errors), "synthetic source-count declaration mistaken for exhaustive conversion proof");
    for (const auto member : std::array<std::uint64_t Map::ConversionProof::*, 4>{{&Map::ConversionProof::source_cells, &Map::ConversionProof::roots_checked, &Map::ConversionProof::bounds_checked, &Map::ConversionProof::routes_checked}})
    {
      auto false_count = compact;
      ++(false_count.proof.*member);
      reseal(false_count);
      check(!Map::VerifyConversion(source, false_count, errors), "re-sealed false source/proof count accepted");
    }
    const auto path = directory / "capacity-declaration.root";
    {
      TFile output(path.c_str(), "CREATE");
      check(Map::Write(output.mkdir("Runtime"), declared, limits, errors), "cannot write source-count admission fixture");
      check(PHGarfieldFlushAndCloseRootOutput(output, errors), errors.str());
    }
    {
      TFile input(path.c_str(), "READ");
      Map::Data parsed;
      check(!Map::Read(input.GetDirectory("Runtime"), g, expected(compact), parsed, defaults, errors, true) && parsed.nodes.empty(), "default reader accepted source proof over600M");
      check(Map::Read(input.GetDirectory("Runtime"), g, expected(compact), parsed, limits, errors, true) && parsed.nodes.capacity() < 1000, "explicit reader cap failed or eagerly allocated declaration");
    }
  }
}  // namespace

int main(int argc, char** argv)
{
  (void) argc;
  char pattern[] = "/tmp/phgarfield-compact-test-XXXXXX";
  const auto* tmp = mkdtemp(pattern);
  if (!tmp) return 1;
  const std::filesystem::path directory = tmp;
  try
  {
    for (bool outer : {false, true})
    {
      const auto g = geometry();
      Source::Data s;
      Source::BuildSummary summary;
      std::ostringstream log, error;
      const auto boundary = [&](const std::array<double, 3>& x, Source::Point& p, std::ostream&)
      {
        const double remaining = g.readout_abs_z_cm - std::abs(x[2]);
        const double fraction = remaining / (g.readout_abs_z_cm - g.central_membrane_half_width_cm);
        const double endpoint = outer ? g.outer_radius_cm + .5 * Source::EndpointRadialToleranceCm : g.inner_radius_cm - .5 * Source::EndpointRadialToleranceCm;
        p.reaches = 1;
        p.forward_status = static_cast<int>(PHGarfield::ForwardDriftStatus::ReachedReadout);
        p.value = {{(endpoint - x[1]) * fraction, 0., remaining / g.nominal_velocity_cm_ns}};
        return true;
      };
      check(Source::Build(g, Source::Policy{}, 0, Source::RootCount(g), boundary, s, summary, log, error), error.str());
      auto compact = convert(s);
      population(s, compact);
      check(Map::VerifyConversion(s, compact, error), error.str());
    }
    check(Map::Sha256("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "SHA256 empty vector differs");
    check(Map::Sha256("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "SHA256 abc vector differs");
    auto source = training(Source::Policy{}, true, false);
    auto data = convert(source);
    population(source, data);
    capacity(directory);
    check(data.values.size() < source.points.size(), "training-only probe points were retained");
    std::array<bool, 5> classes{};
    for (const auto settings : std::array<std::array<int, 3>, 5>{{{{127, 12, 2048}}, {{1, 0, 2048}}, {{127, 0, 2048}}, {{127, 12, 1}}, {{3, 12, 2048}}}})
    {
      Source::Policy p;
      p.max_nodes_per_root = settings[0];
      p.max_depth = settings[1];
      p.max_queries_per_root = settings[2];
      auto fixture = training(p, settings[2] != 1, true);
      auto compact = convert(fixture);
      population(fixture, compact);
      for (const auto& c : fixture.cells)
        if (c.split_axis == -1) classes[c.leaf_class] = true;
    }
    Source::Policy p;
    p.max_nodes_per_root = 1;
    p.max_depth = 0;
    auto numeric = training(p, true, false);
    auto numeric_data = convert(numeric);
    population(numeric, numeric_data);
    for (const auto& c : numeric.cells)
      if (c.split_axis == -1) classes[c.leaf_class] = true;
    check(std::all_of(classes.begin(), classes.end(), [](bool x)
                      { return x; }),
          "fixtures did not exercise every leaf class");
    auto bad = data;
    bad.root_lookup[0] = Map::Invalid;
    invalid(bad, "missing root");
    const auto internal = std::find_if(data.nodes.begin(), data.nodes.end(), [](const auto& n)
                                       { return (n.tag & 3U) < 3; });
    check(internal != data.nodes.end(), "no anisotropic split fixture");
    const auto index = internal - data.nodes.begin();
    bad = data;
    bad.nodes[index].first = index;
    invalid(bad, "cycle");
    bad = data;
    bad.nodes[index].second = bad.nodes[index].first;
    invalid(bad, "overlap");
    bad = data;
    bad.nodes[index].tag |= 256;
    invalid(bad, "reserved bits");
    bad = data;
    bad.corners[0][0] = Map::Invalid;
    invalid(bad, "dangling corner");
    bad = data;
    bad.values[0][2] = NAN;
    invalid(bad, "nonfinite value");
    bad = data;
    bad.nodes.push_back(data.nodes.back());
    invalid(bad, "orphan node/proof counts");
    bad = data;
    bad.values[0][2] = std::nextafter(bad.values[0][2], INFINITY);
    reseal(bad);
    std::ostringstream errors;
    check(!Map::VerifyConversion(source, bad, errors), "changed corner bits accepted by exhaustive proof");
    auto incomplete_source = source;
    incomplete_source.roots.pop_back();
    auto incomplete_compact = data;
    incomplete_compact.source.content_fingerprint_fnv1a64 = Source::ContentFingerprint(incomplete_source);
    reseal(incomplete_compact);
    check(!Map::VerifyConversion(incomplete_source, incomplete_compact, errors), "incomplete source roots accepted by exhaustive proof");
    std::uint64_t bytes = 0;
    check(!Map::EstimatedBytes(UINT64_MAX, 0, 0, 0, bytes), "size overflow accepted");
    auto limit = Map::Limits{};
    limit.max_nodes = 1;
    check(!Map::Validate(data, limit, errors, true), "node allocation limit ignored");
    const auto root = directory / "runtime.root";
    write(root, data);
    auto owned = read(root, data);
    check(Map::ContentFingerprint(data) == Map::ContentFingerprint(owned), "ROOT roundtrip differs");
    {
      TFile file(root.c_str(), "READ");
      auto expected_bad = expected(data);
      expected_bad.file_sha256 = Map::Sha256("wrong");
      Map::Data rejected;
      check(!Map::Read(file.GetDirectory("Runtime"), data.geometry, expected_bad, rejected, Map::Limits{}, errors, true), "wrong external source accepted");
    }
    const auto authentic_sha = Map::FileSha256(root);
    std::filesystem::remove(root);
    population(source, owned);
    check(!authentic_sha.empty(), "file hash missing");
    auto large_geometry = geometry();
    large_geometry.phi_nodes = 24;
    large_geometry.radial_nodes = 20;
    large_geometry.longitudinal_nodes = 20;
    auto large_source = training(Source::Policy{}, false, false, large_geometry);
    auto large = convert(large_source);
    const auto streamed = directory / "streamed.root";
    write(streamed, large);
    {
      TFile file(streamed.c_str(), "READ");
      for (const char* name : {"Nodes", "Corners", "Values", "Roots"}) persisted(file.GetDirectory("Runtime"), name);
    }
    auto large_copy = read(streamed, large);
    check(Map::ContentFingerprint(large) == Map::ContentFingerprint(large_copy), "multibasket data changed");
    const auto source_payload = 64ULL * large_source.points.size() + 232ULL * large_source.cells.size();
    const auto compact_payload = 12ULL * large.nodes.size() + 32ULL * large.corners.size() + 24ULL * large.values.size() + 4ULL * large.root_lookup.size();
    check(compact_payload * 4 < source_payload, "compact fixture did not remove witness memory");
    const auto declared = directory / "absurd.root";
    write(declared, data);
    {
      TFile file(declared.c_str(), "UPDATE");
      file.cd("Runtime");
      TParameter<Long64_t>("node_count", INT64_MAX).Write("node_count", TObject::kOverwrite);
    }
    {
      TFile file(declared.c_str(), "READ");
      Map::Data rejected;
      check(!Map::Read(file.GetDirectory("Runtime"), data.geometry, expected(data), rejected, Map::Limits{}, errors, true) && rejected.nodes.empty(), "absurd row declaration allocated/accepted");
    }
    // Reproduce missing baskets from a foreign readonly current directory;
    // positive aggregate GetEntry is deliberately not used as proof.
    const auto damaged = directory / "missing-baskets.root";
    write(damaged, large);
    const auto foreign_path = directory / "foreign.root";
    {
      TFile foreign(foreign_path.c_str(), "CREATE");
      TNamed("marker", "foreign readonly fixture").Write();
    }
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
      gErrorIgnoreLevel = saved;
      TFile file(damaged.c_str(), "UPDATE");
      file.cd("Runtime");
      check(tree.Write("Nodes", TObject::kOverwrite) > 0, "cannot write damaged fixture");
      tree.SetDirectory(nullptr);
    }
    {
      TFile file(damaged.c_str(), "READ");
      auto* t = dynamic_cast<TTree*>(file.Get("Runtime/Nodes"));
      auto* branch = t->GetBranch("first");
      check(branch->GetWriteBasket() > 1, "damaged fixture lacked lost baskets");
      Map::Data rejected;
      std::ostringstream err;
      check(!Map::Read(file.GetDirectory("Runtime"), large.geometry, expected(large), rejected, Map::Limits{}, err, true), "missing persisted baskets accepted");
      check(err.str().find("basket") != std::string::npos, "damaged fixture was not rejected by structural basket check");
    }
    const auto source_path = directory / "training.root";
    write_training(source_path, source);
    const auto source_sha = Map::FileSha256(source_path);
    const auto runtime_path = directory / "converted.root", proof_path = directory / "conversion.json";
    const auto converter = std::filesystem::absolute(std::filesystem::path(argv[0]).parent_path() / "PHGarfieldCompactMapConverter");
    check(execute(converter, {"--input", source_path, "--output", runtime_path, "--proof-output", proof_path, "--expect-source-sha256", source_sha, "--expect-source-configuration", "0123456789abcdef", "--expect-source-content", Source::ContentFingerprint(source)}) == 0, "standalone converter failed");
    check(Map::FileSha256(source_path) == source_sha, "converter modified training input");
    std::ifstream proof_file(proof_path);
    const std::string proof_text((std::istreambuf_iterator<char>(proof_file)), {});
    check(proof_text.find(Map::FileSha256(runtime_path)) != std::string::npos && proof_text.find("PASS_CONVERSION_ONLY") != std::string::npos, "external conversion manifest does not pin actual runtime SHA");
    check(execute(converter, {"--input", source_path, "--output", runtime_path, "--proof-output", proof_path, "--expect-source-sha256", source_sha, "--expect-source-configuration", "0123456789abcdef", "--expect-source-content", Source::ContentFingerprint(source)}) != 0, "converter overwrote published output");
    const auto wide_runtime_path = directory / "converted-wide.root", wide_proof_path = directory / "conversion-wide.json";
    check(execute(converter, {"--input", source_path, "--output", wide_runtime_path, "--proof-output", wide_proof_path, "--expect-source-sha256", source_sha,
                              "--expect-source-configuration", "0123456789abcdef", "--expect-source-content", Source::ContentFingerprint(source),
                              "--max-source-points", "1500000000", "--max-runtime-nodes", "250000000", "--max-runtime-values", "1500000000", "--max-runtime-corners", "150000000", "--max-runtime-bytes", "137438953472"}) == 0,
          "explicit converter resource caps failed");
    std::ifstream wide_proof_file(wide_proof_path);
    const std::string wide_proof((std::istreambuf_iterator<char>(wide_proof_file)), {});
    check(wide_proof.find("\"runtime_limits\":") != std::string::npos && wide_proof.find("\"max_source_points\":1500000000") != std::string::npos && wide_proof.find("\"max_bytes\":137438953472") != std::string::npos, "conversion proof omitted resolved explicit limits");
    {
      TFile normal(runtime_path.c_str(), "READ"), wide(wide_runtime_path.c_str(), "READ");
      Map::Data normal_data, wide_data;
      const Map::ExpectedSource expected_source{source_sha, "0123456789abcdef", Source::ContentFingerprint(source)};
      check(Map::Read(normal.GetDirectory("GarfieldFullTransport/Runtime"), source.geometry, expected_source, normal_data, Map::Limits{}, errors, true) &&
                Map::Read(wide.GetDirectory("GarfieldFullTransport/Runtime"), source.geometry, expected_source, wide_data, wide_limits(), errors, true),
            "converter outputs failed independent admission readback");
      check(Map::ContentFingerprint(normal_data) == Map::ContentFingerprint(wide_data) && wide_data.source.canonical_policy == Source::BuildCanonicalPolicy(source), "explicit resource limits changed compact arrays/source policy");
    }
    unsigned int invalid_index = 0;
    for (const auto& [flag, value] : std::array<std::pair<std::string, std::string>, 8>{{{"--max-source-points", "1500000001"}, {"--max-runtime-nodes", "250000001"}, {"--max-runtime-values", "1500000001"}, {"--max-runtime-corners", "150000001"}, {"--max-runtime-bytes", "137438953473"}, {"--max-runtime-nodes", "-1"}, {"--max-source-points", "18446744073709551616"}, {"--max-runtime-bytes", "4096junk"}}})
    {
      const auto label = "invalid-cap-" + std::to_string(invalid_index++);
      const auto out = directory / (label + ".root"), proof = directory / (label + ".json");
      check(execute(converter, {"--input", source_path, "--output", out, "--proof-output", proof, "--expect-source-sha256", source_sha,
                                "--expect-source-configuration", "0123456789abcdef", "--expect-source-content", Source::ContentFingerprint(source), flag, value}) != 0,
            "invalid converter resource cap accepted");
      check(!std::filesystem::exists(out) && !std::filesystem::exists(proof), "invalid converter resource cap published outputs");
    }
    std::cout << "PASS compact conversion/routing/allclasses/corruption/streaming/ownership/standalone fixture; source_points=" << large_source.points.size() << " retained_values=" << large.values.size() << " runtime_payload_bytes=" << compact_payload << '\n';
    std::filesystem::remove_all(directory);
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "FAIL compact fixture preserved at " << directory << ": " << e.what() << '\n';
    return 1;
  }
}
