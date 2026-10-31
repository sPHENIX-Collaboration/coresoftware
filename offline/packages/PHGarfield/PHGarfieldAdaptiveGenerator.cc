#include "PHGarfieldAdaptiveMap.h"
#include "PHGarfieldRootOutput.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <memory>
#include <ostream>
#include <system_error>

namespace
{
  using Map = PHGarfieldAdaptiveMap;
  template <class T>
  void write_parameter(const char* name, T value)
  {
    TParameter<T> object(name, value);
    object.Write();
  }
}  // namespace

bool PHGarfieldAdaptiveMap::Generate(const PHGarfield& garfield,
                                     const PHGarfieldForwardMap::Config& config, const Policy& policy,
                                     std::uint64_t first, std::uint64_t end, BuildSummary& build_summary,
                                     std::ostream& log, std::ostream& error)
{
  build_summary = {};
  if (config.phi_chunk_count != 1 || config.phi_chunk_index != 0 ||
      !config.probe_cell_centers || config.interpolation_tolerance_um != policy.guard_um ||
      garfield.GetElectricFieldMap3D(0).empty() || garfield.GetElectricFieldMap3D(1).empty())
  {
    error << "Adaptive generation requires pinned per-side3D, full stencil probes, matching positive guard and explicit root-range partitioning\n";
    return false;
  }
  const std::filesystem::path target(config.output_filename);
  std::error_code io_error;
  if (!config.overwrite && std::filesystem::exists(target, io_error))
  {
    error << "Adaptive output already exists: " << target << '\n';
    return false;
  }
  if (io_error)
  {
    error << "Cannot inspect adaptive output: " << io_error.message() << '\n';
    return false;
  }
  PHGarfieldForwardMap::Provenance provenance;
  PHGarfieldForwardMap::Summary summary;
  if (!PHGarfieldForwardMap::Prepare(garfield, config, provenance, summary, log, error))
  {
    return false;
  }
  provenance.reach_probe_scheme = ProbeScheme;
  const auto sampler = [&](const std::array<double, 3>& coordinate, Point& point, std::ostream& sampler_error)
  {
    // The exact boundary query is mathematical identity. Encode it explicitly
    // rather than injecting r*cos/sin roundoff into a zero-length trajectory.
    if (std::abs(coordinate[2]) == config.readout_abs_z_cm)
    {
      point.reaches = 1;
      point.forward_status = static_cast<int>(PHGarfield::ForwardDriftStatus::ReachedReadout);
      point.value = {{0., 0., 0.}};
      return true;
    }
    PHGarfield::ForwardDriftResult result;
    const bool ok = garfield.ForwardDrift(coordinate[1] * std::cos(coordinate[0]), coordinate[1] * std::sin(coordinate[0]), coordinate[2],
                                          config.readout_abs_z_cm, config.integration_step_ns, result, true);
    point.forward_status = static_cast<int>(result.status);
    if (ok && result.reached_readout)
    {
      PHGarfieldForwardMap::NodeValue encoded;
      if (!PHGarfieldForwardMap::EncodeNode(coordinate[1], coordinate[0], result, encoded))
      {
        sampler_error << "Adaptive point encoding failed\n";
        return false;
      }
      point.reaches = 1;
      point.value = {{encoded.delta_r_cm, encoded.delta_rphi_cm, encoded.drift_time_ns}};
      return true;
    }
    if (result.status == PHGarfield::ForwardDriftStatus::LeftDriftVolume)
    {
      point.reaches = 0;
      point.value = {{0., 0., 0.}};
      return true;
    }
    sampler_error << "Adaptive direct sample failed: " << PHGarfield::ForwardDriftStatusName(result.status) << '\n';
    return false;
  };
  Data data;
  if (!Build(FromConfig(config), policy, first, end, sampler, data, build_summary, log, error))
  {
    return false;
  }
  summary.maximum_drift_time_ns = data.maximum_drift_time_ns;
  summary.probe_trajectories = build_summary.actual_sampler_queries;
  const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  const std::filesystem::path temporary = target.string() + ".partial." + std::to_string(stamp);
  auto discard = [&]()
  {std::error_code cleanup;std::filesystem::remove(temporary,cleanup); };
  {
    std::unique_ptr<TFile> file(TFile::Open(temporary.c_str(), "CREATE"));
    if (!file || file->IsZombie())
    {
      error << "Cannot create adaptive temporary file " << temporary << '\n';
      discard();
      return false;
    }
    auto* directory = file->mkdir("GarfieldFullTransport");
    auto* metadata = directory ? directory->mkdir("Metadata") : nullptr;
    auto* adaptive = directory ? directory->mkdir("Adaptive") : nullptr;
    const bool complete = first == 0 && end == RootCount(data.geometry);
    if (!PHGarfieldForwardMap::WriteMetadata(metadata, config, provenance, summary, SchemaVersion,
                                             complete ? PHGarfieldForwardMap::SchemaName : "PHG4TpcGarfieldFullTransportLUTShard") ||
        !Write(adaptive, data, error))
    {
      file->Close();
      discard();
      return false;
    }
    metadata->cd();
    TNamed("generation_refinement_order", "largest_cylindrical_unresolved_leaf_volume_first_stable_cell_id_v1").Write();
    write_parameter<Long64_t>("generation_root_first", first);
    write_parameter<Long64_t>("generation_root_end", end);
    write_parameter<Long64_t>("adaptive_sampler_queries", build_summary.actual_sampler_queries);
    write_parameter<Long64_t>("adaptive_reused_points", build_summary.reused_points);
    write_parameter<Long64_t>("adaptive_certified_leaves", build_summary.certified_leaves);
    write_parameter<Long64_t>("adaptive_direct_leaves", build_summary.direct_leaves);
    write_parameter<Long64_t>("adaptive_budget_leaves", build_summary.budget_leaves);
    if (!PHGarfieldFlushAndCloseRootOutput(*file, error))
    {
      discard();
      return false;
    }
  }
  // Read into independent owning data after closing every writer. An incomplete
  // shard is allowed only as an explicitly declared execution range.
  {
    std::unique_ptr<TFile> file(TFile::Open(temporary.c_str(), "READ"));
    Data independent;
    if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered) ||
        !Read(file->GetDirectory("GarfieldFullTransport/Adaptive"), data.geometry, independent, error, false) ||
        ContentFingerprint(independent) != ContentFingerprint(data))
    {
      error << "Independent adaptive ROOT readback failed\n";
      if (file)
      {
        file->Close();
      }
      discard();
      return false;
    }
  }
  // hard-link provides atomic no-clobber publication on the same filesystem.
  if (config.overwrite)
  {
    std::filesystem::rename(temporary, target, io_error);
  }
  else
  {
    std::filesystem::create_hard_link(temporary, target, io_error);
    if (!io_error)
    {
      discard();
    }
  }
  if (io_error)
  {
    error << "Adaptive publication failed: " << io_error.message() << '\n';
    discard();
    return false;
  }
  std::uint64_t estimated_bytes = 0;
  EstimatedBytes(data.points.size(), data.cells.size(), RootCount(data.geometry), estimated_bytes);
  log << "Wrote adaptive root range [" << first << ',' << end << ") points=" << data.points.size() << " cells=" << data.cells.size()
      << " certified_leaves=" << build_summary.certified_leaves << " direct_leaves=" << build_summary.direct_leaves
      << " sampler_queries=" << build_summary.actual_sampler_queries << " estimated_peak_bytes=" << estimated_bytes << '\n';
  return true;
}
