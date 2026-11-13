#include "PHGarfieldAdaptiveMap.h"
#include "PHGarfieldCompactMap.h"
#include "PHGarfieldMapTools.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>

#include <chrono>
#include <cmath>
#include <exception>
#include <filesystem>
#include <memory>
#include <ostream>
#include <system_error>

namespace
{
  template <class T>
  void write_parameter(const char* name, T value)
  {
    TParameter<T> object(name, value);
    object.Write();
  }
}  // namespace

bool PHGarfieldAdaptiveMap::Generate(const PHGarfieldForwardTransport& garfield,
                                     const PHGarfieldForwardMap::Config& config, const Policy& policy,
                                     std::uint64_t first, std::uint64_t end, BuildSummary& build_summary,
                                     std::ostream& log, std::ostream& error)
{
  build_summary = {};
  if (config.interpolation_tolerance_um != policy.guard_um ||
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
      point.forward_status = static_cast<int>(PHGarfieldForwardTransport::ForwardDriftStatus::ReachedReadout);
      point.value = {{0., 0., 0.}};
      return true;
    }
    PHGarfieldForwardTransport::ForwardDriftResult result;
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
    if (result.status == PHGarfieldForwardTransport::ForwardDriftStatus::LeftDriftVolume)
    {
      point.reaches = 0;
      point.value = {{0., 0., 0.}};
      return true;
    }
    sampler_error << "Adaptive direct sample failed: " << PHGarfieldForwardTransport::ForwardDriftStatusName(result.status) << '\n';
    return false;
  };
  Data data;
  if (!Build(FromConfig(config), policy, first, end, sampler, data, build_summary, log, error))
  {
    return false;
  }
  summary.maximum_drift_time_ns = data.maximum_drift_time_ns;
  summary.probe_trajectories = build_summary.actual_sampler_queries;
  const auto point_count = data.points.size();
  const auto cell_count = data.cells.size();
  const bool complete = first == 0 && end == RootCount(data.geometry);
  PHGarfieldCompactMap::Data runtime;
  const auto limits = PHGarfieldMapTools::generation_limits(policy);
  const auto configuration = PHGarfieldForwardMap::Fingerprint(
      PHGarfieldForwardMap::BuildCanonicalConfiguration(config, provenance));
  if (!PHGarfieldCompactMap::Pack(data, configuration, runtime, limits, error, complete))
  {
    return false;
  }
  // Construction samples and refinement diagnostics are not part of the file.
  // Release them before the independent ROOT readback allocates runtime arrays.
  data = {};
  const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  const std::filesystem::path temporary = target.string() + ".partial." + std::to_string(stamp);
  const auto discard = [&]()
  {
    std::error_code cleanup;
    std::filesystem::remove(temporary, cleanup);
  };
  bool verified = false;
  try
  {
    PHGarfieldMapTools::write_compact_output(temporary.string(), runtime, limits, error, [&](TDirectory& metadata)
                                             {
      PHGarfieldMapTools::require(PHGarfieldForwardMap::WriteMetadata(&metadata, config, provenance, summary),
                                   "cannot write physical metadata");
      TNamed("generation_method", "adaptive_direct_compact_v1").Write();
      TNamed("generation_refinement_order", "largest_cylindrical_unresolved_leaf_volume_first_stable_cell_id_v1").Write();
      write_parameter<Long64_t>("generation_root_first", first);
      write_parameter<Long64_t>("generation_root_end", end);
      write_parameter<Long64_t>("generation_shards", 1);
      const std::array<std::uint64_t, 5> counters{{build_summary.actual_sampler_queries, build_summary.reused_points,
          build_summary.certified_leaves, build_summary.direct_leaves, build_summary.budget_leaves}};
      for (std::size_t index = 0; index < counters.size(); ++index)
      {
        write_parameter<Long64_t>(PHGarfieldMapTools::GenerationCounters[index], counters[index]);
      } });
    PHGarfieldMapTools::readback_compact_output(temporary.string(), runtime, limits, error, complete);
    verified = true;
    if (config.overwrite)
    {
      std::filesystem::rename(temporary, target);
    }
    else
    {
      PHGarfieldMapTools::publish_exclusive(temporary.string(), target.string());
    }
  }
  catch (const std::exception& exception)
  {
    error << "Adaptive LUT publication failed: " << exception.what() << '\n';
    if (verified)
    {
      error << "Verified LUT preserved at " << temporary << '\n';
    }
    else
    {
      discard();
    }
    return false;
  }
  log << "Wrote adaptive LUT root range [" << first << ',' << end << ") sampled_points=" << point_count << " cells=" << cell_count
      << " certified_leaves=" << build_summary.certified_leaves << " direct_leaves=" << build_summary.direct_leaves
      << " sampler_queries=" << build_summary.actual_sampler_queries << " stored_values=" << runtime.values.size() << '\n';
  return true;
}
