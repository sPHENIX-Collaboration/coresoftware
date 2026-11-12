#include "TpcConditionsReference.h"
#include "TpcConditionsUtils.h"
#include "TpcPHGarfieldReference.h"
#include "TpcRunAverageConditions.h"

#include <TSystem.h>
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// SDCC CDB snapshot: no event input, GL1 packet, reconstruction or collector build.

namespace TpcConditionsDetail
{
  inline std::string quote(const std::string &text)
  {
    std::string out = "\"";
    for (char c : text)
    {
      if (c == '\"')
      {
        out += '\"';
      }
      out += c;
    }
    return out + '\"';
  }
  inline std::string number(double x)
  {
    if (!std::isfinite(x))
    {
      return "";
    }
    std::ostringstream out;
    out << std::setprecision(15) << x;
    return out.str();
  }
  inline std::string joined(const std::vector<std::string> &values)
  {
    std::string out;
    for (const auto &s : values)
    {
      if (!out.empty())
      {
        out += ';';
      }
      out += s;
    }
    return out;
  }
}  // namespace TpcConditionsDetail

void TpcConditionsReference::export_conditions(int run, const std::string &groupName,
                                               const std::string &outputName,
                                               const std::string &seriesName) const

{
  if (m_cdbTag.empty() || run <= 0)
  {
    throw std::runtime_error("A positive run number and explicit CDB tag are required");
  }
  const char *group = groupName.c_str();
  const char *output = outputName.c_str();
  const char *seriesOutput = seriesName.c_str();
  const char *cdbTag = m_cdbTag.c_str();
  const bool useFallback = m_useFallback;
  using TpcConditionsDetail::par0;
  using TpcConditionsDetail::par1;
  using TpcConditionsDetail::referenceKeff;
  using TpcConditionsDetail::referenceNorth;
  using TpcConditionsDetail::referenceSouth;

  using TpcConditionsDetail::fallbackKeff;
  using TpcConditionsDetail::joined;
  using TpcConditionsDetail::median;
  using TpcConditionsDetail::missing;
  using TpcConditionsDetail::number;
  using TpcConditionsDetail::quote;
  using TpcConditionsDetail::r1Channels;
  using TpcConditionsDetail::summarize;
  auto *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", run);
  rc->set_uint64Flag("TIMESTAMP", run);
  rc->set_StringFlag("CDB_GLOBALTAG", cdbTag);
  auto *cdb = CDBInterface::instance();
  const std::string conditionFile = cdb->getUrl("TPC_CONDITIONS");
  const std::string keffFile = cdb->getUrl("Tpc_PolyClusterizer_kEff");
  std::unique_ptr<CDBTTree> conditions;
  std::unique_ptr<CDBTTree> keff;
  if (!conditionFile.empty())
  {
    conditions = std::make_unique<CDBTTree>(conditionFile);
    conditions->LoadCalibrations();
  }
  if (!keffFile.empty())
  {
    keff = std::make_unique<CDBTTree>(keffFile);
    keff->LoadCalibrations();
  }
  std::map<std::uint64_t, int> channels;
  if (conditions)
  {
    for (const auto &entry : conditions->GetUInt64EntryMap())
    {
      const auto found = entry.second.find("gbco");
      if (found == entry.second.end() || found->second == 0 ||
          found->second == std::numeric_limits<std::uint64_t>::max())
      {
        continue;
      }
      if (!channels.emplace(found->second, entry.first).second)
      {
        throw std::runtime_error("Duplicate conditions BCO");
      }
    }
  }
  std::vector<double> currents[2];
  std::vector<double> pressure;
  std::vector<double> temperature;
  std::ofstream series(seriesOutput);
  if (!series)
  {
    throw std::runtime_error("Cannot write conditions series");
  }
  series << "group,run,side,bco,gem_current_uA,pressure,temperature_K\n";
  for (const auto &entry : channels)
  {
    const double p = conditions->GetFloatValue(entry.second, "gas_pressure");
    double t = conditions->GetFloatValue(entry.second, "gas_temperature");
    if (std::isfinite(t) && t < 200.)
    {
      t += 273.15;
    }
    pressure.push_back(p);
    temperature.push_back(t);
    for (int side = 0; side < 2; ++side)
    {
      std::vector<double> values;
      for (const auto &name : r1Channels(side))
      {
        values.push_back(conditions->GetFloatValue(entry.second, name));
      }
      const double value = median(values);
      currents[side].push_back(value);
      series << quote(group) << ',' << run << ',' << side << ',' << entry.first << ','
             << number(value) << ',' << number(p) << ',' << number(t) << '\n';
    }
  }
  series.close();
  const auto p = summarize(pressure);
  const auto t = summarize(temperature);
  std::ofstream out(output);
  if (!out)
  {
    throw std::runtime_error("Cannot write conditions snapshot");
  }
  out << "group,run,side,gem_current_uA,gem_current_rms_uA,gem_current_min_uA,gem_current_max_uA,"
         "keff,keff_base,keff_applied,keff_kind,keff_source,conditions_status,current_scope,"
         "current_definition,current_channels,conditions_samples,valid_current_samples,bco_first,"
         "bco_last,pressure,temperature_K,cdb_tag,conditions_payload,keff_payload,source,recipe_"
         "reference_keff,recipe_par0,recipe_par1,recipe_reference_current_uA\n";
  for (int side = 0; side < 2; ++side)
  {
    const auto values = summarize(currents[side]);
    double scale = missing();
    std::string source = "missing";
    std::string kind = "missing";
    if (keff)
    {
      const auto &entries = keff->GetSingleFloatEntryMap();
      const auto found = entries.find("Fkeffside" + std::to_string(side));
      if (found != entries.end() && std::isfinite(found->second))
      {
        scale = found->second;
        source = "CDB base coefficient snapshot";
        kind = "cdb_base_not_verified_applied";
      }
      // A present but invalid CDB coefficient must not silently become a fallback.
    }
    else if (useFallback)
    {
      scale = fallbackKeff(values.mean, side, referenceKeff, par0, par1, referenceSouth,
                           referenceNorth);
      if (std::isfinite(scale))
      {
        source = "explicit PHGarfield GEM-current fallback recipe";
        kind = "derived_base_not_verified_applied";
      }
    }
    const bool complete = !channels.empty() && values.n == channels.size();
    std::string status = "missing_or_partial_current";
    if (complete)
    {
      status = std::isfinite(scale) ? "complete_snapshot" : "missing_keff";
    }
    out << quote(group) << ',' << run << ',' << side << ',' << number(values.mean) << ','
        << number(values.rms) << ',' << number(values.low) << ',' << number(values.high) << ','
        << number(scale) << ',' << number(scale) << ",," << quote(kind) << ',' << quote(source)
        << ',' << quote(status)
        << ",full_CDB_run_not_selected_QA_segments,unweighted_mean_of_R1_undamaged_channel_medians,"
        << quote(joined(r1Channels(side))) << ',' << channels.size() << ',' << values.n << ',';
    if (!channels.empty())
    {
      out << channels.begin()->first;
    }
    out << ',';
    if (!channels.empty())
    {
      out << channels.rbegin()->first;
    }
    out << ',' << number(p.mean) << ',' << number(t.mean) << ',' << quote(cdbTag) << ','
        << quote(conditionFile) << ',' << quote(keffFile)
        << ",\"current CDB snapshot; historical applied value not recoverable from ordinary QA "
           "histograms\","
        << number(referenceKeff) << ',' << number(par0) << ',' << number(par1) << ','
        << number(side == 0 ? referenceSouth : referenceNorth) << '\n';
  }
  out.close();
}

void TpcConditionsReference::export_phgarfield(int run, const std::string &groupName,
                                               const std::string &outputName) const

{
  if (m_cdbTag.empty())
  {
    throw std::runtime_error("An explicit CDB tag is required");
  }
  const char *group = groupName.c_str();
  const char *csv = outputName.c_str();
  const char *cdbTag = m_cdbTag.c_str();
  const double referenceField = m_referenceField;
  const double referencePressure = m_referencePressure;
  const double referenceTemperatureK = m_referenceTemperatureK;
  const double t0Bins = m_t0Bins;
  const double timeBinNs = m_timeBinNs;
  const double stepNs = m_stepNs;

  if (run <= 0 || !std::isfinite(referenceField) || referenceField <= 0 ||
      !std::isfinite(referencePressure) || referencePressure <= 0 ||
      !std::isfinite(referenceTemperatureK) || referenceTemperatureK <= 0 ||
      !std::isfinite(t0Bins) || !std::isfinite(timeBinNs) || timeBinNs <= 0 ||
      !std::isfinite(stepNs) || stepNs <= 0)
  {
    throw std::runtime_error("Invalid PHGarfield reference settings");
  }
  std::unique_ptr<Fun4AllServer> server(Fun4AllServer::instance());
  auto *se = server.get();
  auto *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", run);
  rc->set_uint64Flag("TIMESTAMP", run);
  rc->set_StringFlag("CDB_GLOBALTAG", cdbTag);
  se->registerSubsystem(new FlagHandler());
  const std::string geometry = CDBInterface::instance()->getUrl("Tracking_Geometry");
  if (geometry.empty())
  {
    throw std::runtime_error("Tracking_Geometry CDB lookup failed");
  }
  auto *geo = new Fun4AllRunNodeInputManager("TPCPHGGeometry");
  geo->AddFile(geometry);
  se->registerInputManager(geo);
  auto *dummy = new Fun4AllDummyInputManager("TPCPHGConditionsEvent");
  dummy->SetRunNumber(run);
  se->registerInputManager(dummy);
  se->registerSubsystem(new TpcRunAverageConditions());
  auto *qa = new TpcPHGarfieldReference(referenceField, referencePressure, referenceTemperatureK,
                                        t0Bins, timeBinNs, stepNs);
  se->registerSubsystem(qa);
  if (se->run(1) != 0)
  {
    throw std::runtime_error("PHGarfield reference InitRun failed");
  }
  const auto velocity = qa->velocity();
  for (double value : velocity)
  {
    if (!std::isfinite(value) || value <= 0)
    {
      throw std::runtime_error("No valid PHGarfield mean velocity");
    }
  }
  const auto quoted = [](const std::string &text)
  {
    std::string out = "\"";
    for (char c : text)
    {
      if (c == '\"')
      {
        out += '\"';
      }
      out += c;
    }
    return out + '\"';
  };
  std::ofstream out(csv);
  if (!out)
  {
    throw std::runtime_error("Cannot write reference CSV");
  }
  out << std::setprecision(14);
  out << "group,run,side,velocity_cm_us,t0_bins,source,cdb_tag,reference_field_"
         "Vcm,effective_field_Vcm,pressure,temperature_K,time_bin_ns,step_ns,z_"
         "pad_cm,cm_halfwidth_cm,geometry,reference_pressure,reference_"
         "temperature_K\n";
  for (int side = 0; side < 2; ++side)
  {
    out << quoted(group) << ',' << run << ',' << side << ',' << velocity[side] << ',' << t0Bins
        << ",\"PHGarfield full-run-average P/T snapshot equal-FEE baseline; CDB drift "
           "corrections disabled\","
        << quoted(cdbTag) << ',' << referenceField << ',' << qa->effectiveField() << ','
        << qa->pressure() << ',' << qa->temperatureK() << ',' << timeBinNs << ',' << stepNs << ','
        << qa->padPlane() << ',' << qa->cmHalf() << ',' << quoted(geometry) << ','
        << referencePressure << ',' << referenceTemperatureK << '\n';
  }
  out.close();
  if (!out)
  {
    throw std::runtime_error("Failed to finish PHGarfield reference CSV");
  }
  se->End();
}
