#include "TpcPHGarfieldReference.h"
#include <TPolyLine3D.h>
#include <TSystem.h>
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <phgarfield/PHGarfield.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <tpcconditions/TpcConditions.h>
#include "TpcConditionsReference.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>

#include "TpcRunAverageConditions.h"

TpcPHGarfieldReference::TpcPHGarfieldReference(double field, double pressure,
                                               double temperature, double t0,
                                               double clock, double step)
  : SubsysReco("TpcPHGarfieldReference")
  , m_referenceCMFieldVcm(field)
  , m_referencePressure(pressure)
  , m_referenceTemperatureK(temperature)
  , m_tpcTimeBinNs(clock)
  , m_phgT0Bins(t0)
  , m_phgStepNs(step)
{
}
TpcPHGarfieldReference::~TpcPHGarfieldReference() = default;
int TpcPHGarfieldReference::InitRun(PHCompositeNode *topNode)
{
  if (!buildFeeMap(topNode))
  {
    throw std::runtime_error("TPC FEE map/geometry failed for PHGarfield reference");
  }
  if (!initializeDriftCoordinateCalibration(topNode))
  {
    throw std::runtime_error("PHGarfield reference integration failed; see preceding log");
  }
  for (int side = 0; side < 2; ++side)
  {
    m_velocity[side] = 1000. *
                       std::accumulate(m_phgDriftVelocityCmNs[side].begin(),
                                       m_phgDriftVelocityCmNs[side].end(), 0.) /
                       kNFees;
  }
  auto *geometry = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  m_cmHalf = geometry->GetLayerCellGeom(20)->get_CM_halfwidth();
  return Fun4AllReturnCodes::EVENT_OK;
}

bool TpcPHGarfieldReference::buildFeeMap(PHCompositeNode *topNode)
{
  auto *geom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  const std::string url = CDBInterface::instance()->getUrl("TPC_FEE_CHANNEL_MAP");
  if (!geom || url.empty())
  {
    return false;
  }
  CDBTTree map(url);
  map.LoadCalibrations();
  // Same physical FEE -> CDB mapping and sector convention as the TPC unpacker.
  const double pi = std::acos(-1.);
  const int feeMap[26]{4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1,
                       3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};
  const int feeR[26]{2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2,
                     1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};

  m_feeX = {};
  m_feeY = {};
  m_feeCount = {};
  std::map<std::tuple<int, int, unsigned int>, int> owners;
  for (int side = 0; side < 2; ++side)
  {
    owners.clear();
    for (int sec = 0; sec < 12; ++sec)
    {
      for (int fee = 0; fee < 26; ++fee)
      {
        for (int ch = 0; ch < 256; ++ch)
        {
          int offset = 0;
          if (feeR[fee] == 2)
          {
            offset = 6;
          }
          else if (feeR[fee] == 3)
          {
            offset = 14;
          }
          const int mapped = feeMap[fee] + offset;
          const int key = 256 * mapped + ch;
          const int layer = map.GetIntValue(key, "layer");
          if (layer < 7 || layer > 54)
          {
            continue;  // antenna pads are not in TRKR_HITSET
          }
          auto *lg = geom->GetLayerCellGeom(layer);
          if (!lg)
          {
            return false;
          }
          const double phi =
              (side == 1 ? 1. : -1.) * (map.GetDoubleValue(key, "phi") - pi / 2.) + sec * pi / 6.;
          if (!is_good_number(phi))
          {
            return false;
          }
          const unsigned int pad = lg->get_phibin(phi, side);
          const auto lookup = std::make_tuple(sec, layer, pad);
          const int id = 26 * sec + fee;
          auto result = owners.emplace(lookup, id);
          if (!result.second && result.first->second != id)
          {
            return false;
          }
          if (!result.second)
          {
            continue;
          }
          m_feeX[side][id] += lg->get_radius() * std::cos(phi);
          m_feeY[side][id] += lg->get_radius() * std::sin(phi);
          m_feeCount[side][id] += 1.;
        }
      }
    }
  }
  for (int s = 0; s < 2; ++s)
  {
    for (int f = 0; f < kNFees; ++f)
    {
      if (m_feeCount[s][f] <= 0.)
      {
        return false;
      }
      m_feeX[s][f] /= m_feeCount[s][f];
      m_feeY[s][f] /= m_feeCount[s][f];
    }
  }
  return true;
}

bool TpcPHGarfieldReference::initializeDriftCoordinateCalibration(PHCompositeNode *topNode)
{
  for (auto &side_values : m_phgDriftTimeNs)
  {
    side_values.fill(0.0);
  }
  for (auto &side_values : m_phgDriftVelocityCmNs)
  {
    side_values.fill(0.0);
  }

  auto *conditions = findNode::getClass<TpcConditions>(topNode, "TpcConditions");
  if (!conditions || !conditions->get_ConditionsAvailable())
  {
    std::cout << Name() << " TpcConditions unavailable for PHGarfield P/T calibration"
              << std::endl;
    return false;
  }

  auto *geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!geom_container)
  {
    std::cout << Name() << " missing TPCGEOMCONTAINER" << std::endl;
    return false;
  }

  auto *layer_geom = geom_container->GetLayerCellGeom(20);
  if (!layer_geom)
  {
    std::cout << Name() << " missing TPC geometry for layer 20" << std::endl;
    return false;
  }

  m_zPadPlaneCm = layer_geom->get_max_driftlength() + layer_geom->get_CM_halfwidth();
  if (!is_good_number(m_zPadPlaneCm) || m_zPadPlaneCm <= 0.0)
  {
    std::cout << Name() << " invalid pad-plane z = " << m_zPadPlaneCm << std::endl;
    return false;
  }

  const double pressure = conditions->get_Pressure();
  const double temperatureK = temperature_kelvin(conditions->get_Temperature());

  m_pressure = pressure;
  m_temperatureK = temperatureK;
  double effective_field_vcm = m_referenceCMFieldVcm;
  if (m_usePTCorrectedPHGarfield)
  {
    if (!(pressure > 0.0 && temperatureK > 0.0 && m_referencePressure > 0.0 &&
          m_referenceTemperatureK > 0.0))
    {
      std::cout << Name() << " invalid P/T for PHGarfield scaling:"
                << " P=" << pressure << " T=" << temperatureK << std::endl;
      return false;
    }

    // Match the run-dependent reduced field E/N with the reference gas table:
    // E_eff = E_ref * (P_ref/P_run) * (T_run/T_ref).
    effective_field_vcm = m_referenceCMFieldVcm * (m_referencePressure / pressure) *
                          (temperatureK / m_referenceTemperatureK);
  }

  std::cout << Name() << "::initializeDriftCoordinateCalibration"
            << " P=" << pressure << " T=" << temperatureK << " K"
            << " Eref=" << m_referenceCMFieldVcm << " V/cm"
            << " Eeff=" << effective_field_vcm << " V/cm"
            << " zPad=" << m_zPadPlaneCm << " cm"
            << " t0=" << m_phgT0Bins << " bins"
            << " dt=" << m_tpcTimeBinNs << " ns/bin" << std::endl;

  m_driftPHGarfield = std::make_unique<PHGarfield>(Name() + "_DriftPHGarfield");
  m_effectiveField = effective_field_vcm;
  // This private diagnostic must never publish a node that becomes dangling
  // when the temporary integrator is deleted at the end of InitRun.
  m_driftPHGarfield->SetPublishToNodeTree(false);
  m_driftPHGarfield->SetUseCDBDriftCorrections(false);
  m_driftPHGarfield->SetCMVoltageDefault(effective_field_vcm);

  const int garfield_ret = m_driftPHGarfield->InitRun(topNode);
  if (garfield_ret != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << Name() << " PHGarfield::InitRun failed with code " << garfield_ret << std::endl;
    m_driftPHGarfield.reset();
    return false;
  }

  for (int side = 0; side < 2; ++side)
  {
    for (int fee = 0; fee < kNFees; ++fee)
    {
      if (fee % 52 == 0)
      {
        std::cout << Name() << ": integrating side " << side << ", FEE " << fee << "/" << kNFees
                  << std::endl;
      }
      PHGarfield::ReverseDriftStatus status = PHGarfield::ReverseDriftStatus::Running;
      std::unique_ptr<TPolyLine3D> path(m_driftPHGarfield->ReverseDrift(
          m_feeX[side][fee], m_feeY[side][fee], side == 0 ? -m_zPadPlaneCm : m_zPadPlaneCm,
          m_phgStepNs, &status));
      if (!path || path->GetN() < 1 || status != PHGarfield::ReverseDriftStatus::CentralMembrane)
      {
        std::cerr << Name() << ": trajectory did not reach CM: side=" << side << " FEE=" << fee
                  << " status=" << static_cast<int>(status) << std::endl;
        return false;
      }
      // ReverseDrift stores the initial point but not the step that crosses the CM.
      // GetN()*step includes that last step, matching the original scan convention.
      const double time = path->GetN() * m_phgStepNs;
      m_phgDriftTimeNs[side][fee] = time;
      m_phgDriftVelocityCmNs[side][fee] = m_zPadPlaneCm / time;
    }
  }
  m_driftPHGarfield.reset();  // integration only at InitRun, never per hit
  return true;
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
