#include "PHG4GeantinoIonization.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <cstddef>
#include <iostream>

namespace
{
  struct GasFractions
  {
    double neon = 0;
    double argon = 0;
    double cf4 = 0;
    double nitrogen = 0;
    double isobutane = 0;
  };

  // Pure-gas MIP stopping powers in keV/cm. They match the values used by
  // the current TPC and Micromegas hit-reconstruction modules.
  constexpr double neonMipDedx = 1.56;
  constexpr double argonMipDedx = 2.44;
  constexpr double cf4MipDedx = 7.00;
  constexpr double nitrogenMipDedx = 2.127;
  constexpr double isobutaneMipDedx = 5.93;

  // Mean silicon MIP stopping powers in GeV/cm. The MVTX value corresponds
  // to 9.6 keV in 25 microns; the INTT value is the value documented by its
  // hit reconstruction.
  double mvtxMipDedx = 0.00384;
  double inttMipDedx = 0.00387;

  GasFractions tpcGasFractions{0.00, 0.75, 0.20, 0.00, 0.05};
  GasFractions tpotGasFractions{0.00, 0.90, 0.00, 0.00, 0.10};

  double calculateMipDedx(const GasFractions& fractions)
  {
    return fractions.neon * neonMipDedx +
           fractions.argon * argonMipDedx +
           fractions.cf4 * cf4MipDedx +
           fractions.nitrogen * nitrogenMipDedx +
           fractions.isobutane * isobutaneMipDedx;
  }
}  // namespace

PHG4GeantinoIonization::PHG4GeantinoIonization(const std::string& name)
  : SubsysReco(name)
  , m_detectorConfigs{{
        {DetectorId::mvtx, "MVTX", "G4HIT_MVTX", true},
        {DetectorId::intt, "INTT", "G4HIT_INTT", true},
        {DetectorId::tpc, "TPC", "G4HIT_TPC", true},
        {DetectorId::tpot, "MICROMEGAS", "G4HIT_MICROMEGAS", true}}}
{
}

void PHG4GeantinoIonization::set_mvtx_mip_dedx(const double value)
{
  mvtxMipDedx = value;
}

void PHG4GeantinoIonization::set_intt_mip_dedx(const double value)
{
  inttMipDedx = value;
}

void PHG4GeantinoIonization::set_tpc_gas_fractions(
    const double neon,
    const double argon,
    const double cf4,
    const double nitrogen,
    const double isobutane)
{
  tpcGasFractions = {neon, argon, cf4, nitrogen, isobutane};
}

void PHG4GeantinoIonization::set_tpot_gas_fractions(
    const double neon,
    const double argon,
    const double cf4,
    const double nitrogen,
    const double isobutane)
{
  tpotGasFractions = {neon, argon, cf4, nitrogen, isobutane};
}

double PHG4GeantinoIonization::mip_dedx(const DetectorId detector)
{
  switch (detector)
  {
  case DetectorId::mvtx:
    return mvtxMipDedx;
  case DetectorId::intt:
    return inttMipDedx;
  case DetectorId::tpc:
    return 1e-6 * calculateMipDedx(tpcGasFractions);
  case DetectorId::tpot:
    return 1e-6 * calculateMipDedx(tpotGasFractions);
  }

  return 0;
}

int PHG4GeantinoIonization::process_event(PHCompositeNode* topNode)
{
  const auto* truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthInfo)
  {
    std::cout << PHWHERE << " Missing G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (const auto& config : m_detectorConfigs)
  {
    if (!config.enabled)
    {
      continue;
    }

    auto* hits = findNode::getClass<PHG4HitContainer>(topNode, config.hitNodeName);
    if (!hits)
    {
      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << " Missing optional node "
                  << config.hitNodeName << std::endl;
      }
      continue;
    }

    process_detector(hits, truthInfo, config);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4GeantinoIonization::process_detector(
    PHG4HitContainer* hits,
    const PHG4TruthInfoContainer* truthInfo,
    const DetectorConfig& config) const
{
  std::size_t inspected = 0;
  std::size_t modified = 0;
  std::size_t missingParticle = 0;
  std::size_t invalidPath = 0;

  const auto hitRange = hits->getHits();
  for (auto hitIter = hitRange.first; hitIter != hitRange.second; ++hitIter)
  {
    auto* hit = hitIter->second;
    if (!hit)
    {
      continue;
    }

    ++inspected;

    const auto* particle = truthInfo->GetParticle(hit->get_trkid());
    if (!particle)
    {
      ++missingParticle;
      continue;
    }

    if (particle->get_name() != m_particleName)
    {
      continue;
    }

    // Keep the operation idempotent. Current stepping actions store negative
    // edep/eion sentinels for geantinos. A finite nonnegative value means this
    // hit has already been processed.
    const double edep = hit->get_edep();
    const double eion = hit->get_eion();
    if (std::isfinite(edep) && edep >= 0 &&
        std::isfinite(eion) && eion >= 0)
    {
      continue;
    }

    const double dx = hit->get_x(1) - hit->get_x(0);
    const double dy = hit->get_y(1) - hit->get_y(0);
    const double dz = hit->get_z(1) - hit->get_z(0);
    const double pathLength = std::sqrt(dx * dx + dy * dy + dz * dz);

    const double mipDedx = mip_dedx(config.detector);
    if (!std::isfinite(pathLength) || pathLength <= 0 ||
        !std::isfinite(mipDedx) || mipDedx <= 0)
    {
      ++invalidPath;
      continue;
    }

    const double syntheticEnergyDeposit = mipDedx * pathLength;
    const double syntheticIonization = syntheticEnergyDeposit;
    hit->set_edep(syntheticEnergyDeposit);
    hit->set_eion(syntheticIonization);
    ++modified;

    if (Verbosity() > 2)
    {
      std::cout << Name() << " " << config.name
                << " hit " << hitIter->first
                << " track " << hit->get_trkid()
                << " path length " << pathLength << " cm"
                << " synthetic edep " << syntheticEnergyDeposit << " GeV"
                << ", eion " << syntheticIonization << " GeV"
                << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << " " << config.name
              << ": inspected " << inspected
              << ", modified " << modified
              << ", missing particle " << missingParticle
              << ", invalid path " << invalidPath
              << std::endl;
  }
}
