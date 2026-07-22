#include "PHG4GeantinoIonization.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <cstddef>
#include <iostream>

namespace
{
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

  // PHG4MicromegasDetector and PHG4MicromegasHitReco both use a fixed
  // Ar/isobutane 90/10 gas mixture.
  constexpr double tpotMipDedx =
      1e-6 * (0.9 * argonMipDedx + 0.1 * isobutaneMipDedx);
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

int PHG4GeantinoIonization::InitRun(PHCompositeNode* topNode)
{
  if (!m_detectorConfigs[2].enabled)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Read the gas fractions from the TPC geometry parameters, as is done in
  // PHG4TpcElectronDrift. This keeps the synthetic ionization consistent with
  // the geometry built by the macro or loaded from the CDB.
  const auto* tpcParamsContainer =
      findNode::getClass<PHParametersContainer>(topNode, "G4GEO_TPC");
  if (!tpcParamsContainer)
  {
    std::cout << PHWHERE << " Missing G4GEO_TPC" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const PHParameters* tpcParams = tpcParamsContainer->GetParameters(0);
  if (!tpcParams)
  {
    std::cout << PHWHERE << " Missing TPC geometry parameters" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const double neonFraction = tpcParams->get_double_param("Ne_frac");
  const double argonFraction = tpcParams->get_double_param("Ar_frac");
  const double cf4Fraction = tpcParams->get_double_param("CF4_frac");
  const double nitrogenFraction = tpcParams->get_double_param("N2_frac");
  const double isobutaneFraction = tpcParams->get_double_param("isobutane_frac");

  m_tpcMipDedx =
      1e-6 * (neonFraction * neonMipDedx +
              argonFraction * argonMipDedx +
              cf4Fraction * cf4MipDedx +
              nitrogenFraction * nitrogenMipDedx +
              isobutaneFraction * isobutaneMipDedx);

  if (Verbosity() > 0)
  {
    std::cout << Name()
              << " TPC gas fractions (Ne/Ar/CF4/N2/isobutane): "
              << neonFraction << "/" << argonFraction << "/"
              << cf4Fraction << "/" << nitrogenFraction << "/"
              << isobutaneFraction
              << ", MIP dE/dx: " << m_tpcMipDedx << " GeV/cm"
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4GeantinoIonization::set_mvtx_mip_dedx(const double value)
{
  mvtxMipDedx = value;
}

void PHG4GeantinoIonization::set_intt_mip_dedx(const double value)
{
  inttMipDedx = value;
}

double PHG4GeantinoIonization::mip_dedx(const DetectorId detector) const
{
  switch (detector)
  {
  case DetectorId::mvtx:
    return mvtxMipDedx;
  case DetectorId::intt:
    return inttMipDedx;
  case DetectorId::tpc:
    return m_tpcMipDedx;
  case DetectorId::tpot:
    return tpotMipDedx;
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
