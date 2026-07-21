#include "PHG4GeantinoIonization.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>

namespace
{
  // Silicon MIP stopping power used by the INTT digitizer.
  constexpr double siliconMipDedx = 0.003876;  // GeV/cm

  // Mean MIP stopping power for the default TPC gas mixture
  // Ar/CF4/isobutane = 75/20/5.
  constexpr double tpcMipDedx =
      (0.75 * 2.44 + 0.20 * 7.00 + 0.05 * 5.93) * 1e-6;  // GeV/cm

  // Mean MIP stopping power for the default Micromegas gas mixture
  // Ar/isobutane = 90/10.
  constexpr double micromegasMipDedx =
      (0.90 * 2.44 + 0.10 * 5.93) * 1e-6;  // GeV/cm
}  // namespace

PHG4GeantinoIonization::PHG4GeantinoIonization(const std::string& name)
  : SubsysReco(name)
  , m_detectorConfigs{{
        {"MVTX", "G4HIT_MVTX", siliconMipDedx, true},
        {"INTT", "G4HIT_INTT", siliconMipDedx, true},
        {"TPC", "G4HIT_TPC", tpcMipDedx, true},
        {"MICROMEGAS", "G4HIT_MICROMEGAS", micromegasMipDedx, true}}}
{
}

void PHG4GeantinoIonization::set_tpc_gas_fractions(
    const double neon,
    const double argon,
    const double cf4,
    const double nitrogen,
    const double isobutane)
{
  // Keep these values synchronized with PHG4TpcElectronDrift. With
  // eion = <dE/dx> * path length, its electrons-per-GeV conversion gives
  // the corresponding mean number of primary electrons for this mixture.
  constexpr double neonMipDedx = 1.56;       // keV/cm
  constexpr double argonMipDedx = 2.44;      // keV/cm
  constexpr double cf4MipDedx = 7.00;        // keV/cm
  constexpr double nitrogenMipDedx = 2.127;  // keV/cm
  constexpr double isobutaneMipDedx = 5.93;  // keV/cm

  m_detectorConfigs[2].mipDedx =
      (neon * neonMipDedx + argon * argonMipDedx + cf4 * cf4MipDedx +
       nitrogen * nitrogenMipDedx + isobutane * isobutaneMipDedx) *
      1e-6;
}

int PHG4GeantinoIonization::process_event(PHCompositeNode* topNode)
{
  const auto* truthInfo =
      findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
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

    auto* hits =
        findNode::getClass<PHG4HitContainer>(topNode, config.hitNodeName);
    if (!hits)
    {
      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << " Missing optional node "
                  << config.hitNodeName << std::endl;
      }
      continue;
    }

    const auto counters = process_detector(hits, truthInfo, config);
    if (Verbosity() > 0)
    {
      std::cout << Name() << " " << config.name
                << ": inspected " << counters.inspected
                << ", modified " << counters.modified
                << ", missing particle " << counters.missingParticle
                << ", invalid path " << counters.invalidPath
                << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

PHG4GeantinoIonization::DetectorCounters
PHG4GeantinoIonization::process_detector(
    PHG4HitContainer* hits,
    const PHG4TruthInfoContainer* truthInfo,
    const DetectorConfig& config) const
{
  DetectorCounters counters;
  const auto hitRange = hits->getHits();
  for (auto hitIter = hitRange.first; hitIter != hitRange.second; ++hitIter)
  {
    auto* hit = hitIter->second;
    if (!hit)
    {
      continue;
    }

    ++counters.inspected;

    const auto* particle = truthInfo->GetParticle(hit->get_trkid());
    if (!particle)
    {
      ++counters.missingParticle;
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

    if (!std::isfinite(pathLength) || pathLength <= 0 ||
        !std::isfinite(config.mipDedx) || config.mipDedx <= 0)
    {
      ++counters.invalidPath;
      continue;
    }

    const double syntheticIonization = config.mipDedx * pathLength;
    hit->set_edep(syntheticIonization);
    hit->set_eion(syntheticIonization);
    ++counters.modified;

    if (Verbosity() > 2)
    {
      std::cout << Name() << " " << config.name
                << " hit " << hitIter->first
                << " track " << hit->get_trkid()
                << " path length " << pathLength << " cm"
                << " synthetic edep/eion " << syntheticIonization << " GeV"
                << std::endl;
    }
  }

  return counters;
}
