
#include "PHActsKDTreeSeeding.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>


#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/SeedFilterConfig.hpp>
#include <Acts/Seeding/SeedFinderOrthogonalConfig.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>
#include <Acts/Utilities/KDTree.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>
#include <Acts/Seeding/SeedFinderOrthogonal.hpp>

#include <optional>


//____________________________________________________________________________..
PHActsKDTreeSeeding::PHActsKDTreeSeeding(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
PHActsKDTreeSeeding::~PHActsKDTreeSeeding()
{
}

//____________________________________________________________________________..
int PHActsKDTreeSeeding::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsKDTreeSeeding::InitRun(PHCompositeNode*)
{

  Acts::SeedFilterConfig filterCfg;
  filterCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
          Acts::SeedFilter<SimSpacePoint>(filterCfg));

  m_cfg.seedFinderConfig.rMax = m_cfg.rMax;
  m_cfg.seedFinderConfig.deltaRMinTopSP = m_cfg.deltaRMinTopSP;
  m_cfg.seedFinderConfig.deltaRMaxTopSP = m_cfg.deltaRMaxTopSP;
  m_cfg.seedFinderConfig.deltaRMinBottomSP = m_cfg.deltaRMinBottomSP;
  m_cfg.seedFinderConfig.deltaRMaxBottomSP = m_cfg.deltaRMaxBottomSP;
  m_cfg.seedFinderConfig.collisionRegionMin = m_cfg.collisionRegionMin;
  m_cfg.seedFinderConfig.collisionRegionMax = m_cfg.collisionRegionMax;
  m_cfg.seedFinderConfig.zMin = m_cfg.zMin;
  m_cfg.seedFinderConfig.zMax = m_cfg.zMax;
  m_cfg.seedFinderConfig.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_cfg.seedFinderConfig.cotThetaMax = m_cfg.cotThetaMax;
  m_cfg.seedFinderConfig.sigmaScattering = m_cfg.sigmaScattering;
  m_cfg.seedFinderConfig.radLengthPerSeed = m_cfg.radLengthPerSeed;
  m_cfg.seedFinderConfig.minPt = m_cfg.minPt;
  m_cfg.seedFinderConfig.bFieldInZ = m_cfg.bFieldInZ;
  m_cfg.seedFinderConfig.beamPos =
      Acts::Vector2(m_cfg.beamPosX, m_cfg.beamPosY);
  m_cfg.seedFinderConfig.impactMax = m_cfg.impactMax;

  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_cfg.seedFinderConfig.highland =
      13.6 * std::sqrt(m_cfg.seedFinderConfig.radLengthPerSeed) *
      (1 + 0.038 * std::log(m_cfg.seedFinderConfig.radLengthPerSeed));
  float maxScatteringAngle =
      m_cfg.seedFinderConfig.highland / m_cfg.seedFinderConfig.minPt;
  m_cfg.seedFinderConfig.maxScatteringAngle2 =
      maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_cfg.seedFinderConfig.pTPerHelixRadius =
      300. * m_cfg.seedFinderConfig.bFieldInZ;
  m_cfg.seedFinderConfig.minHelixDiameter2 =
      std::pow(m_cfg.seedFinderConfig.minPt * 2 /
                   m_cfg.seedFinderConfig.pTPerHelixRadius,
               2);

  m_cfg.seedFinderConfig.pT2perRadius = std::pow(
      m_cfg.seedFinderConfig.highland / m_cfg.seedFinderConfig.pTPerHelixRadius,
      2);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsKDTreeSeeding::process_event(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsKDTreeSeeding::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
