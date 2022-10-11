// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSKDTREESEEDING_H
#define PHACTSKDTREESEEDING_H

#include "SpacePoint.h"

#include <fun4all/SubsysReco.h>

#include <Acts/Seeding/SeedFilterConfig.hpp>
#include <Acts/Seeding/SeedFinderOrthogonalConfig.hpp>

#include <string>

class PHCompositeNode;

class PHActsKDTreeSeeding : public SubsysReco
{
 public:

  PHActsKDTreeSeeding(const std::string &name = "PHActsKDTreeSeeding");

  ~PHActsKDTreeSeeding() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:

  Acts::SeedFilterConfig m_seedFilterConfig;
  Acts::SeedFinderOrthogonalConfig<SpacePoint> m_seedFinderConfig;
  
  float m_rMax = 200.;
  float m_deltaRMinTopSP = 1.;
  float m_deltaRMaxTopSP = 60.;
  float m_deltaRMinBottomSP = 1.;
  float m_deltaRMaxBottomSP = 60.;
  float m_collisionRegionMin = -250;
  float m_collisionRegionMax = 250.;
  float m_zMin = -2000.;
  float m_zMax = 2000.;
  float m_maxSeedsPerSpM = 1;
  float m_cotThetaMax = 7.40627;  // 2.7 eta
  float m_sigmaScattering = 5;
  float m_radLengthPerSeed = 0.1;
  float m_minPt = 500.;
  float m_bFieldInZ = 0.00199724;
  float m_beamPosX = 0;
  float m_beamPosY = 0;
  float m_impactMax = 3.;
  
};

#endif // PHACTSKDTREESEEDING_H
