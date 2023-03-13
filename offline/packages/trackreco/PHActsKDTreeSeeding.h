// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSKDTREESEEDING_H
#define PHACTSKDTREESEEDING_H

#include <trackbase/SpacePoint.h>

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>

#include <Acts/Seeding/SeedFilterConfig.hpp>
#include <Acts/Seeding/SeedFinderOrthogonalConfig.hpp>

#include <string>

class PHCompositeNode;
class TrkrClusterIterationMapv1;
class ActsGeometry;
class TrkrClusterContainer;
class TrackSeedContainer;
class TrkrCluster;
class PHG4CylinderGeomContainer;

class PHActsKDTreeSeeding : public SubsysReco
{
 public:
  PHActsKDTreeSeeding(const std::string &name = "PHActsKDTreeSeeding");

  ~PHActsKDTreeSeeding() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  
  void useTruthClusters(bool truth) { m_useTruthClusters = truth; }
  void set_cluster_version(int ver) { m_clusterVersion = ver; }
 private:

  Acts::SeedFinderOrthogonalConfig<SpacePoint> configureSeedFinder();
  int getNodes(PHCompositeNode* topNode);
  int createNodes(PHCompositeNode* topNode);
  SeedContainer runSeeder();
  void fillTrackSeedContainer(SeedContainer& seeds);
  std::vector<const SpacePoint*> getMvtxSpacePoints();
  SpacePointPtr makeSpacePoint(const Surface& surf,
			       const TrkrDefs::cluskey key,
			       TrkrCluster* clus);

  /// Projects circle fit to INTT radii to find possible INTT clusters
  /// belonging to MVTX track stub
  void findInttMatches(std::map<TrkrDefs::cluskey, Acts::Vector3>& clusters,
		       TrackSeed& seed);

  void matchInttClusters(std::map<TrkrDefs::cluskey, Acts::Vector3>& clusters,
			 const double xProj[],
			 const double yProj[],
			 const double zProj[]);

  Acts::SeedFilterConfig m_seedFilterConfig;
  Acts::SeedFinderOrthogonalConfig<SpacePoint> m_seedFinderConfig;
  
  /// configured to seed in the MVTX using the middle layer 
  /// as the seed anchor
  /// Defines volume to search for seeds in
  float m_rMax = 200. * Acts::UnitConstants::mm;
  float m_deltaRMinTopSP = 1. * Acts::UnitConstants::mm;
  float m_deltaRMaxTopSP = 20. * Acts::UnitConstants::mm;
  float m_deltaRMinBottomSP = 1. * Acts::UnitConstants::mm;
  float m_deltaRMaxBottomSP = 20. * Acts::UnitConstants::mm;
  float m_collisionRegionMin = -300 * Acts::UnitConstants::mm;
  float m_collisionRegionMax = 300 * Acts::UnitConstants::mm;
  float m_zMin = -300. * Acts::UnitConstants::mm;
  float m_zMax = 300. * Acts::UnitConstants::mm;

  /// max number of seeds a single middle sp can belong to
  float m_maxSeedsPerSpM = 1;
  float m_cotThetaMax = 2.9; 
  float m_sigmaScattering = 5;
  float m_radLengthPerSeed = 0.05;
  float m_minPt = 100.; // MeV
  float m_bFieldInZ = 0.0014; //kTesla
  float m_beamPosX = 0;
  float m_beamPosY = 0;

  /// Maximum transverse PCA allowed
  float m_impactMax = 20. * Acts::UnitConstants::mm;

  /// Middle spacepoint must fall between these two radii
  float m_rMinMiddle = 28. * Acts::UnitConstants::mm;
  float m_rMaxMiddle = 36. * Acts::UnitConstants::mm;

  int m_nIteration = 0;
  std::string m_trackMapName = "SiliconTrackSeedContainer";
  bool m_useTruthClusters = false;
  int m_clusterVersion = 4;
  ClusterErrorPara m_clusErrPara;
  float m_uncfactor = 3.175;
  const static int m_nInttLayers = 4;
  float m_nInttLayerRadii[m_nInttLayers] = {0};
  float m_rPhiSearchWin = 0.1;

  PHG4CylinderGeomContainer* m_geomContainerIntt = nullptr;
  TrkrClusterIterationMapv1* m_iterationMap = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  TrkrClusterContainer* m_clusterMap = nullptr;
  TrackSeedContainer* m_seedContainer = nullptr;
  
};

#endif // PHACTSKDTREESEEDING_H
