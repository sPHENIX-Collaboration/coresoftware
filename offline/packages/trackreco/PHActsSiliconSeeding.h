#ifndef TRACKRECO_PHACTSSILICONSEEDING_H
#define TRACKRECO_PHACTSSILICONSEEDING_H

#include "ActsTrackingGeometry.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Seeding/Seedfinder.hpp>
#include <Acts/Utilities/Units.hpp>

#include <Acts/Seeding/BinFinder.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>

#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>

#include <boost/bimap.hpp>

#include <string>
#include <map>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxVertexMap;
class TrkrCluster;
class TrkrClusterContainer;

using SourceLink = ActsExamples::TrkrClusterSourceLink;
typedef boost::bimap<TrkrDefs::cluskey, unsigned int> CluskeyBimap;

/**
 * A struct for Acts to take cluster information for seeding
 */
struct SpacePoint {
  unsigned int m_hitId;
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  Acts::GeometryIdentifier m_geoId;
  float m_varianceRphi;
  float m_varianceZ;
  
  unsigned int Id() const { return m_hitId; }

  /// These are needed by Acts
  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  float r() const { return m_r; }

};

inline bool operator==(SpacePoint a, SpacePoint b) {
  return (a.m_hitId == b.m_hitId);
}


using SpacePointPtr = std::unique_ptr<SpacePoint>;
using GridSeeds = std::vector<std::vector<Acts::Seed<SpacePoint>>>;


/**
 * This class runs the Acts seeder over the silicon measurements
 * to create track stubs for the rest of the stub matching pattern
 * recognition
 */
class PHActsSiliconSeeding : public SubsysReco
{
 public:
  PHActsSiliconSeeding(const std::string& name = "PHActsSiliconSeeding");
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  /// Set seeding with truth clusters
  void useTruthClusters(bool useTruthClusters)
    { m_useTruthClusters = useTruthClusters; }

 private:

  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  GridSeeds runSeeder();

  /// Configure the seeding parameters for Acts. There
  /// are a number of tunable parameters for the seeder here
  Acts::SeedfinderConfig<SpacePoint> configureSeeder();
  Acts::SpacePointGridConfig configureSPGrid();
  
  /// Take final seeds and fill the SvtxTrackMap
  void makeSvtxTracks(GridSeeds& seedVector);
  
  /// Create a seeding space point out of an Acts::SourceLink
  SpacePointPtr makeSpacePoint(const unsigned int& hitId,
			    const SourceLink& sl);
  
  /// Get all space points for the seeder
  std::vector<const SpacePoint*> getSpacePoints();

  /// Perform circle/line fits with the final seed to get
  /// initial point and momentum estimates for stub matching
  int circleFitSeed(const std::vector<TrkrCluster*>& clusters,
		     double& x, double& y, double& z,
		     double& px, double& py, double& pz);
  void circleFitByTaubin(const std::vector<TrkrCluster*>& clusters,
			 double& R, double& X0, double& Y0);
  void lineFit(const std::vector<TrkrCluster*>& clusters,
	       double& A, double& B);
  void findRoot(const double R, const double X0, const double Y0,
		double& x, double& y);
  int getCharge(const std::vector<TrkrCluster*>& clusters,
		const double circPhi);
  double normPhi2Pi(const double phi);

  std::map<unsigned int, SourceLink> *m_sourceLinks;
  ActsTrackingGeometry *m_tGeometry = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  CluskeyBimap *m_hitIdCluskey;
  TrkrClusterContainer *m_clusterMap = nullptr;
  
  /// Configuration classes for Acts seeding
  Acts::SeedfinderConfig<SpacePoint> m_seedFinderCfg;
  Acts::SpacePointGridConfig m_gridCfg;

  /// Configurable parameters
  /// seed pt has to be in MeV
  float m_minSeedPt = 100;

  /// How many seeds a given hit can be the middle hit of the seed
  int m_maxSeedsPerSpM = 1;

  /// Limiting location of measurements (e.g. detector constraints)
  float m_rMax = 250.;
  float m_rMin = 20.;
  float m_zMax = 500.;
  float m_zMin = -500.;
 
  /// max distance between two measurements in one seed
  float m_deltaRMax = 50;
  
  /// Cot of maximum theta angle. Equivalent to eta=1.1 here
  float m_cotThetaMax = 1.335647;
  
  /// B field value in z direction
  /// bfield for space point grid neds to be in kiloTesla
  float m_bField = 1.4 / 1000.;

  std::shared_ptr<Acts::BinFinder<SpacePoint>> 
    m_bottomBinFinder, m_topBinFinder;

  int m_event = 0;

  /// Whether or not to use truth clusters in hit lookup
  bool m_useTruthClusters = false;

};


#endif
