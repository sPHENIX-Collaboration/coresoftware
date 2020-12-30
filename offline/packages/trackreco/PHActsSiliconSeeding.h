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
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

class PHCompositeNode;
class PHG4CylinderGeomContainer;
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

/// This is needed by the Acts seedfinder 
inline bool operator==(SpacePoint a, SpacePoint b) {
  return (a.m_hitId == b.m_hitId);
}

using SpacePointPtr = std::unique_ptr<SpacePoint>;
using GridSeeds = std::vector<std::vector<Acts::Seed<SpacePoint>>>;

/**
 * This class runs the Acts seeder over the MVTX measurements
 * to create track stubs for the rest of the stub matching pattern
 * recognition. The module also projects the MVTX stubs to the INTT
 * to find possible matches in the INTT to the MVTX triplet.
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

  /// Output some diagnostic histograms
  void seedAnalysis(bool seedAnalysis)
    { m_seedAnalysis = seedAnalysis; }
   
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
  std::vector<const SpacePoint*> getMvtxSpacePoints();

  /// Perform circle/line fits with the final MVTX seed to get
  /// initial point and momentum estimates for stub matching
  int circleFitSeed(std::vector<TrkrCluster*>& clusters,
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

  /// Projects circle fit to INTT radii to find possible INTT clusters
  /// belonging to MVTX track stub
  std::vector<TrkrDefs::cluskey> findInttMatches(
			const std::vector<TrkrCluster*>& clusters,
			const double R,
			const double X0,
			const double Y0,
			const double B,
			const double m);
  std::vector<TrkrDefs::cluskey> matchInttClusters(const double xProj[],
						   const double yProj[],
						   const double zProj[]);
  void circleCircleIntersection(const double layerRadius, 
				const double circRadius,
				const double circX0,
				const double circY0,
				double& xplus,
				double& yplus,
				double& xminus,
				double& yminus);
  void createSvtxTrack(const double x,
		       const double y,
		       const double z,
		       const double px,
		       const double py,
		       const double pz,
		       const int charge,
		       const std::vector<TrkrCluster*> clusters);
  std::map<const unsigned int, std::vector<TrkrCluster*>>
    makePossibleStubs(std::vector<TrkrCluster*> allClusters);
  void createHistograms();
  void writeHistograms();
  double normPhi2Pi(const double phi);

  std::map<unsigned int, SourceLink> *m_sourceLinks;
  ActsTrackingGeometry *m_tGeometry = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  CluskeyBimap *m_hitIdCluskey;
  TrkrClusterContainer *m_clusterMap = nullptr;
  PHG4CylinderGeomContainer *m_geomContainerIntt = nullptr;
  
  /// Configuration classes for Acts seeding
  Acts::SeedfinderConfig<SpacePoint> m_seedFinderCfg;
  Acts::SpacePointGridConfig m_gridCfg;

  /// Configurable parameters
  /// seed pt has to be in MeV
  float m_minSeedPt = 100;

  /// How many seeds a given hit can be the middle hit of the seed
  /// MVTX can only have the middle layer be the middle hit
  int m_maxSeedsPerSpM = 1;

  /// Limiting location of measurements (e.g. detector constraints)
  /// We limit to the MVTX
  float m_rMax = 50.;
  float m_rMin = 23.;
  float m_zMax = 300.;
  float m_zMin = -300.;
 
  /// max distance between two measurements in one seed
  float m_deltaRMax = 15;
  
  /// Cot of maximum theta angle. Equivalent to eta=1.1 here
  float m_cotThetaMax = 1.335647;
  
  /// B field value in z direction
  /// bfield for space point grid neds to be in kiloTesla
  float m_bField = 1.4 / 1000.;

  std::shared_ptr<Acts::BinFinder<SpacePoint>> 
    m_bottomBinFinder, m_topBinFinder;

  int m_event = 0;

  double m_maxSeedPCA = 0.01;
  
  const static unsigned int m_nInttLayers = 4;
  const double m_nInttLayerRadii[m_nInttLayers] = 
    {7.188, 7.732, 9.680,10.262}; /// cm
  
  /// Search window for phi to match intt clusters in cm
  double m_rPhiSearchWin = 0.1;

  /// Whether or not to use truth clusters in hit lookup
  bool m_useTruthClusters = false;

  bool m_seedAnalysis = false;
  TFile *m_file = nullptr;
  TH1 *h_nMvtxHits = nullptr;
  TH1 *h_nInttHits = nullptr;
  TH2 *h_nHits = nullptr;
  TH1 *h_nSeeds = nullptr;
  TH1 *h_nTotSeeds = nullptr;
  TH1 *h_nInputMeas = nullptr;
  TH1 *h_nInputMvtxMeas = nullptr;
  TH1 *h_nInputInttMeas = nullptr;
  TH2 *h_hits = nullptr;
  TH2 *h_zhits = nullptr;
  TH2 *h_projHits = nullptr;
  TH2 *h_zprojHits = nullptr;
  TH2 *h_resids = nullptr;

};


#endif
