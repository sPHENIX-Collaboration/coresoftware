#ifndef TRACKRECO_PHACTSSILICONSEEDING_H
#define TRACKRECO_PHACTSSILICONSEEDING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Seeding/SeedFinder.hpp>

#include <Acts/Seeding/SpacePointGrid.hpp>
#include <Acts/Utilities/GridBinFinder.hpp>

#include <trackbase/SpacePoint.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <map>
#include <string>

class PHCompositeNode;
class PHG4CylinderGeomContainer;
class TrackSeed;
class TrackSeedContainer;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterIterationMap;
class TrkrClusterCrossingAssoc;

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
  PHActsSiliconSeeding(const std::string &name = "PHActsSiliconSeeding");
  ~PHActsSiliconSeeding() override;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void isStreaming()
  {
    m_streaming = true;
  }

  void setStrobeRange(const int low, const int high)
  {
    m_lowStrobeIndex = low;
    m_highStrobeIndex = high;
  }
  void setunc(float unc) { m_uncfactor = unc; }
  /// Set seeding with truth clusters
  void useTruthClusters(bool useTruthClusters)
  {
    m_useTruthClusters = useTruthClusters;
  }

  /// Output some diagnostic histograms
  void seedAnalysis(bool seedAnalysis)
  {
    m_seedAnalysis = seedAnalysis;
  }

  void setinttRPhiSearchWindow(const float win)
  {
    m_inttrPhiSearchWin = win;
  }
  void setinttZSearchWindow(const float &win)
  {
    m_inttzSearchWin = win;
  }
  void setmvtxRPhiSearchWindow(const float win)
  {
    m_mvtxrPhiSearchWin = win;
  }
  void setmvtxZSearchWindow(const float &win)
  {
    m_mvtxzSearchWin = win;
  }
  /// For each MVTX+INTT seed, take the best INTT hits and form
  /// 1 silicon seed per MVTX seed
  void cleanSeeds(bool cleanSeeds)
  {
    m_cleanSeeds = cleanSeeds;
  }

  void rMax(const float rMax)
  {
    m_rMax = rMax;
  }
  void rMin(const float rMin)
  {
    m_rMin = rMin;
  }
  void zMax(const float zMax)
  {
    m_zMax = zMax;
  }
  void zMin(const float zMin)
  {
    m_zMin = zMin;
  }
  void deltaRMax(const float deltaRMax)
  {
    m_deltaRMax = deltaRMax;
  }
  void cotThetaMax(const float cotThetaMax)
  {
    m_cotThetaMax = cotThetaMax;
  }
  void gridFactor(const float gridFactor)
  {
    m_gridFactor = gridFactor;
  }
  void sigmaScattering(const float sigma)
  {
    m_sigmaScattering = sigma;
  }
  void maxPtScattering(const float pt)
  {
    m_maxPtScattering = pt;
  }
  void sigmaError(const float sigma)
  {
    m_sigmaError = sigma;
  }
  void zalign(const float z)
  {
    m_zalign = z;
  }
  void ralign(const float r)
  {
    m_ralign = r;
  }
  void tolerance(const float tolerance)
  {
    m_tolerance = tolerance;
  }
  void helixcut(const float cut)
  {
    m_helixcut = cut;
  }
  void bfield(const float field)
  {
    m_bField = field;
  }
  void minpt(const float pt)
  {
    m_minSeedPt = pt;
  }

  /// A function to run the seeder with large (true)
  /// or small (false) grid spacing
  void largeGridSpacing(const bool spacing);

  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void iteration(int iter) { m_nIteration = iter; }
  void searchInIntt() { m_searchInIntt = true; }
  void strobeWindowLowSearch(const int width) { m_strobeLowWindow = width; }
  void strobeWindowHighSearch(const int width) { m_strobeHighWindow = width; }
 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  
  int m_strobeLowWindow = -1;
  int m_strobeHighWindow = 2;

  void runSeeder();

  /// Configure the seeding parameters for Acts. There
  /// are a number of tunable parameters for the seeder here
  void configureSeeder();
  void configureSPGrid();
  Acts::SeedFilterConfig configureSeedFilter() const;

  /// Take final seeds and fill the TrackSeedContainer
  void makeSvtxTracks(const GridSeeds &seedVector);

  /// Take final seeds and fill the TrackSeedContainer
  void makeSvtxTracksWithTime(const GridSeeds &seedVector, const int &strobe);

  /// Create a seeding space point out of an Acts::SourceLink
  SpacePointPtr makeSpacePoint(
      const Surface &surf,
      const TrkrDefs::cluskey,
      //    const TrkrCluster* clus);
      TrkrCluster *clus);

  /// Get all space points for the seeder
  std::vector<const SpacePoint *> getSiliconSpacePoints(Acts::Extent &rRangeSPExtent,
                                                        const int strobe);
  void printSeedConfigs(Acts::SeedFilterConfig &sfconfig);
  bool isTimingMismatched(TrackSeed& seed) const;
  
      /// Projects circle fit to radii to find possible MVTX/INTT clusters
      /// belonging to track stub
      std::vector<TrkrDefs::cluskey>
      findMatches(
          std::vector<Acts::Vector3> &clusters,
          std::vector<TrkrDefs::cluskey> &keys,
          TrackSeed &seed);

  std::vector<std::vector<TrkrDefs::cluskey>> findMatchesWithTime(
      std::map<TrkrDefs::cluskey, Acts::Vector3> &positions,
      const int &strobe);
  std::vector<std::vector<TrkrDefs::cluskey>> iterateLayers(const int &startLayer,
                                                            const int &endLayer, const int &strobe,
                                                            const std::vector<TrkrDefs::cluskey> &keys,
                                                            const std::vector<Acts::Vector3> &positions);
  std::vector<TrkrDefs::cluskey> matchInttClusters(std::vector<Acts::Vector3> &clusters,
                                                   TrackSeed &seed,
                                                   const double xProj[],
                                                   const double yProj[],
                                                   const double zProj[]);
  short int getCrossingIntt(TrackSeed &si_track);
  std::vector<short int> getInttCrossings(TrackSeed &si_track);

  void createHistograms();
  void writeHistograms();
  double normPhi2Pi(const double phi);
  void clearTreeVariables();

  TrkrClusterCrossingAssoc *_cluster_crossing_map = nullptr;
  TTree *m_tree = nullptr;
  int m_seedid = std::numeric_limits<int>::quiet_NaN();
  std::vector<float> m_mvtxgx = {};
  std::vector<float> m_mvtxgy = {};
  std::vector<float> m_mvtxgz = {};
  std::vector<float> m_mvtxgr = {};
  float m_projgx = std::numeric_limits<float>::quiet_NaN();
  float m_projgy = std::numeric_limits<float>::quiet_NaN();
  float m_projgz = std::numeric_limits<float>::quiet_NaN();
  float m_projgr = std::numeric_limits<float>::quiet_NaN();
  float m_projlx = std::numeric_limits<float>::quiet_NaN();
  float m_projlz = std::numeric_limits<float>::quiet_NaN();
  float m_clusgx = std::numeric_limits<float>::quiet_NaN();
  float m_clusgy = std::numeric_limits<float>::quiet_NaN();
  float m_clusgz = std::numeric_limits<float>::quiet_NaN();
  float m_clusgr = std::numeric_limits<float>::quiet_NaN();
  float m_cluslx = std::numeric_limits<float>::quiet_NaN();
  float m_cluslz = std::numeric_limits<float>::quiet_NaN();

  ActsGeometry *m_tGeometry = nullptr;
  TrackSeedContainer *m_seedContainer = nullptr;
  TrkrClusterContainer *m_clusterMap = nullptr;
  PHG4CylinderGeomContainer *m_geomContainerIntt = nullptr;
  PHG4CylinderGeomContainer *m_geomContainerMvtx = nullptr;
  int m_lowStrobeIndex = 0;
  int m_highStrobeIndex = 1;
  /// Configuration classes for Acts seeding
  Acts::SeedFinderConfig<SpacePoint> m_seedFinderCfg;
  Acts::CylindricalSpacePointGridConfig m_gridCfg;
  Acts::CylindricalSpacePointGridOptions m_gridOptions;
  Acts::SeedFinderOptions m_seedFinderOptions;

  /// boolean whether or not we are going to match the intt clusters
  /// per strobe with crossing information and take all possible matches
  bool m_streaming = false;

  // default to 10 mus
  float m_strobeWidth = 10;
  /// boolean whether or not to include the intt in the acts search windows
  bool m_searchInIntt = false;

  /// Configurable parameters
  /// seed pt has to be in MeV
  float m_minSeedPt = 100 * Acts::UnitConstants::MeV;
  float m_uncfactor = 3.18;

  /// How many seeds a given hit can be the middle hit of the seed
  /// MVTX can only have the middle layer be the middle hit
  int m_maxSeedsPerSpM = 1;

  /// Limiting location of measurements (e.g. detector constraints)
  float m_rMax = 200. * Acts::UnitConstants::mm;
  float m_rMin = 15. * Acts::UnitConstants::mm;
  float m_zMax = 500. * Acts::UnitConstants::mm;
  float m_zMin = -500. * Acts::UnitConstants::mm;

  /// misalignment parameters
  float m_helixcut = 1;
  float m_tolerance = 1.1 * Acts::UnitConstants::mm;
  float m_ralign = 0;
  float m_zalign = 0;
  float m_maxPtScattering = 10;
  float m_sigmaScattering = 5.;
  float m_sigmaError = 5;

  /// Value tuned to provide as large of phi bins as possible.
  /// Increases the secondary finding efficiency
  float m_gridFactor = 2.3809;

  /// max distance between two measurements in one seed
  float m_deltaRMax = 15 * Acts::UnitConstants::mm;
  float m_deltaRMin = 1. * Acts::UnitConstants::mm;
  /// Cot of maximum theta angle
  float m_cotThetaMax = 2.9;

  /// Maximum impact parameter allowed in mm
  float m_impactMax = 20 * Acts::UnitConstants::mm;

  /// Only used in seeding with specified z bin edges, which
  /// is more configuration than we need
  int m_numPhiNeighbors = 1;

  /// B field value in z direction
  /// bfield for space point grid neds to be in kiloTesla
  float m_bField = 1.4 * Acts::UnitConstants::T;
  std::vector<std::pair<int, int>> zBinNeighborsTop;
  std::vector<std::pair<int, int>> zBinNeighborsBottom;
  int nphineighbors = 1;
  std::unique_ptr<const Acts::GridBinFinder<2ul>> m_bottomBinFinder;
  std::unique_ptr<const Acts::GridBinFinder<2ul>> m_topBinFinder;

  int m_event = 0;

  /// Maximum allowed transverse PCA for seed, cm
  double m_maxSeedPCA = 2.;

  /// Search window for phi to match intt clusters in cm
  double m_inttrPhiSearchWin = 0.1;
  float m_inttzSearchWin = 2.0;  // default to one strip width
  double m_mvtxrPhiSearchWin = 0.2;
  float m_mvtxzSearchWin = 0.5;
  /// Whether or not to use truth clusters in hit lookup
  bool m_useTruthClusters = false;

  bool m_cleanSeeds = false;

  int m_nBadUpdates = 0;
  int m_nBadInitialFits = 0;
  TrkrClusterIterationMap *_iteration_map = nullptr;
  int m_nIteration = 0;
  std::string _track_map_name = "SiliconTrackSeedContainer";

  bool m_seedAnalysis = false;
  TFile *m_file = nullptr;
  TH2 *h_nInttProj = nullptr;
  TH1 *h_nMvtxHits = nullptr;
  TH1 *h_nInttHits = nullptr;
  TH1 *h_nMatchedClusters = nullptr;
  TH2 *h_nHits = nullptr;
  TH1 *h_nSeeds = nullptr;
  TH1 *h_nActsSeeds = nullptr;
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
