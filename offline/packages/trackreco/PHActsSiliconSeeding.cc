#include "PHActsSiliconSeeding.h"

#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase/TrkrCluster.h>            
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
 
#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>

PHActsSiliconSeeding::PHActsSiliconSeeding(const std::string& name)
  : SubsysReco(name)
{}

int PHActsSiliconSeeding::Init(PHCompositeNode */*topNode*/)
{
  m_seedFinderCfg = configureSeeder();
  m_gridCfg = configureSPGrid();
  Acts::SeedFilterConfig sfCfg = configureSeedFilter();

  // vector containing the map of z bins in the top and bottom layers
  std::vector<std::pair<int, int> > zBinNeighborsTop;
  std::vector<std::pair<int, int> > zBinNeighborsBottom;
  int nphineighbors = 1;
  m_bottomBinFinder = 
    std::make_shared<Acts::BinFinder<SpacePoint>>(Acts::BinFinder<SpacePoint>(zBinNeighborsBottom, nphineighbors));
  m_topBinFinder = 
    std::make_shared<Acts::BinFinder<SpacePoint>>(Acts::BinFinder<SpacePoint>(zBinNeighborsTop, nphineighbors));

  m_seedFinderCfg.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
     Acts::SeedFilter<SpacePoint>(sfCfg));

  if(m_seedAnalysis)
    {    
      m_file = new TFile(std::string(Name() + ".root").c_str(),"recreate");
      createHistograms();
    }  
  
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::InitRun(PHCompositeNode *topNode)
{

  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }
  
  auto beginend = m_geomContainerIntt->get_begin_end();
  int i = 0;
  for(auto iter = beginend.first; iter != beginend.second; ++iter)
    {
      m_nInttLayerRadii[i] = iter->second->get_radius();
      i++;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::process_event(PHCompositeNode *topNode)
{
  if(_n_iteration>0){
    _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!_iteration_map){
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  auto eventTimer = std::make_unique<PHTimer>("eventTimer");
  eventTimer->stop();
  eventTimer->restart();

  if(Verbosity() > 0)
    std::cout << "Processing PHActsSiliconSeeding event "
	      << m_event << std::endl;

  std::vector<const SpacePoint*> spVec;
  auto seedVector = runSeeder(spVec);

  eventTimer->stop();
  auto seederTime = eventTimer->get_accumulated_time();
  eventTimer->restart();
  
  makeSvtxTracks(seedVector);

  eventTimer->stop();
  auto circleFitTime = eventTimer->get_accumulated_time();

  for(auto sp : spVec)
    { delete sp; }
  spVec.clear();

  if(Verbosity() > 0)
    std::cout << "Finished PHActsSiliconSeeding process_event"
	      << std::endl;

  if(Verbosity() > 0)
    {
      std::cout << "PHActsSiliconSeeding Acts seed time "
		<< seederTime << std::endl;
      std::cout << "PHActsSiliconSeeding circle fit time "
		<< circleFitTime << std::endl;
      std::cout << "PHActsSiliconSeeding total event time " 
		<< circleFitTime + seederTime  << std::endl;
    }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::End(PHCompositeNode */*topNode*/)
{
  if(m_seedAnalysis)
    {
      writeHistograms();
    }

  if(Verbosity() > 1)
    {
      std::cout << "There were " << m_nBadInitialFits 
		<< " bad initial circle fits" << std::endl;
      std::cout << "There were " << m_nBadUpdates 
		<< " bad second circle fits" << std::endl;
    }
     
  return Fun4AllReturnCodes::EVENT_OK;
}

GridSeeds PHActsSiliconSeeding::runSeeder(std::vector<const SpacePoint*>& spVec)
{
  
  Acts::Seedfinder<SpacePoint> seedFinder(m_seedFinderCfg);
  
  /// Covariance converter functor needed by seed finder
  auto covConverter = 
    [=](const SpacePoint& sp, float, float, float)
    -> std::pair<Acts::Vector3, Acts::Vector2> { 
       Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
       Acts::Vector2 cov{sp.m_varianceR, sp.m_varianceZ};
       return std::make_pair(position, cov);
  };

  Acts::Extent rRangeSPExtent;

  spVec = getMvtxSpacePoints(rRangeSPExtent);

  if(m_seedAnalysis)
    { h_nInputMeas->Fill(spVec.size()); }
  
  auto grid = 
    Acts::SpacePointGridCreator::createGrid<SpacePoint>(m_gridCfg);

  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(),
						 spVec.end(),
						 covConverter,
						 m_bottomBinFinder,
						 m_topBinFinder,
						 std::move(grid),
						 m_seedFinderCfg);

  GridSeeds seedVector;
  auto groupIt = spGroup.begin();
  auto endGroup = spGroup.end();
  SeedContainer seeds;
  seeds.clear();
  decltype(seedFinder)::State state;

  for(; !(groupIt == endGroup); ++groupIt)
    {
    
      seedFinder.createSeedsForGroup(state, std::back_inserter(seeds),
				     groupIt.bottom(),
				     groupIt.middle(),
				     groupIt.top(),
				     rRangeSPExtent);

    }
  
  seedVector.push_back(seeds);

  return seedVector;
}

void PHActsSiliconSeeding::makeSvtxTracks(GridSeeds& seedVector)
{
  int numSeeds = 0;
  int numGoodSeeds = 0;
  
  /// Loop over grid volumes
  for(auto& seeds : seedVector)
    {
      /// Loop over actual seeds in this grid volume
      for(auto& seed : seeds)
	{
	  if(Verbosity() > 1) {
	    std::cout << "Seed " << numSeeds << " has "
		      << seed.sp().size() << " measurements " 
		      << std::endl;
	  }

	  numSeeds++;

	  std::vector<TrkrCluster*> clusters;
	  std::vector<Acts::Vector3> globalPositions;
	  ActsTransformations transformer;
	  for(auto& spacePoint : seed.sp())
	    {
	      auto cluskey = spacePoint->m_clusKey;
	      clusters.push_back(m_clusterMap->findCluster(cluskey));

	      globalPositions.push_back(transformer.getGlobalPosition(
				        m_clusterMap->findCluster(cluskey),
					m_surfMaps, m_tGeometry));			      

	      if(Verbosity() > 1) {
		std::cout << "Adding cluster with x,y "
			  << spacePoint->x() <<", " << spacePoint->y()
			  << " mm in detector " 
			  << TrkrDefs::getTrkrId(cluskey)
			  << std::endl;
	      }
	    }
	 
	  double x = NAN, y = NAN, z = seed.z() / Acts::UnitConstants::cm;
	  double px, py, pz;
	  
	  auto fitTimer = std::make_unique<PHTimer>("trackfitTimer");
	  fitTimer->stop();
	  fitTimer->restart();

	  /// Performs circle fit and extrapolates to INTT layers to
	  /// to get additional clusters in this track seed
	  int charge = circleFitSeed(clusters, globalPositions, 
				     x, y, z,
				     px, py, pz);
	  fitTimer->stop();
	  auto circlefittime = fitTimer->get_accumulated_time();
	  fitTimer->restart();
	  /// Bad seed, if x is nan so are y and z
	  if(std::isnan(x))
	    {
	      m_nBadInitialFits++;
	      continue;
	    }

	  numGoodSeeds++;
	  
	  createSvtxTrack(x, y, seed.z() / Acts::UnitConstants::cm,
			  px, py, pz, charge,
			  clusters, globalPositions);
	  fitTimer->stop();
	  auto svtxtracktime = fitTimer->get_accumulated_time();
	  if(Verbosity() > 0)
	    {
	      std::cout << "Circle fit time " << circlefittime << " and svtx time "
		    << svtxtracktime << std::endl;
	    }
	}
    }

  if(m_seedAnalysis)
    {
      h_nSeeds->Fill(numGoodSeeds);
      h_nActsSeeds->Fill(numSeeds);
    }
  if(Verbosity() > 1)
    {
      std::cout << "Total number of seeds found in " 
		<< seedVector.size() << " volume regions gives " 
		<< numSeeds << " Acts seeds " << std::endl;
    }

  return;
  
}
 
void PHActsSiliconSeeding::createSvtxTrack(const double x,
					   const double y,
					   const double z,
					   const double px,
					   const double py,
					   const double pz,
					   const int charge,
					   std::vector<TrkrCluster*>& clusters,
					   std::vector<Acts::Vector3>& clusGlobPos)
{
  auto fitTimer = std::make_unique<PHTimer>("trackfitTimer");
  fitTimer->stop();
  fitTimer->restart();

  auto stubs = makePossibleStubs(clusters, clusGlobPos);

  fitTimer->stop();
  auto possibleStubs = fitTimer->get_accumulated_time();
  fitTimer->restart();

  int numSeedsPerActsSeed = 0;
  
  double trackX = x;
  double trackY = y;
  double trackZ = z;
  double trackPx = px;
  double trackPy = py;
  double trackPz = pz;
  double trackCharge = charge;
  double trackPhi = atan2(py,px);
  double trackEta = atanh(pz / sqrt(px * px + py * py + pz * pz));

  /// Make a track for every stub that was constructed
  /// We use the same xyz and pxpypz given by the mvtx circle
  /// fit since that is the "anchor" for the stub

  if(stubs.size() > 1 and  m_cleanSeeds)
    {
      stubs = identifyBestSeed(stubs);
    }

  for(const auto& [stub, stubClusterPairs] : stubs)
    {
      std::vector<TrkrCluster*> stubClusters = stubClusterPairs.first;
      std::vector<Acts::Vector3> stubClusterPositions = stubClusterPairs.second;

      int nMvtx = 0;
      int nIntt = 0;
      numSeedsPerActsSeed++;
       
      auto svtxTrack = std::make_unique<SvtxTrack_v2>(); 

      svtxTrack->set_id(m_trackMap->size());
      
      for(const auto clus : stubClusters) 
	{
	  const auto cluskey = clus->getClusKey();
	  svtxTrack->insert_cluster_key(cluskey);
	  if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::mvtxId)
	    { nMvtx++; }
	  else if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::inttId)
	    { nIntt++; }	 

	  if(Verbosity() > 2)
	    { std::cout << "Cluskey adding : " << cluskey << std::endl; }
	}
 
      /// Get a less rough estimate of R, and thus, p
      double R, X0, Y0;
      circleFitByTaubin(stubClusterPositions, R, X0, Y0);
      
      /// 0.3 conversion factor, 1.4=B field, 
      /// 100 convert R from cm to m
      float pt = 0.3 * 1.4 * R / 100.;
  
      trackPx = pt * cos(trackPhi);
      trackPy = pt * sin(trackPhi);
      trackPz = pt * sinh(trackEta);
      
      /// Diagnostic
      if(m_seedAnalysis)
	{
	  h_nInttHits->Fill(nIntt);
	  h_nMvtxHits->Fill(nMvtx);
	  h_nHits->Fill(nMvtx, nIntt);
	}

      if(Verbosity() > 1) {
	std::cout << "Setting silicon seed with track id " 
		  << m_trackMap->size()
		  << " and (x,y,z) = " 
		  << trackX << ", " << trackY << ", " << trackZ
		  << std::endl << " and (px,py,pz) " << trackPx 
		  << ", " << trackPy << ", " << trackPz << std::endl
		  << " with charge " << trackCharge << std::endl;
      }

      svtxTrack->set_x(trackX);
      svtxTrack->set_y(trackY);
      svtxTrack->set_z(trackZ);
      svtxTrack->set_px(trackPx);
      svtxTrack->set_py(trackPy);
      svtxTrack->set_pz(trackPz);
      svtxTrack->set_charge(trackCharge);
      
      m_trackMap->insert(svtxTrack.get());
  
    }
  
  fitTimer->stop();
  auto makeTracks = fitTimer->get_accumulated_time();
  fitTimer->restart();


  if(m_seedAnalysis)
    { h_nTotSeeds->Fill(numSeedsPerActsSeed); }

  if(Verbosity() > 0)
    {
      std::cout << "make stub time " << possibleStubs << " create track time "
		<< makeTracks << std::endl;
    }

  if(Verbosity() > 1)
    {
      std::cout << "Found " << numSeedsPerActsSeed << " seeds for one Acts seed"
		<< std::endl;
    }
}

std::map<const unsigned int, std::pair<std::vector<TrkrCluster*>,
				       std::vector<Acts::Vector3>>>
     PHActsSiliconSeeding::identifyBestSeed(
	  std::map<const unsigned int, 
	  std::pair<std::vector<TrkrCluster*>,
	            std::vector<Acts::Vector3>>> allSeeds)
{
  std::map<const unsigned int, std::pair<std::vector<TrkrCluster*>,
					 std::vector<Acts::Vector3>>> returnStub;

  double firstLayerBestResidual = std::numeric_limits<double>::max();
  double secondLayerBestResidual = std::numeric_limits<double>::max();
  TrkrDefs::cluskey firstlayerkey = std::numeric_limits<unsigned long long>::max();
  TrkrDefs::cluskey secondlayerkey = std::numeric_limits<unsigned long long>::max();
  Acts::Vector3 firstlayerpos, secondlayerpos;
  
  std::vector<Acts::Vector3> mvtxClusters;
  std::vector<TrkrCluster*> mvtxTrkrClusters;
  double R=NAN, X0=NAN, Y0=NAN;

  /// Just need to fit the mvtx once since all seeds have the same mvtx hits
  for(const auto& [key, clusterPairs] : allSeeds)
    {
      auto clusterVec = clusterPairs.first;
      auto clusterPosVec = clusterPairs.second;

      /// Circle fit the mvtx triplet only
      for(int i=0; i<clusterVec.size(); i++)
	{
	  if(TrkrDefs::getTrkrId(clusterVec.at(i)->getClusKey()) == TrkrDefs::mvtxId)
	    { 
	      mvtxClusters.push_back(clusterPosVec.at(i)); 
	      mvtxTrkrClusters.push_back(clusterVec.at(i));
	    }
	}

      circleFitByTaubin(mvtxClusters, R, X0, Y0);
         
      break;
    }
      
  /// If mvtx seed wasn't found, just return
  if(mvtxClusters.size() == 0)
    { return returnStub; }
  
  for(const auto& [key, clusterPairs] : allSeeds)
    {
      auto clusterVec = clusterPairs.first;
      auto clusterPosVec = clusterPairs.second;
      
      for(int i=0; i<clusterVec.size(); i++)
	{
	  if(TrkrDefs::getTrkrId(clusterVec.at(i)->getClusKey()) != TrkrDefs::inttId)
	    { continue; }

	  Acts::Vector3 globalPos = clusterPosVec.at(i);

	  double residual = sqrt( pow( globalPos(0) - X0, 2) +
				  pow( globalPos(1) - Y0, 2)) - R;
	  
	  if(Verbosity() > 2) {
	    std::cout << "Residual for cluster " << clusterVec.at(i)->getClusKey()
		      << " and position " << globalPos(0) << ", " 
		      << globalPos(1) << " is " << residual << std::endl;
	  }

	  double r = sqrt( pow(globalPos(0), 2) +
	                   pow(globalPos(1), 2));
	
	  if(r < m_nInttLayerRadii[2] and residual < firstLayerBestResidual)
	    {
	      firstLayerBestResidual = residual;
	      firstlayerkey = clusterVec.at(i)->getClusKey();
	      firstlayerpos = clusterPosVec.at(i);
	    }
	  if( r > m_nInttLayerRadii[2] and residual < secondLayerBestResidual)
	    {
	      secondLayerBestResidual = residual;
	      secondlayerkey = clusterVec.at(i)->getClusKey();
	      secondlayerpos = clusterPosVec.at(i);
	    }    
	}
    }

  // Form the set of clusters that was identified as the smallest residuals
  std::vector<TrkrCluster*> bestClusters = mvtxTrkrClusters;
  std::vector<Acts::Vector3> bestClusterPos = mvtxClusters;
  if(firstlayerkey < std::numeric_limits<unsigned long long>::max())
    { 
      bestClusters.push_back(m_clusterMap->findCluster(firstlayerkey)); 
      bestClusterPos.push_back(firstlayerpos);
    }
  if(secondlayerkey < std::numeric_limits<unsigned long long>::max())
    { 
      bestClusters.push_back(m_clusterMap->findCluster(secondlayerkey)); 
      bestClusterPos.push_back(secondlayerpos);
    }

  if(Verbosity() > 2)
    {
      std::cout << "Best cluster set is " << std::endl;
      for(auto& cluster : bestClusters)
	std::cout << cluster->getClusKey() << std::endl;
    }

  /// Just make it the 0th entry since it is the only one
  returnStub.insert(std::make_pair(0, std::make_pair(bestClusters,bestClusterPos)));

  return returnStub;
}

std::map<const unsigned int, std::pair<std::vector<TrkrCluster*>, std::vector<Acts::Vector3>>> 
PHActsSiliconSeeding::makePossibleStubs(std::vector<TrkrCluster*>& allClusters,
					std::vector<Acts::Vector3>& clusGlobPos)
{

  std::vector<TrkrCluster*> mvtxClusters;
  std::vector<TrkrCluster*> inttFirstLayerClusters;
  std::vector<TrkrCluster*> inttSecondLayerClusters;
  std::vector<Acts::Vector3> mvtxClusPos;
  std::vector<Acts::Vector3> inttFirstLayerClusPos;
  std::vector<Acts::Vector3> inttSecondLayerClusPos;

  std::map<const unsigned int, std::pair<std::vector<TrkrCluster*>,std::vector<Acts::Vector3>>> stubs;
  unsigned int combo = 0;

  for(int i = 0; i < allClusters.size(); i++)
    {
      const auto cluskey = allClusters.at(i)->getClusKey();
      const auto globalPos = clusGlobPos.at(i);
     
      if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::mvtxId)
	{ 
	  mvtxClusters.push_back(allClusters.at(i)); 
	  mvtxClusPos.push_back(globalPos);
	}
      else
	{
	  const double r = sqrt(pow(globalPos(0), 2) +
				pow(globalPos(1), 2));
	  if(r < 8.) 
	    { 
	      inttFirstLayerClusters.push_back(allClusters.at(i)); 
	      inttFirstLayerClusPos.push_back(globalPos);
	    }
	  else
	    { 
	      inttSecondLayerClusters.push_back(allClusters.at(i)); 
	      inttSecondLayerClusPos.push_back(globalPos);
	    }
	}
    }
  
  /// Add the INTT measurements to the track stub
  std::vector<TrkrCluster*> dumVec = mvtxClusters;
  std::vector<Acts::Vector3> dumclusVec = mvtxClusPos;
  for(int i=0; i<inttFirstLayerClusters.size(); i++)
    {
      dumVec.push_back(inttFirstLayerClusters.at(i));
      dumclusVec.push_back(inttFirstLayerClusPos.at(i));
    }
  for(int i=0; i<inttSecondLayerClusters.size(); i++) {
    dumVec.push_back(inttSecondLayerClusters.at(i));
    dumclusVec.push_back(inttSecondLayerClusPos.at(i));
  }

  stubs.insert(std::make_pair(combo, std::make_pair(dumVec, dumclusVec)));
  
  return stubs;
}

int PHActsSiliconSeeding::circleFitSeed(std::vector<TrkrCluster*>& clusters,
					std::vector<Acts::Vector3>& clusGlobPos,
					double& x, double& y, double& z,
					double& px, double& py, double& pz)
{
  if(Verbosity() > 2) {
    for(const auto clus : clusters)
      std::cout << "Evaluating cluster : " << clus->getClusKey()
		<< std::endl;
  }
  

  /// Circle radius at x,y center
  /// Note - units are sPHENIX cm since we are using TrkrClusters
  double R, X0, Y0;
  circleFitByTaubin(clusGlobPos,R, X0, Y0);
  
  if(Verbosity() > 2)
    std::cout << "Circle R, X0, Y0 : " << R << ", " << X0
	      << ", " << Y0 << std::endl;

  findRoot(R, X0, Y0, x, y);
  
  /// If the xy position is O(100s) microns, the initial vertex 
  /// finder will throw an eigen stepper error trying to propagate 
  /// from the PCA. These  are likely bad seeds anyway since the 
  /// MVTX has position resolution O(5) microns. Units are cm
  
  if(fabs(x) > m_maxSeedPCA or fabs(y) > m_maxSeedPCA)
    {
      if(Verbosity() > 1)
	std::cout << "x,y circle fit : " << x << ", " 
		  << y << std::endl;

      x = NAN;
      y = NAN;
      /// Return statement doesn't matter as x = nan will be caught
      return 1;
    }

  int charge = getCharge(clusGlobPos, atan2(Y0,X0));
  
  /// Now determine the line tangent to the circle at this point to get phi
  /// The slope of the line connecting the circle center and PCA is 
  /// m = (y0-y)/(x0-x). So the perpendicular slope (i.e. phi) is then -1/m
  /// For some reason the phi value comes back from atan2 off by 
  /// a factor of pi for positive charged tracks, hence the check
  
  double phi = atan2(-1 * (X0-x), Y0-y);
  if(charge > 0)
    {
      phi += M_PI;
      if(phi > M_PI) 
	phi -= 2. * M_PI;
    }
 
  if(Verbosity() > 2)
    std::cout << "track seed phi : " << phi <<  std::endl;

  double m, B;
  
  /// m is slope as a function of radius, B is z intercept (vertex)
  
  lineFit(clusGlobPos, m, B);
  z = B;
  
  double theta = atan(1./m);

  /// normalize to 0 < theta < pi
  if(theta < 0)
    theta += M_PI;

  if(Verbosity() > 2)
    std::cout << "Track seed theta: " << theta << std::endl;
 
  /// 0.3 conversion factor, 1.4=B field, 100 convert R from cm to m
  /// Get a very rough estimate of p
  float pt = 0.3 * 1.4 * R / 100.;
  float eta = -log(tan(theta/2.));
  float p = pt * cosh(eta);

  /// The only thing that is really needed for the propagation
  /// is the direction
  px = p * sin(theta) * cos(phi);
  py = p * sin(theta) * sin(phi);
  pz = p * cos(theta);
  
  if(Verbosity() > 2)
    {
      std::cout << "Momentum vector estimate: (" << px <<" , " 
		<< py << ", " << pz << ") " << std::endl;
    }
  
  auto fitTimer = std::make_unique<PHTimer>("inttMatchTimer");
  fitTimer->stop();
  fitTimer->restart();

  /// Project to INTT and find matches
  auto additionalClusters = findInttMatches(clusGlobPos, R, X0, Y0, z, m);
  
  /// Add possible matches to cluster list to be parsed when
  /// Svtx tracks are made
  for(auto& cluskey : additionalClusters)
    { clusters.push_back(m_clusterMap->findCluster(cluskey)); }
  
  fitTimer->stop();
  auto addClusters = fitTimer->get_accumulated_time();

  if(Verbosity() > 0)
    {
      std::cout << "find intt clusters time " << addClusters << std::endl;
    }
  
  return charge;

}

std::vector<TrkrDefs::cluskey> PHActsSiliconSeeding::findInttMatches(
		               std::vector<Acts::Vector3>& clusters,
			       const double R,
			       const double X0,
			       const double Y0,
			       const double B,
			       const double m)
{

  double xProj[m_nInttLayers];
  double yProj[m_nInttLayers];
  double zProj[m_nInttLayers];
  
  /// Diagnostic 
  if(m_seedAnalysis)
    {
      for(const auto& glob : clusters)
	{
	  h_hits->Fill(glob(0), glob(1));
	  h_zhits->Fill(glob(2),
			sqrt(pow(glob(0),2) + pow(glob(1),2)));
	  h_projHits->Fill(glob(0), glob(1));
	  h_zprojHits->Fill(glob(2),
			    sqrt(pow(glob(0),2) + pow(glob(1),2)));
	}
    }
  
  

  /// Project the seed to the INTT to find matches
  for(int layer = 0; layer < m_nInttLayers; ++layer)
    {
      double xplus = 0;
      double yplus = 0;
      double xminus = 0;
      double yminus = 0;
      circleCircleIntersection(m_nInttLayerRadii[layer],
			       R, X0, Y0, xplus, yplus,
			       xminus, yminus);
      
      /// If there are no real solutions to the intersection, skip
      if(std::isnan(xplus))
	{
	  if(Verbosity() > 2)
	    {
	      std::cout << "Circle intersection calc failed, skipping" 
			<< std::endl;
	      std::cout << "layer radius " << m_nInttLayerRadii[layer] 
			<< " and circ rad " << R << " with center " << X0 
			<< ", " << Y0 << std::endl;
	    }
	  continue;
	}
      
      /// Figure out which solution is correct based on the position 
      /// of the last layer in the mvtx seed
      const auto lastclusglob = clusters.at(clusters.size() - 1);
      const double lastClusPhi = atan2(lastclusglob(1), lastclusglob(0));
      const double plusPhi = atan2(yplus, xplus);
      const double minusPhi = atan2(yminus, xminus);
      
      if(fabs(lastClusPhi - plusPhi) < fabs(lastClusPhi - minusPhi))
	{
	  xProj[layer] = xplus;
	  yProj[layer] = yplus;
	}
      else
	{
	  xProj[layer] = xminus;
	  yProj[layer] = yminus;
	}
      
      zProj[layer] = m * m_nInttLayerRadii[layer] + B;
      
      if(m_seedAnalysis) {
	h_projHits->Fill(xProj[layer], yProj[layer]);
	h_zprojHits->Fill(zProj[layer], sqrt(pow(xProj[layer],2) + 
					     pow(yProj[layer],2)));
      }

      if(Verbosity() > 2)
	{
	  std::cout << "Projected point is : " << xProj[layer] << ", "
		    << yProj[layer] << ", " << zProj[layer] << std::endl;
	}
    }
  
  return matchInttClusters(clusters, xProj, yProj, zProj);
}

std::vector<TrkrDefs::cluskey> PHActsSiliconSeeding::matchInttClusters(
			       std::vector<Acts::Vector3>& clusters,
			       const double xProj[],
			       const double yProj[],
			       const double zProj[])
{
  std::vector<TrkrDefs::cluskey> matchedClusters;
  ActsTransformations transform;

  for(int inttlayer = 0; inttlayer < m_nInttLayers; inttlayer++)
    {
      auto hitsetrange = m_hitsets->getHitSets(TrkrDefs::TrkrId::inttId, inttlayer+3);
      const double projR = sqrt(pow(xProj[inttlayer], 2) + 
				pow(yProj[inttlayer], 2));
      const double projPhi = atan2(yProj[inttlayer], xProj[inttlayer]);
      const double projRphi = projR * projPhi;

      for (auto hitsetitr = hitsetrange.first;
	   hitsetitr != hitsetrange.second;
	   ++hitsetitr)
	{
	  const int ladderzindex = InttDefs::getLadderZId(hitsetitr->first);
	  const int ladderphiindex = InttDefs::getLadderPhiId(hitsetitr->first);
	  double ladderLocation[3] = {0.,0.,0.};

	  // Add three to skip the mvtx layers for comparison
	  // to projections
	  auto layerGeom = dynamic_cast<CylinderGeomIntt*>
	    (m_geomContainerIntt->GetLayerGeom(inttlayer+3));
	  
	  layerGeom->find_segment_center(ladderzindex, ladderphiindex, ladderLocation);
	  const double ladderphi = atan2(ladderLocation[1], ladderLocation[0]) + layerGeom->get_strip_phi_tilt();
	  const auto stripZSpacing = layerGeom->get_strip_z_spacing();
	  float dphi = ladderphi - projPhi;
	  if(dphi > M_PI)
	    { dphi -= 2. * M_PI; }
	  else if (dphi < -1 * M_PI)
	    { dphi += 2. * M_PI; }
	  
	  /// Check that the projection is within some reasonable amount of the segment
	  /// to reject e.g. looking at segments in the opposite hemisphere. This is about
	  /// the size of one intt segment (256 * 80 micron strips in a segment)
	  if(fabs(dphi) > 0.2)
	    { continue; }

	  TVector3 projectionLocal(0,0,0);
	  TVector3 projectionGlobal(xProj[inttlayer],yProj[inttlayer],zProj[inttlayer]);
	  projectionLocal = layerGeom->get_local_from_world_coords(ladderzindex, 
								   ladderphiindex,
								   projectionGlobal);

	  auto range = m_clusterMap->getClusters(hitsetitr->first);	
	  for(auto clusIter = range.first; clusIter != range.second; ++clusIter )
	    {
	      const auto cluskey = clusIter->first;
	      const auto cluster = clusIter->second;
	      /// Diagnostic
	      if(m_seedAnalysis)
		{ 
		  const auto globalP = transform.getGlobalPosition(cluster, m_surfMaps,
								   m_tGeometry);
		  h_nInttProj->Fill(projectionLocal[1] - cluster->getLocalX(),
				    projectionLocal[2] - cluster->getLocalY()); 
		  h_hits->Fill(globalP(0), globalP(1));
		  h_zhits->Fill(globalP(2),
				sqrt(pow(globalP(0),2)+pow(globalP(1),2)));
		  
		  h_resids->Fill(zProj[inttlayer] - cluster->getLocalY(),
				 projRphi - cluster->getLocalX());
		}
	     
	      /// Z strip spacing is the entire strip, so because we use fabs
	      /// we divide by two
	      if(fabs(projectionLocal[1] - cluster->getLocalX()) < m_rPhiSearchWin and
		 fabs(projectionLocal[2] - cluster->getLocalY()) < stripZSpacing / 2.)
		{
	
		  matchedClusters.push_back(cluskey);
		  /// Cache INTT global positions with seed
		  const auto globalPos = transform.getGlobalPosition(cluster, m_surfMaps,
								     m_tGeometry);
		  clusters.push_back(globalPos);

		 		  	      
		  if(Verbosity() > 4)
		    {
		      std::cout << "Matched INTT cluster with cluskey " << cluskey 
				<< " and position " << globalPos.transpose() 
				<< std::endl << " with projections rphi "
				<< projRphi << " and inttclus rphi " << cluster->getLocalX()
				<< " and proj z " << zProj[inttlayer] << " and inttclus z "
				<< cluster->getLocalY() << " in layer " << inttlayer 
				<< " with search windows " << m_rPhiSearchWin 
				<< " in rphi and strip z spacing " << stripZSpacing 
				<< std::endl;
		    }
		}
	    }
	}  
    }
  
  if(m_seedAnalysis) {
    h_nMatchedClusters->Fill(matchedClusters.size());
  }

  return matchedClusters;
}
void PHActsSiliconSeeding::circleCircleIntersection(const double layerRadius,
						    const double circRadius,
						    const double circX0,
						    const double circY0,
						    double& xplus,
						    double& yplus,
						    double& xminus,
						    double& yminus)
{
  /// Solutions to the circle intersection are (xplus, yplus) and 
  /// (xminus, yminus). The intersection of the two circles occurs when
  /// (x-x1)^2 + (y-y1)^2 = r1^2,  / (x-x2)^2 + (y-y2)^2 = r2^2
  /// Here we assume that circle 1 is an sPHENIX layer centered on x1=y1=0, 
  /// and circle 2 is arbitrary such that they are described by
  ///  x^2 +y^2 = r1^2,   (x-x0)^2 + (y-y0)^2 = r2^2
  /// expand the equations and subtract to eliminate the x^2 and y^2 terms, 
  /// gives the radial line connecting the intersection points
  /// iy = - (2*x2*x - D) / 2*y2, 
  /// then substitute for y in equation of circle 1

  double D = layerRadius*layerRadius - circRadius*circRadius + circX0*circX0 + circY0*circY0;
  double a = 1.0 + (circX0*circX0) / (circY0*circY0);
  double b = - D * circX0/( circY0*circY0);
  double c = D*D / (4.0*circY0*circY0) - layerRadius*layerRadius;

  xplus = (-b + sqrt(b*b - 4.0* a * c) ) / (2.0 * a);
  xminus = (-b - sqrt(b*b - 4.0* a * c) ) / (2.0 * a);

  // both values of x are valid
  // but for each of those values, there are two possible y values on circle 1
  // but only one of those falls on the radical line:

  yplus = - (2*circX0*xplus - D) / (2.0*circY0); 
  yminus = -(2*circX0*xminus - D) / (2.0*circY0);
}

void PHActsSiliconSeeding::findRoot(const double R, const double X0,
				    const double Y0, double& x,
				    double& y)
{
  /**
   * We need to determine the closest point on the circle to the origin
   * since we can't assume that the track originates from the origin
   * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to 
   * minimize d = sqrt((0-x)^2+(0-y)^2), the distance between the 
   * origin and some (currently, unknown) point on the circle x,y.
   * 
   * Solving the circle eqn for x and substituting into d gives an eqn for
   * y. Taking the derivative and setting equal to 0 gives the following 
   * two solutions. We take the smaller solution as the correct one, as 
   * usually one solution is wildly incorrect (e.g. 1000 cm)
   */
  
  double miny = (sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double miny2 = (-sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double minx = sqrt(pow(R, 2) - pow(miny - Y0, 2)) + X0;
  double minx2 = -sqrt(pow(R, 2) - pow(miny2 - Y0, 2)) + X0;
  
  if(Verbosity() > 1)
    std::cout << "minx1 and x2 : " << minx << ", " << minx2 << std::endl
	      << "miny1 and y2 : " << miny << ", " << miny2 << std::endl;

  /// Figure out which of the two roots is actually closer to the origin
  if(fabs(minx) < fabs(minx2))
    x = minx;
  else
    x = minx2;

  if(fabs(miny) < fabs(miny2))
    y = miny;
  else
    y = miny2;
  
  if(Verbosity() > 1)
    {
      std::cout << "Minimum x and y positions " << x << ",  " 
		<< y << std::endl;
    }

}

int PHActsSiliconSeeding::getCharge(const std::vector<Acts::Vector3>& globalPos,
				    const double circPhi)
{

  /**
   * If the circle center phi is positioned clockwise to the seed phi, 
   * the seed is positively charged. If the circle center phi is positioned
   * counter clockwise, the seed is negatively charged
   */

  int charge = 0;
  
  /// Get a crude estimate of the seed phi by taking the average of the
  /// measurements
  double trackPhi = 0;

  for(const auto& pos : globalPos)
    {
      double clusPhi = atan2(pos(1), pos(0));

      /// if it is close to the periodic boundary normalize to 
      /// two pi to avoid -pi and pi issues
      if(fabs(fabs(clusPhi) - M_PI) < 0.2)
	clusPhi = normPhi2Pi(clusPhi);
      trackPhi += clusPhi;
    }

  trackPhi /= globalPos.size();

  /// normalize back
  if(trackPhi > M_PI)
    trackPhi -= 2. * M_PI;

  float quadrants[5] = {-M_PI,-M_PI / 2., 0, M_PI/2., M_PI};
  int quadrant = -1;
  for(int i=0; i<4; i++)
    {
      if(trackPhi > quadrants[i] && trackPhi <= quadrants[i+1])
	{
	  quadrant = i;
	  break;
	}
    }

  if(quadrant == -1)
    std::cout << "quadrant was not set... shouldn't be possible"
	      << std::endl;

  if(quadrant == 1 or quadrant == 2)
    {
      if(circPhi > trackPhi)
	charge = -1;
      else
	charge = 1;
    }
  else
    {
      /// Shift the periodic boundary to make quadrants 0 and 3 away
      /// from boundary
      double normTrackPhi = normPhi2Pi(trackPhi);
      double normCircPhi = normPhi2Pi(circPhi);
  
      if(normCircPhi > normTrackPhi)
	charge = -1;
      else
	charge = 1;
    }

  if(Verbosity() > 1)
    std::cout << "Track seed charge determined to be " 
	      << charge << " in quadrant " << quadrant << std::endl;

  return charge;

}

void PHActsSiliconSeeding::lineFit(const std::vector<Acts::Vector3>& globPos, 
				   double &A, double &B)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  double xsum = 0,x2sum = 0,ysum = 0,xysum = 0;    
  for(const auto& pos : globPos)
    {
      double z = pos(2);
      double r = sqrt(pow(pos(0),2) + pow(pos(1), 2));
      
      xsum=xsum+r;               // calculate sigma(xi)
      ysum=ysum+z;               // calculate sigma(yi)
      x2sum=x2sum+pow(r,2);      // calculate sigma(x^2i)
      xysum=xysum+r*z;           // calculate sigma(xi*yi)
    
      if(Verbosity() > 4)
	{
	  double r = sqrt(pow(pos(0),2) + pow(pos(1), 2));          
	  std::cout << " r " << r << " z " << pos(2)
		    << std::endl; 
	}    
    }
  
  /// calculate slope
  A = (globPos.size()*xysum-xsum*ysum) / (globPos.size()*x2sum-xsum*xsum);

  /// calculate intercept
  B = (x2sum*ysum-xsum*xysum) / (x2sum*globPos.size()-xsum*xsum);
  
  if(Verbosity() > 4)
    {
      for(int i =0; i <m_nInttLayers; i++)
	{
	  std::cout << "intt z_fit layer " << i << " is " 
		    << A * m_nInttLayerRadii[i] + B << std::endl;
	}
    }
  
  return;
}   

void PHActsSiliconSeeding::circleFitByTaubin(const std::vector<Acts::Vector3>& globalPositions,
					     double& R, double& X0, double& Y0)
{
  /**  
   *   Circle fit to a given set of data points (in 2D)
   *   This is an algebraic fit, due to Taubin, based on the journal article
   *   G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
   *               Space Curves Defined By Implicit Equations, With 
   *               Applications To Edge And Range Image Segmentation",
   *               IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
   *  It works well whether data points are sampled along an entire circle 
   *  or along a small arc. 
   *  It still has a small bias and its statistical accuracy is slightly lower 
   *  than that of the geometric fit (minimizing geometric distances),
   *  It provides a very good initial guess for a subsequent geometric fit. 
   *    Nikolai Chernov  (September 2012)
   */
  
  int iter, IterMAX=99;
  
  double Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
  double A0, A1, A2, A22, A3, A33;
  double x, y;
  double DET, Xcenter, Ycenter;
  
  // Compute x- and y- sample means   
  double meanX = 0;
  double meanY = 0;
  double weight = 0;

  for(auto& globalPos : globalPositions)
    {
      meanX += globalPos(0);
      meanY += globalPos(1);
      weight++;
    }

  meanX /= weight;
  meanY /= weight;

  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

  for(auto& pos : globalPositions)
    {
      
      double Xi = pos(0) - meanX;
      double Yi = pos(1) - meanY;
      double Zi = Xi * Xi + Yi * Yi;

      Mxy += Xi*Yi;
      Mxx += Xi*Xi;
      Myy += Yi*Yi;
      Mxz += Xi*Zi;
      Myz += Yi*Zi;
      Mzz += Zi*Zi;
    }

  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;

  Mz = Mxx + Myy;
  Cov_xy = Mxx * Myy - Mxy * Mxy;
  Var_z = Mzz - Mz * Mz;
  A3 = 4 * Mz;
  A2 = -3 * Mz * Mz - Mzz;
  A1 = Var_z * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
  A0 = Mxz * (Mxz * Myy - Myz * Mxy) + 
    Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;

  for (x=0., y=A0, iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
      double Dy = A1 + x * (A22 + A33 * x);
      double xnew = x - y / Dy;
      if ((xnew == x)||(!std::isfinite(xnew))) break;
      double ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
      if (fabs(ynew)>=fabs(y))  break;
      x = xnew;  y = ynew;
    }
  
  //  computing parameters of the fitting circle
  
  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
  
  //  assembling the output
  
  X0 = Xcenter + meanX;
  Y0 = Ycenter + meanY;
  R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);

}

SpacePointPtr PHActsSiliconSeeding::makeSpacePoint(const Surface& surf,
						   const TrkrCluster* clus)
{
  Acts::Vector2 localPos(clus->getLocalX() * Acts::UnitConstants::cm, 
			 clus->getLocalY() * Acts::UnitConstants::cm);
  Acts::Vector3 globalPos(0,0,0);
  Acts::Vector3 mom(1,1,1);

  globalPos = surf->localToGlobal(m_tGeometry->geoContext,
				  localPos, mom);

  Acts::SymMatrix2 localCov = Acts::SymMatrix2::Zero();
  localCov(0,0) = clus->getActsLocalError(0,0) * Acts::UnitConstants::cm2;
  localCov(1,1) = clus->getActsLocalError(1,1) * Acts::UnitConstants::cm2;
  
  float x = globalPos.x();
  float y = globalPos.y();
  float z = globalPos.z();
  float r = std::sqrt(x * x + y * y);

  /// The space point requires only the variance of the transverse and
  /// longitudinal position. Reduce computations by transforming the
  /// covariance directly from local to r/z.
  ///
  /// compute Jacobian from global coordinates to r/z
  ///
  ///         r = sqrt(x² + y²)
  /// dr/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  ///             = 2 * {x,y} / r
  ///       dz/dz = 1 

  Acts::RotationMatrix3 rotLocalToGlobal =
    surf->referenceFrame(m_tGeometry->geoContext, globalPos, mom);
  auto scale = 2 / std::hypot(x,y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
    jacXyzToRhoZ *
    rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  SpacePointPtr spPtr(new SpacePoint{clus->getClusKey(), x, y, z, r, 
	surf->geometryId(), var[0], var[1]});

  if(Verbosity() > 2)
    std::cout << "Space point has " 
	      << x << ", " << y << ", " << z << " with local coords "
	      << localPos.transpose() 
	      << " with rphi/z variances " << localCov(0,0) 
	      << ", " << localCov(1,1) << " and rotated variances "
	      << var[0] << ", " << var[1] 
	      << " and cluster key "
	      << clus->getClusKey() << " and geo id "
	      << surf->geometryId() << std::endl;
  
  return spPtr;

}

std::vector<const SpacePoint*> PHActsSiliconSeeding::getMvtxSpacePoints(Acts::Extent& rRangeSPExtent)
{
  std::vector<const SpacePoint*> spVec;
  unsigned int numSiliconHits = 0;
 
  auto hitsetrange = m_hitsets->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (auto hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    {
      auto range = m_clusterMap->getClusters(hitsetitr->first);
      for( auto clusIter = range.first; clusIter != range.second; ++clusIter )
	{
	  const auto cluskey = clusIter->first;
	  const auto cluster = clusIter->second;
	  if(_iteration_map != NULL && _n_iteration >0 ){
	    //	    std::cout << "map exists entries: " << _iteration_map->size() << std::endl;
	    if(_iteration_map->getIteration(cluskey)>0){
	      continue; // skip hits used in a previous iteration
	    }
	  }

	  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
	  const auto surface = getSurface(hitsetkey);
	  if(!surface)
	    continue;

	  auto sp = makeSpacePoint(surface, cluster).release();
	  spVec.push_back(sp);
	  rRangeSPExtent.check({sp->x(), sp->y(), sp->z()});
	  numSiliconHits++;
	}
    }

  if(m_seedAnalysis)
    { h_nInputMvtxMeas->Fill(numSiliconHits); }

  if(Verbosity() > 1)
    std::cout << "Total number of silicon hits to seed find with is "
	      << numSiliconHits << std::endl;


  return spVec;
}
Surface PHActsSiliconSeeding::getSurface(TrkrDefs::hitsetkey hitsetkey)
{
  /// Only seed with the MVTX, so there is a 1-1 mapping between hitsetkey
  /// and acts surface
  auto surfMap = m_surfMaps->siliconSurfaceMap;
  auto iter = surfMap.find(hitsetkey);
  if(iter != surfMap.end())
    {
      return iter->second;
    }
  
  /// If it can't be found, return nullptr
  return nullptr;

}


Acts::SpacePointGridConfig PHActsSiliconSeeding::configureSPGrid()
{
  Acts::SpacePointGridConfig config;

  config.bFieldInZ = m_bField;
  config.minPt = m_minSeedPt / m_gridFactor;
  config.rMax = m_rMax;
  config.zMax = m_zMax;
  config.zMin = m_zMin;
  config.deltaRMax = m_deltaRMax;
  config.cotThetaMax = m_cotThetaMax;
  config.impactMax = m_impactMax;
  config.numPhiNeighbors = m_numPhiNeighbors;


  return config;
}

Acts::SeedFilterConfig PHActsSiliconSeeding::configureSeedFilter()
{
  Acts::SeedFilterConfig config;
  config.maxSeedsPerSpM = m_maxSeedsPerSpM;
  return config;
}

Acts::SeedfinderConfig<SpacePoint> PHActsSiliconSeeding::configureSeeder()
{
  Acts::SeedfinderConfig<SpacePoint> config;
  
  /// Limiting location of measurements (e.g. detector constraints)
  config.rMax = m_rMax;
  config.rMin = m_rMin;
  config.zMin = m_zMin;
  config.zMax = m_zMax;

  /// Min/max distance between two measurements in one seed
  config.deltaRMin = m_deltaRMin;
  config.deltaRMax = m_deltaRMax;

  /// Limiting collision region in z
  config.collisionRegionMin = -300. * Acts::UnitConstants::mm;
  config.collisionRegionMax = 300. * Acts::UnitConstants::mm;
  config.sigmaScattering = 5.;
  config.maxSeedsPerSpM = m_maxSeedsPerSpM;
  config.cotThetaMax = m_cotThetaMax;
  config.minPt = m_minSeedPt;
  config.bFieldInZ = m_bField;

  /// Average radiation length traversed per seed
  config.radLengthPerSeed = 0.05;

  /// Maximum impact parameter must be smaller than rMin
  config.impactMax = m_impactMax;

  return config;
}

int PHActsSiliconSeeding::getNodes(PHCompositeNode *topNode)
{
  m_surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  if(!m_surfMaps)
    {
      std::cout << PHWHERE << "Acts surface maps not on node tree, can't continue."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_geomContainerIntt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if(!m_geomContainerIntt)
    {
      std::cout << PHWHERE << "CYLINDERGEOM_INTT node not found on node tree"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "No ActsTrackingGeometry on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  if(m_useTruthClusters)
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
							    "TRKR_CLUSTER_TRUTH");
  else
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
							     "TRKR_CLUSTER");

  if(!m_clusterMap)
    {
      std::cout << PHWHERE << "No cluster container on the node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  m_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!m_hitsets)
    {
      std::cout << PHWHERE << "No hitset container on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::createNodes(PHCompositeNode *topNode)
{
 
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsSiliconSeeding::createNodes");
  }
  
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode,_track_map_name);
  if(!m_trackMap)
    {
      m_trackMap = new SvtxTrackMap_v1;
      PHIODataNode<PHObject> *trackNode = 
	new PHIODataNode<PHObject>(m_trackMap,_track_map_name,"PHObject");
      svtxNode->addNode(trackNode);

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsSiliconSeeding::writeHistograms()
{
  m_file->cd();
  h_nInttProj->Write();
  h_nMatchedClusters->Write();
  h_nMvtxHits->Write();
  h_nSeeds->Write();
  h_nActsSeeds->Write();
  h_nTotSeeds->Write();
  h_nInttHits->Write();
  h_nInputMeas->Write();
  h_nHits->Write();
  h_nInputMvtxMeas->Write();
  h_nInputInttMeas->Write();
  h_hits->Write();
  h_zhits->Write();
  h_projHits->Write();
  h_zprojHits->Write();
  h_resids->Write();
  m_file->Write();      
  m_file->Close(); 
}

void PHActsSiliconSeeding::createHistograms()
{
  h_nMatchedClusters = new TH1F("nMatchedClusters",";N_{matches}", 50,0,50);
  h_nInttProj = new TH2F("nInttProj",";l_{0}^{proj}-l_{0}^{clus} [cm]; l_{1}^{proj}-l_{1}^{clus} [cm]",
			 10000,-10,10,10000,-50,50);
  h_nMvtxHits = new TH1I("nMvtxHits",";N_{MVTX}",6,0,6);
  h_nInttHits = new TH1I("nInttHits",";N_{INTT}",80,0,80);
  h_nHits = new TH2I("nHits",";N_{MVTX};N_{INTT}",10,0,10,
		     80,0,80);
  h_nActsSeeds = new TH1I("nActsSeeds",";N_{Seeds}",400,0,400);
  h_nSeeds = new TH1I("nActsGoodSeeds",";N_{Seeds}",400,0,400);
  h_nTotSeeds = new TH1I("nTotSeeds",";N_{Seeds}",500,0,500);
  h_nInputMeas = new TH1I("nInputMeas",";N_{Meas}",2000,0,2000);
  h_nInputMvtxMeas = new TH1I("nInputMvtxMeas",";N_{meas}^{mvtx}",
			      150,0,150);
  h_nInputInttMeas = new TH1I("nInputInttMeas",";N_{meas}^{intt}",
			      150,0,150);
  h_hits = new TH2F("hits",";x [cm]; y [cm]",1000,-20,20,
		    1000,-20,20);
  h_zhits = new TH2F("zhits",";z [cm]; r [cm]",1000,-30,30,
		     1000,-30,30);
  h_projHits = new TH2F("projhits",";x [cm]; y [cm]",1000,-20,20,
			1000,-20,20);
  h_zprojHits = new TH2F("zprojhits",";z [cm]; r [cm]",1000,-30,30,
			 1000,-30,30);
  h_resids = new TH2F("resids",";z_{resid} [cm]; rphi_{resid} [cm]",
		      100,-1,1,100,-1,1);
}

double PHActsSiliconSeeding::normPhi2Pi(const double phi)
{
  double returnPhi = phi;
  if(returnPhi < 0)
    returnPhi += 2 * M_PI;
  return returnPhi;
}


void PHActsSiliconSeeding::largeGridSpacing(const bool spacing)
{
  if(!spacing)
    {
      m_gridFactor = 1.;
      m_rMax = 50.;
      m_cotThetaMax = 1.335647;
      m_maxSeedPCA = 0.1;
    }

}
