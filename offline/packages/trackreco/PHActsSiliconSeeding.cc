#include "PHActsSiliconSeeding.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase/TrkrCluster.h>            
#include <trackbase/TrkrClusterContainer.h>

#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>

PHActsSiliconSeeding::PHActsSiliconSeeding(const std::string& name)
  : SubsysReco(name)
  , m_sourceLinks(nullptr)
  , m_hitIdCluskey(nullptr)
{}

int PHActsSiliconSeeding::Init(PHCompositeNode *topNode)
{
  m_seedFinderCfg = configureSeeder();
  m_gridCfg = configureSPGrid();
  
  m_bottomBinFinder = 
    std::make_shared<Acts::BinFinder<SpacePoint>>(Acts::BinFinder<SpacePoint>());
  m_topBinFinder = 
    std::make_shared<Acts::BinFinder<SpacePoint>>(Acts::BinFinder<SpacePoint>());

  Acts::SeedFilterConfig sfCfg;
  sfCfg.maxSeedsPerSpM = m_maxSeedsPerSpM;

  m_seedFinderCfg.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
     Acts::SeedFilter<SpacePoint>(sfCfg));

  if(m_seedAnalysis)
    {
      
      m_file = new TFile("seedingOutfile.root","recreate");
    }
  createHistograms();
    
  
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::InitRun(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::process_event(PHCompositeNode *topNode)
{

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
    delete sp;
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

int PHActsSiliconSeeding::End(PHCompositeNode *topNode)
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
  auto covConverter = [=](const SpacePoint& sp, float, float, float)
    -> Acts::Vector2D { return {sp.m_varianceRphi, sp.m_varianceZ};
  };

   spVec = getMvtxSpacePoints();

  h_nInputMeas->Fill(spVec.size());

  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid = 
    Acts::SpacePointGridCreator::createGrid<SpacePoint>(m_gridCfg);

  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(),
						 spVec.end(),
						 covConverter,
						 m_bottomBinFinder,
						 m_topBinFinder,
						 std::move(grid),
						 m_seedFinderCfg);

  /// This is a vector of seeds inside of a vector which represents 
  /// volume grids of the detector area. The seeds can be accessed
  /// by iterating over each area, and then collecting the seeds in 
  /// that area
  GridSeeds seedVector;
  auto groupIt = spGroup.begin();
  auto endGroup = spGroup.end();

  for(; !(groupIt == endGroup); ++groupIt)
    {
      seedVector.push_back(seedFinder.createSeedsForGroup(groupIt.bottom(),
							  groupIt.middle(),
							  groupIt.top()));
    }
  
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
	  if(Verbosity() > 1)
	    std::cout << "Seed " << numSeeds << " has "
		      << seed.sp().size() << " measurements " 
		      << std::endl;

	  numSeeds++;

	  std::vector<TrkrCluster*> clusters;
	  for(auto& spacePoint : seed.sp())
	    {
	      auto cluskey = m_hitIdCluskey->right.find(spacePoint->m_hitId)->second;
	      clusters.push_back(m_clusterMap->findCluster(cluskey));

	      if(Verbosity() > 1)
		std::cout << "Adding cluster with x,y "
			  << spacePoint->x() <<", " << spacePoint->y()
			  << " mm in detector " 
			  << TrkrDefs::getTrkrId(cluskey)
			  << std::endl;
	    }
	 
	  double x = NAN, y = NAN, z = seed.z() / Acts::UnitConstants::cm;
	  double px, py, pz;
	  
	  /// Performs circle fit and extrapolates to INTT layers to
	  /// to get additional clusters in this track seed
	  int charge = circleFitSeed(clusters, x, y, z,
				     px, py, pz);

	  /// Bad seed, if x is nan so are y and z
	  if(std::isnan(x))
	    {
	      m_nBadInitialFits++;
	      continue;
	    }

	  numGoodSeeds++;
	  
	  createSvtxTrack(x, y, seed.z() / Acts::UnitConstants::cm,
			  px, py, pz, charge,
			  clusters);
	}
    }

  h_nSeeds->Fill(numGoodSeeds);
  h_nActsSeeds->Fill(numSeeds);

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
					   const std::vector<TrkrCluster*> clusters)
{

  auto stubs = makePossibleStubs(clusters);
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
  for(auto [stub, stubClusters] : stubs)
    {
      int nMvtx = 0;
      int nIntt = 0;
      numSeedsPerActsSeed++;
      
      #if __cplusplus < 201402L
      auto svtxTrack = boost::make_unique<SvtxTrack_v1>();
      #else
      auto svtxTrack = std::make_unique<SvtxTrack_v1>();
      #endif

      svtxTrack->set_id(m_trackMap->size());
      
      for(const auto clus : stubClusters) 
	{
	  const auto cluskey = clus->getClusKey();
	  svtxTrack->insert_cluster_key(cluskey);
	  if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::mvtxId)
	    nMvtx++;
	  else if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::inttId)
	    nIntt++;	 
	}
 
      /// Get a less rough estimate of R, and thus, p
      double R, X0, Y0;
      circleFitByTaubin(stubClusters, R, X0, Y0);
      
      /// 0.3 conversion factor, 1.4=B field, 
      /// 100 convert R from cm to m
      float pt = 0.3 * 1.4 * R / 100.;
  
      trackPx = pt * cos(trackPhi);
      trackPy = pt * sin(trackPhi);
      trackPz = pt * sinh(trackEta);
      
      /// Diagnostic
      h_nInttHits->Fill(nIntt);
      h_nMvtxHits->Fill(nMvtx);
      h_nHits->Fill(nMvtx, nIntt);
      
      if(Verbosity() > 1)
	std::cout << "Setting silicon seed with track id " 
		  << m_trackMap->size()
		  << " and (x,y,z) = " 
		  << trackX << ", " << trackY << ", " << trackZ
		  << std::endl << " and (px,py,pz) " << trackPx 
		  << ", " << trackPy << ", " << trackPz << std::endl
		  << " with charge " << trackCharge << std::endl;
      
      svtxTrack->set_x(trackX);
      svtxTrack->set_y(trackY);
      svtxTrack->set_z(trackZ);
      svtxTrack->set_px(trackPx);
      svtxTrack->set_py(trackPy);
      svtxTrack->set_pz(trackPz);
      svtxTrack->set_charge(trackCharge);
      
      m_trackMap->insert(svtxTrack.release());
  
    }

  h_nTotSeeds->Fill(numSeedsPerActsSeed);
  if(Verbosity() > 1)
    std::cout << "Found " << numSeedsPerActsSeed << " seeds for one Acts seed"
	      << std::endl;

}

std::map<const unsigned int, std::vector<TrkrCluster*>> 
     PHActsSiliconSeeding::makePossibleStubs(std::vector<TrkrCluster*> 
					     allClusters)
{

  std::vector<TrkrCluster*> mvtxClusters;
  std::vector<TrkrCluster*> inttFirstLayerClusters;
  std::vector<TrkrCluster*> inttSecondLayerClusters;
  std::map<const unsigned int, std::vector<TrkrCluster*>> stubs;
  unsigned int combo = 0;

  for(const auto clus : allClusters)
    {
      const auto cluskey = clus->getClusKey();
    
      if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::mvtxId)
	mvtxClusters.push_back(clus);
      else
	{
	  const double r = sqrt(pow(clus->getX(), 2) +
				pow(clus->getY(), 2));
	  if(r < 8.) 
	    inttFirstLayerClusters.push_back(clus);
	  else
	    inttSecondLayerClusters.push_back(clus);
	}
    }
  
  /// if we have no INTT matches, we can just return one track stub
  /// with the mvtx hits
  if(inttFirstLayerClusters.size() == 0 and 
     inttSecondLayerClusters.size() == 0)
    {
      stubs.insert(std::make_pair(combo, mvtxClusters));
    }

  /// If we have only one matched INTT hit, 
  /// make a single stub and return
  else if(inttFirstLayerClusters.size() == 1 and 
     inttSecondLayerClusters.size() == 0)
    {
      std::vector<TrkrCluster*> dumVec = mvtxClusters;
      dumVec.push_back(inttFirstLayerClusters.at(0));
      stubs.insert(std::make_pair(combo, dumVec));
    }
  else if(inttFirstLayerClusters.size() == 0 and 
     inttSecondLayerClusters.size() == 1)
    {
      std::vector<TrkrCluster*> dumVec = mvtxClusters;
      dumVec.push_back(inttSecondLayerClusters.at(0));
      stubs.insert(std::make_pair(combo, dumVec));
    }
  else
    {
      /// Otherwise we have 2 or more INTT matched hits
      /// Find combinations of INTT clusters that were within
      /// the nominal matching windows if each INTT layer group
      /// had at least one cluster
      
      for(const auto firstLayerClus : inttFirstLayerClusters)
	{
	  for(const auto secondLayerClus : inttSecondLayerClusters)
	    {
	      /// Start with the mvtx clusters
	      std::vector<TrkrCluster*> dumVec = mvtxClusters;
	      dumVec.push_back(firstLayerClus);
	      dumVec.push_back(secondLayerClus);
	      stubs.insert(std::make_pair(combo, dumVec));
	      combo++;
	    }
	}
    }

  return stubs;
}

int PHActsSiliconSeeding::circleFitSeed(std::vector<TrkrCluster*>& clusters,
					double& x, double& y, double& z,
					double& px, double& py, double& pz)
{
  if(Verbosity() > 2)
    for(const auto clus : clusters)
      std::cout << "Evaluating cluster : " << clus->getClusKey()
		<< std::endl;

  /// Circle radius at x,y center
  /// Note - units are sPHENIX cm since we are using TrkrClusters
  double R, X0, Y0;
  circleFitByTaubin(clusters, R, X0, Y0);
  
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

  int charge = getCharge(clusters, atan2(Y0,X0));
  
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
  
  lineFit(clusters, m, B);
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
    std::cout << "Momentum vector estimate: (" << px <<" , " 
	      << py << ", " << pz << ") " << std::endl;
    
  /// Project to INTT and find matches
  auto additionalClusters = findInttMatches(clusters, R, X0, Y0, z, m);
  
  /// Add possible matches to cluster list to be parsed when
  /// Svtx tracks are made
  for(auto cluskey : additionalClusters)
    clusters.push_back(m_clusterMap->findCluster(cluskey));
  
  return charge;

}

std::vector<TrkrDefs::cluskey> PHActsSiliconSeeding::findInttMatches(
				      const std::vector<TrkrCluster*>& clusters,
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
  for(auto clus : clusters)
    {
      h_hits->Fill(clus->getX(), clus->getY());
      h_zhits->Fill(clus->getZ(),
		    sqrt(pow(clus->getX(),2) + pow(clus->getY(),2)));
      h_projHits->Fill(clus->getX(), clus->getY());
      h_zprojHits->Fill(clus->getZ(),
			sqrt(pow(clus->getX(),2) + pow(clus->getY(),2)));
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
      const unsigned int lastClus = clusters.size() - 1;
      const double lastClusPhi = atan2(clusters.at(lastClus)->getY(),
				       clusters.at(lastClus)->getX());
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

      h_projHits->Fill(xProj[layer], yProj[layer]);
      h_zprojHits->Fill(zProj[layer], sqrt(pow(xProj[layer],2) + 
					   pow(yProj[layer],2)));
   
      if(Verbosity() > 2)
	{
	  std::cout << "Projected point is : " << xProj[layer] << ", "
		    << yProj[layer] << ", " << zProj[layer] << std::endl;
	}
    }

  auto additionalClusters = matchInttClusters(xProj, yProj, zProj);

  return additionalClusters;
}

std::vector<TrkrDefs::cluskey> PHActsSiliconSeeding::matchInttClusters(
						     const double xProj[],
						     const double yProj[],
						     const double zProj[])
{
  std::vector<TrkrDefs::cluskey> matchedClusters;
  
  TrkrClusterContainer::ConstRange inttClusRange = 
    m_clusterMap->getClusters(TrkrDefs::inttId);
  
  for(TrkrClusterContainer::ConstIterator clusIter = inttClusRange.first;
      clusIter != inttClusRange.second; ++clusIter)
    {
      const auto cluskey = clusIter->first;
      const auto cluster = clusIter->second;
      
      /// Subtract three to subtract off the mvtx layers for comparison
      /// to projections
      const auto projLayer = TrkrDefs::getLayer(cluskey) - 3;

      const auto sphenixLayer = TrkrDefs::getLayer(cluskey);
  
      auto layerGeom = dynamic_cast<CylinderGeomIntt*>
	(m_geomContainerIntt->GetLayerGeom(sphenixLayer));

      const auto stripZSpacing = layerGeom->get_strip_z_spacing();

      const double inttClusZ = cluster->getZ();
      const double inttClusR = sqrt(pow(cluster->getX(), 2) + 
				    pow(cluster->getY(), 2) );
      const double inttClusRphi = inttClusR * atan2(cluster->getY(),
						    cluster->getX());
      const double projR = sqrt(pow(xProj[projLayer], 2) + 
				pow(yProj[projLayer], 2));
      const double projRphi = projR * atan2(yProj[projLayer], xProj[projLayer]);

      h_hits->Fill(cluster->getX(), cluster->getY());
      h_zhits->Fill(cluster->getZ(),
		    inttClusR);

      if(Verbosity() > 4)
	std::cout << "Checking INTT cluster with position " << cluster->getX()
		  << ", " << cluster->getY() << ", " << cluster->getZ()
		  << std::endl << " with projections rphi "
		  << projRphi << " and inttclus rphi " << inttClusRphi
		  << " and proj z " << zProj[projLayer] << " and inttclus z "
		  << inttClusZ << " in layer " << projLayer 
		  << " with search windows " << m_rPhiSearchWin 
		  << " in rphi and strip z spacing " << stripZSpacing 
		  << std::endl;

      h_resids->Fill(zProj[projLayer] - inttClusZ,
		     projRphi - inttClusRphi);

      /// Z strip spacing is the entire strip, so because we use fabs
      /// we divide by two
      if(fabs(projRphi - inttClusRphi) < m_rPhiSearchWin and
	 fabs(zProj[projLayer] - inttClusZ) < stripZSpacing / 2.)
	{
	  matchedClusters.push_back(cluskey);

	  if(Verbosity() > 2)
	    {
	      std::cout << "Found matching projection with cluskey " 
			<< cluskey << std::endl;
	    }

	}

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

int PHActsSiliconSeeding::getCharge(const std::vector<TrkrCluster*>& clusters,
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
  for(auto clus : clusters)
    {
      double clusPhi = atan2(clus->getY(), clus->getX());

      /// if it is close to the periodic boundary normalize to 
      /// two pi to avoid -pi and pi issues
      if(fabs(fabs(clusPhi) - M_PI) < 0.2)
	clusPhi = normPhi2Pi(clusPhi);
      trackPhi += clusPhi;
    }

  trackPhi /= clusters.size();

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

void PHActsSiliconSeeding::lineFit(const std::vector<TrkrCluster*>& clusters, 
				   double &A, double &B)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  
  double xsum = 0,x2sum = 0,ysum = 0,xysum = 0;    
  for(auto& cluster : clusters)
    {
      double z = cluster->getZ();
      double r = sqrt(pow(cluster->getX(),2) + pow(cluster->getY(), 2));
      
      xsum=xsum+r;               // calculate sigma(xi)
      ysum=ysum+z;               // calculate sigma(yi)
      x2sum=x2sum+pow(r,2);      // calculate sigma(x^2i)
      xysum=xysum+r*z;           // calculate sigma(xi*yi)
    }
  
  /// calculate slope
  A = (clusters.size()*xysum-xsum*ysum) / (clusters.size()*x2sum-xsum*xsum);

  /// calculate intercept
  B = (x2sum*ysum-xsum*xysum) / (x2sum*clusters.size()-xsum*xsum);
  
  if(Verbosity() > 4)
    {
      for (auto& cluster : clusters)
	{
	  double r = sqrt(pow(cluster->getX(),2) + pow(cluster->getY(), 2));
	  /// To calculate y(fitted) at given x points
	  double z_fit = A * r + B;               
	  std::cout << " r " << r << " z " << cluster->getZ() 
		    << " z_fit " << z_fit << std::endl; 
	} 
      for(int i =0; i <m_nInttLayers; i++)
	{
	  std::cout << "intt z_fit layer " << i << " is " 
		    << A * m_nInttLayerRadii[i] + B << std::endl;
	}
    }
  
  return;
}   

void PHActsSiliconSeeding::circleFitByTaubin(const std::vector<TrkrCluster*>& clusters,
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
  
  for(auto clus : clusters)
    {
      meanX += clus->getX();
      meanY += clus->getY();
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

  for(auto clus : clusters)
    {
      double Xi = clus->getX() - meanX;
      double Yi = clus->getY() - meanY;
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

SpacePointPtr PHActsSiliconSeeding::makeSpacePoint(const unsigned int& hitId,
						   const SourceLink& sl)
{
  Acts::Vector2D localPos(sl.location()(0), sl.location()(1));
  Acts::Vector3D globalPos(0,0,0);
  Acts::Vector3D mom(1,1,1);

  globalPos = sl.referenceSurface().localToGlobal(m_tGeometry->geoContext,
						  localPos, mom);

  auto cov = sl.covariance();
  float x = globalPos.x();
  float y = globalPos.y();
  float z = globalPos.z();
  float r = std::sqrt(x * x + y * y);
  float varianceRphi = cov(0,0);
  float varianceZ = cov(1,1);
  
  SpacePointPtr spPtr(new SpacePoint{sl.hitID(), x, y, z, r, 
	sl.referenceSurface().geometryId(), varianceRphi, varianceZ});

  if(Verbosity() > 2)
    std::cout << "Space point has " 
	      << x << ", " << y << ", " << z
	      << " with variances " << varianceRphi 
	      << ", " << varianceZ 
	      << " and hit id "
	      << sl.hitID() << " and geo id "
	      << sl.referenceSurface().geometryId() << std::endl;
  
  return spPtr;

}

std::vector<const SpacePoint*> PHActsSiliconSeeding::getMvtxSpacePoints()
{
  std::vector<const SpacePoint*> spVec;
  unsigned int numSiliconHits = 0;

  for(const auto &[hitId, sl] : *m_sourceLinks)
    {
      /// collect only source links in MVTX
      auto volume = sl.referenceSurface().geometryId().volume();

      /// If we run without MMs, volumes are 7, 9, 11 for mvtx, intt, tpc
      /// If we run with MMs, volumes are 10, 12, 14, 16 for mvtx, intt, tpc, mm
      if(volume == 7 or volume == 10)
	{
     	  auto sp = makeSpacePoint(hitId, sl).release();
	  spVec.push_back(sp);
	  numSiliconHits++;
	}
    }
  
  h_nInputMvtxMeas->Fill(numSiliconHits);

  if(Verbosity() > 1)
    std::cout << "Total number of silicon hits to seed find with is "
	      << numSiliconHits << std::endl;

  return spVec;
}

Acts::SpacePointGridConfig PHActsSiliconSeeding::configureSPGrid()
{
  Acts::SpacePointGridConfig config;

  config.bFieldInZ = m_bField;
  config.minPt = m_minSeedPt;
  config.rMax = m_rMax;
  config.zMax = m_zMax;
  config.zMin = m_zMin;
  config.deltaRMax = m_deltaRMax;
  config.cotThetaMax = m_cotThetaMax;

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
  config.deltaRMin = 1.;
  config.deltaRMax = m_deltaRMax;

  /// Limiting collision region in z
  config.collisionRegionMin = -100.;
  config.collisionRegionMax = 100.;
  config.sigmaScattering = 50.;
  config.maxSeedsPerSpM = m_maxSeedsPerSpM;
  config.cotThetaMax = m_cotThetaMax;
  config.minPt = m_minSeedPt;
  config.bFieldInZ = m_bField;

  /// Average radiation length traversed per seed
  config.radLengthPerSeed = 0.05;

  /// Maximum impact parameter must be smaller than rMin
  config.impactMax = 20;

  return config;
}

int PHActsSiliconSeeding::getNodes(PHCompositeNode *topNode)
{
  m_geomContainerIntt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if(!m_geomContainerIntt)
    {
      std::cout << PHWHERE << "CYLINDERGEOM_INTT node not found on node tree"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");
  if(!m_sourceLinks)
    {
      std::cout << PHWHERE << "TrkrClusterSourceLinks node not on node tree. Bailing."
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

  m_hitIdCluskey = findNode::getClass<CluskeyBimap>(topNode,
						    "HitIDClusIDActsMap");
  if(!m_hitIdCluskey)
    {
      std::cout << PHWHERE << "No hit id clus id source link map on node tree. Bailing."
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

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxSiliconTrackMap");
  if(!m_trackMap)
    {
      m_trackMap = new SvtxTrackMap_v1;
      PHIODataNode<PHObject> *trackNode = 
	new PHIODataNode<PHObject>(m_trackMap,"SvtxSiliconTrackMap","PHObject");
      dstNode->addNode(trackNode);

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsSiliconSeeding::writeHistograms()
{
  m_file->cd();
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
