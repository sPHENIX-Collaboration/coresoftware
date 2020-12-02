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

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>
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

  if(Verbosity() > 0)
    std::cout << "Processing PHActsSiliconSeeding event "
	      << m_event << std::endl;

  Acts::Seedfinder<SpacePoint> seedFinder(m_seedFinderCfg);
  
  /// Covariance converter tool needed by seed finder
  auto covConverter = [=](const SpacePoint& sp, float, float, float)
    -> Acts::Vector2D { return {sp.m_varianceRphi, sp.m_varianceZ};
  };

  auto spVec = getSpacePoints();

  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid = 
    Acts::SpacePointGridCreator::createGrid<SpacePoint>(m_gridCfg);

  auto spGroup = Acts::BinnedSPGroup<SpacePoint>( spVec.begin(),
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

  makeSvtxTracks(seedVector);

  if(Verbosity()> 0)
    std::cout << "Finished PHActsSiliconSeeding process_event"
	      << std::endl;

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
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
	  auto svtxTrack = std::make_unique<SvtxTrack_v1>();
	  svtxTrack->set_id(m_trackMap->size());
	  
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
		std::cout << "Adding cluster with radius " 
			  << sqrt(spacePoint->x() * spacePoint->x()
				  + spacePoint->y() * spacePoint->y())
			  << " mm " << std::endl;

	      svtxTrack->insert_cluster_key(cluskey);
	    }
	  
	  double x = NAN,y, z;;
	  double px, py, pz;
	  circleFitSeed(clusters, x, y, z,
			px, py, pz);

	  /// Bad seed, if x is nan so are y and z
	  if(std::isnan(x))
	    continue;

	  if(Verbosity() > 1)
	    std::cout <<"Setting silicon seed with (x,y,z) = " 
		      << x << ", " << y << ", " << seed.z() / 10.
		      << " and (px,py,pz) " << px << ", " << py
		      << ", " << pz << std::endl;
	  numGoodSeeds++;
	  
	  /// x and y were calculated in sPHENIX units
	  svtxTrack->set_x(x);
	  svtxTrack->set_y(y);
	  svtxTrack->set_z(seed.z() / Acts::UnitConstants::cm);
	  svtxTrack->set_px(px);
	  svtxTrack->set_py(py);
	  svtxTrack->set_pz(pz);

	  m_trackMap->insert(svtxTrack.release());
	}
    }

  if(Verbosity() > 1)
    {
      std::cout << "Total number of seeds found in " 
		<< seedVector.size() << " volume regions gives " 
		<< numSeeds << " seeds " << std::endl;
      std::cout << "Number of good seeds added to map : " << numGoodSeeds
		<< std::endl;
    }
  

}


void PHActsSiliconSeeding::circleFitSeed(const std::vector<TrkrCluster*> clusters,
					 double& x, double& y, double& z,
					 double& px, double& py, double& pz)
{
  /// Circle radius at x,y center
  /// Note - units are sPHENIX cm since we are using TrkrClusters
  double R, X0, Y0;
  circleFitByTaubin(clusters, R, X0, Y0);
  
  if(Verbosity() > 2)
    std::cout << "Circle R, X0, Y0 : " << R << ", " << X0
	      << ", " << Y0 << std::endl;

  /**
   * We need to determine the closest point on the circle to the origin
   * since we can't assume that the track originates from the origin
   * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to 
   * minimize d = sqrt((0-x)^2+(0-y)^2), the distance between the 
   * origin and some (currently, unknown) point on the circle x,y.
   * 
   * Solving the circle eqn for x and substituting into d gives an eqn for
   * y. Taking the derivative and setting equal to 0 gives the following 
   * two solutions. Depending on the charge of the track determines the 
   * correct solution
   */
  
  double miny = (sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double miny2 = (-sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double minx = sqrt(pow(R, 2) - pow(miny - Y0, 2)) + X0;
  double minx2 = sqrt(pow(R, 2) - pow(miny2 - Y0, 2)) + X0;
  
  if(Verbosity() > 1)
    std::cout << "minx1 and x2 : " << minx << ", " << minx2 << std::endl
	      << "miny1 and y2 : " << miny << ", " << miny2 << std::endl;

  if(fabs(minx) < fabs(minx2))
    x = minx;
  else
    x = minx2;

  /// determine which y solution is smaller
  if(fabs(miny) < fabs(miny2))
    y = miny;
  else
    y = miny2;
  
  /// If the x or y initial position was found to be greater than 10 cm
  /// it is a bad seed
  if(fabs(x) > 10. or fabs(y) > 10.)
    {
      x = NAN;
      return;
    }

  if(Verbosity() > 1)
    {
      if(!std::isnan(x) && !std::isnan(y))
	std::cout << "Minimum x and y positions " << x << ",  " 
		  << y << std::endl;
    }

  /// Now determine the line tangent to the circle at this point to get phi
  double phi = atan( -1./( (y - Y0) / (x - X0)));
  if(phi > M_PI) phi -= 2. * M_PI;
  if(phi < -M_PI) phi += 2. * M_PI;
  
  if(Verbosity() > 1)
    std::cout << "Track seed phi : " << phi << std::endl;

  double m, B;
  
  /// m is slope as a function of radius, B is z intercept (vertex)
  lineFit(clusters, m, B);

  z = B;

  double theta = atan(1./m);

  if(Verbosity() > 1)
    std::cout << "Track seed theta: " << theta << std::endl;
 
  /// 0.035 is the unit conversion from cmT to GeV
  double p = R * 1.4 * 0.035;
  
  if(Verbosity() > 1)
    std::cout << "track momentum estimate is " << p << std::endl;

  px = p * sin(theta) * cos(phi);
  py = p * sin(theta) * sin(phi);
  pz = p * cos(theta);
  
  if(Verbosity() > 1)
    std::cout << "Momentum estimate: (" << px <<" , " << py 
	      << ", " << pz << ") " << std::endl;
  
  /// normalize to unit vector because we just need the direction
  px /= p;
  py /= p;
  pz /= p;
  
  return;

}
void PHActsSiliconSeeding::lineFit(std::vector<TrkrCluster*> clusters, 
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
  
  if(Verbosity() > 10)
    {
      for (auto& cluster : clusters)
	{
	  double r = sqrt(pow(cluster->getX(),2) + pow(cluster->getY(), 2));
	  /// To calculate y(fitted) at given x points
	  double z_fit = A * r + B;               
	  std::cout << " r " << r << " z " << cluster->getZ() 
		    << " z_fit " << z_fit << std::endl; 
	} 
    }
  
  return;
}   

void PHActsSiliconSeeding::circleFitByTaubin(const std::vector<TrkrCluster*> clusters,
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

std::vector<const SpacePoint*> PHActsSiliconSeeding::getSpacePoints()
{
  std::vector<const SpacePoint*> spVec;
  unsigned int numSiliconHits = 0;
  for(auto &[hitId, sl] : *m_sourceLinks)
    {
      /// collect only source links in silicon
      auto volume = sl.referenceSurface().geometryId().volume();
      
      /// If we run without MMs, volumes are 7, 9, 11 for mvtx, intt, tpc
      /// If we run with MMs, volumes are 10, 12, 14, 16 for mvtx, intt, tpc, mm
      if(volume == 11 or volume > 12)
	continue;
      
      auto sp = makeSpacePoint(hitId, sl).release();
      spVec.push_back(sp);
      numSiliconHits++;
    }
  
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

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode,
						  "SvtxVertexMap");
  if(!m_vertexMap)
    {
      std::cout << PHWHERE << "No SvtxVertexMap on node tree. Bailing."
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
    throw std::runtime_error("Failed to find DST node in PHActsTracks::createNodes");
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
