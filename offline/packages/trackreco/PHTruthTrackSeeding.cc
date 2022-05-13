#include "PHTruthTrackSeeding.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxEvaluator.h>
#include <g4eval/SvtxTrackEval.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHRandomSeed.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssoc.h> 
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTra...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, Svtx...
#include <trackbase_historic/SvtxTrack_FastSim_v3.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                            // for gsl_rng_alloc

#include <cassert>
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, std::endl 
#include <map>       // for multimap, map<>::c...
#include <memory>
#include <set>
#include <utility>  // for pair

class PHCompositeNode;

namespace
{ template< class T> inline constexpr T square( const T& x ) { return x*x; } }

PHTruthTrackSeeding::PHTruthTrackSeeding(const std::string& name)
  : PHTrackSeeding(name)
  , _clustereval(nullptr)
{
  // initialize random generator
  const uint seed = PHRandomSeed();
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );
}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode)
{
  std::cout << "Enter PHTruthTrackSeeding:: Setup" << std::endl ;

  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  _clustereval = new  SvtxClusterEval(topNode);
  _clustereval->do_caching(true);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process(PHCompositeNode* topNode)
{
  _clustereval->next_event(topNode);
  _track_map->Clear();

  using TrkClustersMap = std::map<int, std::set<TrkrCluster*> >;
  TrkClustersMap m_trackID_clusters;

  std::vector<TrkrDefs::cluskey> ClusterKeyList; 

  PHG4TruthInfoContainer::ConstRange range = m_g4truth_container->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter){
    ClusterKeyList.clear();
    PHG4Particle* g4particle = iter->second;

    if (g4particle==NULL){
      std::cout <<__PRETTY_FUNCTION__<<" - validity check failed: missing truth particle" << std::endl;
      exit(1);
    }

    const float gtrackID = g4particle->get_track_id();
    const int vertexID = g4particle->get_vtx_id();
    // monentum cut-off
    if (_min_momentum>0){
      const double monentum2 =
	g4particle->get_px() * g4particle->get_px()+
	g4particle->get_py() * g4particle->get_py()+
	g4particle->get_pz() * g4particle->get_pz();
      
      if (monentum2 < _min_momentum * _min_momentum){
	if (Verbosity() >= 3){
	  std::cout <<__PRETTY_FUNCTION__<<" ignore low momentum particle "<< gtrackID << std::endl;
	  g4particle->identify();
	}
	continue;
      }
    }

    for(unsigned int layer = _min_layer;layer < _max_layer;layer++){
      TrkrDefs::cluskey cluskey = _clustereval->best_cluster_by_nhit(gtrackID, layer);
      if(cluskey!=0) 
	ClusterKeyList.push_back(cluskey);
    }
    if(ClusterKeyList.size()< _min_clusters_per_track)
      continue;

    auto svtx_track = std::make_unique<SvtxTrack_FastSim_v3>();
    svtx_track->set_id(_track_map->size());
    svtx_track->set_truth_track_id(gtrackID);
    ///g4 vertex id starts at 1, svtx vertex map starts at 0
    // This is te truth particle vertex, and does not necessarily give the ID in the vertex map
    svtx_track->set_vertex_id(vertexID-1);
    
    // set track charge
    /*
     * having the proper track charge is necessary for the ACTS fit to work properly
     * with normal tracking, it is set by the track seeding. Here we get it from the G4Particle
     * unfortunately, particle charge is not stored into PHG4Particle.
     * we need to retrieve it from TParticlePDG instead
     */
    {
      const auto particle = TDatabasePDG::Instance()->GetParticle(g4particle->get_pid());
      if( particle ) svtx_track->set_charge(particle->Charge());
    }

    /* 
     * Smear the truth values out by +/-5% 
     * so that the seed momentum and 
     * position aren't completely exact
     */
    const auto random = gsl_ran_flat(m_rng.get(), 0.95, 1.05);
    svtx_track->set_px(g4particle->get_px() * random);
    svtx_track->set_py(g4particle->get_py() * random);
    svtx_track->set_pz(g4particle->get_pz() * random);

    // assign the track position using the truth vertex for this track
    auto g4vertex = m_g4truth_container->GetVtx(vertexID);
    svtx_track->set_x(g4vertex->get_x() * random);
    svtx_track->set_y(g4vertex->get_y() * random);
    svtx_track->set_z(g4vertex->get_z() * random);
    
    // associate all cluster keys to the track
    for( const auto& cluskey : ClusterKeyList){
      svtx_track->insert_cluster_key(cluskey);
    }

    // set intt crossing
    if(_min_layer < 7)
      {
	// silicon tracklet
	/* inspired from PHtruthSiliconAssociation */
	const auto intt_crossings = getInttCrossings(svtx_track.get());
	if(intt_crossings.empty()) 
	  {
	    if(Verbosity() > 1)  std::cout << "PHTruthTrackSeeding::Process - Silicon track " << svtx_track->get_id() << " has no INTT clusters" << std::endl;
	    continue ;
	  } else if( intt_crossings.size() > 1 ) {
	  
	  if(Verbosity() > 1) 
	    { std::cout << "PHTruthTrackSeeding::Process - INTT crossings not all the same for track " << svtx_track->get_id() << " crossing_keep - dropping this match " << std::endl; }
	  
	} else {
	  
	  const auto& crossing = *intt_crossings.begin();
	  svtx_track->set_crossing(crossing);
	  if(Verbosity() > 1)
	    std::cout << "PHTruthTrackSeeding::Process - Combined track " << svtx_track->get_id()  << " bunch crossing " << crossing << std::endl;           
	}
      }  // end if _min_layer
    else
      {
	// no INTT layers, crossing is unknown
	svtx_track->set_crossing(SHRT_MAX);	
      }

    // track fit
    if(m_helicalTrackFit)
    {
      double x, y, z, px, py, pz;
      circleFitSeed(ClusterKeyList, x, y, z,
        px, py, pz, svtx_track->get_charge());
      svtx_track->set_x(x);
      svtx_track->set_y(y);
      svtx_track->set_z(z);
      svtx_track->set_px(px);
      svtx_track->set_py(py);
      svtx_track->set_pz(pz);
    }
    
    svtx_track->set_ndf(ClusterKeyList.size()*3-5);
    svtx_track->set_chisq(1.5*ClusterKeyList.size()*3-5);
    
    _track_map->insert(svtx_track.get());
  }

  if (Verbosity() >= 5)
  {
    std::cout << "Loop over TrackMap " << _track_map_name << " entries " << std::endl ;
    for( const auto& [key,svtx_track]:*_track_map )
    {
    
      std::cout << "Track ID: " << svtx_track->get_id() << ", Dummy Track pT: "
           << svtx_track->get_pt() << ", Truth Track/Particle ID: "
           << svtx_track->get_truth_track_id() 
	   << " (X,Y,Z) " << svtx_track->get_x() << ", " << svtx_track->get_y() << ", " << svtx_track->get_z()  
	   << std::endl ;
      std::cout << " nhits: " << svtx_track->size_cluster_keys()<< std::endl ;
      //Print associated clusters;
      ActsTransformations transformer;
      for (SvtxTrack::ConstClusterKeyIter iter_clus =
               svtx_track->begin_cluster_keys();
           iter_clus != svtx_track->end_cluster_keys(); ++iter_clus)
      {
        TrkrDefs::cluskey cluster_key = *iter_clus;
        std::cout << "Key: "  << cluster_key<< std::endl ;
        TrkrCluster* cluster = m_clusterMap->findCluster(cluster_key);
        Acts::Vector3 global = transformer.getGlobalPosition(
          cluster_key, cluster,
          surfmaps,
          tgeometry);
        float radius = std::sqrt(square(global(0)) + square(global(1)));
        std::cout << "       cluster ID: "
             << cluster_key << ", cluster radius: " << radius
             << std::endl ;
      }
    }
  }

  //==================================

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::GetNodes(PHCompositeNode* topNode)
{

  tgeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  
  if((!tgeometry or !surfmaps) and m_helicalTrackFit) 
    {
      std::cerr << PHWHERE << "Error, can' find needed Acts nodes " << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  
  if(!m_clusterMap)
    {
      std::cerr << PHWHERE << "Error: Can't find node TRKR_CLUSTER" << std::endl ;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!m_cluster_crossing_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << std::endl ;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  m_g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_g4truth_container)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << std::endl ;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    std::cout << PHWHERE << "Failed to find TRKR_HITTRUTHASSOC node, quit!" << std::endl ;
    exit(1);
  }

  using nodePair = std::pair<std::string, PHG4HitContainer*&>;
  std::initializer_list<nodePair> nodes =
  {
    { "G4HIT_TPC", phg4hits_tpc },
    { "G4HIT_INTT", phg4hits_intt },
    { "G4HIT_MVTX", phg4hits_mvtx },
    { "G4HIT_MICROMEGAS", phg4hits_micromegas }
  };
  
  for( auto&& node: nodes )
  {
    if( !( node.second = findNode::getClass<PHG4HitContainer>( topNode, node.first ) ) )
    { std::cerr << PHWHERE << " PHG4HitContainer " << node.first << " not found" << std::endl; }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::End()
{ return 0; }

//_____________________________________________________________________________________________
std::set<short int> PHTruthTrackSeeding::getInttCrossings(SvtxTrack *si_track) const
{
  std::set<short int> intt_crossings;

  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  // loop over associated clusters to get keys for silicon cluster
  
  for( auto iter = si_track->begin_cluster_keys(); iter != si_track->end_cluster_keys(); ++iter)
  {
    
    const TrkrDefs::cluskey& cluster_key = *iter;
    const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
    if(trkrid == TrkrDefs::inttId)
    {
      
      // get layre from cluster key
      const unsigned int layer = TrkrDefs::getLayer(cluster_key);
      
      // get the bunch crossings for all hits in this cluster
      const auto crossings = m_cluster_crossing_map->getCrossings(cluster_key);
      for(auto iter = crossings.first; iter != crossings.second; ++iter)
      {
        const auto& [key, crossing] = *iter;
        if( Verbosity() )
        { std::cout << "PHTruthTrackSeeding::getInttCrossings - si Track " << si_track->get_id() << " cluster " << key << " layer " << layer << " crossing " << crossing  << std::endl; }
        intt_crossings.insert(crossing);
      }
    }
  }
  
  return intt_crossings;
}

void PHTruthTrackSeeding::circleFitSeed(std::vector<TrkrDefs::cluskey> clusters,
					double& x, double& y, double& z,
					double& px, double& py, double& pz, 
					int charge) const 
{
  double R, X0, Y0;
  auto globalPositions = circleFitByTaubin(clusters, R, X0, Y0);

  findRoot(R, X0, Y0, x, y);

  double phi = atan2(-1 * (X0-x), Y0-y);

  if(charge > 0)
    {
      phi += M_PI;
      if(phi > M_PI) 
	phi -= 2. * M_PI;
    }
 
  double m, B;
  lineFit(globalPositions, m, B);
  z = B;
  
  double theta = atan(1./m);
  if(theta < 0)
    theta =+ M_PI;
  
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
  
}

void PHTruthTrackSeeding::lineFit(std::vector<Acts::Vector3>& clusterPositions,
				  double& A, double& B) const 
{
  
  double xsum = 0,x2sum = 0,ysum = 0,xysum = 0;    
  for(const auto& pos : clusterPositions)
    {
      double z = pos(2);
      double r = sqrt(pow(pos(0),2) + pow(pos(1), 2));
      
      xsum=xsum+r;               // calculate sigma(xi)
      ysum=ysum+z;               // calculate sigma(yi)
      x2sum=x2sum+pow(r,2);      // calculate sigma(x^2i)
      xysum=xysum+r*z;           // calculate sigma(xi*yi)
    }
  
  /// calculate slope
  A = (clusterPositions.size()*xysum-xsum*ysum) / (clusterPositions.size()*x2sum-xsum*xsum);

  /// calculate intercept
  B = (x2sum*ysum-xsum*xysum) / (x2sum*clusterPositions.size()-xsum*xsum);
  
}
void PHTruthTrackSeeding::findRoot(const double& R, const double& X0, const double& Y0, double& x, double& y) const 
{
  
  double miny = (sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double miny2 = (-sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double minx = sqrt(pow(R, 2) - pow(miny - Y0, 2)) + X0;
  double minx2 = -sqrt(pow(R, 2) - pow(miny2 - Y0, 2)) + X0;
  
  /// Figure out which of the two roots is actually closer to the origin
  if(fabs(minx) < fabs(minx2))
    x = minx;
  else
    x = minx2;

  if(fabs(miny) < fabs(miny2))
    y = miny;
  else
    y = miny2;
  
}

std::vector<Acts::Vector3> PHTruthTrackSeeding::circleFitByTaubin(std::vector<TrkrDefs::cluskey>& clusters,
								   double& R, double& X0, double& Y0) const
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
  ActsTransformations transformer;
  std::vector<Acts::Vector3> globalPositions;
  for(auto& cluskey : clusters)
    {
      const auto clus = m_clusterMap->findCluster(cluskey);
      const auto glob = transformer.getGlobalPosition(cluskey, clus,surfmaps,tgeometry);
      globalPositions.push_back(glob);
      meanX += glob(0);
      meanY += glob(1);
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

  return globalPositions;
}
