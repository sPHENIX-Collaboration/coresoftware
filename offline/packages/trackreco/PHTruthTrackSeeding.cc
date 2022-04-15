#include "PHTruthTrackSeeding.h"

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTra...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, Svtx...
#include <trackbase_historic/SvtxTrack_FastSim_v3.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>


#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/ActsTransformations.h>

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

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl
#include <map>       // for multimap, map<>::c...
#include <memory>
#include <utility>  // for pair
#include <cassert>
#include <set>

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

class PHCompositeNode;

using namespace std;

PHTruthTrackSeeding::PHTruthTrackSeeding(const std::string& name)
  : PHTrackSeeding(name)
  , _clustereval(nullptr)
{}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode)
{
  cout << "Enter PHTruthTrackSeeding:: Setup" << endl;

  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  _clustereval = new  SvtxClusterEval(topNode);
  _clustereval->do_caching(true);
  // _clustereval.set_strict(strict);
  //  _clustereval.set_verbosity(verbosity);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process(PHCompositeNode* topNode)
{
  _clustereval->next_event(topNode);
  _track_map->Clear();

  typedef std::map<int, std::set<TrkrCluster*> > TrkClustersMap;
  TrkClustersMap m_trackID_clusters;

  vector<TrkrDefs::cluskey> ClusterKeyList; 

  PHG4TruthInfoContainer::ConstRange range = _g4truth_container->GetPrimaryParticleRange();
  //  Float_t gntracks = (Float_t) truthinfo->GetNumPrimaryVertexParticles();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter){
    ClusterKeyList.clear();
    PHG4Particle* g4particle = iter->second;

    if (g4particle==NULL){
      cout <<__PRETTY_FUNCTION__<<" - validity check failed: missing truth particle" <<endl;
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
	  cout <<__PRETTY_FUNCTION__<<" ignore low momentum particle "<< gtrackID <<endl;
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

    // Smear the truth values out by 5% so that the seed momentum and
    // position aren't completely exact
    
    double random = ((double) rand() / (RAND_MAX)) * 0.05;
    // make it negative sometimes
    if(rand() % 2)
      random *= -1;
    svtx_track->set_px(g4particle->get_px() * (1 + random));
    svtx_track->set_py(g4particle->get_py() * (1 + random));
    svtx_track->set_pz(g4particle->get_pz() * (1 + random));

    // assign the track position using the truth vertex for this track
    auto g4vertex = _g4truth_container->GetVtx(vertexID);
    svtx_track->set_x(g4vertex->get_x() * (1 + random));
    svtx_track->set_y(g4vertex->get_y() * (1 + random));
    svtx_track->set_z(g4vertex->get_z() * (1 + random));
      
    for (TrkrDefs::cluskey cluskey : ClusterKeyList){
      svtx_track->insert_cluster_key(cluskey);
    }

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
    cout << "Loop over TrackMap " << _track_map_name << " entries " << endl;
    for (SvtxTrackMap::Iter iter = _track_map->begin();
         iter != _track_map->end(); ++iter)
    {
      SvtxTrack* svtx_track = iter->second;

      //svtx_track->identify();
      //continue;

      cout << "Track ID: " << svtx_track->get_id() << ", Dummy Track pT: "
           << svtx_track->get_pt() << ", Truth Track/Particle ID: "
           << svtx_track->get_truth_track_id() 
	   << " (X,Y,Z) " << svtx_track->get_x() << ", " << svtx_track->get_y() << ", " << svtx_track->get_z()  
	   << endl;
      cout << " nhits: " << svtx_track->size_cluster_keys()<< endl;
      //Print associated clusters;
      ActsTransformations transformer;
      for (SvtxTrack::ConstClusterKeyIter iter_clus =
               svtx_track->begin_cluster_keys();
           iter_clus != svtx_track->end_cluster_keys(); ++iter_clus)
      {
        TrkrDefs::cluskey cluster_key = *iter_clus;
	cout << "Key: "  << cluster_key<< endl;
        TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
	Acts::Vector3 global = transformer.getGlobalPosition(cluster,
							      surfmaps,
							      tgeometry);
        float radius = sqrt(global(0) * global(0) + global(1) * global(1));
        cout << "       cluster ID: "
             << cluster->getClusKey() << ", cluster radius: " << radius
             << endl;
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

  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
							  "TRKR_CLUSTER");
  
  if(!m_clusterMap)
    {
      cerr << PHWHERE << "Error: Can't find node TRKR_CLUSTER" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  _g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_g4truth_container)
  {
    cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    cout << PHWHERE << "Failed to find TRKR_HITTRUTHASSOC node, quit!" << endl;
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
{
  return 0;
}



void PHTruthTrackSeeding::circleFitSeed(std::vector<TrkrDefs::cluskey> clusters,
					double& x, double& y, double& z,
					double& px, double& py, double& pz, 
					int charge)
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
				  double& A, double& B)
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
void PHTruthTrackSeeding::findRoot(const double& R, const double& X0, const double& Y0, double& x, double& y)
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
  ActsTransformations transformer;
  std::vector<Acts::Vector3> globalPositions;
  for(auto& cluskey : clusters)
    {
      auto clus = m_clusterMap->findCluster(cluskey);
      auto glob = transformer.getGlobalPosition(clus,surfmaps,tgeometry);
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
