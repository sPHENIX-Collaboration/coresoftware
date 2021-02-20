#include "PHTpcTrackSeedVertexAssoc.h"

#include "AssocInfoContainer.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <TF1.h>

#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>
using namespace std;

//____________________________________________________________________________..
PHTpcTrackSeedVertexAssoc::PHTpcTrackSeedVertexAssoc(const std::string &name):
 PHTrackPropagating(name)
 , _track_map_name_silicon("SvtxSiliconTrackMap")
{
  //cout << "PHTpcTrackSeedVertexAssoc::PHTpcTrackSeedVertexAssoc(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHTpcTrackSeedVertexAssoc::~PHTpcTrackSeedVertexAssoc()
{

}

//____________________________________________________________________________..
int PHTpcTrackSeedVertexAssoc::Setup(PHCompositeNode *topNode)
{
  // put these in the output file
  cout << PHWHERE << " p0 " << _par0 << " p1 " << _par1 << " p2 " 
       << _par2 << " Search windows: phi " << _phi_search_win << " eta " 
       << _eta_search_win << endl;

  // corrects the PHTpcTracker phi bias
  fdphi = new TF1("f1", "[0] + [1]/x^[2]");
  fdphi->SetParameter(0, _par0);
  fdphi->SetParameter(1, _par1);
  fdphi->SetParameter(2, _par2);

  // corrects the space charge distortion phi bias
  if(!_is_ca_seeder)
    {
      // PHTpcTracker correction is opposite in sign
      // and different in magnitude - why?
      _parsc0 *= -1.0 * 0.7;
      _parsc1 *= -1.0 * 0.7;
    }
  fscdphi = new TF1("f2","[0] + [1]*x^2");
  fscdphi->SetParameter(0, _parsc0 * _collision_rate / _reference_collision_rate);
  fscdphi->SetParameter(1, _parsc1 * _collision_rate / _reference_collision_rate);

  int ret = PHTrackPropagating::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHTpcTrackSeedVertexAssoc::Process()
{
  // _track_map contains the TPC seed track stubs
  // We want to associate these TPC track seeds with a collision vertex
  // Then we add the collision vertex position as the track seed position

  // All we need is to project the TPC clusters in Z to the beam line.


  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size()  << endl;
  /*
 // We remember the original size of the TPC track map here
  const unsigned int original_track_map_lastkey = _track_map->end()->first;
  */

  // loop over the TPC track seeds
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      /*
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first >= original_track_map_lastkey)  break;
      */      
      _tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet_tpc->get_phi()
	    << endl;
	}

      // get the tpc track seed cluster positions in z and r

      // Get the outermost TPC clusters for this tracklet
      std::map<unsigned int, TrkrCluster*> tpc_clusters_map;
      std::vector<TrkrCluster*> clusters;
      
      for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_tpc->begin_cluster_keys();
	   key_iter != _tracklet_tpc->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  if(layer < _min_tpc_layer) continue;
	  if(layer >= _max_tpc_layer) continue;

	  // get the cluster
	  TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);

	  tpc_clusters_map.insert(std::make_pair(layer, tpc_clus));
	  clusters.push_back(tpc_clus);

	  if(Verbosity() > 10) 
	    std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		      << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " outer_clusters.size() " << tpc_clusters_map.size() << std::endl;
	}


      // need at least 3 clusters to fit a circle
      if(tpc_clusters_map.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough outer clusters " << std::endl; 
	  continue;  // skip to the next TPC tracklet
	}

      /*
      // fit a circle to the clusters
      double R, X0, Y0;
      CircleFitByTaubin(clusters, R, X0, Y0);
      if(Verbosity() > 10) std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // toss tracks for which the fitted circle could not have come from the vertex
      if(R < 40.0) continue;
      */

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit_clusters(clusters, A, B);
      if(Verbosity() > 10) std::cout << " Fitted line has A " << A << " B " << B << std::endl;

      // Project this TPC tracklet  to the beam line and store the projections
      //bool skip_tracklet = false;
      // z projection is unique
      _z_proj = B;
      
      // Find the nearest collision vertex
      // if find multiple close matches, duplicate the track?

      int trackVertexId = 9999;
      double dz = 9999.;	  
      for(SvtxVertexMap::Iter viter = _vertex_map->begin();
	  viter != _vertex_map->end();
	  ++viter)
	{
	  /*
	  const double trackX = track->get_x();
	  const double trackY = track->get_y();
	  const double trackZ = track->get_z();
	  double dx = 9999.;
	  double dy = 9999.;
	  */

	  auto vertexKey = viter->first;
	  auto vertex = viter->second;
	  if(Verbosity() > 100)
	    vertex->identify();

	  /*	  
	  const double vertexX = vertex->get_x();
	  const double vertexY = vertex->get_y();
	  */
	  const double vertexZ = vertex->get_z();
	  
	  if( 
	     //fabs(trackX - vertexX) < dx &&
	     //fabs(trackY - vertexY) < dy &&
	      fabs(_z_proj - vertexZ) < dz )
	    {
	      //dx = fabs(trackX - vertexX);
	      //dy = fabs(trackY - vertexY);
	      dz = fabs(_z_proj - vertexZ);
	      trackVertexId = vertexKey;
	    }	  
	}  // end loop over collision vertices

      _tracklet_tpc->set_vertex_id(trackVertexId);
      auto vertex = _vertex_map->find(trackVertexId)->second;
      vertex->insert_track(phtrk_iter->first);

      // set the track position to the vertex position
      _tracklet_tpc->set_x(vertex->get_x());
      _tracklet_tpc->set_y(vertex->get_y());
      _tracklet_tpc->set_z(vertex->get_z());

      if(Verbosity() > 1)
	{
	  std::cout << "TPC seed track " << phtrk_iter->first << " matched to vertex " << trackVertexId << endl; 
	} 

      // Repeat the z fit including the vertex position, get theta
      std::vector<std::pair<double, double>> points;

      double r_vertex = sqrt(vertex->getX()*vertex->getX() + vertex->getY()*vertex->getY());
      double z_vertex = vertex->getZ();
      points.push_back(make_pair(r_vertex, z_vertex));

      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  double z = clusters[i]->getZ();
	  double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
	  
	  points.push_back(make_pair(r,z));
	}
      
      line_fit(points, A, B);
      if(Verbosity() > 10) std::cout << " Fitted line including vertex has A " << A << " B " << B << std::endl;      

      // extract the track theta, update pz of track?


      // add circle fit including vertex as point


      // Update track pT magnitude from circle fit?


      // extract the track phi (tangent to circle at r_vertex, maybe at y = y_vertex), update px, py of track?


      // update track on node tree, done

      
    }  // end loop over TPC track seeds
  
  if(Verbosity() > 0)  
    cout << " Final track map size " << _track_map->size() << endl;
  
  if (Verbosity() >= 1)
    cout << "PHTpcTrackSeedVertexAssoc::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcTrackSeedVertexAssoc::End()
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHTpcTrackSeedVertexAssoc::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------
  /*
  _track_map_silicon = findNode::getClass<SvtxTrackMap>(topNode,  "SvtxSiliconTrackMap");
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxSiliconTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */  
  return Fun4AllReturnCodes::EVENT_OK;
}


void  PHTpcTrackSeedVertexAssoc::line_fit(std::vector<std::pair<double,double>> points, double &a, double &b)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  
   double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
   for (unsigned int i=0; i<points.size(); ++i)
    {
      //double z = clusters[i]->getZ();
      //double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
      double r = points[i].first;
      double z = points[i].second;

      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+pow(r,2);                //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }
   a=(points.size()*xysum-xsum*ysum)/(points.size()*x2sum-xsum*xsum);            //calculate slope
   b=(x2sum*ysum-xsum*xysum)/(x2sum*points.size()-xsum*xsum);            //calculate intercept

   if(Verbosity() > 10)
     {
       for (unsigned int i=0;i<points.size(); ++i)
	 {
	   double r = points[i].first;
	   double z_fit = a * r + b;                    //to calculate z(fitted) at given r points
	   std::cout << " r " << r << " z " << points[i].second << " z_fit " << z_fit << std::endl; 
	 } 
     }

    return;
}   

void  PHTpcTrackSeedVertexAssoc::line_fit_clusters(std::vector<TrkrCluster*> clusters, double &a, double &b)
{
  std::vector<std::pair<double,double>> points;
  
   for (unsigned int i=0; i<clusters.size(); ++i)
     {
       double z = clusters[i]->getZ();
       double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));

       points.push_back(make_pair(r,z));
     }

   line_fit(points, a, b);

    return;
}
