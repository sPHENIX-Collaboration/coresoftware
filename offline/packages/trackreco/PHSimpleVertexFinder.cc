#include "PHSimpleVertexFinder.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair
#include <iomanip>

#include <vector>
#include <cassert>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <functional>

#include <Eigen/Dense>

//____________________________________________________________________________..
PHSimpleVertexFinder::PHSimpleVertexFinder(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
PHSimpleVertexFinder::~PHSimpleVertexFinder()
{

}

//____________________________________________________________________________..
int PHSimpleVertexFinder::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  return ret;
}

//____________________________________________________________________________..
int PHSimpleVertexFinder::process_event(PHCompositeNode */*topNode*/)
{

  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size()  << std::endl;

  // reset maps
  _vertex_track_map.clear();
  _track_pair_map.clear();
  _track_pair_pca_map.clear();
  _vertex_position_map.clear();
  _vertex_covariance_map.clear();
  _vertex_set.clear();
  
  // Find all instances where two tracks have a dca of < _dcacut,  and capture the pair details
  // Fills _track_pair_map and _track_pair_pca_map
  checkDCAs();

  /// If we didn't find any matches, try again with a slightly larger DCA cut
  if(_track_pair_map.size() == 0)
    {
      _dcacut *= 1.5;
      checkDCAs();
    }

  // get all connected pairs of tracks by looping over the track_pair map
  std::vector<std::set<unsigned int>> connected_tracks = findConnectedTracks();

  // make vertices - each set of connected tracks is a vertex
  for(unsigned int ivtx = 0; ivtx < connected_tracks.size(); ++ivtx)
    {
      if(Verbosity() > 1) std::cout << "process vertex " << ivtx << std::endl;

      for(auto it : connected_tracks[ivtx])
	{
	  unsigned int id = it;
	  _vertex_track_map.insert(std::make_pair(ivtx, id));
	  if(Verbosity() > 2)  std::cout << " adding track " << id << " to vertex " << ivtx << std::endl;	  
	}      
    }

  // make a list of vertices
  for(auto it : _vertex_track_map)
    {
      if(Verbosity() > 1) std::cout << " vertex " << it.first << " track " << it.second << std::endl;      
      _vertex_set.insert(it.first);
    }
  
  // this finds average vertex positions after removal of outlying track pairs 
  removeOutlierTrackPairs();

  // average covariance for accepted tracks
  for(auto it : _vertex_set)
    {
      matrix_t avgCov = matrix_t::Zero();
      double cov_wt = 0.0;
      
      auto ret = _vertex_track_map.equal_range(it);
      for (auto cit=ret.first; cit!=ret.second; ++cit)
	{
	  unsigned int trid = cit->second;
	  matrix_t cov;
	  auto track = _track_map->get(trid);
	  for(int i = 0; i < 3; ++i)
	    for(int j = 0; j < 3; ++j)
	      cov(i,j) = track->get_error(i,j);
	  
	  avgCov += cov;
	  cov_wt++;
	}
      
      avgCov /= sqrt(cov_wt);
      if(Verbosity() > 2)
	{
	  std::cout << "Average covariance for vertex " << it << " is:" << std::endl;
	  std::cout << std::setprecision(8) << avgCov << std::endl; 
	}
      _vertex_covariance_map.insert(std::make_pair(it, avgCov));
    }
  
  // Write the vertices to the vertex map on the node tree
  //==============================================
  if(_vertex_track_map.size() > 0)
    _svtx_vertex_map->clear();

  for(auto it : _vertex_set)
    {
      auto svtxVertex = std::make_unique<SvtxVertex_v1>();

      svtxVertex->set_chisq(0.0);
      svtxVertex->set_ndof(0);
      svtxVertex->set_t0(0);
      svtxVertex->set_id(it);

      auto ret = _vertex_track_map.equal_range(it);
      for (auto cit=ret.first; cit!=ret.second; ++cit)
      {
	unsigned int trid = cit->second;

	if(Verbosity() > 1) std::cout << "   vertex " << it << " insert track " << trid << std::endl; 
	svtxVertex->insert_track(trid);
	_track_map->get(trid)->set_vertex_id(it);
      }

      Eigen::Vector3d pos = _vertex_position_map.find(it)->second;
      svtxVertex->set_x(pos.x());  
      svtxVertex->set_y(pos.y());
      svtxVertex->set_z(pos.z());
      if(Verbosity() > 1) std::cout << "   vertex " << it << " insert pos.x " << pos.x() << " pos.y " << pos.y() << " pos.z " << pos.z() << std::endl; 

      auto vtxCov = _vertex_covariance_map.find(it)->second;
      for(int i = 0; i < 3; ++i) 
	{
	  for(int j = 0; j < 3; ++j)
	    {
	      svtxVertex->set_error(i, j, vtxCov(i,j)); 
	    }
	}
      
      _svtx_vertex_map->insert(svtxVertex.release());      
    }
  
  /// Iterate through the tracks and assign the closest vtx id to 
  /// the track position for propagating back to the vtx. Catches any
  /// tracks that were missed or were not  compatible with any of the
  /// identified vertices
  //=================================================
  for(const auto& [trackkey, track] : *_track_map)
    {
      auto vtxid = track->get_vertex_id();

      /// If there is a vertex already assigned, keep going
      if(_svtx_vertex_map->get(vtxid))
	{ continue; }
      
      float maxdz = std::numeric_limits<float>::max();
      unsigned int newvtxid = std::numeric_limits<unsigned int>::max();
      
      for(const auto& [vtxkey, vertex] : *_svtx_vertex_map)
	{
	  float dz = track->get_z() - vertex->get_z();
	  if(fabs(dz) < maxdz)
	    {
	      maxdz = dz;
	      newvtxid = vtxkey;
	    }
	}
      
      track->set_vertex_id(newvtxid);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int PHSimpleVertexFinder::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleVertexFinder::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  /// Check that it is there
  if (!dstNode)
  {
    std::cerr << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in PHActsInitialVertexFinder::createNodes");
  }

  /// Get the tracking subnode
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  _svtx_vertex_map = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  
  if(!_svtx_vertex_map)
    {
      _svtx_vertex_map = new SvtxVertexMap_v1;
      PHIODataNode<PHObject>* vertexNode = new PHIODataNode<PHObject>( 
		   _svtx_vertex_map, "SvtxVertexMap","PHObject");

      svtxNode->addNode(vertexNode);

    }

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHSimpleVertexFinder::GetNodes(PHCompositeNode* topNode)
{

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

 

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHSimpleVertexFinder::checkDCAs()
{
  // Loop over tracks and check for close DCA match with all other tracks
  for(auto tr1_it = _track_map->begin(); tr1_it != _track_map->end(); ++tr1_it)
    {
      auto id1 = tr1_it->first;
      auto tr1 = tr1_it->second;
      if(tr1->get_quality() > _qual_cut) continue;
      if(_require_mvtx)
	{
	  unsigned int nmvtx = 0;
	  TrackSeed *siliconseed = tr1->get_silicon_seed();
	  if(!siliconseed)
	    { continue; }

	  for(auto clusit = siliconseed->begin_cluster_keys(); clusit != siliconseed->end_cluster_keys(); ++clusit)
	    {
	      if(TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId )
		{
		  nmvtx++;
		}
	      if(nmvtx >= _nmvtx_required) break;
	    }
	  if(nmvtx < _nmvtx_required) continue;
	  if(Verbosity() > 3) std::cout << " tr1 has nmvtx at least " << nmvtx << std::endl;
	}
      
      // look for close DCA matches with all other such tracks
      for(auto tr2_it = std::next(tr1_it); tr2_it != _track_map->end(); ++tr2_it)
	{
	  auto id2 = tr2_it->first;
	  auto tr2 = tr2_it->second;
	  if(tr2->get_quality() > _qual_cut) continue;
	  if(_require_mvtx)
	    {
	      unsigned int nmvtx = 0;
	      TrackSeed *siliconseed = tr2->get_silicon_seed();
	      if(!siliconseed)
		{ continue; }
	     
	      for(auto clusit = siliconseed->begin_cluster_keys(); clusit != siliconseed->end_cluster_keys(); ++clusit)
		{
		  if(TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId)
		    {
		      nmvtx++;
		    }
		  if(nmvtx >= _nmvtx_required) break;
		}
	      if(nmvtx < _nmvtx_required) continue;
	      if(Verbosity() > 3)  std::cout << " tr2 has nmvtx at least " << nmvtx << std::endl;
	    }
	  
	  // find DCA of these two tracks
	  if(Verbosity() > 3) std::cout << "Check DCA for tracks " << id1 << " and  " << id2 << std::endl;
	  
	  findDcaTwoTracks(tr1, tr2);

	}
    }


}

void PHSimpleVertexFinder::findDcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2)
{
  if(tr1->get_pt() < _track_pt_cut) return;
  if(tr2->get_pt() < _track_pt_cut) return;

  unsigned int id1 = tr1->get_id();
  unsigned int id2 = tr2->get_id();

  // get the line equations for the tracks
  
  Eigen::Vector3d a1(tr1->get_x(), tr1->get_y(), tr1->get_z());
  Eigen::Vector3d b1(tr1->get_px() / tr1->get_p(), tr1->get_py() / tr1->get_p(), tr1->get_pz() / tr1->get_p());
  Eigen::Vector3d a2(tr2->get_x(), tr2->get_y(), tr2->get_z());
  Eigen::Vector3d b2(tr2->get_px() / tr2->get_p(), tr2->get_py() / tr2->get_p(), tr2->get_pz() / tr2->get_p());
  
  Eigen::Vector3d PCA1(0,0,0);
  Eigen::Vector3d PCA2(0,0,0);  
  double dca = dcaTwoLines(a1, b1, a2,  b2, PCA1, PCA2);

  // check dca cut is satisfied, and that PCA is close to beam line
  if( fabs(dca) < _dcacut && (fabs(PCA1.x()) < _beamline_xy_cut && fabs(PCA1.y()) < _beamline_xy_cut) )
    {
      if(Verbosity() > 3)
	{
	  std::cout << " good match for tracks " << tr1->get_id() << " and " << tr2->get_id() << " with pT " << tr1->get_pt()  << " and " << tr2->get_pt() << std::endl;
	  std::cout << "    a1.x " << a1.x() << " a1.y " << a1.y() << " a1.z " << a1.z() << std::endl;
	  std::cout << "    a2.x  " << a2.x()  << " a2.y " << a2.y() << " a2.z " << a2.z() << std::endl;
	  std::cout << "    PCA1.x() " << PCA1.x() << " PCA1.y " << PCA1.y() << " PCA1.z " << PCA1.z() << std::endl;
	  std::cout << "    PCA2.x() " << PCA2.x() << " PCA2.y " << PCA2.y() << " PCA2.z " << PCA2.z() << std::endl;      
	  std::cout << "    dca " << dca << std::endl;
	}  

      // capture the results for successful matches
      _track_pair_map.insert(std::make_pair(id1,std::make_pair(id2, dca)));
      _track_pair_pca_map.insert( std::make_pair(id1, std::make_pair(id2, std::make_pair(PCA1, PCA2))) );
    }

  return;
}

double PHSimpleVertexFinder::dcaTwoLines(const Eigen::Vector3d &a1,const Eigen::Vector3d &b1,
					 const Eigen::Vector3d &a2,const Eigen::Vector3d &b2,
					 Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2)
{
  // The shortest distance between two skew lines described by
  //  a1 + c * b1
  //  a2 + d * b2
  // where a1, a2, are vectors representing points on the lines, b1, b2 are direction vectors, and c and d are scalars
  // is:
  // dca = (b1 x b2) .(a2-a1) / |b1 x b2|

  // bcrossb/mag_bcrossb is a unit vector perpendicular to both direction vectors b1 and b2
  auto bcrossb = b1.cross(b2);
  auto mag_bcrossb = bcrossb.norm();
  // a2-a1 is the vector joining any arbitrary points on the two lines
  auto aminusa = a2-a1;

  // The DCA of these two lines is the projection of a2-a1 along the direction of the perpendicular to both 
  // remember that a2-a1 is longer than (or equal to) the dca by definition
  double dca = 999;
  if( mag_bcrossb != 0)
    dca = bcrossb.dot(aminusa) / mag_bcrossb;
  else
    return dca;   // same track, skip combination
  
  // get the points at which the normal to the lines intersect the lines, where the lines are perpendicular
  // Assume the shortest distance occurs at points A on line 1, and B on line 2, call the line AB
  //    AB = a2+d*b2 - (a1+c*b1)
  // we need to find c and d where AB is perpendicular to both lines. so AB.b1 = 0 and AB.b2 = 0
  //    ( (a2 -a1) + d*b2 -c*b1 ).b1 = 0
  //    ( (a2 -a1) + d*b2 -c*b1 ).b2 = 0
  // so we have two simultaneous equations in 2 unknowns
  //    (a2-a1).b1 + d*b2.b1 -c*b1.b1 = 0 => d*b2.b1 = c*b1.b1 - (a2-a1).b1 => d = (1/b2.b1) * (c*b1.b1 - (a2-a1).b1) 
  //    (a2-a1).b2 + d*b2.b2 - c*b1.b2 = 0 => c*b1.b2 =  (a2-a1).b2 + [(1/b2.b1) * (c*b1*b1 -(a2-a1).b1)}*b2.b2
  //    c*b1.b2 - (1/b2.b1) * c*b1.b1*b2.b2 = (a2-a1).b2 - (1/b2.b1)*(a2-a1).b1*b2.b2
  //    c *(b1.b2 - (1/b2.b1)*b1.b1*b2.b2)  = (a2-a1).b2 - (1/b2.b1)*(a2-a1).b1*b2.b2
  // call this: c*X = Y
  // plug back into the d equation
  //     d = c*b1.b1 / b2.b1 - (a2-a1).b1 / b2.b1
  // and call the d equation: d = c * F - G 

  double X =  b1.dot(b2) - b1.dot(b1) * b2.dot(b2) / b2.dot(b1);
  double Y =  (a2.dot(b2) - a1.dot(b2)) - (a2.dot(b1) - a1.dot(b1)) * b2.dot(b2) / b2.dot(b1) ;
  double c = Y/X;

  double F = b1.dot(b1) / b2.dot(b1); 
  double G = - (a2.dot(b1) - a1.dot(b1)) / b2.dot(b1);
  double d = c * F + G;

  // then the points of closest approach are:
  PCA1 = a1+c*b1;
  PCA2 = a2+d*b2;

  return dca;


}

std::vector<std::set<unsigned int>> PHSimpleVertexFinder::findConnectedTracks()
{
 std::vector<std::set<unsigned int>> connected_tracks;
  std::set<unsigned int> connected;
  std::set<unsigned int> used;
  for(auto it : _track_pair_map)
    {
      unsigned int id1 = it.first;
      unsigned int id2 = it.second.first;

      if( (used.find(id1) != used.end()) && (used.find(id2) != used.end()) )
	{
	  if(Verbosity() > 3) std::cout << " tracks " << id1 << " and " << id2 << " are both in used , skip them" << std::endl;
	  continue;
	}
      else if( (used.find(id1) == used.end()) && (used.find(id2) == used.end()))
	{
	  if(Verbosity() > 3) std::cout << " tracks " << id1 << " and " << id2 << " are both not in used , start a new connected set" << std::endl;
	  // close out and start a new connections set
	  if(connected.size() > 0)
	    {
	      connected_tracks.push_back(connected);
	      connected.clear();
	      if(Verbosity() > 3) std::cout << "           closing out set " << std::endl;
	    }
	}

      // get everything connected to id1 and id2
      connected.insert(id1);
      used.insert(id1);
      connected.insert(id2);
      used.insert(id2);
      for(auto cit :  _track_pair_map)
	{
	  unsigned int id3 = cit.first;
	  unsigned int id4 = cit.second.first;
	  if( (connected.find(id3) != connected.end()) || (connected.find(id4) != connected.end()) )
	    {
	      if(Verbosity() > 3) std::cout << " found connection to " << id3 << " and " << id4 << std::endl;
	      connected.insert(id3);
	      used.insert(id3);
	      connected.insert(id4);
	      used.insert(id4);
	    }
	}
    }
  
  // close out the last set
  if(connected.size() > 0)
    {
      connected_tracks.push_back(connected);
      connected.clear();
      if(Verbosity() > 3) std::cout << "           closing out last set " << std::endl;
    }
    
  if(Verbosity() > 3)   std::cout << "connected_tracks size " << connected_tracks.size() << std::endl;

  return connected_tracks;
}

void PHSimpleVertexFinder::removeOutlierTrackPairs()
{
  //  Note: std::multimap<unsigned int, std::pair<unsigned int, std::pair<Eigen::Vector3d,  Eigen::Vector3d>>>  _track_pair_pca_map
 
  for(auto it : _vertex_set) 
    {
      unsigned int vtxid = it;
      if(Verbosity() > 1) std::cout << "calculate average position for vertex " << vtxid << std::endl; 

      // we need the median values of the x and y positions 
      std::vector<double> vx;
      std::vector<double> vy;
      std::vector<double> vz;

      double pca_median_x = 0.;
      double pca_median_y = 0.;
      double pca_median_z = 0.;

      Eigen::Vector3d new_pca_avge(0.,0.,0.);
      double new_wt = 0.0;
      
      auto ret = _vertex_track_map.equal_range(vtxid);
      
      // Start by getting the positions for this vertex into vectors for the median calculation
      for (auto cit=ret.first; cit!=ret.second; ++cit)
	{	  
	  unsigned int tr1id = cit->second;
	  if(Verbosity() > 2) std::cout << "   vectors: get entries for track " << tr1id << " for vertex " << vtxid << std::endl; 
	  
	  // find all pairs for this vertex with tr1id
	  auto pca_range = _track_pair_pca_map.equal_range(tr1id);
	  for (auto pit=pca_range.first; pit!=pca_range.second; ++pit)
	    {
	      unsigned int tr2id = pit->second.first;
	      
	      Eigen::Vector3d PCA1 = pit->second.second.first;
	      Eigen::Vector3d PCA2 = pit->second.second.second;
	      
	      if(Verbosity() > 2)
		std::cout << " vectors: tr1id " << tr1id << " tr2id " << tr2id
			  << " PCA1 " << PCA1.x() << "  " << PCA1.y() << "  " << PCA1.z()
			  << " PCA2 " << PCA2.x() << "  " << PCA2.y() << "  " << PCA2.z()
			  << std::endl;
	      
	      vx.push_back(PCA1.x());
	      vx.push_back(PCA2.x());
	      vy.push_back(PCA1.y());
	      vy.push_back(PCA2.y());
	      vz.push_back(PCA1.z());
	      vz.push_back(PCA2.z());
	    }
	}
      
      // Get the medians for this vertex
      // Using the median as a reference for rejecting outliers only makes sense for more than 2 tracks
      if(vx.size() < 3)
	{
	  new_pca_avge.x() = getAverage(vx);
	  new_pca_avge.y() = getAverage(vy);
	  new_pca_avge.z() = getAverage(vz);
	  _vertex_position_map.insert(std::make_pair(vtxid, new_pca_avge));
	  if(Verbosity() > 1) 
	    std::cout << " Vertex has only 2 tracks, use average for PCA: " << new_pca_avge.x() << "  " << new_pca_avge.y() << "  " << new_pca_avge.z() << std::endl; 

	  // done with this vertex
	  continue;
	}
      
      pca_median_x = getMedian(vx);
      pca_median_y = getMedian(vy);
      pca_median_z = getMedian(vz);
      if(Verbosity() > 1) std::cout << "Median values: x " << pca_median_x << " y " << pca_median_y << " z : " << pca_median_z << std::endl;
      
      // Make the average vertex position with outlier rejection wrt the median
      for (auto cit=ret.first; cit!=ret.second; ++cit)
	{	  
	  unsigned int tr1id = cit->second;
	  if(Verbosity() > 2) std::cout << "   average: get entries for track " << tr1id << " for vertex " << vtxid << std::endl; 
	  
	  // find all pairs for this vertex with tr1id
	  auto pca_range = _track_pair_pca_map.equal_range(tr1id);
	  for (auto pit=pca_range.first; pit!=pca_range.second; ++pit)
	    {
	      unsigned int tr2id = pit->second.first;
	      
	      Eigen::Vector3d PCA1 = pit->second.second.first;
	      Eigen::Vector3d PCA2 = pit->second.second.second;
	      
	      if(
		 fabs(PCA1.x() - pca_median_x) < _outlier_cut &&
		 fabs(PCA1.y() - pca_median_y) < _outlier_cut &&
		 fabs(PCA2.x() - pca_median_x) < _outlier_cut &&
		 fabs(PCA2.y() - pca_median_y) < _outlier_cut 
		 )
		{
		  // good track pair, add to new average
		  
		  new_pca_avge += PCA1;
		  new_wt++;
		  new_pca_avge += PCA2;
		  new_wt++;	  
		}
	      else
		{
		  if(Verbosity() > 1) std::cout << "Reject pair with tr1id " << tr1id << " tr2id " << tr2id << std::endl;
		}
	    }
	}
      if(new_wt > 0.0)      
	new_pca_avge = new_pca_avge / new_wt;
      else
	{
	  // There were no pairs that survived the track cuts, use the median values
	  new_pca_avge.x() = pca_median_x;
	  new_pca_avge.y() = pca_median_y;
	  new_pca_avge.z() = pca_median_z;
	}

      _vertex_position_map.insert(std::make_pair(vtxid, new_pca_avge));

    }
  
  return;
 }
  
double PHSimpleVertexFinder::getMedian(std::vector<double> &v)
{
  double median = 0.0;

  if( (v.size() % 2) == 0)
    {
      // even number of entries
      // we want the average of the middle two numbers, v.size()/2 and v.size()/2-1
      auto m1 = v.begin() + v.size()/2;
      std::nth_element(v.begin(), m1, v.end());
      double median1 =  v[v.size()/2]; 

      auto m2 = v.begin() + v.size()/2 - 1;
      std::nth_element(v.begin(), m2, v.end());
      double median2 =  v[v.size()/2 - 1]; 

      median = (median1 + median2) / 2.0; 
      if(Verbosity() > 2) std::cout << "The vector size is " << v.size() 
				    << " element m is " << v.size() / 2  << " = " << v[v.size()/2] 
				    << " element m-1 is " << v.size() / 2 -1 << " = " << v[v.size()/2-1] 
				    <<  std::endl;
    } 
  else
    {
      // odd number of entries
      auto m = v.begin() + v.size()/2;
      std::nth_element(v.begin(), m, v.end());
      median =  v[v.size()/2];
      if(Verbosity() > 2) std::cout << "The vector size is " << v.size() << " element m is " << v.size() / 2 << " = " << v[v.size()/2] <<  std::endl;
    }

    return median ;
}
double PHSimpleVertexFinder::getAverage(std::vector<double> &v)
{
  double avge = 0.0;
  double wt = 0.0;
  for(auto it : v)
    {
      avge += it;
      wt++;
    }

  avge /= wt;
  if(Verbosity() > 2)
    {
      std::cout << " average = " << avge << std::endl;
    }
  
  return avge;
}


