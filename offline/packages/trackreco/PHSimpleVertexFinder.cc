#include "PHSimpleVertexFinder.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TpcSeedTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex_v1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

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
  _track_pca_map.clear();
  _vertex_position_map.clear();
  _connected.clear();
  _vertex_set.clear();
  
  // define maps
  std::set<unsigned int> track_used_list;
  
  // make a list of tracks that did not make the keep list
  for(auto tr1_it = _track_map->begin(); tr1_it != _track_map->end(); ++tr1_it)
    {
      auto id1 = tr1_it->first;
      
      auto set1_it = track_used_list.find(id1);
      if(set1_it != track_used_list.end())
	{
	  if(Verbosity() > 1)
	    {
	      std::cout << "   track 1 id " << id1 << " already used in a vertex, skip it" << std::endl;
	    }
	  continue;
	}
      
      // this track is not used in a vertex yet, look for close DCA matches with all other such tracks
      for(auto tr2_it = std::next(tr1_it); tr2_it != _track_map->end(); ++tr2_it)
	{
	  auto id2 = tr2_it->first;
	  
	  auto set2_it = track_used_list.find(id2);
	  if(set2_it != track_used_list.end())
	    {
	      if(Verbosity() > 1)
		{
		  std::cout << "   track 2 id " << id2 << " already used in a vertex, skip it" << std::endl;
		}
	      continue;
	    }
	  
	  // Check DCA of these two tracks
	  if(Verbosity() > 1) std::cout << "Check DCA for tracks " << id1 << " and " << id2 << std::endl;
	  
	  auto tr1 = tr1_it->second;
	  auto tr2 = tr2_it->second;
	  double dca = findDcaTwoTracks(tr1, tr2);
	  if(Verbosity() > 1) std::cout << " final returned dca = " << dca << std::endl; 
	  
	  if(fabs(dca) < dcacut)
	    {
	      if(Verbosity() > 1) std::cout << "      good match for tracks " << id1 << " and " << id2 << std::endl;	      
	      track_used_list.insert(id2);
	      track_used_list.insert(id1);
	    }
	  
	}
      
    }

  // get all connected pairs of tracks, make vertices      
  unsigned int vtxid = 0;
  for(auto it : _connected)
    {
      unsigned int id1 = it;
      _vertex_track_map.insert(std::make_pair(vtxid, id1));
      if(Verbosity() > 2)  std::cout << " adding track " << id1 << " to vertex " << vtxid << std::endl;	  
      auto ret = _track_pair_map.equal_range(id1);

      for (auto cit=ret.first; cit!=ret.second; ++cit)
	{
	  unsigned int id2 = cit->second.first;
	  _vertex_track_map.insert(std::make_pair(vtxid, id2));
	  if(Verbosity() > 2) std::cout << " adding track " << id2 << " to vertex " << vtxid << std::endl;	  
	}
      vtxid++;
    }
  
  for(auto it : _vertex_track_map)
    {
      if(Verbosity() > 1) std::cout << " vertex " << it.first << " track " << it.second << std::endl;      
      _vertex_set.insert(it.first);
    }

  // Calculate the vertex positions
  for(auto it : _vertex_set) 
    {
      unsigned int vtxid = it;
      if(Verbosity() > 2) std::cout << "calculate average position for vertex " << vtxid << std::endl; 

      Eigen::Vector3d pca_avge(0,0,0);
      double wt = 0.0;

      auto ret = _vertex_track_map.equal_range(vtxid);
      for (auto cit=ret.first; cit!=ret.second; ++cit)
      {
	unsigned int trid = cit->second;
	if(Verbosity() > 2) std::cout << "   get entries for track " << trid << " for vertex " << vtxid << std::endl; 

	auto pca_range = _track_pca_map.equal_range(trid);
	for (auto pit=pca_range.first; pit!=pca_range.second; ++pit)
	  {
	    unsigned int tr2id = pit->second.first;
	    Eigen::Vector3d pca = pit->second.second;
	    if(Verbosity() > 2) std::cout << "       vertex " << vtxid << " tr1 " << trid << " tr2 " << tr2id << " pca.x " << pca.x() << " pca.y " << pca.y() << " pca.z " << pca.z()  << std::endl; 

	    pca_avge += pca;
	    wt++;
	  }
      }

      pca_avge /= wt;
      if(Verbosity() > 1) std::cout << "vertex " << vtxid << " average pca.x " << pca_avge.x() << " pca.y " << pca_avge.y() << " pca.z " << pca_avge.z() << std::endl;

      // capture the position
      _vertex_position_map.insert(std::make_pair(vtxid, pca_avge));
    }       

  // Write the vertices to the vertex map on the node tree
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
      }

      Eigen::Vector3d pos = _vertex_position_map.find(it)->second;
      svtxVertex->set_x(pos.x());  
      svtxVertex->set_y(pos.y());
      svtxVertex->set_z(pos.z());
      if(Verbosity() > 1) std::cout << "   vertex " << it << " insert pos " << pos << std::endl; 

      for(int i = 0; i < 3; ++i) 
	{
	  for(int j = 0; j < 3; ++j)
	    {
	      svtxVertex->set_error(i, j, 0.0); 
	    }
	}
      
      _svtx_vertex_map->insert(svtxVertex.release());      
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}



int PHSimpleVertexFinder::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHSimpleVertexFinder::GetNodes(PHCompositeNode* topNode)
{

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _svtx_vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_svtx_vertex_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxVertexMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

double PHSimpleVertexFinder::findDcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2)
{
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
  
  if( fabs(dca) < dcacut)
    {
      if(Verbosity() > 1)
	{
	  std::cout << " good match for tracks with z locations " << a1.z()  << " and " << a2.z() << std::endl;
	  std::cout << "    a1.x " << a1.x() << " a1.y " << a1.y() << " a1.z " << a1.z() << std::endl;
	  std::cout << "    a2.x  " << a2.x()  << " a2.y " << a2.y() << " a2.z " << a2.z() << std::endl;
	  std::cout << "    PCA1.x() " << PCA1.x() << " PCA1.y " << PCA1.y() << " PCA1.z " << PCA1.z() << std::endl;
	  std::cout << "    PCA2.x() " << PCA2.x() << " PCA2.y " << PCA2.y() << " PCA2.z " << PCA2.z() << std::endl;      
	  std::cout << "    dca " << dca << std::endl;
	}  

      // capture the results for successful matches
      // we want to add this pair to a vertex
      // base the vertex on the first track ID

      _connected.insert(id1);
      _track_pair_map.insert(std::make_pair(id1,std::make_pair(id2, dca)));
      _track_pca_map.insert(std::make_pair(id1, std::make_pair(id2, PCA1)));
      _track_pca_map.insert(std::make_pair(id2, std::make_pair(id1, PCA2)));
    }



  return dca;
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

  //std::cout << " PCA1.x() " << PCA1.x() << " PCA1.y " << PCA1.y() << " PCA1.z " << PCA1.z() << std::endl;
  //std::cout << " PCA2.x() " << PCA2.x() << " PCA2.y " << PCA2.y() << " PCA2.z " << PCA2.z() << std::endl;
 
  return dca;


}
