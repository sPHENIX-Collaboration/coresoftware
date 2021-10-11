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
  _revised_vertex_track_map.clear();
  _vertex_cov_map.clear();
  _track_pair_map.clear();
  _track_pair_pca_map.clear();
  _vertex_position_map.clear();
  _revised_vertex_position_map.clear();
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

  // Get the initial estimate of the vertex positions from the connected track pairs PCA values
  // This fills   _vertex_position_map
  calculateInitialVertexPositions();

  // Find and reject  pairs with outlier PCA values
  findOutlierTrackPairs();
  
  // Calculate the PCA of the tracks to the initial vertex estimates, reject outliers, recalculate the vertex positions
  // Fills  _revised_vertex_track_map,  _revised_vertex_position_map,  _vertex_cov_map
  calculateRevisedVertexPositions();

  // Write the revised vertices to the vertex map on the node tree
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

      //auto ret = _revised_vertex_track_map.equal_range(it);
      auto ret = _vertex_track_map.equal_range(it);
      for (auto cit=ret.first; cit!=ret.second; ++cit)
      {
	unsigned int trid = cit->second;

	if(Verbosity() > 1) std::cout << "   vertex " << it << " insert track " << trid << std::endl; 
	svtxVertex->insert_track(trid);
	_track_map->get(trid)->set_vertex_id(it);
      }

      //Eigen::Vector3d pos = _revised_vertex_position_map.find(it)->second;
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
	  for(auto clusit = tr1->begin_cluster_keys(); clusit != tr1->end_cluster_keys(); ++clusit)
	    {
	      if(TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId )
		{
		  nmvtx++;
		}
	      if(nmvtx > 2) break;
	    }
	  if(nmvtx <3) continue;
	  if(Verbosity() > 1) std::cout << " tr1 has nmvtx " << nmvtx << std::endl;
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
	      for(auto clusit = tr2->begin_cluster_keys(); clusit != tr2->end_cluster_keys(); ++clusit)
		{
		  if(TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId)
		    {
		      nmvtx++;
		    }
		  if(nmvtx > 2) break;
		}
	      if(nmvtx <3) continue;
	      if(Verbosity() > 1)  std::cout << " tr2 has nmvtx " << nmvtx << std::endl;
	    }
	  
	  // find DCA of these two tracks
	  if(Verbosity() > 2) std::cout << "Check DCA for tracks " << id1 << " and  " << id2 << std::endl;
	  
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
      if(Verbosity() > 1)
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
	  if(Verbosity() > 2) std::cout << " tracks " << id1 << " and " << id2 << " are both in used , skip them" << std::endl;
	  continue;
	}
      else if( (used.find(id1) == used.end()) && (used.find(id2) == used.end()))
	{
	  if(Verbosity() > 2) std::cout << " tracks " << id1 << " and " << id2 << " are both not in used , start a new connected set" << std::endl;
	  // close out and start a new connections set
	  if(connected.size() > 0)
	    {
	      connected_tracks.push_back(connected);
	      connected.clear();
	      if(Verbosity() > 2) std::cout << "           closing out set " << std::endl;
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
	      if(Verbosity() > 2) std::cout << " found connection to " << id3 << " and " << id4 << std::endl;
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
      if(Verbosity() > 2) std::cout << "           closing out last set " << std::endl;
    }
    
  if(Verbosity() > 1)   std::cout << "connected_tracks size " << connected_tracks.size() << std::endl;

  return connected_tracks;
}

void PHSimpleVertexFinder::calculateInitialVertexPositions()
{
  // Calculate the vertex positions from the pair  PCA's
  // This is the initial estimate, but it is pretty good
  //=======================================
  for(auto it : _vertex_set) 
    {
      unsigned int vtxid = it;
      if(Verbosity() > 2) std::cout << "calculate average position for vertex " << vtxid << std::endl; 

      Eigen::Vector3d pca_avge(0,0,0);
      double wt = 0.0;
      
      auto ret = _vertex_track_map.equal_range(vtxid);
      for (auto cit=ret.first; cit!=ret.second; ++cit)
      {
	unsigned int tr1id = cit->second;
	if(Verbosity() > 2) std::cout << "   get entries for track " << tr1id << " for vertex " << vtxid << std::endl; 

	// get the covariance matrix for the track position
	matrix_t cov1 = matrix_t::Zero();
	  auto track1 = _track_map->get(tr1id);
	  for(int i = 0; i < 3; ++i)
	    for(int j = 0; j < 3; ++j)
	      {cov1(i,j) = track1->get_error(i,j);}
	  //matrix_t cov1_wt = cov1.fullPivLu().inverse();

	// find all pairs with trid
	auto pca_range = _track_pair_pca_map.equal_range(tr1id);
	for (auto pit=pca_range.first; pit!=pca_range.second; ++pit)
	  {
	    unsigned int tr2id = pit->second.first;
	    matrix_t cov2 = matrix_t::Zero();
	    auto track2 = _track_map->get(tr2id);
	    for(int i = 0; i < 3; ++i)
	      for(int j = 0; j < 3; ++j)
		{cov2(i,j) =  track2->get_error(i,j);}
	    //matrix_t cov2_wt = cov2.fullPivLu().inverse();

	    Eigen::Vector3d pca1 = pit->second.second.first;
	    pca_avge += pca1;
	    wt += 1.0;

	    Eigen::Vector3d pca2 = pit->second.second.second;
	    pca_avge += pca2;
	    wt += 1.0;
	    
	    if(Verbosity() > 2) std::cout << "       vertex " << vtxid << " tr1 " << tr1id << " pt " << track1->get_pt() << " tr2 " << tr2id << " pt " << track2->get_pt() << " pca1.x " << pca1.x() << " pca1.y " << pca1.y() << " pca1.z " << pca1.z()  
					  << " pca2.x " << pca2.x() << " pca2.y " << pca2.y() << " pca2.z " << pca2.z()  << std::endl; 
	  }
      }

      pca_avge = pca_avge / wt;
 
      if(Verbosity() > 1) std::cout << "vertex " << vtxid << " average pca.x " << pca_avge.x() << " pca.y " << pca_avge.y() << " pca.z " << pca_avge.z() << std::endl;

      _vertex_position_map.insert(std::make_pair(vtxid, pca_avge));
    }
}


void PHSimpleVertexFinder::calculateRevisedVertexPositions()
{
  // start with the initial vertex list
  for(auto it : _vertex_set) 
    {
      unsigned int vtxid = it;
      Eigen::Vector3d pca_avge = _vertex_position_map.find(vtxid)->second;

     if(Verbosity() > 2) std::cout << "calculate revised average position for vertex " << vtxid << std::endl; 

      // We have a good initial estimate of the vertex position
      // loop over all of the tracks in the vertex and get the PCA of each track to our initial vertex estimate
      // Use those PCA's to recalculate the vertex position
      //===========================================================================
      Eigen::Vector3d new_pca_avge(0,0,0);
      double new_wt = 0.0;
      matrix_t cov_wt = matrix_t::Zero();

      auto ret = _vertex_track_map.equal_range(vtxid);
      for (auto cit=ret.first; cit!=ret.second; ++cit)
      {
	unsigned int trid = cit->second;
	auto tr = _track_map->get(trid);
	std::cout << " track " << trid << " with pt " << tr->get_pt() << std::endl;

	// We have the line equation	
	Eigen::Vector3d a(tr->get_x(), tr->get_y(), tr->get_z());
	Eigen::Vector3d b(tr->get_px() / tr->get_p(), tr->get_py() / tr->get_p(), tr->get_pz() / tr->get_p());

	// calculate the DCA of this line to our vertex estimate
	Eigen::Vector3d Aa = pca_avge - a;
	Eigen::Vector3d perp = b.cross(Aa);
	double dca = perp.norm() / b.norm();
	if(Verbosity() > 1) std::cout << " point line dca " << dca << " trid " << trid << std::endl;

	// get the point of closest approach. The point is vector A. The closest point on the line is B
	// The line is a + c*b, where c is a scalar, we want to solve for c when the ,ine is at B
	//    AB = A - (a + c*b) = A - a - c*b
	//  When AB is perpendicular to the line:
	//    AB.b = 0 =  (A-a).b - c*b.b
	//    so: c = (A-a).b / b.b,   and then: B = a + c*b

	double c = Aa.dot(b) / b.dot(b);
	Eigen::Vector3d B = a + c*b;
	if(Verbosity() > 1) std::cout <<  " trid " << trid << " point line PCA " << B << std::endl;

	// reject any bad outliers
	//double cut_mag = 2.0;  // 2 seems to be optimum
	double cut_mag = 10.0;  // 2 seems to be optimum
	if( (abs(B.x()-pca_avge.x()) < cut_mag * _dcacut) && (abs(B.y()-pca_avge.y()) < cut_mag * _dcacut) ) 
	  {
	    // track is OK, remake the average
	    _revised_vertex_track_map.insert(std::make_pair(vtxid, trid));

	    if(_use_track_covariance)
	      {
		// get the covariance matrix for the track position
		matrix_t cov = matrix_t::Zero();
		for(int i = 0; i < 3; ++i)
		  for(int j = 0; j < 3; ++j)
		    {cov(i,j) = tr->get_error(i,j);
		      std::cout << std::setprecision(8) << " cov(" << i << "," << j << ") = " << cov(i,j) << std::endl;
		    }
		matrix_t cov_inv = cov.fullPivLu().inverse();

		std::cout << " Covariance: adding to new_pca_avge: " << cov_inv * B << std::endl;
		
		new_pca_avge += cov_inv * B;
		cov_wt += cov_inv;
	      }
	    else
	      {
		std::cout << " Equal weight: adding to new_pca_avge: " << B << std::endl;
		new_pca_avge += B;
		new_wt++;
	      }	    
	  }
	else
	  {
	    if(Verbosity() > 1) std::cout << "rejecting track " << trid << " from vertex " << vtxid << std::endl;
	  }
      }
      
      if(_use_track_covariance)
	{
	  if(Verbosity() > 1)  std::cout << " Covariance: new_pca_avge cum " << new_pca_avge << " cov_wt cum " << cov_wt << std::endl;
	  new_pca_avge = cov_wt.fullPivLu().inverse() * new_pca_avge;
	}
      else
	{
	  if(Verbosity() > 1)  std::cout << " Equal weight: new_pca_avge cum " << new_pca_avge << " new_wt " << new_wt << std::endl;
	  new_pca_avge = new_pca_avge / new_wt;
	}

      if(Verbosity() > 1)  
	{
	  std::cout << " pca_avge " << std::setprecision(8) << pca_avge.x() << "  " << pca_avge.y() << "  " << pca_avge.z() << std::endl; 
	  std::cout << " new_pca_avge " << new_pca_avge.x() << "  " << new_pca_avge.y() << "  " << new_pca_avge.z() << std::endl; 
	}

      // capture the position
      _revised_vertex_position_map.insert(std::make_pair(vtxid, new_pca_avge));
    }       

  // average covariance for accepted tracks
  for(auto it : _vertex_set)
    {
      matrix_t avgCov = matrix_t::Zero();
      double cov_wt = 0.0;

      auto ret = _revised_vertex_track_map.equal_range(it);
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
      _vertex_cov_map.insert(std::make_pair(it, avgCov));
    }
}

void PHSimpleVertexFinder::findOutlierTrackPairs()
{
  // We want to use the average PCA vector as a reference, and loop over the pair map to find outlying pairs

  for( auto it = _track_pair_pca_map.begin(); it != _track_pair_pca_map.end(); ++it)
    {
      //  Note: std::multimap<unsigned int, std::pair<unsigned int, std::pair<Eigen::Vector3d,  Eigen::Vector3d>>>  _track_pair_pca_map

      unsigned int tr1id = it->first;
      unsigned int tr2id = it->second.first;
      Eigen::Vector3d PCA1 = it->second.second.first;
      Eigen::Vector3d PCA2 = it->second.second.second;

      std::cout << " tr1id " << tr1id << " tr2id " << tr2id
		<< " PCA1 " << PCA1.x() << "  " << PCA1.y() << "  " << PCA1.z()
		<< " PCA2 " << PCA2.x() << "  " << PCA2.y() << "  " << PCA2.z()
		<< std::endl;


    }
  return;
}
