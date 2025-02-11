#include "PHSimpleVertexFinder.h"

/// Tracking includes

#include <trackbase/TrackVertexCrossingAssoc_v1.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>            // for cluskey, getLayer, TrkrId
#include <trackbase/TrackFitUtils.h>
#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap_v2.h>

#include <globalvertex/SvtxVertexMap_v1.h>
#include <globalvertex/SvtxVertex_v2.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <cmath>  // for sqrt, fabs, atan2, cos
#include <iomanip>
#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for map
#include <set>       // for _Rb_tree_const_iterator
#include <utility>   // for pair, make_pair

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>

//____________________________________________________________________________..
PHSimpleVertexFinder::PHSimpleVertexFinder(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int PHSimpleVertexFinder::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  return ret;
}

//____________________________________________________________________________..
int PHSimpleVertexFinder::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;
  }

  _active_dcacut = _base_dcacut;

  if (_vertex_track_map.size() > 0)
  {
    _svtx_vertex_map->clear();
  }

  // Write to a new map on the node tree that contains (crossing, trackid) pairs for all tracks
  // Later, will add to it a map  containing (crossing, vertexid)

  std::set<short int> crossings;
  for (const auto &[trackkey, track] : *_track_map)
  {
    auto crossing = track->get_crossing();
    auto siseed = track->get_silicon_seed();

    if (Verbosity() > 0)
      {
	std::cout << "track id " << trackkey << " crossing " << crossing 
		  << " x y z " << track->get_x() << "  " << track->get_y() << "  " << track->get_z() << std::endl;
      }

    // crossing zero contains unmatched TPC tracks
    // Here we skip those crossing = zero tracks that do not have silicon seeds
    if( (crossing == 0) & !siseed)
      {
	continue;
      }
    
    crossings.insert(crossing);
    _track_vertex_crossing_map->addTrackAssoc(crossing, trackkey);    
  }
  
  unsigned int vertex_id = 0;

  for (auto cross : crossings)
  {
    // reset maps for each crossing
    _vertex_track_map.clear();
    _track_pair_map.clear();
    _track_pair_pca_map.clear();
    _vertex_position_map.clear();
    _vertex_covariance_map.clear();
    _vertex_set.clear();

    if (Verbosity() > 0)
    {
      std::cout << "process tracks for beam crossing " << cross << std::endl;
    }

    // get the subset of tracks for this crossing
    auto crossing_track_index = _track_vertex_crossing_map->getTracks(cross);
    SvtxTrackMap *crossing_tracks = new SvtxTrackMap_v2;
    for (auto iter = crossing_track_index.first; iter != crossing_track_index.second; ++iter)
    {
      unsigned int trackkey = (*iter).second;
      SvtxTrack *track = _track_map->get(trackkey);
      if(!track)
      {
        continue;
      }
      crossing_tracks->insertWithKey(track, trackkey);
    }

    // Find all instances where two tracks have a dca of < _dcacut,  and capture the pair details
    // Fills _track_pair_map and _track_pair_pca_map
    if(_zero_field)
      {
	checkDCAsZF(crossing_tracks);
      }
    else
      {
	checkDCAs(crossing_tracks);
      }

    /// If we didn't find any matches, try again with a slightly larger DCA cut
    if (_track_pair_map.size() == 0)
    {
      _active_dcacut = 3.0 * _base_dcacut;
      if(_zero_field)
	{
	  checkDCAsZF(crossing_tracks);
	}
      else
	{
	  checkDCAs(crossing_tracks);
	}
    }
    
    if (Verbosity() > 0)
    {
      std::cout << "crossing " << cross << " track pair map size " << _track_pair_map.size() << std::endl;
    }

    // get all connected pairs of tracks by looping over the track_pair map
    std::vector<std::set<unsigned int>> connected_tracks = findConnectedTracks();

    // make vertices - each set of connected tracks is a vertex
    for (unsigned int ivtx = 0; ivtx < connected_tracks.size(); ++ivtx)
    {
      if (Verbosity() > 0)
      {
        std::cout << "process vertex " << ivtx + vertex_id << std::endl;
      }

      for (auto it : connected_tracks[ivtx])
      {
        unsigned int id = it;
        _vertex_track_map.insert(std::make_pair(ivtx, id));
        if (Verbosity() > 0)
        {
          std::cout << "  adding track " << id << " to vertex " << ivtx + vertex_id << std::endl;
        }
      }
    }

    // make a list of vertices
    for (auto it : _vertex_track_map)
    {
      if (Verbosity() > 1)
      {
        std::cout << " vertex " << it.first + vertex_id << " track " << it.second << std::endl;
      }
      _vertex_set.insert(it.first);
    }

    // this finds average vertex positions after removal of outlying track pairs
    removeOutlierTrackPairs();

    // average covariance for accepted tracks
    for (auto it : _vertex_set)
    {
      matrix_t avgCov = matrix_t::Zero();
      double cov_wt = 0.0;

      auto ret = _vertex_track_map.equal_range(it);
      for (auto cit = ret.first; cit != ret.second; ++cit)
      {
        unsigned int trid = cit->second;
        matrix_t cov;
        auto track = _track_map->get(trid);
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            cov(i, j) = track->get_error(i, j);
          }
        }

        avgCov += cov;
        cov_wt++;
      }

      avgCov /= sqrt(cov_wt);
      if (Verbosity() > 2)
      {
        std::cout << "Average covariance for vertex " << it << " is:" << std::endl;
        std::cout << std::setprecision(8) << avgCov << std::endl;
      }
      _vertex_covariance_map.insert(std::make_pair(it, avgCov));
    }

    // Write the vertices to the vertex map on the node tree
    //==============================================

    for (auto it : _vertex_set)
    {
      unsigned int thisid = it + vertex_id;  // the address of the vertex in the event

      auto svtxVertex = std::make_unique<SvtxVertex_v2>();

      svtxVertex->set_chisq(0.0);
      svtxVertex->set_ndof(0);
      svtxVertex->set_t0(0);
      svtxVertex->set_id(thisid);
      svtxVertex->set_beam_crossing(cross);

      auto ret = _vertex_track_map.equal_range(it);
      for (auto cit = ret.first; cit != ret.second; ++cit)
      {
        unsigned int trid = cit->second;

        if (Verbosity() > 1)
        {
          std::cout << "   vertex " << thisid << " insert track " << trid << std::endl;
        }
        svtxVertex->insert_track(trid);
        _track_map->get(trid)->set_vertex_id(thisid);
      }

      Eigen::Vector3d pos = _vertex_position_map.find(it)->second;
      svtxVertex->set_x(pos.x());
      svtxVertex->set_y(pos.y());
      svtxVertex->set_z(pos.z());
      if (Verbosity() > 1)
      {
        std::cout << "   vertex " << thisid << " insert pos.x " << pos.x() << " pos.y " << pos.y() << " pos.z " << pos.z() << std::endl;
      }

      auto vtxCov = _vertex_covariance_map.find(it)->second;
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          svtxVertex->set_error(i, j, vtxCov(i, j));
        }
      }

      _svtx_vertex_map->insert(svtxVertex.release());
    }

    vertex_id += _vertex_set.size();

    /// Iterate through the tracks and assign the closest vtx id to
    /// the track position for propagating back to the vtx. Catches any
    /// tracks that were missed or were not  compatible with any of the
    /// identified vertices
    //=================================================
    for (const auto &[trackkey, track] : *crossing_tracks)
    {
      auto thistrack = _track_map->get(trackkey);  // get the original, not the copy
      auto vtxid = thistrack->get_vertex_id();
      if (Verbosity() > 1)
      {
        std::cout << "        track " << trackkey << " track vtxid " << vtxid << std::endl;
      }

      /// If there is a vertex already assigned, keep going
      if (vtxid != std::numeric_limits<unsigned int>::max())
      {
        continue;
      }

      float maxdz = std::numeric_limits<float>::max();
      unsigned int newvtxid = std::numeric_limits<unsigned int>::max();

      for (auto it : _vertex_set)
      {
        unsigned int thisid = it + vertex_id - _vertex_set.size();

        if (Verbosity() > 1)
        {
          std::cout << "                test vertex " << thisid << std::endl;
        }

        auto thisvertex = _svtx_vertex_map->get(thisid);
        float dz = thistrack->get_z() - thisvertex->get_z();
        if (std::fabs(dz) < maxdz)
        {
          maxdz = dz;
          newvtxid = thisid;
        }
      }

      // this updates the track, but does not add it to the vertex primary track list
      if (newvtxid != std::numeric_limits<unsigned int>::max())
      {
        thistrack->set_vertex_id(newvtxid);
        if (Verbosity() > 1)
        {
          std::cout << "                assign vertex " << newvtxid << " to additional track " << trackkey << std::endl;
        }
      }
    }

    delete crossing_tracks;

  }  // end loop over crossings

  // update the crossing vertex map with the results
  for (const auto &[vtxkey, vertex] : *_svtx_vertex_map)
  {
    short int crossing = vertex->get_beam_crossing();
    _track_vertex_crossing_map->addVertexAssoc(crossing, vtxkey);

    if (Verbosity() > 1)
    {
      std::cout << "Vertex ID: " << vtxkey << " vertex crossing " << crossing << " list of tracks: " << std::endl;
      for (auto trackiter = vertex->begin_tracks(); trackiter != vertex->end_tracks(); ++trackiter)
      {
        SvtxTrack *track = _track_map->get(*trackiter);
        if (!track)
        {
          continue;
        }

        auto siseed = track->get_silicon_seed();
        short int intt_crossing = siseed->get_crossing();

        // the track crossing may be from the INTT clusters or from geometric matching if there are no INTT clusters
        short int track_crossing = track->get_crossing();
        std::cout << " vtxid " << vtxkey << " vertex crossing " << crossing
                  << " track crossing " << track_crossing
                  << " intt crossing " << intt_crossing
                  << " trackID " << *trackiter
                  << " track Z " << track->get_z()
                  << " X " << track->get_x()
                  << " Y " << track->get_y()
                  << " quality " << track->get_quality()
                  << " pt " << track->get_pt()
                  << std::endl;
        if (Verbosity() > 3)
        {
          siseed->identify();
        }
      }
    }
  }

  if (Verbosity() > 2)
  {
    _track_vertex_crossing_map->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleVertexFinder::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleVertexFinder::CreateNodes(PHCompositeNode *topNode)
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

  _svtx_vertex_map = findNode::getClass<SvtxVertexMap>(topNode, _vertex_map_name);

  if (!_svtx_vertex_map)
  {
    _svtx_vertex_map = new SvtxVertexMap_v1;
    PHIODataNode<PHObject> *vertexNode = new PHIODataNode<PHObject>(
        _svtx_vertex_map, _vertex_map_name, "PHObject");

    svtxNode->addNode(vertexNode);
  }

  _track_vertex_crossing_map = findNode::getClass<TrackVertexCrossingAssoc>(topNode, "TrackVertexCrossingAssocMap");
  if (!_track_vertex_crossing_map)
  {
    _track_vertex_crossing_map = new TrackVertexCrossingAssoc_v1;
    PHIODataNode<PHObject> *trackvertexNode = new PHIODataNode<PHObject>(
        _track_vertex_crossing_map, "TrackVertexCrossingAssocMap", "PHObject");

    svtxNode->addNode(trackvertexNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHSimpleVertexFinder::GetNodes(PHCompositeNode *topNode)
{
  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << _track_map_name << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER container" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tGeometry)
  {
    std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHSimpleVertexFinder::checkDCAs(SvtxTrackMap *track_map)
{
  // Loop over tracks and check for close DCA match with all other tracks
  for (auto tr1_it = track_map->begin(); tr1_it != track_map->end(); ++tr1_it)
  {
    auto id1 = tr1_it->first;
    auto tr1 = tr1_it->second;
    if (tr1->get_quality() > _qual_cut)
    {
      continue;
    }
    if (_require_mvtx)
    {
      unsigned int nmvtx = 0;
      TrackSeed *siliconseed = tr1->get_silicon_seed();
      if (!siliconseed)
      {
        continue;
      }

      for (auto clusit = siliconseed->begin_cluster_keys(); clusit != siliconseed->end_cluster_keys(); ++clusit)
      {
        if (TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId)
        {
          nmvtx++;
        }
        if (nmvtx >= _nmvtx_required)
        {
          break;
        }
      }
      if (nmvtx < _nmvtx_required)
      {
        continue;
      }
      if (Verbosity() > 3)
      {
        std::cout << " tr1 id " << id1 << " has nmvtx at least " << nmvtx << std::endl;
      }
    }

    // look for close DCA matches with all other such tracks
    for (auto tr2_it = std::next(tr1_it); tr2_it != track_map->end(); ++tr2_it)
    {
      auto id2 = tr2_it->first;
      auto tr2 = tr2_it->second;
      if (tr2->get_quality() > _qual_cut)
      {
        continue;
      }
      if (_require_mvtx)
      {
        unsigned int nmvtx = 0;
        TrackSeed *siliconseed = tr2->get_silicon_seed();
        if (!siliconseed)
        {
          continue;
        }

        for (auto clusit = siliconseed->begin_cluster_keys(); clusit != siliconseed->end_cluster_keys(); ++clusit)
        {
          if (TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId)
          {
            nmvtx++;
          }
          if (nmvtx >= _nmvtx_required)
          {
            break;
          }
        }
        if (nmvtx < _nmvtx_required)
        {
          continue;
        }
        if (Verbosity() > 3)
        {
          std::cout << " tr2 id " << id2 << " has nmvtx at least " << nmvtx << std::endl;
        }
      }

      // find DCA of these two tracks
      if (Verbosity() > 3)
      {
        std::cout << "Check DCA for tracks " << id1 << " and  " << id2 << std::endl;
      }

      findDcaTwoTracks(tr1, tr2);
    }
  }
}

void PHSimpleVertexFinder::checkDCAsZF(SvtxTrackMap *track_map)
{
  // ZF tracks do not have an Acts fit, and the seeding does not give
  // reliable track parameters - refit clusters with straight lines
  // No distortion corrections applied in TPC at present

  std::vector<unsigned int> cumulative_trackid_vec;
  std::vector<unsigned int> cumulative_nmvtx_vec;
  std::vector<unsigned int> cumulative_nintt_vec;
  std::vector<std::vector<Acts::Vector3>> cumulative_global_vec;
  std::vector<std::vector<TrkrDefs::cluskey>> cumulative_cluskey_vec;
  std::vector<std::vector<float>> cumulative_fitpars_vec;

  for (auto & tr1_it : *track_map)
  {
    auto id1 = tr1_it.first;
    auto tr1 = tr1_it.second;

    //    tr1->identify();
 
    TrackSeed *siliconseed = tr1->get_silicon_seed();
    if (_require_mvtx)
      {
	if (!siliconseed)
	  {
	    continue;
	  }
      }
    TrackSeed *tpcseed = tr1->get_tpc_seed();

    std::vector<Acts::Vector3> global_vec;
    std::vector<TrkrDefs::cluskey> cluskey_vec;
    
    // Get a vector of cluster keys from the silicon seed and TPC seed
    if(siliconseed)
      {
	getTrackletClusterList(siliconseed, cluskey_vec);
	if(Verbosity() > 0) 
	  { 
	    std::cout << "  after silicon: silicon cluskey_vec size " << cluskey_vec.size() << std::endl; 
	    for(unsigned long i : cluskey_vec)
	      {  
		std::cout << i << std::endl;
	      }
	  }
      }
    if(tpcseed)
      {

	getTrackletClusterList(tpcseed, cluskey_vec);
	if(Verbosity() > 0) 
	  { 
	    std::cout << "  after tpc: cluskey_vec size " << cluskey_vec.size() << std::endl; 
	    for(unsigned long i : cluskey_vec)
	      {  
		std::cout << i << std::endl;
	      }
	  }
      }

    unsigned int nmvtx = 0;
    unsigned int nintt = 0;
    for (auto& key : cluskey_vec)
      {
        if (TrkrDefs::getTrkrId(key) == TrkrDefs::mvtxId)
	  {
	    nmvtx++;
	  }
	if(TrkrDefs::getTrkrId(key) == TrkrDefs::inttId)
	  {
	    nintt++;
	  }
      }
    
    // store cluster global positions in a vector
    TrackFitUtils::getTrackletClusters(_tGeometry, _cluster_map, global_vec, cluskey_vec);
    
    std::vector<float> fitpars = TrackFitUtils::fitClustersZeroField(global_vec, cluskey_vec, true);
    
    if (Verbosity() > 1)
      {
	if(fitpars.size() == 4)
	  {
	    std::cout << " Track " << id1 << " dy/dx " << fitpars[0] << " y intercept " << fitpars[1] 
		      << " dx/dz " << fitpars[2] << " Z0 " << fitpars[3] << std::endl;
	  }
	else
	  {
	    std::cout << " Track " << id1 << " ZF line fits failed, fitpars is empty" << std::endl;
	  }
      }

    cumulative_trackid_vec.push_back(id1);   
    cumulative_nmvtx_vec.push_back(nmvtx);   
    cumulative_nintt_vec.push_back(nintt);   
    cumulative_global_vec.push_back(global_vec);
    cumulative_cluskey_vec.push_back(cluskey_vec);
    cumulative_fitpars_vec.push_back(fitpars);
  }

  for(unsigned int i1 = 0; i1 < cumulative_trackid_vec.size(); ++i1)
    {
      if(cumulative_fitpars_vec[i1].size() == 0) { continue; }

      for(unsigned int i2 = i1; i2 < cumulative_trackid_vec.size(); ++i2)
	{
	  if(cumulative_fitpars_vec[i2].size() == 0) { continue; }

	  //  For straight line: fitpars[4] = { xyslope, y0, xzslope, z0 }
	  Eigen::Vector3d a1(0.0, cumulative_fitpars_vec[i1][1],cumulative_fitpars_vec[i1][3]);      // point on track 1 at x = 0
	  Eigen::Vector3d a2(0.0, cumulative_fitpars_vec[i2][1],cumulative_fitpars_vec[i2][3]);      // point on track 2 at x = 0
	  // direction vectors made from dy/dx = xyslope and dz/dx = xzslope
	  Eigen::Vector3d b1(1.0, cumulative_fitpars_vec[i1][0],cumulative_fitpars_vec[i1][2]);      // direction vector of track 1
	  Eigen::Vector3d b2(1.0, cumulative_fitpars_vec[i2][0],cumulative_fitpars_vec[i2][2]);      // direction vector of track 2
	  	  
	  Eigen::Vector3d PCA1(0, 0, 0);
	  Eigen::Vector3d PCA2(0, 0, 0);
	  double dca = dcaTwoLines(a1, b1, a2, b2, PCA1, PCA2);


	  // check dca cut is satisfied, and that PCA is close to beam line
	  if (fabs(dca) < _active_dcacut && (fabs(PCA1.x()) < _beamline_xy_cut && fabs(PCA1.y()) < _beamline_xy_cut))
	    {
	      int id1 = cumulative_trackid_vec[i1];
	      int id2 = cumulative_trackid_vec[i2];

	      if (Verbosity() > 3)
		{
		  std::cout << " good match for tracks " << id1 << " and " << id2 << std::endl;
		  std::cout << "    a1.x " << a1.x() << " a1.y " << a1.y() << " a1.z " << a1.z() << std::endl;
		  std::cout << "    a2.x  " << a2.x() << " a2.y " << a2.y() << " a2.z " << a2.z() << std::endl;
		  std::cout << "    PCA1.x() " << PCA1.x() << " PCA1.y " << PCA1.y() << " PCA1.z " << PCA1.z() << std::endl;
		  std::cout << "    PCA2.x() " << PCA2.x() << " PCA2.y " << PCA2.y() << " PCA2.z " << PCA2.z() << std::endl;
		  std::cout << "    dca " << dca << std::endl;
		}
	      
	      // capture the results for successful matches
	      _track_pair_map.insert(std::make_pair(id1, std::make_pair(id2, dca)));
	      _track_pair_pca_map.insert(std::make_pair(id1, std::make_pair(id2, std::make_pair(PCA1, PCA2))));
	    }	  	  
	}
    }

  return; 
}

void PHSimpleVertexFinder::getTrackletClusterList(TrackSeed* tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  for (auto clusIter = tracklet->begin_cluster_keys();
       clusIter != tracklet->end_cluster_keys();
       ++clusIter)
  {
    auto key = *clusIter;
    auto cluster = _cluster_map->findCluster(key);
    if (!cluster)
    {
      std::cout << PHWHERE << "Failed to get cluster with key " << key << std::endl;
      continue;
    }

    /// Make a safety check for clusters that couldn't be attached to a surface
    auto surf = _tGeometry->maps().getSurface(key, cluster);
    if (!surf)
    {
      continue;
    }

    // drop some bad layers in the TPC completely
    unsigned int layer = TrkrDefs::getLayer(key);
    if (layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39)
    {
      continue;
    }

    // drop INTT clusters for now  -- TEMPORARY!
    if (layer > 2 && layer < 7)
    {
      continue;
    }


    cluskey_vec.push_back(key);

  }  // end loop over clusters for this track
}

void PHSimpleVertexFinder::checkDCAs()
{
  // Loop over tracks and check for close DCA match with all other tracks
  for (auto tr1_it = _track_map->begin(); tr1_it != _track_map->end(); ++tr1_it)
  {
    auto id1 = tr1_it->first;
    auto tr1 = tr1_it->second;
    if (tr1->get_quality() > _qual_cut)
    {
      continue;
    }
    if (_require_mvtx)
    {
      unsigned int nmvtx = 0;
      TrackSeed *siliconseed = tr1->get_silicon_seed();
      if (!siliconseed)
      {
        continue;
      }

      for (auto clusit = siliconseed->begin_cluster_keys(); clusit != siliconseed->end_cluster_keys(); ++clusit)
      {
        if (TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId)
        {
          nmvtx++;
        }
        if (nmvtx >= _nmvtx_required)
        {
          break;
        }
      }
      if (nmvtx < _nmvtx_required)
      {
        continue;
      }
      if (Verbosity() > 3)
      {
        std::cout << " tr1 id " << id1 << " has nmvtx at least " << nmvtx << std::endl;
      }
    }

    // look for close DCA matches with all other such tracks
    for (auto tr2_it = std::next(tr1_it); tr2_it != _track_map->end(); ++tr2_it)
    {
      auto id2 = tr2_it->first;
      auto tr2 = tr2_it->second;
      if (tr2->get_quality() > _qual_cut)
      {
        continue;
      }
      if (_require_mvtx)
      {
        unsigned int nmvtx = 0;
        TrackSeed *siliconseed = tr2->get_silicon_seed();
        if (!siliconseed)
        {
          continue;
        }

        for (auto clusit = siliconseed->begin_cluster_keys(); clusit != siliconseed->end_cluster_keys(); ++clusit)
        {
          if (TrkrDefs::getTrkrId(*clusit) == TrkrDefs::mvtxId)
          {
            nmvtx++;
          }
          if (nmvtx >= _nmvtx_required)
          {
            break;
          }
        }
        if (nmvtx < _nmvtx_required)
        {
          continue;
        }
        if (Verbosity() > 3)
        {
          std::cout << " tr2 id " << id2 << " has nmvtx at least " << nmvtx << std::endl;
        }
      }

      // find DCA of these two tracks
      if (Verbosity() > 3)
      {
        std::cout << "Check DCA for tracks " << id1 << " and  " << id2 << std::endl;
      }

      findDcaTwoTracks(tr1, tr2);
    }
  }
}

void PHSimpleVertexFinder::findDcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2)
{
  if (tr1->get_pt() < _track_pt_cut)
  {
    return;
  }
  if (tr2->get_pt() < _track_pt_cut)
  {
    return;
  }

  unsigned int id1 = tr1->get_id();
  unsigned int id2 = tr2->get_id();

  // get the line equations for the tracks

  Eigen::Vector3d a1(tr1->get_x(), tr1->get_y(), tr1->get_z());
  Eigen::Vector3d b1(tr1->get_px() / tr1->get_p(), tr1->get_py() / tr1->get_p(), tr1->get_pz() / tr1->get_p());
  Eigen::Vector3d a2(tr2->get_x(), tr2->get_y(), tr2->get_z());
  Eigen::Vector3d b2(tr2->get_px() / tr2->get_p(), tr2->get_py() / tr2->get_p(), tr2->get_pz() / tr2->get_p());

  Eigen::Vector3d PCA1(0, 0, 0);
  Eigen::Vector3d PCA2(0, 0, 0);
  double dca = dcaTwoLines(a1, b1, a2, b2, PCA1, PCA2);

  if (Verbosity() > 3)
  {
    std::cout << " pair dca is " << dca << " _active_dcacut is " << _active_dcacut
              << " PCA1.x " << PCA1.x() << " PCA1.y " << PCA1.y()
              << " PCA2.x " << PCA2.x() << " PCA2.y " << PCA2.y() << std::endl;
  }

  // check dca cut is satisfied, and that PCA is close to beam line
  if (fabs(dca) < _active_dcacut && (fabs(PCA1.x()) < _beamline_xy_cut && fabs(PCA1.y()) < _beamline_xy_cut))
  {
    if (Verbosity() > 3)
    {
      std::cout << " good match for tracks " << tr1->get_id() << " and " << tr2->get_id() << " with pT " << tr1->get_pt() << " and " << tr2->get_pt() << std::endl;
      std::cout << "    a1.x " << a1.x() << " a1.y " << a1.y() << " a1.z " << a1.z() << std::endl;
      std::cout << "    a2.x  " << a2.x() << " a2.y " << a2.y() << " a2.z " << a2.z() << std::endl;
      std::cout << "    PCA1.x() " << PCA1.x() << " PCA1.y " << PCA1.y() << " PCA1.z " << PCA1.z() << std::endl;
      std::cout << "    PCA2.x() " << PCA2.x() << " PCA2.y " << PCA2.y() << " PCA2.z " << PCA2.z() << std::endl;
      std::cout << "    dca " << dca << std::endl;
    }

    // capture the results for successful matches
    _track_pair_map.insert(std::make_pair(id1, std::make_pair(id2, dca)));
    _track_pair_pca_map.insert(std::make_pair(id1, std::make_pair(id2, std::make_pair(PCA1, PCA2))));
  }

  return;
}

double PHSimpleVertexFinder::dcaTwoLines(const Eigen::Vector3d &a1, const Eigen::Vector3d &b1,
                                         const Eigen::Vector3d &a2, const Eigen::Vector3d &b2,
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
  auto aminusa = a2 - a1;

  // The DCA of these two lines is the projection of a2-a1 along the direction of the perpendicular to both
  // remember that a2-a1 is longer than (or equal to) the dca by definition
  double dca = 999;
  if (mag_bcrossb != 0)
  {
    dca = bcrossb.dot(aminusa) / mag_bcrossb;
  }
  else
  {
    return dca;  // same track, skip combination
  }

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

  double X = b1.dot(b2) - b1.dot(b1) * b2.dot(b2) / b2.dot(b1);
  double Y = (a2.dot(b2) - a1.dot(b2)) - (a2.dot(b1) - a1.dot(b1)) * b2.dot(b2) / b2.dot(b1);
  double c = Y / X;

  double F = b1.dot(b1) / b2.dot(b1);
  double G = -(a2.dot(b1) - a1.dot(b1)) / b2.dot(b1);
  double d = c * F + G;

  // then the points of closest approach are:
  PCA1 = a1 + c * b1;
  PCA2 = a2 + d * b2;

  return dca;
}

std::vector<std::set<unsigned int>> PHSimpleVertexFinder::findConnectedTracks()
{
  std::vector<std::set<unsigned int>> connected_tracks;
  std::set<unsigned int> connected;
  std::set<unsigned int> used;
  for (auto it : _track_pair_map)
  {
    unsigned int id1 = it.first;
    unsigned int id2 = it.second.first;

    if ((used.find(id1) != used.end()) && (used.find(id2) != used.end()))
    {
      if (Verbosity() > 3)
      {
        std::cout << " tracks " << id1 << " and " << id2 << " are both in used , skip them" << std::endl;
      }
      continue;
    }
    else if ((used.find(id1) == used.end()) && (used.find(id2) == used.end()))
    {
      if (Verbosity() > 3)
      {
        std::cout << " tracks " << id1 << " and " << id2 << " are both not in used , start a new connected set" << std::endl;
      }
      // close out and start a new connections set
      if (connected.size() > 0)
      {
        connected_tracks.push_back(connected);
        connected.clear();
        if (Verbosity() > 3)
        {
          std::cout << "           closing out set " << std::endl;
        }
      }
    }

    // get everything connected to id1 and id2
    connected.insert(id1);
    used.insert(id1);
    connected.insert(id2);
    used.insert(id2);
    for (auto cit : _track_pair_map)
    {
      unsigned int id3 = cit.first;
      unsigned int id4 = cit.second.first;
      if ((connected.find(id3) != connected.end()) || (connected.find(id4) != connected.end()))
      {
        if (Verbosity() > 3)
        {
          std::cout << " found connection to " << id3 << " and " << id4 << std::endl;
        }
        connected.insert(id3);
        used.insert(id3);
        connected.insert(id4);
        used.insert(id4);
      }
    }
  }

  // close out the last set
  if (connected.size() > 0)
  {
    connected_tracks.push_back(connected);
    connected.clear();
    if (Verbosity() > 3)
    {
      std::cout << "           closing out last set " << std::endl;
    }
  }

  if (Verbosity() > 3)
  {
    std::cout << "connected_tracks size " << connected_tracks.size() << std::endl;
  }

  return connected_tracks;
}

void PHSimpleVertexFinder::removeOutlierTrackPairs()
{
  //  Note: std::multimap<unsigned int, std::pair<unsigned int, std::pair<Eigen::Vector3d,  Eigen::Vector3d>>>  _track_pair_pca_map

  for (auto it : _vertex_set)
  {
    unsigned int vtxid = it;
    if (Verbosity() > 1)
    {
      std::cout << "calculate average position for vertex " << vtxid << std::endl;
    }

    // we need the median values of the x and y positions
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;

    double pca_median_x = 0.;
    double pca_median_y = 0.;
    double pca_median_z = 0.;

    Eigen::Vector3d new_pca_avge(0., 0., 0.);
    double new_wt = 0.0;

    auto ret = _vertex_track_map.equal_range(vtxid);

    // Start by getting the positions for this vertex into vectors for the median calculation
    for (auto cit = ret.first; cit != ret.second; ++cit)
    {
      unsigned int tr1id = cit->second;
      if (Verbosity() > 2)
      {
        std::cout << "   vectors: get entries for track " << tr1id << " for vertex " << vtxid << std::endl;
      }

      // find all pairs for this vertex with tr1id
      auto pca_range = _track_pair_pca_map.equal_range(tr1id);
      for (auto pit = pca_range.first; pit != pca_range.second; ++pit)
      {
        unsigned int tr2id = pit->second.first;

        Eigen::Vector3d PCA1 = pit->second.second.first;
        Eigen::Vector3d PCA2 = pit->second.second.second;

        if (Verbosity() > 2)
        {
          std::cout << " vectors: tr1id " << tr1id << " tr2id " << tr2id
                    << " PCA1 " << PCA1.x() << "  " << PCA1.y() << "  " << PCA1.z()
                    << " PCA2 " << PCA2.x() << "  " << PCA2.y() << "  " << PCA2.z()
                    << std::endl;
        }

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
    if (vx.size() < 3)
    {
      new_pca_avge.x() = getAverage(vx);
      new_pca_avge.y() = getAverage(vy);
      new_pca_avge.z() = getAverage(vz);
      _vertex_position_map.insert(std::make_pair(vtxid, new_pca_avge));
      if (Verbosity() > 1)
      {
        std::cout << " Vertex has only 2 tracks, use average for PCA: " << new_pca_avge.x() << "  " << new_pca_avge.y() << "  " << new_pca_avge.z() << std::endl;
      }

      // done with this vertex
      continue;
    }

    pca_median_x = getMedian(vx);
    pca_median_y = getMedian(vy);
    pca_median_z = getMedian(vz);
    if (Verbosity() > 1)
    {
      std::cout << "Median values: x " << pca_median_x << " y " << pca_median_y << " z : " << pca_median_z << std::endl;
    }

    // Make the average vertex position with outlier rejection wrt the median
    for (auto cit = ret.first; cit != ret.second; ++cit)
    {
      unsigned int tr1id = cit->second;
      if (Verbosity() > 2)
      {
        std::cout << "   average: get entries for track " << tr1id << " for vertex " << vtxid << std::endl;
      }

      // find all pairs for this vertex with tr1id
      auto pca_range = _track_pair_pca_map.equal_range(tr1id);
      for (auto pit = pca_range.first; pit != pca_range.second; ++pit)
      {
        unsigned int tr2id = pit->second.first;

        Eigen::Vector3d PCA1 = pit->second.second.first;
        Eigen::Vector3d PCA2 = pit->second.second.second;

        if (
            fabs(PCA1.x() - pca_median_x) < _outlier_cut &&
            fabs(PCA1.y() - pca_median_y) < _outlier_cut &&
            fabs(PCA2.x() - pca_median_x) < _outlier_cut &&
            fabs(PCA2.y() - pca_median_y) < _outlier_cut)
        {
          // good track pair, add to new average

          new_pca_avge += PCA1;
          new_wt++;
          new_pca_avge += PCA2;
          new_wt++;
        }
        else
        {
          if (Verbosity() > 1)
          {
            std::cout << "Reject pair with tr1id " << tr1id << " tr2id " << tr2id << std::endl;
          }
        }
      }
    }
    if (new_wt > 0.0)
    {
      new_pca_avge = new_pca_avge / new_wt;
    }
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

  if ((v.size() % 2) == 0)
  {
    // even number of entries
    // we want the average of the middle two numbers, v.size()/2 and v.size()/2-1
    auto m1 = v.begin() + v.size() / 2;
    std::nth_element(v.begin(), m1, v.end());
    double median1 = v[v.size() / 2];

    auto m2 = v.begin() + v.size() / 2 - 1;
    std::nth_element(v.begin(), m2, v.end());
    double median2 = v[v.size() / 2 - 1];

    median = (median1 + median2) / 2.0;
    if (Verbosity() > 2)
    {
      std::cout << "The vector size is " << v.size()
                << " element m is " << v.size() / 2 << " = " << v[v.size() / 2]
                << " element m-1 is " << v.size() / 2 - 1 << " = " << v[v.size() / 2 - 1]
                << std::endl;
    }
  }
  else
  {
    // odd number of entries
    auto m = v.begin() + v.size() / 2;
    std::nth_element(v.begin(), m, v.end());
    median = v[v.size() / 2];
    if (Verbosity() > 2)
    {
      std::cout << "The vector size is " << v.size() << " element m is " << v.size() / 2 << " = " << v[v.size() / 2] << std::endl;
    }
  }

  return median;
}
double PHSimpleVertexFinder::getAverage(std::vector<double> &v)
{
  double avge = 0.0;
  double wt = 0.0;
  for (auto it : v)
  {
    avge += it;
    wt++;
  }

  avge /= wt;
  if (Verbosity() > 2)
  {
    std::cout << " average = " << avge << std::endl;
  }

  return avge;
}
