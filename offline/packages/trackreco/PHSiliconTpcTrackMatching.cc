#include "PHSiliconTpcTrackMatching.h"

#include "AssocInfoContainer.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TpcSeedTrackMapv1.h>     
#include <trackbase/TrkrClusterContainerv3.h>   
#include <trackbase/TrkrClusterv3.h>   
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/ActsTransformations.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>


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
PHSiliconTpcTrackMatching::PHSiliconTpcTrackMatching(const std::string &name):
  SubsysReco(name)
 , _track_map_name_silicon("SvtxSiliconTrackMap")
{
  //cout << "PHSiliconTpcTrackMatching::PHSiliconTpcTrackMatching(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHSiliconTpcTrackMatching::~PHSiliconTpcTrackMatching()
{

}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::InitRun(PHCompositeNode *topNode)
{
  // put these in the output file
  cout << PHWHERE "_is_ca_seeder " << _is_ca_seeder << " Search windows: phi " << _phi_search_win << " eta " 
       << _eta_search_win << " pp_mode " << _pp_mode << endl;

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
 
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::process_event(PHCompositeNode*)
{
  // _track_map contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // We will add the silicon clusters to the TPC tracks already on the node tree
  // We will have to expand the number of tracks whenever we find multiple matches to the silicon

  _seed_track_map->Reset();
  _z_mismatch_map.clear();

  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size() << " Silicon track map size " << _track_map_silicon->size() << endl;

  if(_track_map->size() == 0)
    return Fun4AllReturnCodes::EVENT_OK;
 
  if(Verbosity() > 2)
    {
      // list silicon tracks
      for (auto phtrk_iter_si = _track_map_silicon->begin();
	   phtrk_iter_si != _track_map_silicon->end(); 
	   ++phtrk_iter_si)
	{
	  _tracklet_si = phtrk_iter_si->second;	  
	  
	  double si_phi = atan2(_tracklet_si->get_py(), _tracklet_si->get_px());
	  double si_eta = _tracklet_si->get_eta();
	  
	  cout << " Si track " << _tracklet_si->get_id()  << " si_phi " << si_phi  << " si_eta " << si_eta << endl;
	}  
    }
  
  // Find all matches of tpc and si tracklets in eta and phi, x and y
  //     If _pp_mode is not set, a match in z is also required - gives same behavior as old code
  // In any case, multiple matches are handled by duplicating the tpc tracklet into a new track
  //     "_seed_track_map" records (original id,duplicate id) so that the track cleaner can choose the one with the best Acts fit later
  std::multimap<unsigned int, unsigned int> tpc_matches;
  std::set<unsigned int> tpc_matched_set;
  findEtaPhiMatches(tpc_matched_set, tpc_matches);
  
  // We have a complete list of all eta/phi matched tracks in the map "tpc_matches"
  // In all cases, the crossing value is only rough at this point (it is only a tag)
  std::multimap<int, std::pair<unsigned int, unsigned int>> crossing_matches;
  std::map<unsigned int, int> tpc_crossing_map;
  std::set<int> crossing_set;
  if(_pp_mode)
    {
      // This section is to correct the TPC z positions of tracks for all bunch crossings.  
      //=================================================================================

      // All tracks treated as if we do not know the bunch crossing
      // The crossing estimate here is crude, for now
      tagMatchCrossing(tpc_matches, crossing_set, crossing_matches, tpc_crossing_map);

      // Sort candidates by the silicon tracklet Z position by putting them in si_sorted map
      //      -- captures crossing, tpc_id and si_id
      std::multimap<double, std::pair<unsigned int, unsigned int>> si_sorted_map;
      for( auto ncross : crossing_set)
	{
	  if(Verbosity() > 1) std::cout << " ncross = " << ncross << std::endl;
	  
	  auto ret  = crossing_matches. equal_range(ncross);            
	  for(auto it = ret.first; it != ret.second; ++it)
	    {
	      if(Verbosity() > 1) 
		std::cout << "  crossing " << it->first << " tpc_id " << it->second.first << " si_id " << it->second.second 
			  << " pT " << _track_map->get(it->second.first)->get_pt() 
			  << " eta " << _track_map->get(it->second.first)->get_eta() 
			  << " tpc_z " << _track_map->get(it->second.first)->get_z() 
			  << " si_z " << _track_map_silicon->get(it->second.second)->get_z() 
			  << std::endl;
	      
	      double z_si = _track_map_silicon->get(it->second.second)->get_z();
	      si_sorted_map.insert(std::make_pair(z_si, std::make_pair(it->second.first, it->second.second)));
	    }	
	}

      // make a list of silicon vertices, and a multimap with all associated tpc/si pairs 
      std::multimap<unsigned int, std::pair<unsigned int, unsigned int>> vertex_map;
      std::vector<double> vertex_list;
      getSiVertexList(si_sorted_map, vertex_list, vertex_map);
      
      // find the crossing number for every track from the mismatch with the associated silicon vertex z position
      std::map<unsigned int, double> vertex_crossings_map;
      getCrossingNumber(vertex_list, vertex_map, vertex_crossings_map);

      // Remove tracks from vertex_map where the vertex crossing is badly inconsistent with the initial crossing estimate
      cleanVertexMap( vertex_crossings_map, vertex_map, tpc_crossing_map );
      
      // correct the TPC cluster z values for the bunch crossing offset
      correctTpcClusterZ(vertex_crossings_map, vertex_map);

      // add silicon clusters to the surviving tracks
      addSiliconClusters(vertex_map);

      //===================================================
    }
  else
    {
      // only crossing zero has been added to the tpc_matches map, just add silicon clusters
      addSiliconClusters(tpc_matches);
    }

  // loop over all tracks and copy the silicon clusters to the corrected cluster map
  if(_corrected_cluster_map)
    copySiliconClustersToCorrectedMap();
  
  if(Verbosity() > 0)  
    cout << " Final track map size " << _track_map->size() 
	 << " seed-track map size " << _seed_track_map->size() << endl;
  
  if (Verbosity() > 0)
    cout << "PHSiliconTpcTrackMatching::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
  return Fun4AllReturnCodes::EVENT_OK;
 }

int PHSiliconTpcTrackMatching::End(PHCompositeNode* )
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHSiliconTpcTrackMatching::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

   _assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
  if (!_assoc_container)
  {
    cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_silicon = findNode::getClass<SvtxTrackMap>(topNode, _silicon_track_map_name);
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxSiliconTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

   _seed_track_map  = findNode::getClass<TpcSeedTrackMap>(topNode, _tpcseed_track_map_name);
  if(!_seed_track_map)
    {
      std::cout << "Creating node TpcSeedTrackMap" << std::endl;

      /// Get the DST Node
      PHNodeIterator iter(topNode);
      PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  
      /// Check that it is there
      if (!dstNode)
	{
	  std::cerr << "DST Node missing, quitting" << std::endl;
	  throw std::runtime_error("failed to find DST node in PHActsSourceLinks::createNodes");
	}

      /// Get the tracking subnode
      PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
      
      /// Check that it is there
      if (!svtxNode)
	{
	  svtxNode = new PHCompositeNode("SVTX");
	  dstNode->addNode(svtxNode);
	}
      
      _seed_track_map = new TpcSeedTrackMapv1();
      PHIODataNode<PHObject> *node
	= new PHIODataNode<PHObject>(_seed_track_map, _tpcseed_track_map_name);
      svtxNode->addNode(node);
    }
  
  _corrected_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"CORRECTED_TRKR_CLUSTER");
  if(_corrected_cluster_map)
    {
      std::cout << " Found CORRECTED_TRKR_CLUSTER node " << std::endl;
    }
  
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!_cluster_map)
      {
	std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
	return Fun4AllReturnCodes::ABORTEVENT;
      }

   _tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
   if(!_surfmaps)
     {
       std::cout << PHWHERE << "Error, can't find acts surface maps" << std::endl;
       return Fun4AllReturnCodes::ABORTEVENT;
     }
 
     return Fun4AllReturnCodes::EVENT_OK;
} 

double PHSiliconTpcTrackMatching::getBunchCrossing(unsigned int trid, double z_mismatch )
{
  // the bunch crossing separation is 110 nsec
  double vdrift = 8.00;  // cm /microsecond
  double z_bunch_separation = 0.107 * vdrift;  // 107 ns bunch crossing interval

  // The sign of z_mismatch will depend on which side of the TPC the tracklet is in
  SvtxTrack *track = _track_map->get(trid);

  double crossings = z_mismatch / z_bunch_separation;

  // Check the sign of z for the first cluster in the track
  // loop over associated clusters to get hits for TPC only, add to new track copy
  double clus_z = 0.0;
  ActsTransformations transformer;
  for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
       iter != track->end_cluster_keys();
       ++iter)
    {
      TrkrDefs::cluskey cluster_key = *iter;
      unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
      if(trkrid == TrkrDefs::tpcId)
	{

	  TrkrCluster *tpc_clus;
	  if(_corrected_cluster_map)
	    tpc_clus =  _corrected_cluster_map->findCluster(cluster_key);
	  else
	    tpc_clus =  _cluster_map->findCluster(cluster_key);

	  auto global = transformer.getGlobalPosition(tpc_clus,
						       _surfmaps,
						       _tGeometry);
	  clus_z = global[2];
	  break;   // we only need the first one  
	}
    }

  // if true z > 0, the z offset will be negative from a positive t0, so z_tpc - z_si will be negative for a positive crossing offset
  //       -- reverse the sign of crossings to get a positive crossing for a positive t0
  // if true z < 0, the z offset will be positive from a positive t0, z_tpc - z_si will be positive for a positive crossing offset
  //       -- the sign of crossings is OK
  // have to correct clus_z for the z_mismatch to see if it was really positive or negative
  if(clus_z -z_mismatch > 0)
    crossings *= -1.0;
  
  if(Verbosity() > 1) std::cout << "             trackid " << trid << " clus_z " << clus_z << " z_mismatch " << z_mismatch << " crossings " << crossings << std::endl;

  return crossings;
}

double PHSiliconTpcTrackMatching::getMedian(std::vector<double> &v)
{
  if(v.size() == 0) return NAN;

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

void PHSiliconTpcTrackMatching::addSiliconClusters( std::multimap<unsigned int, unsigned int> &tpc_matches )
{

  for(auto it = tpc_matches.begin(); it != tpc_matches.end(); ++it)
    {
      unsigned int tpcid = it->first;
      SvtxTrack *tpc_track = _track_map->get(tpcid);
      if(Verbosity() > 1) std::cout << "  tpcid " << tpcid << " original z " << tpc_track->get_z() << std::endl;
      
      // add the silicon cluster keys to the track
      unsigned int si_id = it->second;
      SvtxTrack *si_track = _track_map_silicon->get(si_id);
      if(Verbosity() > 1) std::cout << "  si track id " << si_id << std::endl;
      for (SvtxTrack::ConstClusterKeyIter si_iter = si_track->begin_cluster_keys();
	   si_iter != si_track->end_cluster_keys();
	   ++si_iter)
	{
	  TrkrDefs::cluskey si_cluster_key = *si_iter;
	  
	  if(Verbosity() > 1) 
	    cout << "   inserting si cluster key " << si_cluster_key << " into existing TPC track " << tpc_track->get_id() << endl;
	  
	  tpc_track->insert_cluster_key(si_cluster_key);
	  _assoc_container->SetClusterTrackAssoc(si_cluster_key, tpc_track->get_id());
	  

	}

      // update the track position to the si one
      tpc_track->set_x(si_track->get_x());
      tpc_track->set_y(si_track->get_y());
      tpc_track->set_z(si_track->get_z());
      
      if(Verbosity() > 2)
	std::cout << " TPC seed track ID " << tpc_track->get_id() << " si track id " << si_track->get_id()
		  << " new nclus " << tpc_track->size_cluster_keys() << std::endl;

      if(Verbosity() > 2)
	tpc_track->identify(); 
    }

  return;
}	  

// not used
void PHSiliconTpcTrackMatching::addSiliconClusters( std::multimap<int, std::pair<unsigned int, unsigned int>> &crossing_matches )
{

  for(auto it = crossing_matches.begin(); it != crossing_matches.end(); ++it)
    {
      unsigned int tpcid = it->second.first;
      SvtxTrack *tpc_track = _track_map->get(tpcid);
      if(Verbosity() > 1) std::cout << "  tpcid " << tpcid << " original z " << tpc_track->get_z() << std::endl;
      
      // add the silicon cluster keys to the track
      unsigned int si_id = it->second.second;
      SvtxTrack *si_track = _track_map_silicon->get(si_id);
      if(Verbosity() > 1) std::cout << "  si track id " << si_id << std::endl;
      for (SvtxTrack::ConstClusterKeyIter si_iter = si_track->begin_cluster_keys();
	   si_iter != si_track->end_cluster_keys();
	   ++si_iter)
	{
	  TrkrDefs::cluskey si_cluster_key = *si_iter;
	  
	  if(Verbosity() > 1) 
	    cout << "   inserting si cluster key " << si_cluster_key << " into existing TPC track " << tpc_track->get_id() << endl;
	  
	  tpc_track->insert_cluster_key(si_cluster_key);
	  _assoc_container->SetClusterTrackAssoc(si_cluster_key, tpc_track->get_id());
	  
	}

      // update the track position to the si one
      tpc_track->set_x(si_track->get_x());
      tpc_track->set_y(si_track->get_y());
      tpc_track->set_z(si_track->get_z());
      
      if(Verbosity() > 2)
	std::cout << " TPC seed track ID " << tpc_track->get_id() << " si track id " << si_track->get_id()
		  << " new nclus " << tpc_track->size_cluster_keys() << std::endl;

      if(Verbosity() > 2)
	tpc_track->identify();      
    }
  
  return;
}	  

void PHSiliconTpcTrackMatching::addSiliconClusters(  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>> &vertex_map)
{

  for(auto it = vertex_map.begin(); it != vertex_map.end(); ++it)
    {
      unsigned int tpcid = it->second.first;
      SvtxTrack *tpc_track = _track_map->get(tpcid);
      if(Verbosity() > 1) std::cout << "  tpcid " << tpcid << " original z " << tpc_track->get_z() << std::endl;
      
      // add the silicon cluster keys to the track
      unsigned int si_id = it->second.second;
      SvtxTrack *si_track = _track_map_silicon->get(si_id);
      if(Verbosity() > 1) std::cout << "  si track id " << si_id << std::endl;
      for (SvtxTrack::ConstClusterKeyIter si_iter = si_track->begin_cluster_keys();
	   si_iter != si_track->end_cluster_keys();
	   ++si_iter)
	{
	  TrkrDefs::cluskey si_cluster_key = *si_iter;
	  
	  if(Verbosity() > 1) 
	    cout << "   inserting si cluster key " << si_cluster_key << " into existing TPC track " << tpc_track->get_id() << endl;
	  
	  tpc_track->insert_cluster_key(si_cluster_key);
	  _assoc_container->SetClusterTrackAssoc(si_cluster_key, tpc_track->get_id());
	  
	}
  
      // update the track position to the si one
      tpc_track->set_x(si_track->get_x());
      tpc_track->set_y(si_track->get_y());
      tpc_track->set_z(si_track->get_z());
      
      if(Verbosity() > 2)
	{
	  std::cout << " TPC seed track ID " << tpc_track->get_id() << " si track id " << si_track->get_id()
		    << " new nclus " << tpc_track->size_cluster_keys() << std::endl;
	}
      
      if(Verbosity() > 2)  
	tpc_track->identify();
    }
  
  return;
}	  
      
void PHSiliconTpcTrackMatching::correctTpcClusterZ(  
						   std::map<unsigned int, double> &vertex_crossings_map,
						   std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map )
{
  if(vertex_crossings_map.size() == 0) return;

  // Correct the z positions of the clusters
  double vdrift = 8.00;  // cm /microsecond
  double z_bunch_separation = 0.107 * vdrift;  // 107 ns bunch crossing interval
  ActsTransformations transformer;
  for(auto [ivert, crossing] : vertex_crossings_map)
    {
      // the direction of the z shift depends on 
      // a) the sign of crossings (ie. t0 earlier or later),   b) the side of the TPC
      //      positive (corrected) z values: z values decrease for positive crossings, correction is positive = crossings * z_bunch_separation
      //      negative z values: z values increase for positive crossings, correction is negative = - crossings * z_bunch_separation

      if(Verbosity() > 1) std::cout << "      ivert " << ivert << " crossing " << crossing << std::endl;		  

      if(crossing == 0) continue;
   
      // loop over all tracks associated with this vertex
      // remember that any tracks not associated with one of the vertices are unusable
      auto  ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  unsigned int tpcid = it->second.first;
	  SvtxTrack *tpc_track = _track_map->get(tpcid);
	  if(Verbosity() > 1) std::cout << "  tpcid " << tpcid << " original z " << tpc_track->get_z() << std::endl;

	  // loop over associated clusters to get hits for TPC only, shift the z positions by z_shift
	  for (SvtxTrack::ConstClusterKeyIter iter = tpc_track->begin_cluster_keys();
	       iter != tpc_track->end_cluster_keys();
	       ++iter)
	    {
	      TrkrDefs::cluskey cluster_key = *iter;
	      unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
	      if(trkrid == TrkrDefs::tpcId)
		{
		  // get the cluster z
		  TrkrCluster *tpc_clus;
		  if(_corrected_cluster_map)
		    tpc_clus =  _corrected_cluster_map->findCluster(cluster_key);
		  else
		    tpc_clus =  _cluster_map->findCluster(cluster_key);

		  if(Verbosity() > 2) std::cout << "      original local cluster z " << tpc_clus->getLocalY() << std::endl;
		  auto global = transformer.getGlobalPosition(tpc_clus,
							      _surfmaps,
							      _tGeometry);
		  double clus_z = global[2];
		  if(Verbosity() > 2) std::cout << "  global: x " << global[0] << " y " << global[1] << " z " << global[2] << std::endl;

		  double corrected_clus_z;
		  if(clus_z > 0)
		    corrected_clus_z = clus_z + crossing * z_bunch_separation;
		  else
		    corrected_clus_z = clus_z - crossing * z_bunch_separation;
		  
		  auto surface =  transformer.getSurface(tpc_clus,  _surfmaps);
		  Acts::Vector3D center = surface->center(_tGeometry->geoContext)  / Acts::UnitConstants::cm;
		  double surfZCenter = center[2];
		  double local_z = corrected_clus_z - surfZCenter; 
		  tpc_clus->setLocalY(local_z);
		  if(Verbosity() > 2) std::cout << "      original clus_z " << clus_z << " corrected_clus_z " << corrected_clus_z << " local z " << local_z << std::endl;		  


		  if(Verbosity() > 2)
		    {
		      TrkrCluster *tpc_clus_check;
		      if(_corrected_cluster_map)
			tpc_clus_check =  _corrected_cluster_map->findCluster(cluster_key);
		      else
			tpc_clus_check =  _cluster_map->findCluster(cluster_key);
		      
		      auto global_check = transformer.getGlobalPosition(tpc_clus_check,
									_surfmaps,
									_tGeometry);
		      std::cout << "  global check: x " << global_check[0] << " y " << global_check[1] << " z " << global_check[2] << std::endl;
		    }
		}
	    }
	}
    }
  return;
}

void PHSiliconTpcTrackMatching::cleanVertexMap(  
						   std::map<unsigned int, double> &vertex_crossings_map,
						   std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map,
						   std::map<unsigned int, int> &tpc_crossing_map )
{
  if(vertex_crossings_map.size() == 0) return;

  // Sanity check: is the initial rough crossing estimate consistent with the crossing number for this si vertex?

  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>> bad_map;
  for(auto [ivert, crossing] : vertex_crossings_map)
    {
      if(Verbosity() > 1) std::cout << " CleanVertexMap:   ivert " << ivert << " crossing " << crossing << std::endl;		  
   
      // loop over all tracks associated with this vertex
      // remember that any tracks not associated with one of the vertices are unusable, they are ignored
      auto  ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  unsigned int tpcid = it->second.first;
	  unsigned int si_id = it->second.second;

	  // need the initial crossing estimate for comparison - stored in tpc_ crossing_map
	  auto sit = tpc_crossing_map.find(tpcid);	  
	  if(abs(sit->second - crossing) > 2)
	    {
	      // Fails the sanity check, mark for removal from the vertex map
	      if(Verbosity() > 1) 
		std::cout << "      Crossing mismatch: ivert " << ivert << " tpc id " << tpcid  << " vert crossing " << crossing << " track crossing " << sit->second << std::endl;

	      bad_map.insert(std::make_pair(ivert, std::make_pair(tpcid, si_id)));
	    }
	}
    }

  // If we delete an entry from vertex map, the TPC track is not associated with a bunch crossing, 
  //    -- the cluster z values are not changed, and the silicon clusters are not added. 
  // The track remains in the track map as a TPC-only track, by default assumed to be in bunch crossing zero
  // If multiple entries are deleted from vertex map, the TPC-only track copies are left in the track map, and also in _seed-track_map
  // The track cleaner will remove all but one of them, based on chisq from the Acts fitter.

  // remove bad entries from vertex_map so the wrong silicon is not associated
  for(auto [ivert, id_pair] : bad_map)
    {
      unsigned int tpc_id = id_pair.first;
      unsigned int si_id = id_pair.second;

      // Have to iterate over vertex_map and examine each pair to find the one matching bad_map
      // this logic works because we call the equal range on vertex_map for every id_pair
      // so we only delete one entry per equal range call
      auto ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  if(it->second.first == tpc_id && it->second.second == si_id)
	    {
	      if(Verbosity() > 1) std::cout << "   erasing entry for ivert " << ivert << " tpc_id " << tpc_id << " si_id " << si_id << std::endl;
	      vertex_map.erase(it);
	      break;  // the iterator is no longer valid
	    }
	}
    } 
      
  return;
}

void PHSiliconTpcTrackMatching::getCrossingNumber(  
						  std::vector<double> &vertex_list,
						  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map, 
						  std::map<unsigned int, double> &vertex_crossings_map)
{
  if(vertex_list.size() == 0) return;

  for(unsigned int ivert = 0; ivert<vertex_list.size(); ++ivert)
    {     
      std::vector<double> crossing_vec;

      double vert_z = vertex_list[ivert];
      if(Verbosity() > 1) std::cout << "Vertex " << ivert << " vertex z " << vert_z << std::endl;
      
      auto  ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  unsigned int tpcid = it->second.first;
	  double tpc_z =  _track_map->get(tpcid)->get_z();
	  double z_mismatch = tpc_z - vert_z;

	  double crossings = getBunchCrossing(tpcid, z_mismatch);
	  crossing_vec.push_back(crossings);

	  if(Verbosity() > 1) std::cout << "   tpc_id " << tpcid << " tpc pT " << _track_map->get(tpcid)->get_pt() 
					<< " tpc eta " << _track_map->get(tpcid)->get_eta() << " si_id " << it->second.second 
					<< " tpc z " << tpc_z  << " crossings " << crossings << std::endl;
	}

      if(crossing_vec.size() == 0) continue;
      double crossing_median = getMedian(crossing_vec);
  
      double crossing_avge = 0.0;
      double crossing_wt = 0.0;
      for(auto cross : crossing_vec)
	{
	  if(fabs(cross-crossing_median) < 1.0)
	    {
	      crossing_avge += cross;
	      crossing_wt++;	      
	    }
	}
      crossing_avge /= crossing_wt;      
      double crossing_avge_rounded  =  round(crossing_avge);
      vertex_crossings_map.insert(std::make_pair(ivert, crossing_avge_rounded));

      if(Verbosity() > 1) std::cout << "     crossing_median " << crossing_median << " crossing average = " << crossing_avge << " crossing_wt " << crossing_wt 
				    << " crossing integer " << crossing_avge_rounded << std::endl;
    }
  return;
}

void PHSiliconTpcTrackMatching::getSiVertexList(  
						std::multimap<double, std::pair<unsigned int, unsigned int>> &si_sorted_map,
						  std::vector<double> &vertex_list,
						std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map)
{
  if(si_sorted_map.size() == 0) return;

  // process si_sorted_map to get the vertex list, and the average silicon z value for each vertex
  double zkeep = 999;
  double zavge = 0.0;
  double zwt = 0.0;
  unsigned int nvert = 0;
  for(auto it = si_sorted_map.begin(); it != si_sorted_map.end(); ++it)
    {
      double z_si = it->first;
      std::pair<unsigned int, unsigned int> id_pair = it->second;
      if(Verbosity() > 1) std::cout << " z_si " << z_si << " tpc_id " << it->second.first << " si_id " << it->second.second << std::endl; 
     
      if( (zkeep == 999) || (fabs(z_si - zkeep) < _si_vertex_dzmax) )
	{
	  vertex_map.insert(std::make_pair(nvert, id_pair));
	  zavge+= z_si;
	  zwt++;
	  zkeep = z_si;
	}
      else
	{
	  zavge /= zwt;
	  vertex_list.push_back(zavge);

	  // start a new vertex
	  nvert++;
	  zavge = z_si;
	  zwt = 1.0;
	  zkeep = z_si;
	  vertex_map.insert(std::make_pair(nvert, id_pair));	    
	}
    }
  // get the last one
  zavge /= zwt;
  vertex_list.push_back(zavge);    

  return;
}

void PHSiliconTpcTrackMatching::findEtaPhiMatches(  
						  std::set<unsigned int> &tpc_matched_set,
						  std::multimap<unsigned int, unsigned int> &tpc_matches )
{
  // loop over the original TPC tracks
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      // we may add tracks to the map, so we stop at the last original track
      
      _tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() > 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet_tpc->get_phi()
	    << endl;
	}

      // This will always end up as a final track, no matter what - add it to the seed-track-map
      if(Verbosity() > 1) 
	std::cout << " TPC seed ID " << _tracklet_tpc->get_id() << " original nclus " << _tracklet_tpc->size_cluster_keys() << std::endl;
      _seed_track_map->addAssoc(_tracklet_tpc->get_id(), _tracklet_tpc->get_id());

      double tpc_phi = atan2(_tracklet_tpc->get_py(), _tracklet_tpc->get_px());
      double tpc_eta = _tracklet_tpc->get_eta();
      double tpc_pt = sqrt( pow(_tracklet_tpc->get_px(),2) + pow(_tracklet_tpc->get_py(),2) );

      // phi correction for PHTpcTracker tracklets is charge dependent
      double sign_phi_correction = _tracklet_tpc->get_charge();

      /// Correct the correction for the field direction
      /// Kludge to get the phi matching correct based on the field
      /// direction
      if(_field.find("2d") != std::string::npos)
	{
	  sign_phi_correction *= -1;
	  if(_fieldDir > 0)
	    sign_phi_correction *= -1;
	}

      // this factor will increase the window size at low pT
      // otherwise the matching efficiency drops off at low pT
      // it would be better if this was a smooth function
      double mag = 1.0;
      if(tpc_pt < 6.0) mag = 2;
      if(tpc_pt < 3.0)  mag = 4.0;

      if(Verbosity() > 3)
	{
	  cout << "Original TPC tracklet:" << endl;
	  _tracklet_tpc->identify();
	}

      // This is no longer applicable I believe, and now wrong - look into it
      //=============================
      // correct the TPC tracklet phi for the space charge offset, if this is the calib pass
      // this is done just to let us tighten up the matching window
      if(_sc_calib_flag)
	{
	  tpc_phi -= fscdphi->Eval(tpc_eta);
	}
      // the distortion correction can push tpc_phi outside +/- M_PI
      if(tpc_phi < - M_PI) tpc_phi += 2.0*M_PI;
      if(tpc_phi > M_PI) tpc_phi -= 2.0*M_PI;
      //=============================

      double tpc_x = _tracklet_tpc->get_x();
      double tpc_y = _tracklet_tpc->get_y();
      double tpc_z = _tracklet_tpc->get_z();

      // Now search the silicon track list for a match in eta and phi
      for (auto phtrk_iter_si = _track_map_silicon->begin();
	   phtrk_iter_si != _track_map_silicon->end(); 
	   ++phtrk_iter_si)
	{
	  _tracklet_si = phtrk_iter_si->second;	  

	  double si_phi = atan2(_tracklet_si->get_py(), _tracklet_si->get_px());
	  double si_eta = _tracklet_si->get_eta();
	  double si_x = _tracklet_si->get_x();
	  double si_y = _tracklet_si->get_y();
	  double si_z = _tracklet_si->get_z();

	  if(Verbosity() > 2)
	    {
	      cout << " testing for a match for TPC track " << _tracklet_tpc->get_id() << " with pT " << _tracklet_tpc->get_pt() 
		   << " with Si track " << _tracklet_si->get_id() << endl;	  
	      cout << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " <<   tpc_phi-si_phi << " phi search " << _phi_search_win*mag  << " tpc_eta " << tpc_eta 
		   << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << " eta search " << _eta_search_win*mag << endl;
	      std::cout << "      tpc x " << tpc_x << " si x " << si_x << " tpc y " << tpc_y << " si y " << si_y << " tpc_z " << tpc_z  << " si z " << si_z << std::endl;
	      std::cout << "      x search " << _x_search_win*mag << " y search " << _y_search_win*mag << " z search " << _z_search_win*mag  << std::endl;
	    }

	  bool eta_match = false;
	  bool phi_match = false;
	  bool position_match = false;
	  if(  fabs(tpc_eta - si_eta) < _eta_search_win * mag) eta_match = true;

	  // PHTpcTracker has a bias in the tracklet phi that depends on charge sign, PHCASeeding does not
	  if(_is_ca_seeder)
	    {
	      if(  fabs(tpc_phi - si_phi)  < _phi_search_win * mag) phi_match = true;
	    }
	  else
	    {
	      // PHTpcTracker
	      //double si_pt = sqrt( pow(_tracklet_si->get_px(),2) + pow(_tracklet_si->get_py(),2) );
	      double phi_search_win_lo = fdphi->Eval(tpc_pt) * sign_phi_correction -  _phi_search_win * mag;
	      double phi_search_win_hi = fdphi->Eval(tpc_pt) * sign_phi_correction +  _phi_search_win * mag;

	      if(Verbosity() > 10) 
		cout << " phi_search_win_lo " << phi_search_win_lo << " phi_search_win_hi " << phi_search_win_hi << endl;

	      if(  (tpc_phi - si_phi) > phi_search_win_lo && (tpc_phi - si_phi) < phi_search_win_hi) phi_match = true;	      
	    }

	  if(_pp_mode)
	    {
	      if(
		 fabs(tpc_x - si_x) < _x_search_win * mag
		 && fabs(tpc_y - si_y) < _y_search_win * mag 
		 )
		position_match = true;
	    }
	  else
	    {
	      if(
		 fabs(tpc_x - si_x) < _x_search_win * mag
		 && fabs(tpc_y - si_y) < _y_search_win * mag 
		 && fabs(tpc_z - si_z) < _z_search_win * mag 
		 )
		position_match = true;	    
	    }

	  if(eta_match && phi_match && position_match)
	    {
	      // got a match, add to the list
	      // These stubs are matched in eta, phi, x and y already
	      tpc_matches.insert(std::make_pair(_tracklet_tpc->get_id(), _tracklet_si->get_id()));
	      tpc_matched_set.insert(_tracklet_tpc->get_id());

	      if(Verbosity() > 1)  
		{
		  cout << " found a match for TPC track " << _tracklet_tpc->get_id() << " with Si track " << _tracklet_si->get_id() << endl;
		  cout << "          tpc_phi " << tpc_phi << " si_phi " <<  si_phi << " phi_match " << phi_match 
		       << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " eta_match " << eta_match << endl;
		  std::cout << "      tpc x " << tpc_x << " si x " << si_x << " tpc y " << tpc_y << " si y " << si_y << " tpc_z " << tpc_z  << " si z " << si_z << std::endl;
		}

	      // temporary!
	      if(_test_windows)
		cout << " Try_silicon:  pt " << tpc_pt << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi-si_phi  
		     << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << " tpc_x " << tpc_x << " tpc_y " << tpc_y << " tpc_z " << tpc_z 
		     << " dx " << tpc_x - si_x << " dy " << tpc_y - si_y << " dz " << tpc_z - si_z  
		     << endl;
	    }
	}
    }
  
  // We have all of the candidates for association, but we have to deal with cases of multiple matches for a tpc tracklet
  // We cannot move the TPC clusters independently for multiple matches, so we duplicate the tpc tracklet
  std::multimap<unsigned int, unsigned int> additional_tpc_matches;
  std::multimap<unsigned int, unsigned int> remove_tpc_matches;
  std::set<unsigned int> additional_tpc_matched_set;
  for(unsigned int tpcid : tpc_matched_set)
    {
      auto ret = tpc_matches.equal_range(tpcid);
      
      unsigned int size = std::distance(ret.first, ret.second);
      if(Verbosity() > 1) std::cout << " tpcid " << tpcid << " number of matches " << size << std::endl;

      if(size == 1)  continue;
      
      // make copy of all but the first TPC tracklet
      int counter = -1;
      for(auto it = ret.first; it != ret.second; ++it)
	{	  
	  counter++;
	  if(counter == 0)  continue;

	  unsigned int si_id = it->second;	  
	  auto this_track = _track_map->get(tpcid);
	  auto newTrack = std::make_unique<SvtxTrack_v2>();
	  const unsigned int lastTrackKey =  _track_map->empty() ? 0:std::prev(_track_map->end())->first; 
	  if(Verbosity() > 1) cout << "Extra match, add a new track to node tree with key " <<  lastTrackKey + 1 << endl;
	  
	  newTrack->set_id(lastTrackKey+1);
	  
	  newTrack->set_x(this_track->get_x());
	  newTrack->set_y(this_track->get_y());
	  newTrack->set_z(this_track->get_z());
	  
	  newTrack->set_charge(this_track->get_charge());
	  newTrack->set_px(this_track->get_px());
	  newTrack->set_py(this_track->get_py());
	  newTrack->set_pz(this_track->get_pz());
	  
	  for(int i = 0; i < 6; ++i)
	    {
	      for(int j = 0; j < 6; ++j)
		{
		  newTrack->set_error(i,j, this_track->get_error(i,j));
		}
	    }

	  _seed_track_map->addAssoc(this_track->get_id(), newTrack->get_id());
	  
	  // loop over associated clusters to get hits for TPC only, add to new track copy
	  for (SvtxTrack::ConstClusterKeyIter iter = this_track->begin_cluster_keys();
	       iter != this_track->end_cluster_keys();
	       ++iter)
	    {
	      TrkrDefs::cluskey cluster_key = *iter;
	      unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
	      if(trkrid == TrkrDefs::tpcId)
		{
		  newTrack->insert_cluster_key(cluster_key);
		  _assoc_container->SetClusterTrackAssoc(cluster_key, newTrack->get_id());
		}
	    }
	  
	  if(Verbosity() > 1) cout << "  -- inserting new track with id " << newTrack->get_id() << " from TPC tracklet " << tpcid << " into trackmap " << endl;
	  _track_map->insert(newTrack.get());

	  // add a map remove_tpc_matches so we can remove old matches with the same tpc id
	  remove_tpc_matches.insert(std::make_pair(tpcid, si_id));

	  additional_tpc_matches.insert(std::make_pair(newTrack->get_id(), si_id));
	  additional_tpc_matched_set.insert(newTrack->get_id());	  
	}   
    }
  for(auto [tpcid, si_id] : additional_tpc_matches)
    {
      tpc_matched_set.insert(tpcid);
      tpc_matches.insert(std::make_pair(tpcid, si_id));
      if(Verbosity() > 1) std::cout << "add to tpc_matches tpc_id " << tpcid << " si_id " << si_id << std::endl; 
    }
  for(auto [tpcid, si_id] : remove_tpc_matches)
    {
      if(Verbosity() > 1)  std::cout << "remove from tpc_matches tpc_id " << tpcid << " si_id " << si_id << std::endl; 

      // Have to iterate over the multimap and examine each si_id to find the one matching remove_tpc_matches
      auto ret = tpc_matches.equal_range(tpcid);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  if(it->first == tpcid && it->second == si_id)
	    {
	      tpc_matches.erase(it);
	      break;  // the iterator is no longer valid
	    }
	}
    }

  return;
}

void PHSiliconTpcTrackMatching::tagInTimeTracks(  
						std::multimap<unsigned int, unsigned int> &tpc_matches,
						std::set<int> &crossing_set,
						std::multimap<int, std::pair<unsigned int, unsigned int>> &crossing_matches,
						std::map<unsigned int, int> &tpc_crossing_map )
{

  for(auto [tpcid, si_id] : tpc_matches)
    {
      SvtxTrack *tpc_track = _track_map->get(tpcid);
      double tpc_z = tpc_track->get_z();
      double tpc_pt = tpc_track->get_pt();

      SvtxTrack *si_track =_track_map_silicon->get(si_id);
      double si_z = si_track->get_z();

      // this factor will increase the window size at low pT
      // otherwise the matching efficiency drops off at low pT
      // it would be better if this was a smooth function
      double mag = 1.0;
      if(tpc_pt < 6.0) mag = 2;
      if(tpc_pt < 3.0)  mag = 4.0;

      // Check for a match in z
      if(fabs(tpc_z - si_z) < _z_search_win * mag)
	{
	  //track from  triggered event
	  int crossing = 0;
	  crossing_matches.insert(std::make_pair(crossing,std::make_pair(tpc_track->get_id(), si_track->get_id())));
	  crossing_set.insert(crossing);
	  tpc_crossing_map.insert(std::make_pair(tpc_track->get_id(), crossing));
	  if(Verbosity() > 1)
	    std::cout << " triggered: tpc_trackid " << tpc_track->get_id() 
		      << " eta " << tpc_track->get_eta()
		      << " pT " << tpc_track->get_pt() 
		      << " si_trackid " << si_track->get_id() 
		      << " z_tpc " << tpc_track->get_z() 
		      << " z_si " << si_track->get_z() 
		      << " crossing " << crossing 
		      << std::endl;				  
	}
    }
  return;
}

void PHSiliconTpcTrackMatching::tagMatchCrossing(  
						std::multimap<unsigned int, unsigned int> &tpc_matches,
						std::set<int> &crossing_set,
						std::multimap<int, std::pair<unsigned int, unsigned int>> &crossing_matches,
						std::map<unsigned int, int> &tpc_crossing_map )
{

  for(auto [tpcid, si_id] : tpc_matches)
    {
      SvtxTrack *tpc_track = _track_map->get(tpcid);
      double tpc_z = tpc_track->get_z();
      //double tpc_pt = tpc_track->get_pt();

      SvtxTrack *si_track =_track_map_silicon->get(si_id);
      double si_z = si_track->get_z();

      // this is an initial estimate of the bunch crossing based on the z-mismatch for this track
      int crossing = (int) getBunchCrossing(tpc_track->get_id(), tpc_z - si_z);
      crossing_matches.insert(std::make_pair(crossing,std::make_pair(tpc_track->get_id(), si_track->get_id())));
      crossing_set.insert(crossing);
      tpc_crossing_map.insert(std::make_pair(tpc_track->get_id(), crossing));
      if(Verbosity() > 1)
	std::cout << " pileup: tpc_trackid " << tpc_track->get_id() 
		  << " eta " << tpc_track->get_eta()
		  << " pT " << tpc_track->get_pt() 
		  << " si_trackid " << si_track->get_id() 
		  << " z_tpc " << tpc_track->get_z() 
		  << " z_si " << si_track->get_z() 
		  << " crossing " << crossing 
		  << std::endl;		
    }

return;
}

void PHSiliconTpcTrackMatching::copySiliconClustersToCorrectedMap( )
{
  // loop over final track map, copy silicon clusters to corrected cluster map
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      SvtxTrack *track = phtrk_iter->second;

      // loop over associated clusters to get keys for silicon cluster
      for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
	   iter != track->end_cluster_keys();
	   ++iter)
	{
	  TrkrDefs::cluskey cluster_key = *iter;
	  const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
	  if(trkrid == TrkrDefs::mvtxId || trkrid == TrkrDefs::inttId)
	    {
	      TrkrCluster *cluster =  _cluster_map->findCluster(cluster_key);	
	      if( !cluster ) continue;
	      
	      TrkrCluster *newclus = _corrected_cluster_map->findOrAddCluster(cluster_key)->second;
	      newclus->CopyFrom( cluster );
	    }
	}      
    }
}


