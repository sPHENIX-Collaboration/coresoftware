#include "../trackreco/PHTruthTrackSeeding.h"
#include "AssocInfoContainer.h"

#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxVertexMap_v1.h>
#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex_v1.h>
#include <g4hough/SvtxTrackMap_v1.h>
#include <g4hough/SvtxTrack_FastSim.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit_v1.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/PHRandomSeed.h>
#include <phool/phool.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <memory>
#include "../trackreco/AssocInfoContainer.h"


#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl


using namespace std;

PHTruthTrackSeeding::PHTruthTrackSeeding(const std::string& name) :
		PHTrackSeeding(name),
		_g4truth_container(nullptr),
		phg4hits_svtx(nullptr),
		phg4hits_intt(nullptr),
		phg4hits_maps(nullptr),
		hitsmap(nullptr),
		cells_svtx(nullptr),
		cells_intt(nullptr),
		cells_maps(nullptr),
		_seeding_layers({7, 13, 19, 25, 31, 37, 40}),
		_min_clusters_per_track(0)
{}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode) {

	int ret = Fun4AllReturnCodes::ABORTRUN;

	ret = PHTrackSeeding::Setup(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = CreateNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process() {

	typedef std::map< int, std::set<SvtxCluster*> > TrkClustersMap;
	TrkClustersMap m_trackID_clusters;


	// Build TrackID -> Clusters map
	for (SvtxClusterMap::ConstIter cluster_itr = _cluster_map->begin();
			cluster_itr != _cluster_map->end(); ++cluster_itr) {
		SvtxCluster *cluster = cluster_itr->second;

		if(_seeding_layers.size() > 0 and
				(_seeding_layers.find(cluster->get_layer()) == _seeding_layers.end()))
			continue;

		SvtxHit* svtxhit = hitsmap->find(*cluster->begin_hits())->second;
		PHG4Cell* cell = nullptr;

		if(!cell and cells_svtx) cell = cells_svtx->findCell(svtxhit->get_cellid());
		if(!cell and cells_intt) cell = cells_intt->findCell(svtxhit->get_cellid());
		if(!cell and cells_maps) cell = cells_maps->findCell(svtxhit->get_cellid());

		if(!cell){
			if(verbosity >= 1) {
				LogError("!cell");
			}
			continue;
		}

		//cell->identify();

		for(PHG4Cell::EdepConstIterator hits_it = cell->get_g4hits().first;
				hits_it != cell->get_g4hits().second; hits_it++){

			PHG4Hit *phg4hit = nullptr;
			if(!phg4hit and phg4hits_svtx) phg4hit = phg4hits_svtx->findHit(hits_it->first);
			if(!phg4hit and phg4hits_intt) phg4hit = phg4hits_intt->findHit(hits_it->first);
			if(!phg4hit and phg4hits_maps) phg4hit = phg4hits_maps->findHit(hits_it->first);

			if(!phg4hit){
				if(verbosity >= 1) {
					LogError("!phg4hit");
				}
				continue;
			}

			//phg4hit->identify();

			int particle_id = phg4hit->get_trkid();

			TrkClustersMap::iterator it = m_trackID_clusters.find(particle_id);

			if(it != m_trackID_clusters.end()){
				it->second.insert(cluster);
			} else {
				std::set<SvtxCluster*> clusters;
				clusters.insert(cluster);
				m_trackID_clusters.insert(std::pair< int, std::set<SvtxCluster*> >(particle_id,clusters));
			}
		}
	}

	// Build track
	for (TrkClustersMap::const_iterator trk_clusers_itr = m_trackID_clusters.begin();
			trk_clusers_itr!=m_trackID_clusters.end(); ++trk_clusers_itr) {
		if(trk_clusers_itr->second.size() > _min_clusters_per_track) {
			std::unique_ptr<SvtxTrack_FastSim> svtx_track(new SvtxTrack_FastSim());

			//TODO implement the track ID
			svtx_track->set_id(_track_map->size());
			svtx_track->set_truth_track_id(trk_clusers_itr->first);
			//to make through minimum pT cut
			svtx_track->set_px(10.);
			svtx_track->set_py(0.);
			svtx_track->set_pz(0.);
			for(SvtxCluster *cluster : trk_clusers_itr->second) {
				svtx_track->insert_cluster(cluster->get_id());
				_assoc_container->SetClusterTrackAssoc(cluster->get_id(), svtx_track->get_id());
			}
			_track_map->insert(svtx_track.get());
		}
	}

	if (verbosity >= 2) {
		for (SvtxTrackMap::Iter iter = _track_map->begin();
				iter != _track_map->end(); ++iter) {
			SvtxTrack* svtx_track = iter->second;
			svtx_track->identify();
			continue;
			//Print associated clusters;
			for (SvtxTrack::ConstClusterIter iter =
					svtx_track->begin_clusters();
					iter != svtx_track->end_clusters(); ++iter) {
				unsigned int cluster_id = *iter;
				SvtxCluster* cluster = _cluster_map->get(cluster_id);
				float radius = sqrt(
						cluster->get_x() * cluster->get_x()
								+ cluster->get_y() * cluster->get_y());
				cout << "Track ID: " << svtx_track->get_id() << ", Track pT: "
						<< svtx_track->get_pt() << ", Particle ID: "
						<< svtx_track->get_truth_track_id() << ", cluster ID: "
						<< cluster->get_id() << ", cluster radius: " << radius
						<< endl;
			}
		}
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::CreateNodes(PHCompositeNode* topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::GetNodes(PHCompositeNode* topNode) {
	_g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
	if (!_g4truth_container) {
		cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	phg4hits_svtx = findNode::getClass<PHG4HitContainer>(
			topNode, "G4HIT_SVTX");

	phg4hits_intt = findNode::getClass<PHG4HitContainer>(
			topNode, "G4HIT_SILICON_TRACKER");

	phg4hits_maps = findNode::getClass<PHG4HitContainer>(
			topNode, "G4HIT_MAPS");

	if (!phg4hits_svtx and phg4hits_intt and !phg4hits_maps) {
		if (verbosity >= 0) {
			cerr << PHWHERE << " ERROR: No PHG4HitContainer found!" << endl;
		}
		return Fun4AllReturnCodes::ABORTRUN;
	}

	hitsmap = nullptr;
	// get node containing the digitized hits
	hitsmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
	if (!hitsmap) {
		cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	cells_svtx = findNode::getClass<PHG4CellContainer>(
			topNode, "G4CELL_SVTX");

	cells_intt = findNode::getClass<PHG4CellContainer>(
			topNode, "G4CELL_SILICON_TRACKER");

	cells_maps = findNode::getClass<PHG4CellContainer>(
			topNode, "G4CELL_MAPS");

	if (!cells_svtx and !cells_intt and !cells_maps) {
		if (verbosity >= 0) {
			cerr << PHWHERE << " ERROR: No PHG4CellContainer found!" << endl;
		}
		return Fun4AllReturnCodes::ABORTRUN;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}
