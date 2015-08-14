
#include "SvtxClusterEval.h"

#include "SvtxHitEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <cstdlib>
#include <set>
#include <map>
#include <float.h>

using namespace std;

SvtxClusterEval::SvtxClusterEval(PHCompositeNode* topNode)
  : _topNode(topNode),
    _hiteval(topNode),
    _cache_all_truth_hits(),
    _cache_max_truth_hit_by_energy(),
    _cache_all_truth_particles(),
    _cache_max_truth_particle_by_energy(),
    _cache_all_clusters_from_particle(),
    _cache_all_clusters_from_g4hit(),
    _cache_get_energy_contribution() {
}

void SvtxClusterEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_max_truth_hit_by_energy.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_energy.clear();
  _cache_all_clusters_from_particle.clear();
  _cache_all_clusters_from_g4hit.clear();
  _cache_get_energy_contribution.clear();

  _hiteval.next_event(topNode);
  
  _topNode = topNode;  
}

std::set<PHG4Hit*> SvtxClusterEval::all_truth_hits(SvtxCluster* cluster) {

  if (_cache_all_truth_hits.find(cluster) != _cache_all_truth_hits.end()) {
    return _cache_all_truth_hits[cluster];
  }
  
  std::set<PHG4Hit*> truth_hits;
  
  // need things off of the DST...
  SvtxHitMap* hitmap = findNode::getClass<SvtxHitMap>(_topNode,"SvtxHitMap");
  if (!hitmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxHitMap" << endl;
    exit(-1);
  }

  PHG4CylinderCellContainer* g4cells = findNode::getClass<PHG4CylinderCellContainer>(_topNode,"G4CELL_SVTX");
  if (!g4cells) {
    cerr << PHWHERE << " ERROR: Can't find G4CELL_SVTX" << endl;
    exit(-1);
  }
    
  PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(_topNode,"G4HIT_SVTX");
  if (!g4hits) {
    cerr << PHWHERE << " ERROR: Can't find G4HIT_SVTX" << endl;
    exit(-1);
  }
  
  // loop over all hit cells
  for (SvtxCluster::ConstHitIter hiter = cluster->begin_hits();
       hiter != cluster->end_hits();
       ++hiter) {
    SvtxHit* hit = hitmap->get(*hiter);

    // hop from reco hit to g4cell
    PHG4CylinderCell *cell = g4cells->findCylinderCell(hit->get_cellid());
    if (!cell) continue;

    // loop over all the g4hits in this cell
    for (PHG4CylinderCell::EdepConstIterator g4iter = cell->get_g4hits().first;
	 g4iter != cell->get_g4hits().second;
	 ++g4iter) {
      
      PHG4Hit* g4hit = g4hits->findHit(g4iter->first);

      // fill output set
      truth_hits.insert(g4hit);
    }
  }

  _cache_all_truth_hits.insert(make_pair(cluster,truth_hits));
  
  return truth_hits;
}

PHG4Hit* SvtxClusterEval::max_truth_hit_by_energy(SvtxCluster* cluster) {
  
  if (_cache_max_truth_hit_by_energy.find(cluster) !=
      _cache_max_truth_hit_by_energy.end()) {
    return _cache_max_truth_hit_by_energy[cluster];
  }
  
  std::set<PHG4Hit*> hits = all_truth_hits(cluster);
  PHG4Hit* max_hit = NULL;
  float max_e = FLT_MIN;
  for (std::set<PHG4Hit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter) {
    PHG4Hit *hit = *iter;
    if (hit->get_edep() > max_e) {
      max_e = hit->get_edep();
      max_hit = hit;
    }
  }

  _cache_max_truth_hit_by_energy.insert(make_pair(cluster,max_hit));
  
  return max_hit;
}
  
std::set<PHG4Particle*> SvtxClusterEval::all_truth_particles(SvtxCluster* cluster) {

  if (_cache_all_truth_particles.find(cluster) !=
      _cache_all_truth_particles.end()) {
    return _cache_all_truth_particles[cluster];
  }
  
  std::set<PHG4Particle*> truth_particles;
  
  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster);
  if (g4hits.empty()) return truth_particles;
  
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(_topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* hit = *iter;
    PHG4Particle* particle = truthinfo->GetHit( hit->get_trkid() );
    if (!particle) continue;
    truth_particles.insert(particle);
  }

  _cache_all_truth_particles.insert(make_pair(cluster,truth_particles));
  
  return truth_particles;
}

PHG4Particle* SvtxClusterEval::max_truth_particle_by_energy(SvtxCluster* cluster) {

  if (_cache_max_truth_particle_by_energy.find(cluster) !=
      _cache_max_truth_particle_by_energy.end()) {
    return _cache_max_truth_particle_by_energy[cluster];
  }
  
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(_topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  // loop over all particles associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = NULL;
  float max_e = FLT_MIN;
  std::set<PHG4Particle*> particles = all_truth_particles(cluster);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter) {

    PHG4Particle* particle = *iter;
    float e = get_energy_contribution(cluster,particle);
    if (e > max_e) {
      max_e = e;
      max_particle = particle;      
    }
  }

  _cache_max_truth_particle_by_energy.insert(make_pair(cluster,max_particle));
  
  return max_particle;
}

std::set<SvtxCluster*> SvtxClusterEval::all_clusters_from(PHG4Particle* truthparticle) { 

  if (_cache_all_clusters_from_particle.find(truthparticle) !=
      _cache_all_clusters_from_particle.end()) {
    return _cache_all_clusters_from_particle[truthparticle];
  }
  
  // need things off of the DST...
  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(_topNode,"SvtxClusterMap");
  if (!clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }

  std::set<SvtxCluster*> clusters;
  
  // loop over all the clusters
  for (SvtxClusterMap::Iter iter = clustermap->begin();
       iter != clustermap->end();
       ++iter) {

    SvtxCluster* cluster = &iter->second;

    // loop over all truth particles connected to this cluster
    std::set<PHG4Particle*> particles = all_truth_particles(cluster);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
	 jter != particles.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;
      if (candidate->get_track_id() == truthparticle->get_track_id()) {
	clusters.insert(cluster);
      }    
    }
  }

  _cache_all_clusters_from_particle.insert(make_pair(truthparticle,clusters));
  
  return clusters;
}

std::set<SvtxCluster*> SvtxClusterEval::all_clusters_from(PHG4Hit* truthhit) {

  if (_cache_all_clusters_from_g4hit.find(truthhit) !=
      _cache_all_clusters_from_g4hit.end()) {
    return _cache_all_clusters_from_g4hit[truthhit];
  }
  
  // need things off of the DST...
  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(_topNode,"SvtxClusterMap");
  if (!clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }

  std::set<SvtxCluster*> clusters;
  
  // loop over all the clusters
  for (SvtxClusterMap::Iter iter = clustermap->begin();
       iter != clustermap->end();
       ++iter) {

    SvtxCluster* cluster = &iter->second;

    // loop over all truth hits connected to this cluster
    std::set<PHG4Hit*> hits = all_truth_hits(cluster);
    for (std::set<PHG4Hit*>::iterator jter = hits.begin();
	 jter != hits.end();
	 ++jter) {
      PHG4Hit* candidate = *jter;
      if (candidate->get_trkid() == truthhit->get_trkid()) {
	clusters.insert(cluster);
      }    
    }
  }

  _cache_all_clusters_from_g4hit.insert(make_pair(truthhit,clusters));
  
  return clusters;
}
  
// overlap calculations
float SvtxClusterEval::get_energy_contribution(SvtxCluster* cluster, PHG4Particle* particle) {

  if (_cache_get_energy_contribution.find(make_pair(cluster,particle)) !=
      _cache_get_energy_contribution.end()) {
    return _cache_get_energy_contribution[make_pair(cluster,particle)];
  }
  
  float energy = 0.0;
  std::set<PHG4Hit*> hits = all_truth_hits(cluster);
  for (std::set<PHG4Hit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter) {
    PHG4Hit* hit = *iter;
    if (hit->get_trkid() == particle->get_track_id()) {
      energy += hit->get_edep();
    }
  }

  _cache_get_energy_contribution.insert(make_pair(make_pair(cluster,particle),energy));
  
  return energy;
}
