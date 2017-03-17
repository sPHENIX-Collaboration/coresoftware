#include "PHG4SvtxDeadArea.h"

#include "SvtxHitMap.h"
#include "SvtxHit.h"

#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHRandomSeed.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <gsl/gsl_rng.h>

#include <iostream>

using namespace std;

PHG4SvtxDeadArea::PHG4SvtxDeadArea(const string &name) :
  SubsysReco(name),
  _hits(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) 
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  return;
}

PHG4SvtxDeadArea::~PHG4SvtxDeadArea()
{
  gsl_rng_free (RandomGenerator);
}

int PHG4SvtxDeadArea::InitRun(PHCompositeNode* topNode) {

  // get node containing the digitized hits
  _hits = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hits) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // setup our random number generator
  unsigned int seed = PHRandomSeed(); // fixed seed handled in PHRandomSeed()
  gsl_rng_set(RandomGenerator,seed);

  FillCylinderDeadAreaMap(topNode);
  FillMapsLadderDeadAreaMap(topNode);
  
  if (verbosity > 0) {
    cout << "====================== PHG4SvtxDeadArea::InitRun() ========================" << endl;
    cout << " Random number seed = " << seed << endl;    
    for (std::map<int,float>::iterator iter = _eff_by_layer.begin();
	 iter != _eff_by_layer.end();
	 ++iter) {
      cout << " Active area in Layer #" << iter->first << " set to " << iter->second*100.0 << "%" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxDeadArea::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();
 
  std::vector<unsigned int> remove_hits;
  
  for (SvtxHitMap::Iter iter = _hits->begin();
       iter != _hits->end();
       ++iter) {
    SvtxHit* hit = iter->second;

    if (gsl_rng_uniform_pos(RandomGenerator) > get_hit_efficiency(hit->get_layer())) {
      remove_hits.push_back(hit->get_id());
      if(verbosity > 5)
	cout << "removing hit" << hit->get_id() << endl;
    }
  }

  for (unsigned int i = 0; i < remove_hits.size(); ++i) {
    _hits->erase(remove_hits[i]);
  }
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxDeadArea::End(PHCompositeNode* topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4SvtxDeadArea::FillCylinderDeadAreaMap(PHCompositeNode* topNode) {

  PHG4CylinderCellGeomContainer *geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
    

  if (!geom_container) return;
  
  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {
    int layer = layeriter->second->get_layer();
 
    if (_eff_by_layer.find(layer) == _eff_by_layer.end()) {
      _eff_by_layer.insert(std::make_pair(layer,1.0));
    }    
  }
  
  return;
}


void
PHG4SvtxDeadArea::FillMapsLadderDeadAreaMap(PHCompositeNode* topNode) {

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_MAPS");
    
  if (!geom_container) return;

  if(verbosity > 0)
    cout << "Found CylinderGeom_MAPS" << endl;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {
    int layer = layeriter->second->get_layer();

    if (_eff_by_layer.find(layer) == _eff_by_layer.end()) {
      _eff_by_layer.insert(std::make_pair(layer,1.0));
    }    
  }
  
  return;
}
