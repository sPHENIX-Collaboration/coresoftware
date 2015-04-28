#include "PHG4SvtxDeadArea.h"

#include "SvtxHitMap.h"
#include "SvtxHit.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>
#include <fun4all/recoConsts.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <iostream>

using namespace std;

PHG4SvtxDeadArea::PHG4SvtxDeadArea(const char* name) :
  SubsysReco(name),
  _eff_by_layer(),
  _hits(NULL),
  _rand(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) {
}

int PHG4SvtxDeadArea::InitRun(PHCompositeNode* topNode) {

  // get node containing the digitized hits
  _hits = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hits) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // setup our random number generator
  recoConsts *rc = recoConsts::instance();
  int seed = rc->get_IntFlag("RANDOMSEED");
  _rand = new TRandom3(seed);

  FillCylinderDeadAreaMap(topNode);
  FillLadderDeadAreaMap(topNode);
  
  if (verbosity >= 0) {
    cout << "====================== PHG4SvtxDeadArea::InitRun() ========================" << endl;
    cout << " CVS Version: $Id: PHG4SvtxDeadArea.C,v 1.6 2015/04/21 23:47:09 pinkenbu Exp $" << endl;
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
    SvtxHit* hit = &iter->second;

    if (_rand->Rndm() > get_hit_efficiency(hit->get_layer())) {
      remove_hits.push_back(hit->get_id());
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

  PHG4CylinderCellContainer *cells = NULL;
  PHG4CylinderCellGeomContainer *geom_container = NULL;
    
  PHTypedNodeIterator<PHG4CylinderCellContainer> celliter(topNode);
  PHIODataNode<PHG4CylinderCellContainer>* cell_container_node = celliter.find("G4CELL_SVTX");
  if (cell_container_node) {
    cells = (PHG4CylinderCellContainer*) cell_container_node->getData();
  }

  PHTypedNodeIterator<PHG4CylinderCellGeomContainer> geomiter(topNode);
  PHIODataNode<PHG4CylinderCellGeomContainer>* PHG4CylinderCellGeomContainerNode = geomiter.find("CYLINDERCELLGEOM_SVTX");
  if (PHG4CylinderCellGeomContainerNode) {
    geom_container = (PHG4CylinderCellGeomContainer*) PHG4CylinderCellGeomContainerNode->getData();
  }

  if (!geom_container || !cells) return;
  
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

void PHG4SvtxDeadArea::FillLadderDeadAreaMap(PHCompositeNode* topNode) {

  PHG4CylinderCellContainer *cells = NULL;
  PHG4CylinderGeomContainer *geom_container = NULL;
    
  PHTypedNodeIterator<PHG4CylinderCellContainer> celliter(topNode);
  PHIODataNode<PHG4CylinderCellContainer>* cell_container_node = celliter.find("G4CELL_SILICON_TRACKER");
  if (cell_container_node) {
    cells = (PHG4CylinderCellContainer*) cell_container_node->getData();
  }

  PHTypedNodeIterator<PHG4CylinderGeomContainer> geomiter(topNode);
  PHIODataNode<PHG4CylinderGeomContainer>* PHG4CylinderGeomContainerNode = geomiter.find("CYLINDERGEOM_SILICON_TRACKER");
  if (PHG4CylinderGeomContainerNode) {
    geom_container = (PHG4CylinderGeomContainer*) PHG4CylinderGeomContainerNode->getData();
  }
  
  if (!geom_container || !cells) return;

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
