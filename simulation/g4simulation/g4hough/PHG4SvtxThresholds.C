#include "PHG4SvtxThresholds.h"

#include "SvtxHitMap.h"
#include "SvtxHit.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <iostream>

using namespace std;

PHG4SvtxThresholds::PHG4SvtxThresholds(const char* name) :
  SubsysReco(name),
  _fraction_of_mip(0.5),
  _thresholds_by_layer(),
  _use_thickness_mip(),
  _hits(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) {
}

int PHG4SvtxThresholds::InitRun(PHCompositeNode* topNode) {

  // get node containing the digitized hits
  _hits = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hits) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  CalculateCylinderThresholds(topNode);
  CalculateLadderThresholds(topNode);
  
  if (verbosity >= 0) {
    cout << "====================== PHG4SvtxThresholds::InitRun() ======================" << endl;
    cout << " CVS Version: $Id: PHG4SvtxThresholds.C,v 1.7 2015/04/21 23:47:10 pinkenbu Exp $" << endl;
    cout << " Fraction of expected MIP energy = " << _fraction_of_mip << endl;
    for (std::map<int,float>::iterator iter = _thresholds_by_layer.begin();
	 iter != _thresholds_by_layer.end();
	 ++iter) {
      cout << " MIP threshold basis for Layer #" << iter->first << ": ";
      if (get_use_thickness_mip(iter->first)) cout << "Radial Thickness Penetration";
      else cout << "Short-Axis Penetration";
      cout << endl;
    }
    for (std::map<int,float>::iterator iter = _thresholds_by_layer.begin();
	 iter != _thresholds_by_layer.end();
	 ++iter) {
      cout << " Cell Threshold in Layer #" << iter->first << " = " << 1.0e6*iter->second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxThresholds::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();

  std::vector<unsigned int> remove_hits;
  
  for (SvtxHitMap::Iter iter = _hits->begin();
       iter != _hits->end();
       ++iter) {
    SvtxHit* hit = &iter->second;

    if (hit->get_e() < get_threshold_by_layer(hit->get_layer())) {
      remove_hits.push_back(hit->get_id());
    }
  }

  for (unsigned int i = 0; i < remove_hits.size(); ++i) {
    _hits->erase(remove_hits[i]);
  }
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxThresholds::End(PHCompositeNode* topNode) {

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4SvtxThresholds::CalculateCylinderThresholds(PHCompositeNode* topNode) {

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
    float thickness = (layeriter->second)->get_thickness();
    float pitch = (layeriter->second)->get_phistep()*(layeriter->second)->get_radius();
    float length = (layeriter->second)->get_zstep();

    if (get_use_thickness_mip(layer)) {
      // Si MIP energy = 3.876 MeV / cm
      float threshold = _fraction_of_mip*0.003876*thickness;
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    } else {
      float minpath = pitch;
      if (length < minpath) minpath = length;
      if (thickness < minpath) minpath = thickness;
	
      // Si MIP energy = 3.876 MeV / cm
      float threshold = _fraction_of_mip*0.003876*minpath;  
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    }
  }
  
  return;
}

void PHG4SvtxThresholds::CalculateLadderThresholds(PHCompositeNode* topNode) {

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
    float thickness = (layeriter->second)->get_thickness();
    float pitch = (layeriter->second)->get_strip_y_spacing();
    float length = (layeriter->second)->get_strip_z_spacing();

    if (get_use_thickness_mip(layer)) {
      // Si MIP energy = 3.876 MeV / cm
      float threshold = _fraction_of_mip*0.003876*thickness;
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    } else {
      float minpath = pitch;
      if (length < minpath) minpath = length;
      if (thickness < minpath) minpath = thickness;
      
      // Si MIP energy = 3.876 MeV / cm
      float threshold = _fraction_of_mip*0.003876*minpath;  
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    }
  }
  
  return;
}
