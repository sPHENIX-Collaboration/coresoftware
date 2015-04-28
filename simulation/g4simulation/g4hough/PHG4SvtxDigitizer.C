#include "PHG4SvtxDigitizer.h"

#include "SvtxHitMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>

#include <iostream>
#include <cmath>

using namespace std;

PHG4SvtxDigitizer::PHG4SvtxDigitizer(const char* name) :
  SubsysReco(name),
  _max_adc(),
  _energy_scale(),
  _hitmap(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) {
}

int PHG4SvtxDigitizer::InitRun(PHCompositeNode* topNode) {

  //-------------
  // Add Hit Node
  //-------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode 
    = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
    
  // Create the SVX node if required
  PHCompositeNode* svxNode 
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","SVTX"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }
  
  // Create the Hit node if required
  SvtxHitMap *svxhits = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!svxhits) {
    svxhits = new SvtxHitMap();
    PHIODataNode<PHObject> *SvtxHitMapNode =
      new PHIODataNode<PHObject>(svxhits, "SvtxHitMap", "PHObject");
    svxNode->addNode(SvtxHitMapNode);
  }

  CalculateCylinderCellADCScale(topNode);
  CalculateLadderCellADCScale(topNode);
  
  //----------------
  // Report Settings
  //----------------
  
  if (verbosity >= 0) {
    cout << "====================== PHG4SvtxDigitizer::InitRun() =====================" << endl;
    cout << " CVS Version: $Id: PHG4SvtxDigitizer.C,v 1.5 2015/04/21 23:47:10 pinkenbu Exp $" << endl;
    for (std::map<int,unsigned int>::iterator iter = _max_adc.begin();
	 iter != _max_adc.end();
	 ++iter) {
      cout << " Max ADC in Layer #" << iter->first << " = " << iter->second << endl;
    }
    for (std::map<int,float>::iterator iter = _energy_scale.begin();
	 iter != _energy_scale.end();
	 ++iter) {
      cout << " Energy per ADC in Layer #" << iter->first << " = " << 1.0e6*iter->second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxDigitizer::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();

  _hitmap = NULL;
  PHTypedNodeIterator<SvtxHitMap> iter(topNode);
  PHIODataNode<SvtxHitMap> *SvtxHitMapNode = iter.find("SvtxHitMap");
  if (!SvtxHitMapNode) {
    cout << PHWHERE << " ERROR: Can't find SvtxHitMap." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  } else {
    _hitmap = (SvtxHitMap*)SvtxHitMapNode->getData();
  }

  _hitmap->Reset();
  
  DigitizeCylinderCells(topNode);
  DigitizeLadderCells(topNode);

  PrintHits(topNode);
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4SvtxDigitizer::CalculateCylinderCellADCScale(PHCompositeNode *topNode) {

  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

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
   
    float minpath = pitch;
    if (length < minpath) minpath = length;
    if (thickness < minpath) minpath = thickness;
    float mip_e = 0.003876*minpath;  

    if (_max_adc.find(layer) == _max_adc.end()) {
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }    

  return;
}

void PHG4SvtxDigitizer::CalculateLadderCellADCScale(PHCompositeNode *topNode) {

  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

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
   
    float minpath = pitch;
    if (length < minpath) minpath = length;
    if (thickness < minpath) minpath = thickness;
    float mip_e = 0.003876*minpath;  

    if (_max_adc.find(layer) == _max_adc.end()) {
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }    

  return;
}

void PHG4SvtxDigitizer::DigitizeCylinderCells(PHCompositeNode *topNode) {

  //----------
  // Get Nodes
  //----------
 
  PHG4CylinderCellContainer* cells = 0;
  PHTypedNodeIterator<PHG4CylinderCellContainer> celliter(topNode);
  PHIODataNode<PHG4CylinderCellContainer>* cell_container_node = celliter.find("G4CELL_SVTX");
  if (cell_container_node) {
    cells = (PHG4CylinderCellContainer*) cell_container_node->getData();
  }
  if (!cells) return; 
  
  //-------------
  // Digitization
  //-------------

  PHG4CylinderCellContainer::ConstRange cellrange = cells->getCylinderCells();
  for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
      celliter != cellrange.second;
      ++celliter) {
    
    PHG4CylinderCell* cell = celliter->second;
    
    SvtxHit hit;

    hit.set_layer(cell->get_layer());
    hit.set_cellid(cell->get_cell_id());

    unsigned int adc = cell->get_edep() / _energy_scale[hit.get_layer()];
    if (adc > _max_adc[hit.get_layer()]) adc = _max_adc[hit.get_layer()]; 
    float e = _energy_scale[hit.get_layer()] * adc;

    hit.set_adc(adc);
    hit.set_e(e);

    SvtxHit* ptr = _hitmap->insert(hit);      
    if (!ptr->IsValid()) {
      static bool first = true;
      if (first) {
	cout << PHWHERE << "ERROR: Incomplete SvtxHits are being created" << endl;
	ptr->identify();
	first = false;
      }
    }
  }
  
  return;
}

void PHG4SvtxDigitizer::DigitizeLadderCells(PHCompositeNode *topNode) {

  //----------
  // Get Nodes
  //----------
 
  PHG4CylinderCellContainer* cells = 0;
  PHTypedNodeIterator<PHG4CylinderCellContainer> celliter(topNode);
  PHIODataNode<PHG4CylinderCellContainer>* cell_container_node = celliter.find("G4CELL_SILICON_TRACKER");
  if (cell_container_node) {
    cells = (PHG4CylinderCellContainer*) cell_container_node->getData();
  }
  if (!cells) return; 
  
  //-------------
  // Digitization
  //-------------

  vector<PHG4CylinderCell*> cell_list;
  PHG4CylinderCellContainer::ConstRange cellrange = cells->getCylinderCells();
  for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
      celliter != cellrange.second;
      ++celliter) {
    
    PHG4CylinderCell* cell = celliter->second;
    
    SvtxHit hit;

    hit.set_layer(cell->get_layer());
    hit.set_cellid(cell->get_cell_id());

    unsigned int adc = cell->get_edep() / _energy_scale[hit.get_layer()];
    if (adc > _max_adc[hit.get_layer()]) adc = _max_adc[hit.get_layer()]; 
    float e = _energy_scale[hit.get_layer()] * adc;
    
    hit.set_adc(adc);
    hit.set_e(e);
        
    SvtxHit* ptr = _hitmap->insert(hit);      
    if (!ptr->IsValid()) {
      static bool first = true;
      if (first) {
	cout << PHWHERE << "ERROR: Incomplete SvtxHits are being created" << endl;
	ptr->identify();
	first = false;
      }
    }
  }
  
  return;
}

void PHG4SvtxDigitizer::PrintHits(PHCompositeNode *topNode) {

  if (verbosity >= 1) {

    PHTypedNodeIterator<SvtxHitMap> hititer(topNode);
    PHIODataNode<SvtxHitMap> *SvtxHitMapNode = hititer.find("SvtxHitMap");
    if (!SvtxHitMapNode) return;
    
    cout << "================= PHG4SvtxDigitizer::process_event() ====================" << endl;
  
    SvtxHitMap *hitlist = (SvtxHitMap*)SvtxHitMapNode->getData();

    cout << " Found and recorded the following " << hitlist->size() << " hits: " << endl;

    unsigned int ihit = 0;
    for (SvtxHitMap::Iter iter = hitlist->begin();
	 iter != hitlist->end();
	 ++iter) {

      SvtxHit* hit = &iter->second;
      cout << ihit << " of " << hitlist->size() << endl;
      hit->identify();
      ++ihit;
    }
    
    cout << "===========================================================================" << endl;
  }
  
  return;
}
