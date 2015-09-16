#include "PHG4SvtxDigitizer.h"

#include "SvtxHitMap.h"
#include "SvtxHit.h"
#include "SvtxHit_v1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
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

PHG4SvtxDigitizer::PHG4SvtxDigitizer(const string &name) :
  SubsysReco(name),
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
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
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
  
  if (verbosity > 0) {
    cout << "====================== PHG4SvtxDigitizer::InitRun() =====================" << endl;
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

  _hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hitmap) 
    {
      cout << PHWHERE << " ERROR: Can't find SvtxHitMap." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
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

  PHG4CylinderCellContainer *cells = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SVTX");
  PHG4CylinderCellGeomContainer *geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
    

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

  PHG4CylinderCellContainer *cells = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SILICON_TRACKER");
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");
    
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
 
  PHG4CylinderCellContainer* cells = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SVTX");
  if (!cells) return; 
  
  //-------------
  // Digitization
  //-------------

  PHG4CylinderCellContainer::ConstRange cellrange = cells->getCylinderCells();
  for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
      celliter != cellrange.second;
      ++celliter) {
    
    PHG4CylinderCell* cell = celliter->second;
    
    SvtxHit_v1 hit;

    hit.set_layer(cell->get_layer());
    hit.set_cellid(cell->get_cell_id());

    unsigned int adc = cell->get_edep() / _energy_scale[hit.get_layer()];
    if (adc > _max_adc[hit.get_layer()]) adc = _max_adc[hit.get_layer()]; 
    float e = _energy_scale[hit.get_layer()] * adc;

    hit.set_adc(adc);
    hit.set_e(e);

    SvtxHit* ptr = _hitmap->insert(&hit);      
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
 
  PHG4CylinderCellContainer* cells = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SILICON_TRACKER");
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
    
    SvtxHit_v1 hit;

    hit.set_layer(cell->get_layer());
    hit.set_cellid(cell->get_cell_id());

    unsigned int adc = cell->get_edep() / _energy_scale[hit.get_layer()];
    if (adc > _max_adc[hit.get_layer()]) adc = _max_adc[hit.get_layer()]; 
    float e = _energy_scale[hit.get_layer()] * adc;
    
    hit.set_adc(adc);
    hit.set_e(e);
        
    SvtxHit* ptr = _hitmap->insert(&hit);      
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

    SvtxHitMap *hitlist = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
    if (!hitlist) return;
    
    cout << "================= PHG4SvtxDigitizer::process_event() ====================" << endl;
  

    cout << " Found and recorded the following " << hitlist->size() << " hits: " << endl;

    unsigned int ihit = 0;
    for (SvtxHitMap::Iter iter = hitlist->begin();
	 iter != hitlist->end();
	 ++iter) {

      SvtxHit* hit = iter->second;
      cout << ihit << " of " << hitlist->size() << endl;
      hit->identify();
      ++ihit;
    }
    
    cout << "===========================================================================" << endl;
  }
  
  return;
}
