#include "PHG4SvtxThresholds.h"

#include "SvtxHitMap.h"
#include "SvtxHit.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <iostream>

using namespace std;

PHG4SvtxThresholds::PHG4SvtxThresholds(const string &name) :
  SubsysReco(name),
  _fraction_of_mip(),
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
  CalculateMapsLadderThresholds(topNode);

  if (Verbosity() > 0) {
    cout << "====================== PHG4SvtxThresholds::InitRun() ======================" << endl;
    for (std::map<int,float>::iterator iter = _fraction_of_mip.begin();
	 iter != _fraction_of_mip.end();
	 ++iter) {
      cout << " Fraction of expected MIP energy for Layer #" << iter->first << ": ";
      cout << iter->second;
      cout << endl;
    }
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
    SvtxHit* hit = iter->second;
   
    if (hit->get_e() < get_threshold_by_layer(hit->get_layer())) {
      if(Verbosity() > 2)
	cout << "Removing hit with energy " << hit->get_e() <<  " in layer " << hit->get_layer() << " threshold = " << get_threshold_by_layer(hit->get_layer()) << endl;
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

  PHG4CylinderCellGeomContainer *geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
    
  if (!geom_container) return;
  
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
      float threshold = 0.0;      
      if (_fraction_of_mip.find(layer) != _fraction_of_mip.end()) {	
	threshold = _fraction_of_mip[layer]*0.003876*thickness;
	if(Verbosity() >2) 
	  cout << " using thickness: threshold = " << threshold << " thickness = " << thickness << " fraction of mip = " << _fraction_of_mip[layer] << " layer " << layer << endl;
      }	  
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    } else {
      float minpath = pitch;
      if (length < minpath) minpath = length;
      if (thickness < minpath) minpath = thickness;
	
      // Si MIP energy = 3.876 MeV / cm
      float threshold = 0.0;      
      if (_fraction_of_mip.find(layer) != _fraction_of_mip.end()) {	
	threshold = _fraction_of_mip[layer]*0.003876*minpath;  
      }	      
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
      if(Verbosity() >2) 
	cout << " not using thickness: threshold = " << threshold << " thickness = " << thickness << " fraction of mip = " << _fraction_of_mip[layer] << " layer " << layer << endl;
    }
  }
  
  return;
}

void PHG4SvtxThresholds::CalculateLadderThresholds(PHCompositeNode* topNode) {

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");

  if (!geom_container) return;
  
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
      float threshold = 0.0;
      if (_fraction_of_mip.find(layer) != _fraction_of_mip.end()) {
	threshold = _fraction_of_mip[layer]*0.003876*thickness;
      }
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    } else {
      float minpath = pitch;
      if (length < minpath) minpath = length;
      if (thickness < minpath) minpath = thickness;
      
      // Si MIP energy = 3.876 MeV / cm
      float threshold = 0.0;
      if (_fraction_of_mip.find(layer) != _fraction_of_mip.end()) {
	threshold = _fraction_of_mip[layer]*0.003876*minpath;
      }
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    }
  }
  
  return;
}

void PHG4SvtxThresholds::CalculateMapsLadderThresholds(PHCompositeNode* topNode) {

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_MAPS");

  if (!geom_container) return;
  
  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_pixel_thickness();
    float pitch = (layeriter->second)->get_pixel_x();
    float length = (layeriter->second)->get_pixel_z();

    if (get_use_thickness_mip(layer)) {
      // Si MIP energy = 3.876 MeV / cm
      float threshold = 0.0;
      if (_fraction_of_mip.find(layer) != _fraction_of_mip.end()) {
	threshold = _fraction_of_mip[layer]*0.003876*thickness;
	if(Verbosity() >2) 
	  cout << " using thickness: threshold = " << threshold << " thickness = " << thickness << " fraction of mip = " << _fraction_of_mip[layer] << " layer " << layer << endl;
      }
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
    } else {
      float minpath = pitch;
      if (length < minpath) minpath = length;
      if (thickness < minpath) minpath = thickness;
      
      // Si MIP energy = 3.876 MeV / cm
      float threshold = 0.0;
      if (_fraction_of_mip.find(layer) != _fraction_of_mip.end()) {
	threshold = _fraction_of_mip[layer]*0.003876*minpath;
      }
      _thresholds_by_layer.insert(std::make_pair(layer,threshold));
	if(Verbosity() >2) 
	  cout << " not using thickness: threshold = " << threshold << " thickness = " << thickness << " fraction of mip = " << _fraction_of_mip[layer] << " layer " << layer << endl;
    }
  }
  
  return;
}
