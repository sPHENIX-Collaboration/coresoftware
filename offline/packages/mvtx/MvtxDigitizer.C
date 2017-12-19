#include "MvtxDigitizer.h"

#include <tracker/TrackerDefs.h>
#include <tracker/TrackerHit.h>
#include <tracker/TrackerHitv1.h>
#include <tracker/TrackerHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>

#include <iostream>
#include <cmath>
#include <bitset>

using namespace std;

MvtxDigitizer::MvtxDigitizer(const string &name) :
  SubsysReco(name),
  _hitmap(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) {

  if (verbosity > 0)
    cout << PHWHERE << "Creating MvtxDigitizer with name = " << name << endl;
}

int MvtxDigitizer::InitRun(PHCompositeNode* topNode) {

  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVX node if required
  PHCompositeNode* svxNode
    = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "SVTX"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }

  // Create the Hit node if required
  TrackerHitContainer *mvtxhits =
    findNode::getClass<TrackerHitContainer>(dstNode, "TrackerHitContainer");
  if (!mvtxhits) {
    mvtxhits = new TrackerHitContainer();
    PHIODataNode<PHObject> *TrackerHitContainerNode =
      new PHIODataNode<PHObject>(mvtxhits, "TrackerHitContainer", "PHObject");
    svxNode->addNode(TrackerHitContainerNode);
  }

  CalculateMapsLadderThresholds(topNode);

  //----------------
  // Report Settings
  //----------------

  if (verbosity > 0) {
    cout << "====================== MvtxDigitizer::InitRun() =====================" << endl;
    for (std::map<int, float>::iterator iter = _fraction_of_mip.begin();
         iter != _fraction_of_mip.end();
         ++iter) {
      cout << " Fraction of expected MIP energy for Layer #" << iter->first << ": ";
      cout << iter->second;
      cout << endl;
    }
    for (std::map<int, float>::iterator iter = _thresholds_by_layer.begin();
         iter != _thresholds_by_layer.end();
         ++iter) {
      cout << " Cell Threshold in Layer #" << iter->first << " = " << 1.0e6 * iter->second
           << " keV based on Short-Axis Penetration" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxDigitizer::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();

  _hitmap = findNode::getClass<TrackerHitContainer>(topNode, "TrackerHitContainer");
  if (!_hitmap)
  {
    cout << PHWHERE << " ERROR: Can't find TrackerHitContainer." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _hitmap->Reset();

  DigitizeMapsLadderCells(topNode);

  PrintHits(topNode);

  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxDigitizer::CalculateMapsLadderThresholds(PHCompositeNode* topNode)
{

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MAPS");

  if (!geom_container) return;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {

    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_pixel_thickness();
    float pitch = (layeriter->second)->get_pixel_x();
    float length = (layeriter->second)->get_pixel_z();

    float minpath = pitch;
    if (length < minpath) minpath = length;
    if (thickness < minpath) minpath = thickness;

    // Si MIP energy = 3.876 MeV / cm
    float threshold = 0.0;
    if (_fraction_of_mip.find(layer) != _fraction_of_mip.end())
    {
      threshold = _fraction_of_mip[layer] * 0.003876 * minpath;
    }
    _thresholds_by_layer.insert(std::make_pair(layer, threshold));
    if (verbosity > 2)
    {
      cout << PHWHERE
           << " not using thickness:"
           << " layer " << layer
           << " threshold = " << threshold
           << " thickness = " << thickness
           << " fraction of mip = " << _fraction_of_mip[layer]
           << endl;
    }
  }

  return;
}

void MvtxDigitizer::DigitizeMapsLadderCells(PHCompositeNode *topNode)
{

  //----------
  // Get Nodes
  //----------

  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_MAPS");
  if (!cells) 
  {
    cout << PHWHERE << " - Unable to find node G4CELL_MAPS. Aborting" << endl;
    return;
  }

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MAPS");
  if (!geom_container)
  {
    cout << PHWHERE << " - Unable to find node CYLINDERGEOM_MAPS. Aborting" << endl;
    return;
  } 


  //-------------
  // Digitization
  //-------------

  vector<PHG4Cell*> cell_list;
  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for (PHG4CellContainer::ConstIterator celliter = cellrange.first;
       celliter != cellrange.second;
       ++celliter)
  {

    PHG4Cell* cell = celliter->second;

    // check that the cell energy is above threshold
    if ( cell->get_edep() < get_threshold_by_layer(cell->get_layer()) )
      continue;

    // get the layer geometry helper for indexing
    PHG4CylinderGeom_MAPS *geom = (PHG4CylinderGeom_MAPS*) geom_container->GetLayerGeom(cell->get_layer());

    int pixel_number = cell->get_pixel_index();
    // binphi is the cell index in the phi direction in the sensor
    unsigned short binphi = geom->get_pixel_X_from_pixel_number(pixel_number);
    // binz is the cell index in the z direction in the sensor
    unsigned short binz = geom->get_pixel_Z_from_pixel_number(pixel_number);

    // Make a hit
    TrackerHitv1 * hit = new TrackerHitv1();

    TrackerDefs::keytype key =
      TrackerDefs::MVTXBinning::genhitkey(TrackerDefs::TRACKERID::mvtx_id,
                                           cell->get_layer(),
                                           cell->get_stave_index(),
                                           cell->get_chip_index(),
                                           binphi, binz);


    if ( verbosity > 2 )
    {
      cout << PHWHERE << " - new hit: " 
           << endl
           << "       input:  "
           << " detid:" << TrackerDefs::TRACKERID::mvtx_id
           << " layer:" << cell->get_layer()
           << " ladder:" << cell->get_stave_index()
           << " chip:" << cell->get_chip_index()
           << " pixel:" << pixel_number
           << " row:" << binphi
           << " col:" << binz
           << " key:0x" << hex << key << dec
           << endl
           << "       output: "
           << " detid:" << int(TrackerDefs::get_trackerid(key))
           << " layer:" << int(TrackerDefs::get_layer(key))
           << " ladder:" << int(TrackerDefs::MVTXBinning::get_ladder(key))
           << " chip:" << int(TrackerDefs::MVTXBinning::get_chip(key))
           << " pixel:" << pixel_number
           << " row:" << TrackerDefs::MVTXBinning::get_row(key)
           << " col:" << TrackerDefs::MVTXBinning::get_col(key)
           << " key:0x" << hex << key << dec
           << endl;
      TrackerDefs::print_bits(key);
    }

    hit->set_hitid(key);

    // copy g4hits from cell to TrackerHit
    if ( verbosity > 2 )
    {
      cout << PHWHERE << " - adding g4 hits to TrackerHit" << endl;
    }
    PHG4Cell::EdepConstRange g4hitrange = cell->get_g4hits();
    for (PHG4Cell::EdepConstIterator hititr = g4hitrange.first;
         hititr != g4hitrange.second;
         ++hititr)
    {
      if ( verbosity > 2 )
      {
        cout << PHWHERE << "        g4hit e:" << hititr->second << " key:0x" << hex << hititr->first << endl;
      }
      hit->add_edep(hititr->first, hititr->second);
    }

    _hitmap->AddHitSpecifyKey(key, hit);
    if ( verbosity >= 1 )
      cout << PHWHERE << " - Added new TrackerHit with key 0x" << hex << key << endl;
  }

  return;
}

void MvtxDigitizer::PrintHits(PHCompositeNode *topNode) {

  if (verbosity >= 1) {

    cout << "================= MvtxDigitizer::process_event() ====================" << endl;

    TrackerHitContainer *hitlist = findNode::getClass<TrackerHitContainer>(topNode, "TrackerHitContainer");
    if (!hitlist) return;

    TrackerHitContainer::ConstRange mvtxhitrange = hitlist->getHits(TrackerDefs::TRACKERID::mvtx_id);



    cout << " Found and recorded the following hits: " << endl;

    unsigned int ihit = 0;
    for (TrackerHitContainer::ConstIterator iter = mvtxhitrange.first;
         iter != mvtxhitrange.second;
         ++iter)
    {

      TrackerHit* hit = iter->second;
      cout << ihit << endl;
      hit->identify();
      ++ihit;
    }

    cout << "===========================================================================" << endl;
  }

  return;
}
