#include "PHG4MVTXDigitizer.h"

#include <g4main/PHG4Hit.h>

#include <trackbase_historic/SvtxHit.h>
#include <trackbase_historic/SvtxHitMap.h>
#include <trackbase_historic/SvtxHitMap_v1.h>
#include <trackbase_historic/SvtxHit_v1.h>

// Move to new storage containers
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <mvtx/MvtxDefs.h>
#include <mvtx/MvtxHit.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>
#include <g4detectors/PHG4Cellv1.h>
#include <g4detectors/PHG4Cellv2.h>
#include <phool/PHRandomSeed.h>

#include <gsl/gsl_randist.h>

#include <cmath>
#include <iostream>
#include <limits>

using namespace std;

PHG4MVTXDigitizer::PHG4MVTXDigitizer(const string &name)
  : SubsysReco(name)
  , _hitmap(NULL)
  , _timer(PHTimeServer::get()->insert_new(name))
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  cout << Name() << " random seed: " << seed << endl;
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, seed);

  if (Verbosity() > 0)
    cout << "Creating PHG4MVTXDigitizer with name = " << name << endl;
}

int PHG4MVTXDigitizer::InitRun(PHCompositeNode *topNode)
{
  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVX node if required
  PHCompositeNode *svxNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", "SVTX"));
  if (!svxNode)
  {
    cout << PHWHERE << "Creating SVTX node " << endl;
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }

  // Create the Hit node if required
  SvtxHitMap *svxhits = findNode::getClass<SvtxHitMap>(dstNode, "SvtxHitMap");
  if (!svxhits)
  {
    cout << PHWHERE << "Creating SvtxHitMap node " << endl;
    svxhits = new SvtxHitMap_v1();
    PHIODataNode<PHObject> *SvtxHitMapNode =
        new PHIODataNode<PHObject>(svxhits, "SvtxHitMap", "PHObject");
    svxNode->addNode(SvtxHitMapNode);
  }

  CalculateMVTXLadderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== PHG4MVTXDigitizer::InitRun() =====================" << endl;
    for (std::map<int, unsigned int>::iterator iter = _max_adc.begin();
         iter != _max_adc.end();
         ++iter)
    {
      cout << " Max ADC in Layer #" << iter->first << " = " << iter->second << endl;
    }
    for (std::map<int, float>::iterator iter = _energy_scale.begin();
         iter != _energy_scale.end();
         ++iter)
    {
      cout << " Energy per ADC in Layer #" << iter->first << " = " << 1.0e6 * iter->second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4MVTXDigitizer::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();

  _hitmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!_hitmap)
  {
    cout << PHWHERE << " ERROR: Can't find SvtxHitMap." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //Jin: don't clear up node. Fun4all server does that. Extra cleaning usually cause problems
  //  _hitmap->Reset();

  // This code now only does the MVTX
  DigitizeMVTXLadderCells(topNode);

  PrintHits(topNode);

  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4MVTXDigitizer::CalculateMVTXLadderCellADCScale(PHCompositeNode *topNode)
{
  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  if (!geom_container) return;

  if (Verbosity())
    cout << "Found CYLINDERGEOM_MVTX node" << endl;

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
    float mip_e = 0.003876 * minpath;

    if (Verbosity())
      cout << "mip_e = " << mip_e << endl;

    if (_max_adc.find(layer) == _max_adc.end())
    {
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }

  return;
}

void PHG4MVTXDigitizer::DigitizeMVTXLadderCells(PHCompositeNode *topNode)
{
  //----------
  // Get Nodes
  //----------

  // old containers
  //============
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_MVTX");
  if (!cells) return;

  // Digitization

  vector<PHG4Cell *> cell_list;
  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for (PHG4CellContainer::ConstIterator celliter = cellrange.first;
       celliter != cellrange.second;
       ++celliter)
  {
    PHG4Cell *cell = celliter->second;

    SvtxHit_v1 hit;

    hit.set_layer(cell->get_layer());
    hit.set_cellid(cell->get_cellid());

    unsigned int adc = cell->get_edep() / _energy_scale[hit.get_layer()];
    if (adc > _max_adc[hit.get_layer()]) adc = _max_adc[hit.get_layer()];
    float e = _energy_scale[hit.get_layer()] * adc;

    hit.set_adc(adc);
    hit.set_e(e);

    cout << "    OLD: PHG4MVTXDigitizer: found cell in layer " << hit.get_layer() << " with signal " << cell->get_edep() << " and adc " << adc << endl;

    SvtxHit *ptr = _hitmap->insert(&hit);
    if (!ptr->isValid())
    {
      static bool first = true;
      if (first)
      {
        cout << PHWHERE << "ERROR: Incomplete SvtxHits are being created" << endl;
        ptr->identify();
        first = false;
      }
    }
  }
  // end old containers  
  //===============

  // new containers
  //=============
  // Get the TrkrHitSetContainer node
  TrkrHitSetContainer *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!trkrhitsetcontainer)
    {
      cout << "Could not locate TRKR_HITSET node, quit! " << endl;
      exit(1);
    }

  // Digitization

  // We want all hitsets for the MVTX
  TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
       hitset_iter != hitset_range.second;
       ++hitset_iter)
    {
      // we have an itrator to one TrkrHitSet for the mvtx from the trkrHitSetContainer
      // get the hitset key so we can find the layer
      TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      int layer = TrkrDefs::getLayer(hitsetkey);
      cout << "PHG4MVTXDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << endl;

      // get all of the hits from this hitset      
      TrkrHitSet *hitset = hitset_iter->second;
      TrkrHitSet::ConstRange hit_range = hitset->getHits();
      for(TrkrHitSet::ConstIterator hit_iter = hit_range.first;
	  hit_iter != hit_range.second;
	  ++hit_iter)
	{
	  TrkrHit *hit = hit_iter->second;
      
	  // Convert the signal value to an ADC value and write that to the hit	  
	  unsigned int adc =   hit->getEnergy() / _energy_scale[layer];
	  if (adc > _max_adc[layer]) adc = _max_adc[layer];
	  hit->setAdc(adc);

	  cout << "    PHG4MVTXDigitizer: found hit with key: " << hit_iter->first << " and signal " << hit->getEnergy() << " and adc " << adc << endl;
	}
    }

  
  // end new containers  
  //===============
  
  return;
}

void PHG4MVTXDigitizer::PrintHits(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
  {
    SvtxHitMap *hitlist = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
    if (!hitlist) return;

    cout << "================= PHG4MVTXDigitizer::process_event() ====================" << endl;

    cout << " Found and recorded the following " << hitlist->size() << " hits: " << endl;

    unsigned int ihit = 0;
    for (SvtxHitMap::Iter iter = hitlist->begin();
         iter != hitlist->end();
         ++iter)
    {
      SvtxHit *hit = iter->second;
      cout << ihit << " of " << hitlist->size() << endl;
      hit->identify();
      ++ihit;
    }

    cout << "===========================================================================" << endl;
  }

  return;
}
