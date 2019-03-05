#include "PHG4MVTXDigitizer.h"

#include <g4main/PHG4Hit.h>

#include <trackbase_historic/SvtxHit.h>
#include <trackbase_historic/SvtxHitMap.h>
#include <trackbase_historic/SvtxHitMap_v1.h>
#include <trackbase_historic/SvtxHit_v1.h>

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
  , _energy_scale(0.95e-06)   // 99 electrons default
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

    cout << PHWHERE << " MVTX ADC threshold is " <<  _energy_scale << endl;

  }

  //CalculateMVTXLadderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== PHG4MVTXDigitizer::InitRun() =====================" << endl;
    cout << " Energy per ADC in MVTX " << _energy_scale << " keV" << endl;
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

  return;
}

void PHG4MVTXDigitizer::DigitizeMVTXLadderCells(PHCompositeNode *topNode)
{
  //----------
  // Get Nodes
  //----------

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_MVTX");
  if (!cells) return;

  //-------------
  // Digitization
  //-------------

  vector<PHG4Cell *> cell_list;
  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for (PHG4CellContainer::ConstIterator celliter = cellrange.first;
       celliter != cellrange.second;
       ++celliter)
  {
    PHG4Cell *cell = celliter->second;

    unsigned int adc = cell->get_edep() / _energy_scale;
    if(adc == 0) continue;

    float e = cell->get_edep();

    SvtxHit_v1 hit;

    hit.set_layer(cell->get_layer());
    hit.set_cellid(cell->get_cellid());


    hit.set_adc(adc);
    hit.set_e(e);

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
