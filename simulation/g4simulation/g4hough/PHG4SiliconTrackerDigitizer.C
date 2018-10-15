#include "PHG4SiliconTrackerDigitizer.h"

#include "SvtxHitMap.h"
#include "SvtxHitMap_v1.h"
#include "SvtxHit.h"
#include "SvtxHit_v1.h"
#include "SvtxDeadMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>

#include <iostream>
#include <cmath>

using namespace std;

PHG4SiliconTrackerDigitizer::PHG4SiliconTrackerDigitizer(const string &name) :
  SubsysReco(name),
  _hitmap(NULL),
  _timer(PHTimeServer::get()->insert_new(name)),
  m_nCells(0),
  m_nDeadCells(0)
{
}

int PHG4SiliconTrackerDigitizer::InitRun(PHCompositeNode* topNode) {

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
    svxhits = new SvtxHitMap_v1();
    PHIODataNode<PHObject> *SvtxHitMapNode =
      new PHIODataNode<PHObject>(svxhits, "SvtxHitMap", "PHObject");
    svxNode->addNode(SvtxHitMapNode);
  }

  CalculateLadderCellADCScale(topNode);
  
  //----------------
  // Report Settings
  //----------------
  
  if (Verbosity() > 0) {
    cout << "====================== PHG4SiliconTrackerDigitizer::InitRun() =====================" << endl;
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

int PHG4SiliconTrackerDigitizer::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();

  _hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hitmap) 
    {
      cout << PHWHERE << " ERROR: Can't find SvtxHitMap." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _hitmap->Reset();
  
  DigitizeLadderCells(topNode);

  PrintHits(topNode);
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4SiliconTrackerDigitizer::CalculateLadderCellADCScale(PHCompositeNode *topNode) {

  // FPHX 3-bit ADC, thresholds are set in "set_fphx_adc_scale".

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SILICON_TRACKER");
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");
    
  if (!geom_container || !cells) return;
  
  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();
    if (_max_fphx_adc.find(layer)==_max_fphx_adc.end())
      assert(!"Error: _max_fphx_adc is not available.");

    float thickness = (layeriter->second)->get_thickness(); // cm
    float mip_e     = 0.003876 *thickness; // GeV
    _energy_scale.insert(std::make_pair(layer, mip_e));
 } 

  return;
}

void PHG4SiliconTrackerDigitizer::DigitizeLadderCells(PHCompositeNode *topNode) {

  //----------
  // Get Nodes
  //----------
 
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SILICON_TRACKER");
  if (!cells) return;

  const SvtxDeadMap *deadmap = findNode::getClass<SvtxDeadMap>(topNode, "DEADMAP_SILICON_TRACKER");
  if (Verbosity()>=VERBOSITY_MORE)
  {
    if (deadmap)
    {
      cout << "PHG4SiliconTrackerDigitizer::DigitizeLadderCells - Use deadmap ";
      deadmap->identify();
    }
    else
    {
      cout << "PHG4SiliconTrackerDigitizer::DigitizeLadderCells - Can not find deadmap, all channels enabled "<<endl;
    }
  }

  //-------------
  // Digitization
  //-------------

  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for(PHG4CellContainer::ConstIterator celliter = cellrange.first;
      celliter != cellrange.second;
      ++celliter) {
    
    PHG4Cell* cell = celliter->second;

    ++m_nCells;
    if (deadmap)
    {
      if (deadmap->isDeadChannelINTT(
              cell->get_layer(),             //const int layer,
              cell->get_ladder_phi_index(),  //const int ladder_phi,
              cell->get_ladder_z_index(),    //const int ladder_z,
              cell->get_zbin(),              //const int strip_z,
              cell->get_phibin()             //const int strip_phi
              ))
      {
        ++m_nDeadCells;
        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << "PHG4SiliconTrackerDigitizer::DigitizeLadderCells - dead cell at layer "<<cell->get_layer() <<": ";
          cell->identify();
        }
        continue;
      }
    }  //    if (deadmap)

    SvtxHit_v1 hit;

    const int layer = cell->get_layer();

    hit.set_layer(layer);
    hit.set_cellid(cell->get_cellid());

    if (_energy_scale.count(layer)>1)
      assert(!"Error: _energy_scale has two or more keys.");

    const float mip_e = _energy_scale[layer];

    std::vector< std::pair<double, double> > vadcrange = _max_fphx_adc[layer];

    int adc = -1;
    for (unsigned int irange=0; irange<vadcrange.size(); ++irange)
      if (cell->get_edep()>=vadcrange[irange].first*(double)mip_e && cell->get_edep()<vadcrange[irange].second*(double)mip_e)
	adc = (int)irange;

    if (adc<0) // TODO, underflow is temporarily assigned to ADC=0.
      adc = 0;

    float e = 0.0;
    if (adc>=0 && adc<int(vadcrange.size())-1)
      e = 0.5*(vadcrange[adc].second - vadcrange[adc].first)*mip_e;
    else if (adc==int(vadcrange.size())-1) // overflow
      e = vadcrange[adc].first*mip_e;
    else // underflow
      e = 0.5*vadcrange[0].first*mip_e;
    
    hit.set_adc(adc);
    hit.set_e(e);
        
    SvtxHit* ptr = _hitmap->insert(&hit);      
    if (!ptr->isValid()) {
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

//! end of process
int PHG4SiliconTrackerDigitizer::End(PHCompositeNode *topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "PHG4SiliconTrackerDigitizer::End - processed "
        << m_nCells << " cell with "
        << m_nDeadCells << " dead cells masked"
        <<" ("<<100.*m_nDeadCells/m_nCells<<"%)"<< endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4SiliconTrackerDigitizer::PrintHits(PHCompositeNode *topNode) {

  if (Verbosity() >= VERBOSITY_EVEN_MORE) {

    SvtxHitMap *hitlist = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
    if (!hitlist) return;
    
    cout << "================= PHG4SiliconTrackerDigitizer::process_event() ====================" << endl;
  

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
