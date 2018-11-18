/**
 * @file g4mvtx/PHG4MvtxDigitizer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of PHG4MvtxDigitizer
 */
#include "PHG4MvtxDigitizer.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <mvtx/MvtxDefs.h>
#include <mvtx/MvtxHit.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom_MVTX.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <iostream>
#include <cmath>

using namespace std;

PHG4MvtxDigitizer::PHG4MvtxDigitizer(const string &name)
  : SubsysReco(name)
  , m_hitsets(nullptr)
{

  if(verbosity > 0)
    cout << "Creating PHG4MvtxDigitizer with name = " << name << endl;
}

int PHG4MvtxDigitizer::InitRun(PHCompositeNode* topNode) 
{

  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode 
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) 
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);
    
  // Create the TRKR node if required
  PHCompositeNode* trkrNode 
    = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode","TRKR"));
  if (!trkrNode) 
  {
    trkrNode = new PHCompositeNode("TRKR");
    dstNode->addNode(trkrNode);
  }
  
  // Create the Hit node if required
  TrkrHitSetContainer *trkrhitsets = findNode::getClass<TrkrHitSetContainer>(dstNode,"TrkrHitSetContainer");
  if ( !trkrhitsets ) 
  {
    trkrhitsets = new TrkrHitSetContainer();
    
    PHIODataNode<PHObject> *TrkrHitSetContainerNode =
      new PHIODataNode<PHObject>(trkrhitsets, "TrkrHitSetContainer", "PHObject");
    trkrNode->addNode(TrkrHitSetContainerNode);
  }
  
  CalculateADCScale(topNode);
  
  //----------------
  // Report Settings
  //----------------
  
  if (verbosity > 0) 
  {
    cout << "====================== PHG4MvtxDigitizer::InitRun() =====================" << endl;
    for (std::map<int,unsigned int>::iterator iter = m_maxAdc.begin();
	 iter != m_maxAdc.end();
	 ++iter) 
    {
      cout << " Max ADC in Layer #" << iter->first << " = " << iter->second << endl;
    }
    for (std::map<int,float>::iterator iter = m_energyScale.begin();
	 iter != m_energyScale.end();
	 ++iter) 
    {
      cout << " Energy per ADC in Layer #" << iter->first << " = " << 1.0e6*iter->second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4MvtxDigitizer::process_event(PHCompositeNode *topNode) 
{
  m_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode,"TrkrHitSetContainer");
  if (!m_hitsets) 
    {
      cout << PHWHERE << " ERROR: Can't find SvtxHitMap." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  m_hitsets->Reset();
  
  DigitizeCells(topNode);

  PrintHits(topNode);
  
  return Fun4AllReturnCodes::EVENT_OK;
}


void PHG4MvtxDigitizer::CalculateADCScale(PHCompositeNode *topNode) 
{

  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_MVTX");
    
  if (!geom_container) return;

  if(Verbosity())
    cout << "Found CYLINDERGEOM_MVTX node" << endl;
  
  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) 
  {
    
    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_pixel_thickness();
    float pitch = (layeriter->second)->get_pixel_x();
    float length = (layeriter->second)->get_pixel_z();
    
    float minpath = pitch;
    if (length < minpath) 
    {
      minpath = length;
    }
    if (thickness < minpath) 
    {
      minpath = thickness;
    }
    float mip_e = 0.003876*minpath;  
    
    if (Verbosity())
    {
      cout << "mip_e = " << mip_e << endl;
    }

    if (m_maxAdc.find(layer) == m_maxAdc.end()) 
    {
      m_maxAdc[layer] = 255;
      m_energyScale[layer] = mip_e / 64;
    }
  }    
  
  return;
}


void PHG4MvtxDigitizer::DigitizeCells(PHCompositeNode *topNode) 
{

  //----------
  // Get Nodes
  //----------

  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_MVTX");
  if (!geom_container) return;
 
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_MVTX");
  if (!cells) return; 
  
  //-------------
  // Digitization
  //-------------

  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for(PHG4CellContainer::ConstIterator celliter = cellrange.first;
      celliter != cellrange.second;
      ++celliter) 
  {
    
    PHG4Cell* cell = celliter->second;

    // get cell indeces
    int lyr = cell->get_layer();
    int stave = cell->get_stave_index();
    int chip = cell->get_chip_index();
    int pxl = cell->get_pixel_index();

    PHG4CylinderGeom_MVTX *geom = (PHG4CylinderGeom_MVTX*) geom_container->GetLayerGeom(lyr);
    int col = geom->get_pixel_Z_from_pixel_number(pxl);
    int row = geom->get_pixel_X_from_pixel_number(pxl);

    // get the TrkrHitSet for this hit
    // (this is inefficient)
    TrkrDefs::hitsetkey hskey = MvtxDefs::genHitSetKey(lyr,stave,chip);
    TrkrHitSet* hset = (m_hitsets->findOrAddHitSet(hskey))->second;

    // make and add the hit to the hitset
    TrkrDefs::hitkey hkey = MvtxDefs::genHitKey(col, row);
    
    MvtxHit* mhit = new MvtxHit();

    hset->addHitSpecificKey(hkey, mhit);

  }
  
  return;
}

void PHG4MvtxDigitizer::PrintHits(PHCompositeNode *topNode) 
{

  if (verbosity >= 1) 
  {
    cout << "================= PHG4MvtxDigitizer::process_event() ====================" << endl;

    TrkrHitSetContainer::ConstRange hits = m_hitsets->getHitSets(TrkrDefs::TrkrId::mvtxId);
    for ( auto itr = hits.first; itr != hits.second; ++itr)
    {
      cout << " key:" << itr->first << endl;
      (itr->second)->identify();
    }
    
    cout << "===========================================================================" << endl;
  }
  
  return;
}
