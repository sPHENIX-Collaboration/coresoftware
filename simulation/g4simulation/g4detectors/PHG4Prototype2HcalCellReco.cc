#include "PHG4Prototype2HcalCellReco.h"
#include "PHG4ScintillatorSlatContainer.h"
#include "PHG4ScintillatorSlatv1.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>  // std::numeric_limits
#include <sstream>

using namespace std;
// for hcal dimension
// prototype outer hcal, 1st and 2nd prototype inner hcal, 3rd prototype inner hcal is 17 
#define ROWDIM 21 
// 12 scintillator tiles per row (index starting at 1)
#define COLUMNDIM 13
static array<array<PHG4ScintillatorSlat *, COLUMNDIM>, ROWDIM> slatarray = {nullptr};

PHG4Prototype2HcalCellReco::PHG4Prototype2HcalCellReco(const string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_CheckEnergyConservationFlag(0)
  , m_Tmin(NAN)
  , m_Tmax(NAN)
{
  InitializeParameters();
}

int PHG4Prototype2HcalCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    gSystem->Exit(1);
  }
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));

  string paramnodename = "G4CELLPARAM_" + m_Detector;
  string geonodename = "G4CELLGEO_" + m_Detector;

  m_HitNodeName = "G4HIT_" + m_Detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << m_HitNodeName << endl;
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Print("NODETREE");
    gSystem->Exit(1);
  }
  m_CellNodeName = "G4CELL_" + m_Detector;
  PHG4ScintillatorSlatContainer *slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, m_CellNodeName);
  if (!slats)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(m_Detector);
      dstNode->addNode(DetNode);
    }
    slats = new PHG4ScintillatorSlatContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(slats, m_CellNodeName, "PHObject");
    DetNode->addNode(newNode);
  }
  UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", m_Detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_Detector);
    runNode->addNode(RunDetNode);
  }
  SaveToNodeTree(RunDetNode, paramnodename);
  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", m_Detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(m_Detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);

  m_Tmin = get_double_param("tmin");
  m_Tmax = get_double_param("tmax");
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4Prototype2HcalCellReco::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << m_HitNodeName << endl;
    exit(1);
  }
  PHG4ScintillatorSlatContainer *slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, m_CellNodeName);
  if (!slats)
  {
    cout << "could not locate cell node " << m_CellNodeName << endl;
    exit(1);
  }

  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; ++layer)
  {
    PHG4HitContainer::ConstIterator hiter;
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);
    for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      if (hiter->second->get_t(0) > m_Tmax) continue;
      if (hiter->second->get_t(1) < m_Tmin) continue;
      short icolumn = hiter->second->get_scint_id();
      short irow = hiter->second->get_row();
      if (irow >= ROWDIM || irow < 0)
      {
        cout << "row " << irow
             << " exceed array size: " << ROWDIM
             << " adjust ROWDIM and recompile" << endl;
        gSystem->Exit(1);
      }

      if (icolumn >= COLUMNDIM || icolumn < 0)
      {
        cout << "column: " << icolumn
             << " exceed array size: " << COLUMNDIM
             << " adjust COLUMNDIM and recompile" << endl;
        gSystem->Exit(1);
      }

      if (!slatarray[irow][icolumn])
      {
        slatarray[irow][icolumn] = new PHG4ScintillatorSlatv1();
      }
      slatarray[irow][icolumn]->add_edep(hiter->second->get_edep(),
                                         hiter->second->get_eion(),
                                         hiter->second->get_light_yield());
      slatarray[irow][icolumn]->add_hit_key(hiter->first);
      // cout << "row: " << hiter->second->get_row()
      // 	   << ", column: " << hiter->second->get_scint_id() << endl;
      // checking ADC timing integration window cut
    }  // end loop over g4hits
    int nslathits = 0;
    for (int irow = 0; irow < ROWDIM; irow++)
    {
      for (int icolumn = 0; icolumn < COLUMNDIM; icolumn++)
      {
        if (slatarray[irow][icolumn])
        {
          PHG4ScintillatorSlatDefs::keytype key = PHG4ScintillatorSlatDefs::genkey(irow, icolumn);
          slats->AddScintillatorSlat(key, slatarray[irow][icolumn]);
          slatarray[irow][icolumn] = NULL;
          nslathits++;
        }
      }
    }
    if (Verbosity() > 0)
    {
      cout << Name() << ": found " << nslathits << " slats with energy deposition" << endl;
    }
  }

  if (m_CheckEnergyConservationFlag)
  {
    CheckEnergy(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4Prototype2HcalCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  PHG4ScintillatorSlatContainer *slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, m_CellNodeName);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    sum_energy_g4hit += hiter->second->get_edep();
  }
  PHG4ScintillatorSlatContainer::ConstRange cell_begin_end = slats->getScintillatorSlats();
  PHG4ScintillatorSlatContainer::ConstIterator citer;
  for (citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
  {
    sum_energy_cells += citer->second->get_edep();
  }
  // the fractional eloss for particles traversing eta bins leads to minute rounding errors
  if (fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit > 1e-6)
  {
    cout << "energy mismatch between cells: " << sum_energy_cells
         << " and hits: " << sum_energy_g4hit
         << " diff sum(cells) - sum(hits): " << sum_energy_cells - sum_energy_g4hit
         << endl;
    return -1;
  }
  else
  {
    if (Verbosity() > 0)
    {
      cout << Name() << ": total energy for this event: " << sum_energy_g4hit << " GeV" << endl;
    }
  }
  return 0;
}

void PHG4Prototype2HcalCellReco::SetDefaultParameters()
{
  set_default_double_param("tmax", 60.0);
  set_default_double_param("tmin", 0.0);  
  return;
}

void PHG4Prototype2HcalCellReco::set_timing_window(const double tmi, const double tma)
{
  set_double_param("tmin", tmi);
  set_double_param("tmax", tma);
}
