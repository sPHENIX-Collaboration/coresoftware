#include "G4CellNtuple.h"

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>  // for get_phibin
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>

#include <cassert>
#include <cmath>     // for isfinite
#include <iostream>  // for operator<<
#include <sstream>
#include <utility>  // for pair

using namespace std;

G4CellNtuple::G4CellNtuple(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , nblocks(0)
  , hm(nullptr)
  , _filename(filename)
  , ntup(nullptr)
  , outfile(nullptr)
{
}

G4CellNtuple::~G4CellNtuple()
{
  //  delete ntup;
  delete hm;
}

int G4CellNtuple::Init(PHCompositeNode *)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(_filename.c_str(), "RECREATE");
  ntup = new TNtuple("cellntup", "G4Cells", "detid:layer:phi:eta:edep");
  //  ntup->SetDirectory(0);
  TH1 *h1 = new TH1F("edep1GeV", "edep 0-1GeV", 1000, 0, 1);
  eloss.push_back(h1);
  h1 = new TH1F("edep100GeV", "edep 0-100GeV", 1000, 0, 100);
  eloss.push_back(h1);
  return 0;
}

int G4CellNtuple::process_event(PHCompositeNode *topNode)
{
  ostringstream nodename;
  ostringstream geonodename;
  set<string>::const_iterator iter;
  vector<TH1 *>::const_iterator eiter;
  for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
  {
    int detid = (_detid.find(*iter))->second;
    nodename.str("");
    nodename << "G4CELL_" << *iter;
    geonodename.str("");
    geonodename << "CYLINDERCELLGEOM_" << *iter;
    PHG4TpcCylinderGeomContainer *cellgeos = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, geonodename.str());
    if (!cellgeos)
    {
      cout << "no geometry node " << geonodename.str() << " for " << *iter << endl;
      continue;
    }
    PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, nodename.str());
    if (cells)
    {
      int previouslayer = -1;
      PHG4TpcCylinderGeom *cellgeom = nullptr;
      double esum = 0;
      //          double numcells = cells->size();
      //          ncells[i]->Fill(numcells);
      //	  cout << "number of cells: " << cells->size() << endl;
      PHG4CellContainer::ConstRange cell_range = cells->getCells();
      for (PHG4CellContainer::ConstIterator cell_iter = cell_range.first; cell_iter != cell_range.second; cell_iter++)

      {
        double edep = cell_iter->second->get_edep();
        if (!isfinite(edep))
        {
          cout << "invalid edep: " << edep << endl;
        }
        esum += cell_iter->second->get_edep();
        int phibin = ~0x0;
        int etabin = ~0x0;
        double phi = NAN;
        double eta = NAN;
        int layer = cell_iter->second->get_layer();
        // to search the map fewer times, cache the geom object until the layer changes
        if (layer != previouslayer)
        {
          cellgeom = cellgeos->GetLayerCellGeom(layer);
          previouslayer = layer;
        }
        if (cell_iter->second->has_binning(PHG4CellDefs::etaphibinning))
        {
          phibin = PHG4CellDefs::EtaPhiBinning::get_phibin(cell_iter->second->get_cellid());
          etabin = PHG4CellDefs::EtaPhiBinning::get_etabin(cell_iter->second->get_cellid());
          phi = cellgeom->get_phicenter(phibin);
          eta = cellgeom->get_etacenter(etabin);
        }
        else if (cell_iter->second->has_binning(PHG4CellDefs::sizebinning))
        {
          phibin = PHG4CellDefs::SizeBinning::get_phibin(cell_iter->second->get_cellid());
          etabin = PHG4CellDefs::SizeBinning::get_zbin(cell_iter->second->get_cellid());
          phi = cellgeom->get_phicenter(phibin);
          eta = cellgeom->get_zcenter(etabin);
        }
        assert(cellgeom != nullptr);
        ntup->Fill(detid,
                   cell_iter->second->get_layer(),
                   phi,
                   eta,
                   cell_iter->second->get_edep());
      }
      for (eiter = eloss.begin(); eiter != eloss.end(); ++eiter)
      {
        (*eiter)->Fill(esum);
      }
    }
  }
  return 0;
}

int G4CellNtuple::End(PHCompositeNode */*topNode*/)
{
  outfile->cd();
  ntup->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(_filename, "UPDATE");
  return 0;
}

void G4CellNtuple::AddNode(const std::string &name, const int detid)
{
  _node_postfix.insert(name);
  _detid[name] = detid;
  return;
}
