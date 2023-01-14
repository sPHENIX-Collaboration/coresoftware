#include "PHG4CylinderCellReco.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderGeom.h"
#include "PHG4CylinderGeomContainer.h"

#include "PHG4Cell.h"  // for PHG4Cell, PHG...
#include "PHG4CellContainer.h"
#include "PHG4CellDefs.h"
#include "PHG4Cellv1.h"

#include <phparameter/PHParameterContainerInterface.h>  // for PHParameterCo...
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <cmath>
#include <cstdlib>
#include <cstring>  // for memset
#include <iostream>
#include <iterator>  // for reverse_iterator
#include <sstream>
#include <vector>  // for vector

using namespace std;

PHG4CylinderCellReco::PHG4CylinderCellReco(const string &name)
  : SubsysReco(name)
  , PHParameterContainerInterface(name)
{
  SetDefaultParameters();
  memset(nbins, 0, sizeof(nbins));
}

int PHG4CylinderCellReco::ResetEvent(PHCompositeNode * /*topNode*/)
{
  sum_energy_before_cuts = 0.;
  sum_energy_g4hit = 0.;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4CylinderCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator nodeiter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(nodeiter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  PHCompositeNode *runNode;
  runNode = dynamic_cast<PHCompositeNode *>(nodeiter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    cout << Name() << "DST Node missing, doing nothing." << endl;
    exit(1);
  }
  PHNodeIterator runiter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runiter.findFirst("PHCompositeNode", detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(detector);
    runNode->addNode(RunDetNode);
  }

  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << hitnodename << endl;
    exit(1);
  }
  cellnodename = "G4CELL_" + outdetector;
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", detector));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }
    cells = new PHG4CellContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells, cellnodename, "PHObject");
    DetNode->addNode(newNode);
  }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename);
  if (!geo)
  {
    cout << "Could not locate geometry node " << geonodename << endl;
    exit(1);
  }
  if (Verbosity() > 0)
  {
    geo->identify();
  }
  seggeonodename = "CYLINDERCELLGEOM_" + outdetector;
  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename);
  if (!seggeo)
  {
    seggeo = new PHG4CylinderCellGeomContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename, "PHObject");
    runNode->addNode(newNode);
  }

  GetParamsContainerModify()->set_name(detector);

  UpdateParametersWithMacro();

  map<int, PHG4CylinderGeom *>::const_iterator miter;
  pair<map<int, PHG4CylinderGeom *>::const_iterator, map<int, PHG4CylinderGeom *>::const_iterator> begin_end = geo->get_begin_end();
  map<int, std::pair<double, double> >::iterator sizeiter;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
  {
    PHG4CylinderGeom *layergeom = miter->second;
    int layer = layergeom->get_layer();
    if (!ExistDetid(layer))
    {
      cout << Name() << ": No parameters for detid/layer " << layer
           << ", hits from this detid/layer will not be accumulated into cells" << endl;
      continue;
    }
    implemented_detid.insert(layer);
    set_size(layer, get_double_param(layer, "size_long"), get_double_param(layer, "size_perp"));
    tmin_max.insert(std::make_pair(layer, std::make_pair(get_double_param(layer, "tmin"), get_double_param(layer, "tmax"))));
    m_DeltaTMap.insert(std::make_pair(layer,get_double_param(layer, "delta_t")));
    double circumference = layergeom->get_radius() * 2 * M_PI;
    double length_in_z = layergeom->get_zmax() - layergeom->get_zmin();
    sizeiter = cell_size.find(layer);
    if (sizeiter == cell_size.end())
    {
      cout << "no cell sizes for layer " << layer << endl;
      exit(1);
    }
    // create geo object and fill with variables common to all binning methods
    PHG4CylinderCellGeom *layerseggeo = new PHG4CylinderCellGeom();
    layerseggeo->set_layer(layergeom->get_layer());
    layerseggeo->set_radius(layergeom->get_radius());
    layerseggeo->set_thickness(layergeom->get_thickness());
    if (binning[layer] == PHG4CellDefs::etaphibinning)
    {
      // calculate eta at radius+ thickness (outer radius)
      // length via eta coverage is calculated using the outer radius
      double etamin = PHG4Utils::get_eta(layergeom->get_radius() + layergeom->get_thickness(), layergeom->get_zmin());
      double etamax = PHG4Utils::get_eta(layergeom->get_radius() + layergeom->get_thickness(), layergeom->get_zmax());
      zmin_max[layer] = make_pair(etamin, etamax);
      double etastepsize = (sizeiter->second).first;
      double d_etabins;
      // if the eta cell size is larger than the eta range, make one bin
      if (etastepsize > etamax - etamin)
      {
        d_etabins = 1;
      }
      else
      {
        // it is unlikely that the eta range is a multiple of the eta cell size
        // then fract is 0, if not - add 1 bin which makes the
        // cells a tiny bit smaller but makes them fit
        double fract = modf((etamax - etamin) / etastepsize, &d_etabins);
        if (fract != 0)
        {
          d_etabins++;
        }
      }
      etastepsize = (etamax - etamin) / d_etabins;
      (sizeiter->second).first = etastepsize;
      int etabins = d_etabins;
      double etalow = etamin;
      double etahi = etalow + etastepsize;
      for (int i = 0; i < etabins; i++)
      {
        if (etahi > (etamax + 1.e-6))  // etahi is a tiny bit larger due to numerical uncertainties
        {
          cout << "etahi: " << etahi << ", etamax: " << etamax << endl;
        }
        etahi += etastepsize;
      }

      double phimin = -M_PI;
      double phimax = M_PI;
      double phistepsize = (sizeiter->second).second;
      double d_phibins;
      if (phistepsize >= phimax - phimin)
      {
        d_phibins = 1;
      }
      else
      {
        // it is unlikely that the phi range is a multiple of the phi cell size
        // then fract is 0, if not - add 1 bin which makes the
        // cells a tiny bit smaller but makes them fit
        double fract = modf((phimax - phimin) / phistepsize, &d_phibins);
        if (fract != 0)
        {
          d_phibins++;
        }
      }
      phistepsize = (phimax - phimin) / d_phibins;
      (sizeiter->second).second = phistepsize;
      int phibins = d_phibins;
      double philow = phimin;
      double phihi = philow + phistepsize;
      for (int i = 0; i < phibins; i++)
      {
        if (phihi > phimax)
        {
          cout << "phihi: " << phihi << ", phimax: " << phimax << endl;
        }
        phihi += phistepsize;
      }
      pair<int, int> phi_z_bin = make_pair(phibins, etabins);
      n_phi_z_bins[layer] = phi_z_bin;
      layerseggeo->set_binning(PHG4CellDefs::etaphibinning);
      layerseggeo->set_etabins(etabins);
      layerseggeo->set_etamin(etamin);
      layerseggeo->set_etastep(etastepsize);
      layerseggeo->set_phimin(phimin);
      layerseggeo->set_phibins(phibins);
      layerseggeo->set_phistep(phistepsize);
      phistep[layer] = phistepsize;
      etastep[layer] = etastepsize;
    }
    else if (binning[layer] == PHG4CellDefs::sizebinning)
    {
      zmin_max[layer] = make_pair(layergeom->get_zmin(), layergeom->get_zmax());
      double size_z = (sizeiter->second).second;
      double size_r = (sizeiter->second).first;
      double bins_r;
      // if the size is larger than circumference, make it one bin
      if (size_r >= circumference)
      {
        bins_r = 1;
      }
      else
      {
        // unlikely but if the circumference is a multiple of the cell size
        // use result of division, if not - add 1 bin which makes the
        // cells a tiny bit smaller but makes them fit
        double fract = modf(circumference / size_r, &bins_r);
        if (fract != 0)
        {
          bins_r++;
        }
      }
      nbins[0] = bins_r;
      size_r = circumference / bins_r;
      (sizeiter->second).first = size_r;
      double phistepsize = 2 * M_PI / bins_r;
      double phimin = -M_PI;
      double phimax = phimin + phistepsize;
      phistep[layer] = phistepsize;
      for (int i = 0; i < nbins[0]; i++)
      {
        if (phimax > (M_PI + 1e-9))
        {
          cout << "phimax: " << phimax << ", M_PI: " << M_PI
               << "phimax-M_PI: " << phimax - M_PI << endl;
        }
        phimax += phistepsize;
      }
      // if the size is larger than length, make it one bin
      if (size_z >= length_in_z)
      {
        bins_r = 1;
      }
      else
      {
        // unlikely but if the length is a multiple of the cell size
        // use result of division, if not - add 1 bin which makes the
        // cells a tiny bit smaller but makes them fit
        double fract = modf(length_in_z / size_z, &bins_r);
        if (fract != 0)
        {
          bins_r++;
        }
      }
      nbins[1] = bins_r;
      pair<int, int> phi_z_bin = make_pair(nbins[0], nbins[1]);
      n_phi_z_bins[layer] = phi_z_bin;
      // update our map with the new sizes
      size_z = length_in_z / bins_r;
      (sizeiter->second).second = size_z;
      double zlow = layergeom->get_zmin();
      double zhigh = zlow + size_z;
      ;
      for (int i = 0; i < nbins[1]; i++)
      {
        if (zhigh > (layergeom->get_zmax() + 1e-9))
        {
          cout << "zhigh: " << zhigh << ", zmax "
               << layergeom->get_zmax()
               << ", zhigh-zmax: " << zhigh - layergeom->get_zmax()
               << endl;
        }
        zhigh += size_z;
      }
      layerseggeo->set_binning(PHG4CellDefs::sizebinning);
      layerseggeo->set_zbins(nbins[1]);
      layerseggeo->set_zmin(layergeom->get_zmin());
      layerseggeo->set_zstep(size_z);
      layerseggeo->set_phibins(nbins[0]);
      layerseggeo->set_phistep(phistepsize);
    }
    // add geo object filled by different binning methods
    seggeo->AddLayerCellGeom(layerseggeo);
    if (Verbosity() > 1)
    {
      layerseggeo->identify();
    }
  }

  // print out settings
  if (Verbosity() > 0)
  {
    cout << "===================== PHG4CylinderCellReco::InitRun() =====================" << endl;
    cout << " " << outdetector << " Segmentation Description: " << endl;
    for (std::map<int, int>::iterator iter = binning.begin();
         iter != binning.end(); ++iter)
    {
      int layer = iter->first;

      if (binning[layer] == PHG4CellDefs::etaphibinning)
      {
        // phi & eta bin is usually used to make projective towers
        // so just print the first layer
        cout << " Layer #" << binning.begin()->first << "-" << binning.rbegin()->first << endl;
        cout << "   Nbins (phi,eta): (" << n_phi_z_bins[layer].first << ", " << n_phi_z_bins[layer].second << ")" << endl;
        cout << "   Cell Size (phi,eta): (" << cell_size[layer].first << " rad, " << cell_size[layer].second << " units)" << endl;
        break;
      }
      else if (binning[layer] == PHG4CellDefs::sizebinning)
      {
        cout << " Layer #" << layer << endl;
        cout << "   Nbins (phi,z): (" << n_phi_z_bins[layer].first << ", " << n_phi_z_bins[layer].second << ")" << endl;
        cout << "   Cell Size (phi,z): (" << cell_size[layer].first << " cm, " << cell_size[layer].second << " cm)" << endl;
      }
    }
    cout << "===========================================================================" << endl;
  }
  string nodename = "G4CELLPARAM_" + GetParamsContainer()->Name();
  SaveToNodeTree(RunDetNode, nodename);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4CylinderCellReco::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << hitnodename << endl;
    exit(1);
  }
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
  {
    cout << "could not locate cell node " << cellnodename << endl;
    exit(1);
  }

  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename);
  if (!seggeo)
  {
    cout << "could not locate geo node " << seggeonodename << endl;
    exit(1);
  }

  map<int, std::pair<double, double> >::iterator sizeiter;
  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  //   cout << "number of layers: " << g4hit->num_layers() << endl;
  //   cout << "number of hits: " << g4hit->size() << endl;
  //   for (layer = layer_begin_end.first; layer != layer_begin_end.second; layer++)
  //     {
  //       cout << "layer number: " << *layer << endl;
  //     }
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; layer++)
  {
    // only handle layers/detector ids which have parameters set
    if (implemented_detid.find(*layer) == implemented_detid.end())
    {
      continue;
    }
    PHG4HitContainer::ConstIterator hiter;
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);
    PHG4CylinderCellGeom *geo = seggeo->GetLayerCellGeom(*layer);
    int nphibins = n_phi_z_bins[*layer].first;
    int nzbins = n_phi_z_bins[*layer].second;

    // ------- eta/phi binning ------------------------------------------------------------------------
    if (binning[*layer] == PHG4CellDefs::etaphibinning)
    {
      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
      {
        sum_energy_before_cuts += hiter->second->get_edep();
        // checking ADC timing integration window cut
        if (hiter->second->get_t(0) > tmin_max[*layer].second) continue;
        if (hiter->second->get_t(1) < tmin_max[*layer].first) continue;
	if (hiter->second->get_t(1) - hiter->second->get_t(0) > m_DeltaTMap[*layer]) continue;


        pair<double, double> etaphi[2];
        double phibin[2];
        double etabin[2];
        for (int i = 0; i < 2; i++)
        {
          etaphi[i] = PHG4Utils::get_etaphi(hiter->second->get_x(i), hiter->second->get_y(i), hiter->second->get_z(i));
          etabin[i] = geo->get_etabin(etaphi[i].first);
          phibin[i] = geo->get_phibin(etaphi[i].second);
        }
        // check bin range
        if (phibin[0] < 0 || phibin[0] >= nphibins || phibin[1] < 0 || phibin[1] >= nphibins)
        {
          continue;
        }
        if (etabin[0] < 0 || etabin[0] >= nzbins || etabin[1] < 0 || etabin[1] >= nzbins)
        {
          continue;
        }

        if (etabin[0] < 0)
        {
          if (Verbosity() > 0)
          {
            hiter->second->identify();
          }
          continue;
        }
        sum_energy_g4hit += hiter->second->get_edep();
        int intphibin = phibin[0];
        int intetabin = etabin[0];
        int intphibinout = phibin[1];
        int intetabinout = etabin[1];

        // Determine all fired cells

        double ax = (etaphi[0]).second;  // phi
        double ay = (etaphi[0]).first;   // eta
        double bx = (etaphi[1]).second;
        double by = (etaphi[1]).first;
        if (intphibin > intphibinout)
        {
          int tmp = intphibin;
          intphibin = intphibinout;
          intphibinout = tmp;
        }
        if (intetabin > intetabinout)
        {
          int tmp = intetabin;
          intetabin = intetabinout;
          intetabinout = tmp;
        }

        double trklen = sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
        // if entry and exit hit are the same (seems to happen rarely), trklen = 0
        // which leads to a 0/0 and an NaN in edep later on
        // this code does for particles in the same cell a trklen/trklen (vdedx[ii]/trklen)
        // so setting this to any non zero number will do just fine
        // I just pick -1 here to flag those strange hits in case I want t oanalyze them
        // later on
        if (trklen == 0)
        {
          trklen = -1.;
        }
        vector<int> vphi;
        vector<int> veta;
        vector<double> vdedx;

        if (intphibin == intphibinout && intetabin == intetabinout)  // single cell fired
        {
          if (Verbosity() > 0) cout << "SINGLE CELL FIRED: " << intphibin << " " << intetabin << endl;
          vphi.push_back(intphibin);
          veta.push_back(intetabin);
          vdedx.push_back(trklen);
        }
        else
        {
          double phistep_half = geo->get_phistep() / 2.;
          double etastep_half = geo->get_etastep() / 2.;
          for (int ibp = intphibin; ibp <= intphibinout; ibp++)
          {
            double cx = geo->get_phicenter(ibp) - phistep_half;
            double dx = geo->get_phicenter(ibp) + phistep_half;
            for (int ibz = intetabin; ibz <= intetabinout; ibz++)
            {
              double cy = geo->get_etacenter(ibz) - etastep_half;
              double dy = geo->get_etacenter(ibz) + etastep_half;

              //cout << "##### line: " << ax << " " << ay << " " << bx << " " << by << endl;
              //cout << "####### cell: " << cx << " " << cy << " " << dx << " " << dy << endl;
              pair<bool, double> intersect = PHG4Utils::line_and_rectangle_intersect(ax, ay, bx, by, cx, cy, dx, dy);
              if (intersect.first)
              {
                if (Verbosity() > 0) cout << "CELL FIRED: " << ibp << " " << ibz << " " << intersect.second << endl;
                vphi.push_back(ibp);
                veta.push_back(ibz);
                vdedx.push_back(intersect.second);
              }
            }
          }
        }
        if (Verbosity() > 0) cout << "NUMBER OF FIRED CELLS = " << vphi.size() << endl;

        double tmpsum = 0.;
        for (unsigned int ii = 0; ii < vphi.size(); ii++)
        {
          tmpsum += vdedx[ii];
          vdedx[ii] = vdedx[ii] / trklen;
          if (Verbosity() > 0) cout << "  CELL " << ii << "  dE/dX = " << vdedx[ii] << endl;
        }
        if (Verbosity() > 0) cout << "    TOTAL TRACK LENGTH = " << tmpsum << " " << trklen << endl;

        for (unsigned int i1 = 0; i1 < vphi.size(); i1++)  // loop over all fired cells
        {
          int iphibin = vphi[i1];
          int ietabin = veta[i1];

          // This is the key for cellptmap
          // It is constructed using the phi and z (or eta) bin index values
          // It will be unique for a given phi and z (or eta) bin combination
          unsigned long long tmp = iphibin;
          unsigned long long key = tmp << 32;
          key += ietabin;
          if (Verbosity() > 1)
          {
            cout << " iphibin " << iphibin << " ietabin " << ietabin << " key 0x" << hex << key << dec << endl;
          }
          PHG4Cell *cell = nullptr;
          it = cellptmap.find(key);
          if (it != cellptmap.end())
          {
            cell = it->second;
          }
          else
          {
            PHG4CellDefs::keytype cellkey = PHG4CellDefs::EtaPhiBinning::genkey(*layer, ietabin, iphibin);
            cell = new PHG4Cellv1(cellkey);
            cellptmap[key] = cell;
          }
          if (!isfinite(hiter->second->get_edep() * vdedx[i1]))
          {
            cout << "hit 0x" << hex << hiter->first << dec << " not finite, edep: "
                 << hiter->second->get_edep() << " weight " << vdedx[i1] << endl;
          }
          cell->add_edep(hiter->first, hiter->second->get_edep() * vdedx[i1]);  // add hit with edep to g4hit list
          cell->add_edep(hiter->second->get_edep() * vdedx[i1]);                // add edep to cell
          if (hiter->second->has_property(PHG4Hit::prop_light_yield))
          {
            cell->add_light_yield(hiter->second->get_light_yield() * vdedx[i1]);
          }
          cell->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep() * vdedx[i1]);
          // just a sanity check - we don't want to mess up by having Nan's or Infs in our energy deposition
          if (!isfinite(hiter->second->get_edep() * vdedx[i1]))
          {
            cout << PHWHERE << " invalid energy dep " << hiter->second->get_edep()
                 << " or path length: " << vdedx[i1] << endl;
          }
        }
        vphi.clear();
        veta.clear();

      }  // end loop over g4hits

      int numcells = 0;

      for (it = cellptmap.begin(); it != cellptmap.end(); ++it)
      {
        cells->AddCell(it->second);
        numcells++;
        if (Verbosity() > 1)
        {
          cout << "Adding cell in bin phi: " << PHG4CellDefs::EtaPhiBinning::get_phibin(it->second->get_cellid())
               << " phi: " << geo->get_phicenter(PHG4CellDefs::EtaPhiBinning::get_phibin(it->second->get_cellid())) * 180. / M_PI
               << ", eta bin: " << PHG4CellDefs::EtaPhiBinning::get_etabin(it->second->get_cellid())
               << ", eta: " << geo->get_etacenter(PHG4CellDefs::EtaPhiBinning::get_etabin(it->second->get_cellid()))
               << ", energy dep: " << it->second->get_edep()
               << endl;
        }
      }

      if (Verbosity() > 0)
      {
        cout << Name() << ": found " << numcells << " eta/phi cells with energy deposition" << endl;
      }
    }

    else  // ------ size binning ---------------------------------------------------------------
    {
      sizeiter = cell_size.find(*layer);
      if (sizeiter == cell_size.end())
      {
        cout << "logical screwup!!! no sizes for layer " << *layer << endl;
        exit(1);
      }
      double zstepsize = (sizeiter->second).second;
      double phistepsize = phistep[*layer];

      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
      {
        sum_energy_before_cuts += hiter->second->get_edep();
        // checking ADC timing integration window cut
        if (hiter->second->get_t(0) > tmin_max[*layer].second) continue;
        if (hiter->second->get_t(1) < tmin_max[*layer].first) continue;
	if (hiter->second->get_t(1) - hiter->second->get_t(0) > 100) continue;

        double xinout[2];
        double yinout[2];
        double px[2];
        double py[2];
        double phi[2];
        double z[2];
        double phibin[2];
        double zbin[2];
        if (Verbosity() > 0) cout << "--------- new hit in layer # " << *layer << endl;

        for (int i = 0; i < 2; i++)
        {
          xinout[i] = hiter->second->get_x(i);
          yinout[i] = hiter->second->get_y(i);
          px[i] = hiter->second->get_px(i);
          py[i] = hiter->second->get_py(i);
          phi[i] = atan2(hiter->second->get_y(i), hiter->second->get_x(i));
          z[i] = hiter->second->get_z(i);
          phibin[i] = geo->get_phibin(phi[i]);
          zbin[i] = geo->get_zbin(hiter->second->get_z(i));

          if (Verbosity() > 0) cout << " " << i << "  phibin: " << phibin[i] << ", phi: " << phi[i] << ", stepsize: " << phistepsize << endl;
          if (Verbosity() > 0) cout << " " << i << "  zbin: " << zbin[i] << ", z = " << hiter->second->get_z(i) << ", stepsize: " << zstepsize << " offset: " << zmin_max[*layer].first << endl;
        }
        // check bin range
        if (phibin[0] < 0 || phibin[0] >= nphibins || phibin[1] < 0 || phibin[1] >= nphibins)
        {
          continue;
        }
        if (zbin[0] < 0 || zbin[0] >= nzbins || zbin[1] < 0 || zbin[1] >= nzbins)
        {
          continue;
        }

        if (zbin[0] < 0)
        {
          hiter->second->identify();
          continue;
        }
        sum_energy_g4hit += hiter->second->get_edep();

        int intphibin = phibin[0];
        int intzbin = zbin[0];
        int intphibinout = phibin[1];
        int intzbinout = zbin[1];

        if (Verbosity() > 0)
        {
          cout << "    phi bin range: " << intphibin << " to " << intphibinout << " phi: " << phi[0] << " to " << phi[1] << endl;
          cout << "    Z bin range: " << intzbin << " to " << intzbinout << " Z: " << z[0] << " to " << z[1] << endl;
          cout << "    phi difference: " << (phi[1] - phi[0]) * 1000. << " milliradians." << endl;
          cout << "    phi difference: " << 2.5 * (phi[1] - phi[0]) * 10000. << " microns." << endl;
          cout << "    path length = " << sqrt((xinout[1] - xinout[0]) * (xinout[1] - xinout[0]) + (yinout[1] - yinout[0]) * (yinout[1] - yinout[0])) << endl;
          cout << "       px = " << px[0] << " " << px[1] << endl;
          cout << "       py = " << py[0] << " " << py[1] << endl;
          cout << "       x = " << xinout[0] << " " << xinout[1] << endl;
          cout << "       y = " << yinout[0] << " " << yinout[1] << endl;
        }

        // Determine all fired cells

        double ax = phi[0];
        double ay = z[0];
        double bx = phi[1];
        double by = z[1];
        if (intphibin > intphibinout)
        {
          int tmp = intphibin;
          intphibin = intphibinout;
          intphibinout = tmp;
        }
        if (intzbin > intzbinout)
        {
          int tmp = intzbin;
          intzbin = intzbinout;
          intzbinout = tmp;
        }

        double trklen = sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
        // if entry and exit hit are the same (seems to happen rarely), trklen = 0
        // which leads to a 0/0 and an NaN in edep later on
        // this code does for particles in the same cell a trklen/trklen (vdedx[ii]/trklen)
        // so setting this to any non zero number will do just fine
        // I just pick -1 here to flag those strange hits in case I want to analyze them
        // later on
        if (trklen == 0)
        {
          trklen = -1.;
        }
        vector<int> vphi;
        vector<int> vz;
        vector<double> vdedx;

        if (intphibin == intphibinout && intzbin == intzbinout)  // single cell fired
        {
          if (Verbosity() > 0)
          {
            cout << "SINGLE CELL FIRED: " << intphibin << " " << intzbin << endl;
          }
          vphi.push_back(intphibin);
          vz.push_back(intzbin);
          vdedx.push_back(trklen);
        }
        else
        {
          double phistep_half = geo->get_phistep() / 2.;
          double zstep_half = geo->get_zstep() / 2.;
          for (int ibp = intphibin; ibp <= intphibinout; ibp++)
          {
            double cx = geo->get_phicenter(ibp) - phistep_half;
            double dx = geo->get_phicenter(ibp) + phistep_half;
            for (int ibz = intzbin; ibz <= intzbinout; ibz++)
            {
              double cy = geo->get_zcenter(ibz) - zstep_half;
              double dy = geo->get_zcenter(ibz) + zstep_half;
              //cout << "##### line: " << ax << " " << ay << " " << bx << " " << by << endl;
              //cout << "####### cell: " << cx << " " << cy << " " << dx << " " << dy << endl;
              pair<bool, double> intersect = PHG4Utils::line_and_rectangle_intersect(ax, ay, bx, by, cx, cy, dx, dy);
              if (intersect.first)
              {
                if (Verbosity() > 0) cout << "CELL FIRED: " << ibp << " " << ibz << " " << intersect.second << endl;
                vphi.push_back(ibp);
                vz.push_back(ibz);
                vdedx.push_back(intersect.second);
              }
            }
          }
        }
        if (Verbosity() > 0) cout << "NUMBER OF FIRED CELLS = " << vz.size() << endl;

        double tmpsum = 0.;
        for (unsigned int ii = 0; ii < vz.size(); ii++)
        {
          tmpsum += vdedx[ii];
          vdedx[ii] = vdedx[ii] / trklen;
          if (Verbosity() > 0) cout << "  CELL " << ii << "  dE/dX = " << vdedx[ii] << endl;
        }
        if (Verbosity() > 0) cout << "    TOTAL TRACK LENGTH = " << tmpsum << " " << trklen << endl;

        for (unsigned int i1 = 0; i1 < vphi.size(); i1++)  // loop over all fired cells
        {
          int iphibin = vphi[i1];
          int izbin = vz[i1];

          unsigned long long tmp = iphibin;
          unsigned long long key = tmp << 32;
          key += izbin;
          if (Verbosity() > 1)
          {
            cout << " iphibin " << iphibin << " izbin " << izbin << " key 0x" << hex << key << dec << endl;
          }
          // check to see if there is already an entry for this cell
          PHG4Cell *cell = nullptr;
          it = cellptmap.find(key);

          if (it != cellptmap.end())
          {
            cell = it->second;
            if (Verbosity() > 1)
            {
              cout << "  add energy to existing cell for key = " << cellptmap.find(key)->first << endl;
            }

            if (Verbosity() > 1 && hiter->second->has_property(PHG4Hit::prop_light_yield) && std::isnan(hiter->second->get_light_yield() * vdedx[i1]))
            {
              cout << "    NAN lighy yield with vdedx[i1] = " << vdedx[i1]
                   << " and hiter->second->get_light_yield() = " << hiter->second->get_light_yield() << endl;
            }
          }
          else
          {
            if (Verbosity() > 1)
            {
              cout << "    did not find a previous entry for key = 0x" << hex << key << dec << " create a new one" << endl;
            }
            PHG4CellDefs::keytype cellkey = PHG4CellDefs::SizeBinning::genkey(*layer, izbin, iphibin);
            cell = new PHG4Cellv1(cellkey);
            cellptmap[key] = cell;
          }
          if (!isfinite(hiter->second->get_edep() * vdedx[i1]))
          {
            cout << "hit 0x" << hex << hiter->first << dec << " not finite, edep: "
                 << hiter->second->get_edep() << " weight " << vdedx[i1] << endl;
          }
          cell->add_edep(hiter->first, hiter->second->get_edep() * vdedx[i1]);
          cell->add_edep(hiter->second->get_edep() * vdedx[i1]);  // add edep to cell
          if (hiter->second->has_property(PHG4Hit::prop_light_yield))
          {
            cell->add_light_yield(hiter->second->get_light_yield() * vdedx[i1]);
            if (Verbosity() > 1 && !std::isfinite(hiter->second->get_light_yield() * vdedx[i1]))
            {
              cout << "    NAN lighy yield with vdedx[i1] = " << vdedx[i1]
                   << " and hiter->second->get_light_yield() = " << hiter->second->get_light_yield() << endl;
            }
          }
          cell->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep() * vdedx[i1]);
        }
        vphi.clear();
        vz.clear();

      }  // end loop over hits

      int numcells = 0;

      for (it = cellptmap.begin(); it != cellptmap.end(); ++it)
      {
        cells->AddCell(it->second);
        numcells++;
        if (Verbosity() > 1)
        {
          cout << "Adding cell for key " << hex << it->first << dec << " in bin phi: " << PHG4CellDefs::SizeBinning::get_phibin(it->second->get_cellid())
               << " phi: " << geo->get_phicenter(PHG4CellDefs::SizeBinning::get_phibin(it->second->get_cellid())) * 180. / M_PI
               << ", z bin: " << PHG4CellDefs::SizeBinning::get_zbin(it->second->get_cellid())
               << ", z: " << geo->get_zcenter(PHG4CellDefs::SizeBinning::get_zbin(it->second->get_cellid()))
               << ", energy dep: " << it->second->get_edep()
               << endl;
        }
      }

      if (Verbosity() > 0)
      {
        cout << "found " << numcells << " z/phi cells with energy deposition" << endl;
      }
    }

    //==========================================================
    // now reset the cell map before moving on to the next layer
    if (Verbosity() > 1)
    {
      cout << "cellptmap for layer " << *layer << " has final length " << cellptmap.size();
    }
    while (cellptmap.begin() != cellptmap.end())
    {
      // Assumes that memmory is freed by the cylinder cell container when it is destroyed
      cellptmap.erase(cellptmap.begin());
    }
    if (Verbosity() > 1)
    {
      cout << " reset it to " << cellptmap.size() << endl;
    }
  }
  if (chkenergyconservation)
  {
    CheckEnergy(topNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4CylinderCellReco::cellsize(const int detid, const double sr, const double sz)
{
  if (binning.find(detid) != binning.end())
  {
    cout << "size for layer " << detid << " already set" << endl;
    return;
  }
  binning[detid] = PHG4CellDefs::sizebinning;
  set_double_param(detid, "size_long", sz);
  set_double_param(detid, "size_perp", sr);
}

void PHG4CylinderCellReco::etaphisize(const int detid, const double deltaeta, const double deltaphi)
{
  if (binning.find(detid) != binning.end())
  {
    cout << "size for layer " << detid << " already set" << endl;
    return;
  }
  binning[detid] = PHG4CellDefs::etaphibinning;
  set_double_param(detid, "size_long", deltaeta);
  set_double_param(detid, "size_perp", deltaphi);
  return;
}

void PHG4CylinderCellReco::set_size(const int i, const double sizeA, const double sizeB)
{
  cell_size[i] = std::make_pair(sizeA, sizeB);
  return;
}

void PHG4CylinderCellReco::set_timing_window(const int detid, const double tmin, const double tmax)
{
  set_double_param(detid, "tmin", tmin);
  set_double_param(detid, "tmax", tmax);
  return;
}

int PHG4CylinderCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  double sum_energy_cells = 0.;
  double sum_energy_stored_hits = 0.;
  double sum_energy_stored_showers = 0.;
  PHG4CellContainer::ConstRange cell_begin_end = cells->getCells();
  PHG4CellContainer::ConstIterator citer;
  for (citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
  {
    sum_energy_cells += citer->second->get_edep();
    PHG4Cell::EdepConstRange cellrange = citer->second->get_g4hits();
    for (PHG4Cell::EdepConstIterator iter = cellrange.first; iter != cellrange.second; ++iter)
    {
      sum_energy_stored_hits += iter->second;
    }
    PHG4Cell::ShowerEdepConstRange shwrrange = citer->second->get_g4showers();
    for (PHG4Cell::ShowerEdepConstIterator iter = shwrrange.first; iter != shwrrange.second; ++iter)
    {
      sum_energy_stored_showers += iter->second;
    }
  }
  // the fractional eloss for particles traversing eta bins leads to minute rounding errors
  if (sum_energy_stored_hits > 0)
  {
    if (fabs(sum_energy_cells - sum_energy_stored_hits) / sum_energy_cells > 1e-6)
    {
      cout << "energy mismatch between cell energy " << sum_energy_cells
           << " and stored hit energies " << sum_energy_stored_hits
           << endl;
    }
  }
  if (sum_energy_stored_showers > 0)
  {
    if (fabs(sum_energy_cells - sum_energy_stored_showers) / sum_energy_cells > 1e-6)
    {
      cout << "energy mismatch between cell energy " << sum_energy_cells
           << " and stored shower energies " << sum_energy_stored_showers
           << endl;
    }
  }
  if (fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit > 1e-6)
  {
    cout << "energy mismatch between cells: " << sum_energy_cells
         << " and hits: " << sum_energy_g4hit
         << " diff sum(cells) - sum(hits): " << sum_energy_cells - sum_energy_g4hit
         << " cut val " << fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit
         << endl;
    cout << Name() << ":total energy for this event: " << sum_energy_g4hit << " GeV" << endl;
    cout << Name() << ": sum cell energy: " << sum_energy_cells << " GeV" << endl;
    cout << Name() << ": sum shower energy: " << sum_energy_stored_showers << " GeV" << endl;
    cout << Name() << ": sum stored hit energy: " << sum_energy_stored_hits << " GeV" << endl;
    cout << Name() << ": hit energy before cuts: " << sum_energy_before_cuts << " GeV" << endl;
    return -1;
  }
  else
  {
    if (Verbosity() > 0)
    {
      cout << Name() << ":total energy for this event: " << sum_energy_g4hit << " GeV" << endl;
      cout << Name() << ": sum cell energy: " << sum_energy_cells << " GeV" << endl;
      cout << Name() << ": sum shower energy: " << sum_energy_stored_showers << " GeV" << endl;
      cout << Name() << ": sum stored hit energy: " << sum_energy_stored_hits << " GeV" << endl;
      cout << Name() << ": hit energy before cuts: " << sum_energy_before_cuts << " GeV" << endl;
    }
  }
  return 0;
}

void PHG4CylinderCellReco::Detector(const std::string &d)
{
  detector = d;
  // only set the outdetector nodename if it wasn't set already
  // in case the order in the macro is outdetector(); detector();
  if (outdetector.size() == 0)
  {
    OutputDetector(d);
  }
  return;
}

void PHG4CylinderCellReco::SetDefaultParameters()
{
  set_default_double_param("size_long", NAN);
  set_default_double_param("size_perp", NAN);
  set_default_double_param("tmax", 60.0);
  set_default_double_param("tmin", -20.0);  // collision has a timing spread around the triggered event. Accepting negative time too.
  set_default_double_param("delta_t", 100.);
  return;
}
