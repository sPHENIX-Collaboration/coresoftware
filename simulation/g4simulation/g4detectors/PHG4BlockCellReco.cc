#include "PHG4BlockCellReco.h"

#include "PHG4BlockCellGeom.h"
#include "PHG4BlockCellGeomContainer.h"
#include "PHG4BlockGeom.h"
#include "PHG4BlockGeomContainer.h"
#include "PHG4Cell.h"  // for PHG4Cell, PHG...
#include "PHG4CellContainer.h"
#include "PHG4Cellv1.h"
#include "PHG4CellDefs.h"

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
#include <iostream>
#include <sstream>
#include <vector>  // for vector

namespace
{
  std::vector<PHG4Cell *> cellptarray;
}

PHG4BlockCellReco::PHG4BlockCellReco(const std::string &name)
  : SubsysReco(name)
  , PHParameterContainerInterface(name)
  , sum_energy_g4hit(0)
  , chkenergyconservation(0)
{
  SetDefaultParameters();
}

int PHG4BlockCellReco::ResetEvent(PHCompositeNode * /*topNode*/)
{
  sum_energy_g4hit = 0.;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4BlockCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << Name() << " DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode =
      dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode",
                                                        detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(detector);
    dstNode->addNode(DetNode);
  }

  PHCompositeNode *runNode;
  runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << Name() << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  PHNodeIterator runiter(runNode);
  PHCompositeNode *RunDetNode =
      dynamic_cast<PHCompositeNode *>(runiter.findFirst("PHCompositeNode",
                                                        detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(detector);
    runNode->addNode(RunDetNode);
  }

  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    std::cout << Name() << " Could not locate g4 hit node " << hitnodename << std::endl;
    exit(1);
  }

  cellnodename = "G4CELL_" + detector;
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
  {
    cells = new PHG4CellContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells, cellnodename, "PHObject");
    DetNode->addNode(newNode);
  }

  geonodename = "BLOCKGEOM_" + detector;
  PHG4BlockGeomContainer *geo = findNode::getClass<PHG4BlockGeomContainer>(topNode, geonodename);
  if (!geo)
  {
    std::cout << Name() << " Could not locate geometry node " << geonodename << std::endl;
    exit(1);
  }

  if (Verbosity() > 0)
  {
    geo->identify();
  }

  seggeonodename = "BLOCKCELLGEOM_" + detector;
  PHG4BlockCellGeomContainer *seggeo = findNode::getClass<PHG4BlockCellGeomContainer>(topNode, seggeonodename);
  if (!seggeo)
  {
    seggeo = new PHG4BlockCellGeomContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename, "PHObject");
    runNode->addNode(newNode);
  }
  GetParamsContainerModify()->set_name(detector);

  UpdateParametersWithMacro();

  std::map<int, PHG4BlockGeom *>::const_iterator miter;
  std::pair<std::map<int, PHG4BlockGeom *>::const_iterator, std::map<int, PHG4BlockGeom *>::const_iterator> begin_end = geo->get_begin_end();
  std::map<int, std::pair<double, double> >::iterator sizeiter;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
  {
    PHG4BlockGeom *layergeom = miter->second;
    int layer = layergeom->get_layer();
    if (!ExistDetid(layer))
    {
      std::cout << Name() << ": No parameters for detid/layer " << layer
           << ", hits from this detid/layer will not be accumulated into cells" << std::endl;
      continue;
    }
    implemented_detid.insert(layer);
    double radius = sqrt(pow(layergeom->get_center_x(), 2) + pow(layergeom->get_center_y(), 2));
    double width = layergeom->get_size_x();
    double length_in_z = layergeom->get_size_z();
    double zmin = layergeom->get_center_z() - length_in_z / 2.;
    double zmax = zmin + length_in_z;
    set_size(layer, get_double_param(layer, "deltaeta"), get_double_param(layer, "deltax"), PHG4CellDefs::etaphibinning);
    tmin_max.insert(std::make_pair(layer, std::make_pair(get_double_param(layer, "tmin"), get_double_param(layer, "tmax"))));
    sizeiter = cell_size.find(layer);
    if (sizeiter == cell_size.end())
    {
      std::cout << Name() << "no cell sizes for layer " << layer << std::endl;
      exit(1);
    }

    // create geo object and fill with variables common to all binning methods
    PHG4BlockCellGeom *layerseggeo = new PHG4BlockCellGeom();
    layerseggeo->set_layer(layergeom->get_layer());
    // layerseggeo->set_radius(layergeom->get_radius());
    // layerseggeo->set_thickness(layergeom->get_thickness());

    if (binning[layer] == PHG4CellDefs::etaphibinning)
    {
      // calculate eta at radius+ thickness (outer radius)
      // length via eta coverage is calculated using the outer radius
      double etamin = PHG4Utils::get_eta(radius + 0.5 * layergeom->get_size_y(), zmin);
      double etamax = PHG4Utils::get_eta(radius + 0.5 * layergeom->get_size_y(), zmax);
      zmin_max[layer] = std::make_pair(etamin, etamax);
      double etastepsize = (sizeiter->second).first;
      double d_etabins;
      double fract = modf((etamax - etamin) / etastepsize, &d_etabins);
      if (fract != 0)
      {
        d_etabins++;
      }
      etastepsize = (etamax - etamin) / d_etabins;
      (sizeiter->second).first = etastepsize;
      int etabins = d_etabins;
      double etahi = etamin + etastepsize;
      for (int i = 0; i < etabins; i++)
      {
        if (etahi > (etamax + 1.e-6))  // etahi is a tiny bit larger due to numerical uncertainties
        {
          std::cout << "etahi: " << etahi << ", etamax: " << etamax << std::endl;
        }
        etahi += etastepsize;
      }

      double xmin = -layergeom->get_width() / 2.;
      // double xmax = -xmin;
      double xstepsize = (sizeiter->second).second;
      double d_xbins;
      fract = modf(width / xstepsize, &d_xbins);
      if (fract != 0)
      {
        d_xbins++;
      }

      xstepsize = width / d_xbins;
      (sizeiter->second).second = xstepsize;
      int xbins = d_xbins;
      // double xlow = xmin;
      // double xhi = xlow + xstepsize;

      // for (int i = 0; i < xbins; i++)
      // {
      //   if (xhi > xmax)
      //   {
      //     std::cout << "xhi: " << xhi << ", xmax: " << xmax << std::endl;
      //   }
      //   xlow = xhi;
      //   xhi +=  xstepsize;
      // }

      std::pair<int, int> x_z_bin = std::make_pair(xbins, etabins);
      n_x_z_bins[layer] = x_z_bin;
      layerseggeo->set_binning(PHG4CellDefs::etaphibinning);
      layerseggeo->set_etabins(etabins);
      layerseggeo->set_etamin(etamin);
      layerseggeo->set_etastep(etastepsize);
      layerseggeo->set_xmin(xmin);
      layerseggeo->set_xbins(xbins);
      layerseggeo->set_xstep(xstepsize);
      xstep[layer] = xstepsize;
      etastep[layer] = etastepsize;
    }

    // add geo object filled by different binning methods
    seggeo->AddLayerCellGeom(layerseggeo);
    if (Verbosity() > 1)
    {
      layerseggeo->identify();
    }
  }
  std::string nodename = "G4CELLPARAM_" + GetParamsContainer()->Name();
  SaveToNodeTree(RunDetNode, nodename);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4BlockCellReco::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << hitnodename << std::endl;
    exit(1);
  }
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
  {
    std::cout << "could not locate cell node " << cellnodename << std::endl;
    exit(1);
  }

  PHG4BlockCellGeomContainer *seggeo = findNode::getClass<PHG4BlockCellGeomContainer>(topNode, seggeonodename);
  if (!seggeo)
  {
    std::cout << "could not locate geo node " << seggeonodename << std::endl;
    exit(1);
  }

  PHG4HitContainer::LayerIter layer;
  std::pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  //   std::cout << "number of layers: " << g4hit->num_layers() << std::endl;
  //   std::cout << "number of hits: " << g4hit->size() << std::endl;
  //   for (layer = layer_begin_end.first; layer != layer_begin_end.second; layer++)
  //     {
  //       std::cout << "layer number: " << *layer << std::endl;
  //     }

  for (layer = layer_begin_end.first; layer != layer_begin_end.second; layer++)
  {
    // only handle layers/detector ids which have parameters set
    if (!implemented_detid.contains(*layer))
    {
      continue;
    }
    PHG4HitContainer::ConstIterator hiter;
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);
    PHG4BlockCellGeom *geo = seggeo->GetLayerCellGeom(*layer);
    int nxbins = n_x_z_bins[*layer].first;
    int nzbins = n_x_z_bins[*layer].second;
    unsigned int nbins = nxbins * nzbins;

    if (cellptarray.size() < nbins)
    {
      cellptarray.resize(nbins, nullptr);
    }

    // ------- eta/x binning ------------------------------------------------------------------------
    if (binning[*layer] == PHG4CellDefs::etaphibinning)
    {
      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
      {
        // checking ADC timing integration window cut
        if (hiter->second->get_t(0) > tmin_max[*layer].second)
        {
          continue;
        }
        if (hiter->second->get_t(1) < tmin_max[*layer].first)
        {
          continue;
        }

        std::pair<double, double> etax[2];
        double xbin[2];
        double etabin[2];
        for (int i = 0; i < 2; i++)
        {
          etax[i] = PHG4Utils::get_etaphi(hiter->second->get_x(i), hiter->second->get_y(i), hiter->second->get_z(i));
          etabin[i] = geo->get_etabin(etax[i].first);
          xbin[i] = geo->get_xbin(etax[i].second);
        }

        // check bin range
        if (xbin[0] < 0 || xbin[0] >= nxbins || xbin[1] < 0 || xbin[1] >= nxbins)
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
        int intxbin = xbin[0];
        int intetabin = etabin[0];
        int intxbinout = xbin[1];
        int intetabinout = etabin[1];

        // Determine all fired cells

        double ax = (etax[0]).second;  // x
        double ay = (etax[0]).first;   // eta
        double bx = (etax[1]).second;
        double by = (etax[1]).first;
        if (intxbin > intxbinout)
        {
          int tmp = intxbin;
          intxbin = intxbinout;
          intxbinout = tmp;
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
        std::vector<int> vx;
        std::vector<int> veta;
        std::vector<double> vdedx;

        if (intxbin == intxbinout && intetabin == intetabinout)  // single cell fired
        {
          if (Verbosity() > 0)
          {
            std::cout << "SINGLE CELL FIRED: " << intxbin << " " << intetabin << std::endl;
          }
          vx.push_back(intxbin);
          veta.push_back(intetabin);
          vdedx.push_back(trklen);
        }
        else
        {
          for (int ibp = intxbin; ibp <= intxbinout; ibp++)
          {
            for (int ibz = intetabin; ibz <= intetabinout; ibz++)
            {
              double cx = geo->get_xcenter(ibp) - geo->get_xstep() / 2.;
              double dx = geo->get_xcenter(ibp) + geo->get_xstep() / 2.;
              double cy = geo->get_etacenter(ibz) - geo->get_etastep() / 2.;
              double dy = geo->get_etacenter(ibz) + geo->get_etastep() / 2.;
              // std::cout << "##### line: " << ax << " " << ay << " " << bx << " " << by << std::endl;
              // std::cout << "####### cell: " << cx << " " << cy << " " << dx << " " << dy << std::endl;
              std::pair<bool, double> intersect = PHG4Utils::line_and_rectangle_intersect(ax, ay, bx, by, cx, cy, dx, dy);
              if (intersect.first)
              {
                if (Verbosity() > 0)
                {
                  std::cout << "CELL FIRED: " << ibp << " " << ibz << " " << intersect.second << std::endl;
                }
                vx.push_back(ibp);
                veta.push_back(ibz);
                vdedx.push_back(intersect.second);
              }
            }
          }
        }
        if (Verbosity() > 0)
        {
          std::cout << "NUMBER OF FIRED CELLS = " << vx.size() << std::endl;
        }

        double tmpsum = 0.;
        for (unsigned int ii = 0; ii < vx.size(); ii++)
        {
          tmpsum += vdedx[ii];
          vdedx[ii] = vdedx[ii] / trklen;
          if (Verbosity() > 0)
          {
            std::cout << "  CELL " << ii << "  dE/dX = " << vdedx[ii] << std::endl;
          }
        }
        if (Verbosity() > 0)
        {
          std::cout << "    TOTAL TRACK LENGTH = " << tmpsum << " " << trklen << std::endl;
        }

        for (unsigned int i1 = 0; i1 < vx.size(); i1++)  // loop over all fired cells
        {
          int ixbin = vx[i1];
          int ietabin = veta[i1];
          int ibin = ixbin * nzbins + ietabin;
          if (!cellptarray[ibin])
          {
	    // this is just messed up, the args are swapped but consistently - don't try to fix it
	    //NOLINTNEXTLINE(readability-suspicious-call-argument)
            PHG4CellDefs::keytype key = PHG4CellDefs::EtaXsizeBinning::genkey(*layer, ixbin, ietabin);
            cellptarray[ibin] = new PHG4Cellv1(key);
          }
          cellptarray[ibin]->add_edep(hiter->first, hiter->second->get_edep() * vdedx[i1]);
          cellptarray[ibin]->add_edep(hiter->second->get_edep() * vdedx[i1]);
          if (hiter->second->has_property(PHG4Hit::prop_light_yield))
          {
            cellptarray[ibin]->add_light_yield(hiter->second->get_light_yield() * vdedx[i1]);
          }
          cellptarray[ibin]->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep() * vdedx[i1]);

          // just a sanity check - we don't want to mess up by having Nan's or Infs in our energy deposition
          if (!std::isfinite(hiter->second->get_edep() * vdedx[i1]))
          {
            std::cout << PHWHERE << " invalid energy dep " << hiter->second->get_edep()
                 << " or path length: " << vdedx[i1] << std::endl;
          }
        }

        vx.clear();
        veta.clear();
      }  // end loop over g4hits

      int numcells = 0;
      for (int ix = 0; ix < nxbins; ix++)
      {
        for (int iz = 0; iz < nzbins; iz++)
        {
          int ibin = ix * nzbins + iz;

          if (cellptarray[ibin])
          {
            cells->AddCell(cellptarray[ibin]);
            numcells++;
            if (Verbosity() > 1)
            {
              std::cout << "Adding cell in bin x: " << ix
                   << " x: " << geo->get_xcenter(ix) * 180. / M_PI
                   << ", eta bin: " << iz
                   << ", eta: " << geo->get_etacenter(iz)
                   << ", energy dep: " << cellptarray[ibin]->get_edep()
                   << std::endl;
            }

            cellptarray[ibin] = nullptr;
          }
        }
      }

      if (Verbosity() > 0)
      {
        std::cout << Name() << ": found " << numcells << " eta/x cells with energy deposition" << std::endl;
      }
    }
  }

  if (chkenergyconservation)
  {
    CheckEnergy(topNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4BlockCellReco::etaxsize(const int detid, const double deltaeta, const double deltax)
{
  set_double_param(detid, "deltaeta", deltaeta);
  set_double_param(detid, "deltax", deltax);
  return;
}

void PHG4BlockCellReco::set_timing_window(const int detid, const double tmin, const double tmax)
{
  set_double_param(detid, "tmin", tmin);
  set_double_param(detid, "tmax", tmax);
  return;
}

void PHG4BlockCellReco::set_size(const int i, const double sizeA, const double sizeB, const int what)
{
  if (binning.contains(i))
  {
    std::cout << "size for layer " << i << " already set" << std::endl;
    return;
  }

  binning[i] = what;
  cell_size[i] = std::make_pair(sizeA, sizeB);

  return;
}

//---------------------------------------------------------------

int PHG4BlockCellReco::CheckEnergy(PHCompositeNode *topNode)
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
      std::cout << "energy mismatch between cell energy " << sum_energy_cells
           << " and stored hit energies " << sum_energy_stored_hits
           << std::endl;
    }
  }
  if (sum_energy_stored_showers > 0)
  {
    if (fabs(sum_energy_cells - sum_energy_stored_showers) / sum_energy_cells > 1e-6)
    {
      std::cout << "energy mismatch between cell energy " << sum_energy_cells
           << " and stored shower energies " << sum_energy_stored_showers
           << std::endl;
    }
  }

  if (fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit > 1e-6)
  {
    std::cout << "energy mismatch between cells: " << sum_energy_cells
         << " and hits: " << sum_energy_g4hit
         << " diff sum(cells) - sum(hits): " << sum_energy_cells - sum_energy_g4hit
         << std::endl;
    return -1;
  }

  
      if (Verbosity() > 0)
    {
      std::cout << Name() << ": sum hit energy: " << sum_energy_g4hit << " GeV" << std::endl;
      std::cout << Name() << ": sum cell energy: " << sum_energy_cells << " GeV" << std::endl;
      std::cout << Name() << ": sum shower energy: " << sum_energy_stored_showers << " GeV" << std::endl;
      std::cout << Name() << ": sum stored hit energy: " << sum_energy_stored_hits << " GeV" << std::endl;
    }
 

  return 0;
}

void PHG4BlockCellReco::SetDefaultParameters()
{
  set_default_double_param("deltaeta", std::numeric_limits<double>::quiet_NaN());
  set_default_double_param("deltax", std::numeric_limits<double>::quiet_NaN());
  set_default_double_param("tmax", 60.0);
  set_default_double_param("tmin", 0.0);
  return;
}
