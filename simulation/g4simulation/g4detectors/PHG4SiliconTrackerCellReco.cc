#include "PHG4SiliconTrackerCellReco.h"
#include "PHG4CellContainer.h"
#include "PHG4Cellv1.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomSiLadders.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <boost/format.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

PHG4SiliconTrackerCellReco::PHG4SiliconTrackerCellReco(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_ChkEnergyConservationFlag(0)
  , m_Tmin(NAN)
  , m_Tmax(NAN)
{
  InitializeParameters();
  Detector(name);

  m_HitNodeName = "G4HIT_" + m_Detector;
  m_CellNodeName = "G4CELL_" + m_Detector;
  m_GeoNodeName = "CYLINDERGEOM_" + m_Detector;
  m_LocalOutVec = gsl_vector_alloc(3);
  m_PathVec = gsl_vector_alloc(3);
  m_SegmentVec = gsl_vector_alloc(3);
}

PHG4SiliconTrackerCellReco::~PHG4SiliconTrackerCellReco()
{
  gsl_vector_free(m_LocalOutVec);
  gsl_vector_free(m_PathVec);
  gsl_vector_free(m_SegmentVec);
}

int PHG4SiliconTrackerCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, exiting." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  PHCompositeNode *runNode;
  runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    cout << Name() << "RUN Node missing, exiting." << endl;
    gSystem->Exit(1);
    exit(1);
  }
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode)
  {
    cout << Name() << "PAR Node missing, exiting." << endl;
    gSystem->Exit(1);
    exit(1);
  }
  string paramnodename = "G4CELLPARAM_" + m_Detector;

  PHNodeIterator runiter(runNode);
  PHCompositeNode *RunDetNode =
      dynamic_cast<PHCompositeNode *>(runiter.findFirst("PHCompositeNode",
                                                        m_Detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_Detector);
    runNode->addNode(RunDetNode);
  }

  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << m_HitNodeName << std::endl;
    exit(1);
  }

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, m_CellNodeName);
  if (!cells)
  {
    PHNodeIterator dstiter(dstNode);

    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));

    if (!DetNode)
    {
      DetNode = new PHCompositeNode(m_Detector);
      dstNode->addNode(DetNode);
    }

    cells = new PHG4CellContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells, m_CellNodeName, "PHObject");
    DetNode->addNode(newNode);
  }

  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, m_GeoNodeName);
  if (!geo)
  {
    std::cout << "Could not locate geometry node " << m_GeoNodeName << std::endl;
    exit(1);
  }

  if (Verbosity() > 0)
  {
    geo->identify();
  }

  UpdateParametersWithMacro();
  SaveToNodeTree(RunDetNode, paramnodename);
  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", m_Detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(m_Detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, m_GeoNodeName);
  m_Tmin = get_double_param("tmin");
  m_Tmax = get_double_param("tmax");

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SiliconTrackerCellReco::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << m_HitNodeName << std::endl;
    exit(1);
  }

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, m_CellNodeName);
  if (!cells)
  {
    std::cout << "could not locate cell node " << m_CellNodeName << std::endl;
    exit(1);
  }
  cells->Reset();

  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, m_GeoNodeName);
  if (!geo)
  {
    std::cout << "Could not locate geometry node " << m_GeoNodeName << std::endl;
    exit(1);
  }

  // loop over all of the layers in the hit container
  // we need the geometry object for this layer
  if (Verbosity() > 2) cout << " PHG4SiliconTrackerCellReco: Loop over hits" << endl;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  for (PHG4HitContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    const int sphxlayer = hiter->second->get_detid();
    PHG4CylinderGeomSiLadders *layergeom = dynamic_cast<PHG4CylinderGeomSiLadders *> (geo->GetLayerGeom(sphxlayer));

    // checking ADC timing integration window cut
    // uses default values for now
    // these should depend on layer radius
    if (hiter->second->get_t(0) > m_Tmax)
      continue;
    if (hiter->second->get_t(1) < m_Tmin)
      continue;

    // I made this (small) diffusion up for now, we will get actual values for the INTT later
    double diffusion_width = 5.0e-04;  // diffusion radius 5 microns, in cm

    const int ladder_z_index = hiter->second->get_ladder_z_index();
    const int ladder_phi_index = hiter->second->get_ladder_phi_index();

    // What we have is a hit in the sensor. We have not yet assigned the strip(s) that were hit, we do that here
    //========================================================================

    int strip_y_index_in, strip_z_index_in, strip_y_index_out, strip_z_index_out;
    layergeom->find_strip_index_values(ladder_z_index, hiter->second->get_local_y(0), hiter->second->get_local_z(0), strip_y_index_in, strip_z_index_in);
    layergeom->find_strip_index_values(ladder_z_index, hiter->second->get_local_y(1), hiter->second->get_local_z(1), strip_y_index_out, strip_z_index_out);

    if (Verbosity() > 5)
    {
      // check to see if we get back the positions from these strip index values
      double check_location[3] = {-1, -1, -1};
      layergeom->find_strip_center_localcoords(ladder_z_index, strip_y_index_in, strip_z_index_in, check_location);
      cout << " G4 entry location = " << hiter->second->get_local_x(0) << "  " << hiter->second->get_local_y(0) << "  " << hiter->second->get_local_z(0) << endl;
      cout << " Check entry location = " << check_location[0] << "  " << check_location[1] << "  " << check_location[2] << endl;
      layergeom->find_strip_center_localcoords(ladder_z_index, strip_y_index_out, strip_z_index_out, check_location);
      cout << " G4 exit location = " << hiter->second->get_local_x(1) << " " << hiter->second->get_local_y(1) << "  " << hiter->second->get_local_z(1) << endl;
      cout << " Check exit location = " << check_location[0] << "  " << check_location[1] << "  " << check_location[2] << endl;
    }

    // Now we find how many strips were crossed by this track, and divide the energy between them
    int minstrip_z = strip_z_index_in;
    int maxstrip_z = strip_z_index_out;
    if (minstrip_z > maxstrip_z) swap(minstrip_z, maxstrip_z);

    int minstrip_y = strip_y_index_in;
    int maxstrip_y = strip_y_index_out;
    if (minstrip_y > maxstrip_y) swap(minstrip_y, maxstrip_y);

    // Use an algorithm similar to the one for the MVTX pixels, since it facilitates adding charge diffusion
    // for now we assume small charge diffusion

    vector<int> vybin;
    vector<int> vzbin;
    //vector<double> vlen;
    vector<pair<double, double> > venergy;

    //====================================================
    // Beginning of charge sharing implementation
    //    Find tracklet line inside sensor
    //    Divide tracklet line into n segments (vary n until answer stabilizes)
    //    Find centroid of each segment
    //    Diffuse charge at each centroid
    //    Apportion charge between neighboring pixels
    //    Add the pixel energy contributions from different track segments together
    //====================================================

    // skip this hit if it involves an unreasonable  number of pixels
    // this skips it if either the xbin or ybin range traversed is greater than 8 (for 8 adding two pixels at each end makes the range 12)
    if (maxstrip_y - minstrip_y > 12 || maxstrip_z - minstrip_z > 12)
    {
      continue;
    }
    // this hit is skipped above if this dimensioning would be exceeded
    double stripenergy[12][12] = {};  // init to 0
    double stripeion[12][12] = {};    // init to 0

    int nsegments = 10;
    // Loop over track segments and diffuse charge at each segment location, collect energy in pixels
    // Get the entry point of the hit in sensor local coordinates
    gsl_vector_set(m_PathVec, 0, hiter->second->get_local_x(0));
    gsl_vector_set(m_PathVec, 1, hiter->second->get_local_y(0));
    gsl_vector_set(m_PathVec, 2, hiter->second->get_local_z(0));
    gsl_vector_set(m_LocalOutVec, 0, hiter->second->get_local_x(1));
    gsl_vector_set(m_LocalOutVec, 1, hiter->second->get_local_y(1));
    gsl_vector_set(m_LocalOutVec, 2, hiter->second->get_local_z(1));
    gsl_vector_sub(m_PathVec, m_LocalOutVec);
    for (int i = 0; i < nsegments; i++)
    {
      // Find the tracklet segment location
      // If there are n segments of equal length, we want 2*n intervals
      // The 1st segment is centered at interval 1, the 2nd at interval 3, the nth at interval 2n -1
      double interval = 2 * (double) i + 1;
      double frac = interval / (double) (2 * nsegments);
      gsl_vector_memcpy(m_SegmentVec, m_PathVec);
      gsl_vector_scale(m_SegmentVec, frac);
      gsl_vector_add(m_SegmentVec, m_LocalOutVec);
      // Caculate the charge diffusion over this drift distance
      // increases from diffusion width_min to diffusion_width_max
      double diffusion_radius = diffusion_width;

      if (Verbosity() > 5)
        cout << " segment " << i
             << " interval " << interval
             << " frac " << frac
             << " local_in.X " << hiter->second->get_local_x(0)
             << " local_in.Z " << hiter->second->get_local_z(0)
             << " local_in.Y " << hiter->second->get_local_y(0)
             << " pathvec.X " << gsl_vector_get(m_PathVec, 0)
             << " pathvec.Z " << gsl_vector_get(m_PathVec, 2)
             << " pathvec.Y " << gsl_vector_get(m_PathVec, 1)
             << " segvec.X " << gsl_vector_get(m_SegmentVec, 0)
             << " segvec.Z " << gsl_vector_get(m_SegmentVec, 2)
             << " segvec.Y " << gsl_vector_get(m_SegmentVec, 1) << endl
             << " diffusion_radius " << diffusion_radius
             << endl;

      // Now find the area of overlap of the diffusion circle with each pixel and apportion the energy
      for (int iz = minstrip_z; iz <= maxstrip_z; iz++)
      {
        for (int iy = minstrip_y; iy <= maxstrip_y; iy++)
        {
          // Find the pixel corners for this pixel number
          double location[3] = {-1, -1, -1};
          layergeom->find_strip_center_localcoords(ladder_z_index, iy, iz, location);
          // note that (y1,z1) is the top left corner, (y2,z2) is the bottom right corner of the pixel - circle_rectangle_intersection expects this ordering
          double y1 = location[1] - layergeom->get_strip_y_spacing() / 2.0;
          double y2 = location[1] + layergeom->get_strip_y_spacing() / 2.0;
          double z1 = location[2] + layergeom->get_strip_z_spacing() / 2.0;
          double z2 = location[2] - layergeom->get_strip_z_spacing() / 2.0;

          // here m_SegmentVec.1 (Y) and m_SegmentVec.2 (Z) are the center of the circle, and diffusion_radius is the circle radius
          // circle_rectangle_intersection returns the overlap area of the circle and the pixel. It is very fast if there is no overlap.
          double striparea_frac = circle_rectangle_intersection(y1, z1, y2, z2, gsl_vector_get(m_SegmentVec, 1), gsl_vector_get(m_SegmentVec, 2), diffusion_radius) / (M_PI * (diffusion_radius * diffusion_radius));
          // assume that the energy is deposited uniformly along the tracklet length, so that this segment gets the fraction 1/nsegments of the energy
          stripenergy[iy - minstrip_y][iz - minstrip_z] += striparea_frac * hiter->second->get_edep() / (float) nsegments;
          if (hiter->second->has_property(PHG4Hit::prop_eion))
          {
            stripeion[iy - minstrip_y][iz - minstrip_z] += striparea_frac * hiter->second->get_eion() / (float) nsegments;
          }
          if (Verbosity() > 5)
          {
            cout << "    strip y index " << iy << " strip z index  " << iz
                 << " strip area fraction of circle " << striparea_frac << " accumulated pixel energy " << stripenergy[iy - minstrip_y][iz - minstrip_z]
                 << endl;
          }
        }
      }
    }  // end loop over segments
    // now we have the energy deposited in each pixel, summed over all tracklet segments. We make a vector of all pixels with non-zero energy deposited
    for (int iz = minstrip_z; iz <= maxstrip_z; iz++)
    {
      for (int iy = minstrip_y; iy <= maxstrip_y; iy++)
      {
        if (stripenergy[iy - minstrip_y][iz - minstrip_z] > 0.0)
        {
          vybin.push_back(iy);
          vzbin.push_back(iz);
          pair<double, double> tmppair = make_pair(stripenergy[iy - minstrip_y][iz - minstrip_z], stripeion[iy - minstrip_y][iz - minstrip_z]);
          venergy.push_back(tmppair);
          if (Verbosity() > 1)
            cout << " Added ybin " << iy << " zbin " << iz << " to vectors with energy " << stripenergy[iy - minstrip_y][iz - minstrip_z] << endl;
        }
      }
    }

    //===================================
    // End of charge sharing implementation
    //===================================

    // Add the strips fired by this hit to the cell list
    //===============================

    for (unsigned int i1 = 0; i1 < vybin.size(); i1++)  // loop over all fired cells
    {
      // this string must be unique - it needs the layer too, or in high multiplicity events it will add g4 hits in different layers with the same key together
      std::string key = (boost::format("%d-%d-%d-%d-%d") % sphxlayer % ladder_z_index % ladder_phi_index % vzbin[i1] % vybin[i1]).str();
      PHG4Cell *cell = nullptr;
      map<string, PHG4Cell *>::iterator it;

      it = m_CellList.find(key);
      // If there is an existing cell to add this hit to, find it
      if (it != m_CellList.end())
      {
        cell = it->second;
        if (Verbosity() > 2)
        {
          cout << " found existing cell with key " << key << endl;
        }
      }

      // There is not an existing cell to add this hit to, start a new cell
      if (!cell)
      {
        if (Verbosity() > 2)
        {
          cout << " did not find existing cell with key " << key << " start a new one" << endl;
        }
        unsigned int index = m_CellList.size();
        index++;
        PHG4CellDefs::keytype cellkey = PHG4CellDefs::MapsBinning::genkey(sphxlayer, index);
        cell = new PHG4Cellv1(cellkey);
        m_CellList[key] = cell;
        // This encodes the z and phi position of the sensor
        //          m_CellList[key]->set_sensor_index((boost::format("%d_%d") %ladder_z_index %ladder_phi_index).str());

        cell->set_ladder_z_index(ladder_z_index);
        cell->set_ladder_phi_index(ladder_phi_index);

        // The z and phi position of the hit strip within the sensor
        cell->set_zbin(vzbin[i1]);
        cell->set_phibin(vybin[i1]);
      }

      // One way or another we have a cell pointer - add this hit to the cell
      cell->add_edep(venergy[i1].first);
      cell->add_edep(hiter->first, venergy[i1].first);  // adds the g4hit association to the cell
      cell->add_eion(venergy[i1].second);
    }

  }  // end loop over g4hits

  int numcells = 0;
  for (std::map<std::string, PHG4Cell *>::const_iterator mapiter = m_CellList.begin(); mapiter != m_CellList.end(); ++mapiter)
  {
    cells->AddCell(mapiter->second);
    numcells++;

    if (Verbosity() > 0)
    {
      std::cout << "Adding cell for "
                << " layer " << mapiter->second->get_layer()
                << " ladder z index: " << mapiter->second->get_ladder_z_index()
                << ", ladder phi index: " << mapiter->second->get_ladder_phi_index()
                << ", srip z: " << mapiter->second->get_zbin()
                << ", strip y: " << mapiter->second->get_phibin()
                << ", energy dep: " << mapiter->second->get_edep()
                << std::endl;
    }
  }
  m_CellList.clear();

  if (Verbosity() > 0)
    std::cout << Name() << ": found " << numcells << " silicon strips with energy deposition" << std::endl;

  if (m_ChkEnergyConservationFlag)
  {
    CheckEnergy(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SiliconTrackerCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, m_CellNodeName);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;

  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    sum_energy_g4hit += hiter->second->get_edep();
  }

  PHG4CellContainer::ConstRange cell_begin_end = cells->getCells();
  PHG4CellContainer::ConstIterator citer;
  for (citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
  {
    sum_energy_cells += citer->second->get_edep();
  }
  // the fractional eloss for particles traversing eta bins leads to minute rounding errors
  if (fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit > 1e-6)
  {
    std::cout << "energy mismatch between cells: " << sum_energy_cells
              << " and hits: " << sum_energy_g4hit
              << " diff sum(cells) - sum(hits): " << sum_energy_cells - sum_energy_g4hit
              << std::endl;

    return -1;
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": total energy for this event: " << sum_energy_g4hit << " GeV" << std::endl;
    }
  }
  return 0;
}

double PHG4SiliconTrackerCellReco::circle_rectangle_intersection(double x1, double y1, double x2, double y2, double mx, double my, double r) const
{
  // Find the area of overlap of a circle and rectangle
  // Calls sA, which uses an analytic formula to determine the integral of the circle between limits set by the corners of the rectangle

  // move the rectangle to the frame where the circle is at (0,0)
  x1 -= mx;
  x2 -= mx;
  y1 -= my;
  y2 -= my;

  if (Verbosity() > 7)
  {
    cout << " mx " << mx << " my " << my << " r " << r << " x1 " << x1 << " x2 " << x2 << " y1 " << y1 << " y2 " << y2 << endl;
    cout << " sA21 " << sA(r, x2, y1)
         << " sA11 " << sA(r, x1, y1)
         << " sA22 " << sA(r, x2, y2)
         << " sA12 " << sA(r, x1, y2)
         << endl;
  }

  return sA(r, x2, y1) - sA(r, x1, y1) - sA(r, x2, y2) + sA(r, x1, y2);
}

double PHG4SiliconTrackerCellReco::sA(double r, double x, double y) const
{
  // Uses analytic formula for the integral of a circle between limits set by the corner of a rectangle
  // It is called repeatedly to find the overlap area between the circle and rectangle
  // I found this code implementing the integral on a web forum called "ars technica",
  // https://arstechnica.com/civis/viewtopic.php?t=306492
  // posted by "memp"

  double a;

  if (x < 0)
  {
    return -sA(r, -x, y);
  }

  if (y < 0)
  {
    return -sA(r, x, -y);
  }

  if (x > r)
  {
    x = r;
  }

  if (y > r)
  {
    y = r;
  }

  if (x * x + y * y > r * r)
  {
    a = r * r * asin(x / r) + x * sqrt(r * r - x * x) + r * r * asin(y / r) + y * sqrt(r * r - y * y) - r * r * M_PI_2;

    a *= 0.5;
  }
  else
  {
    a = x * y;
  }

  return a;
}

void PHG4SiliconTrackerCellReco::SetDefaultParameters()
{
  // if we ever need separate timing windows, don't patch around here!
  // use PHParameterContainerInterface which
  // provides for multiple layers/detector types
  set_default_double_param("tmax", 80.0);   // FVTX NIM paper Fig 32
  set_default_double_param("tmin", -20.0);  // FVTX NIM paper Fig 32, collision has a timing spread around the triggered event. Accepting negative time too.
  return;
}
