#include "PHG4InttHitReco.h"

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeom.h>  // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>  // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <phparameter/PHParameterInterface.h>  // for PHParameterInterface

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

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>      // for _Rb_tree_const_it...
#include <memory>   // for allocator_traits<...
#include <utility>  // for pair, swap, make_...
#include <vector>   // for vector

PHG4InttHitReco::PHG4InttHitReco(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_Detector("INTT")
  , m_Tmin(NAN)
  , m_Tmax(NAN)
  , m_crossingPeriod(NAN)
{
  InitializeParameters();

  m_HitNodeName = "G4HIT_" + m_Detector;
  m_CellNodeName = "G4CELL_" + m_Detector;
  m_GeoNodeName = "CYLINDERGEOM_" + m_Detector;
  m_LocalOutVec = gsl_vector_alloc(3);
  m_PathVec = gsl_vector_alloc(3);
  m_SegmentVec = gsl_vector_alloc(3);
}

PHG4InttHitReco::~PHG4InttHitReco()
{
  gsl_vector_free(m_LocalOutVec);
  gsl_vector_free(m_PathVec);
  gsl_vector_free(m_SegmentVec);
}

int PHG4InttHitReco::InitRun(PHCompositeNode *topNode)
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
    std::cout << Name() << "RUN Node missing, exiting." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode)
  {
    std::cout << Name() << "PAR Node missing, exiting." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  std::string paramnodename = "G4CELLPARAM_" + m_Detector;

  PHNodeIterator runiter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runiter.findFirst("PHCompositeNode", m_Detector));
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

  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    hitsetcontainer = new TrkrHitSetContainerv1;
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    DetNode->addNode(newNode);
  }

  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    hittruthassoc = new TrkrHitTruthAssocv1;
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
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
  m_crossingPeriod = get_double_param("beam_crossing_period");

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InttHitReco::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << m_HitNodeName << std::endl;
    exit(1);
  }

  // Get the TrkrHitSetContainer node
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    std::cout << "Could not locate TRKR_HITSET node, quit! " << std::endl;
    exit(1);
  }

  // Get the TrkrHitTruthAssoc node
  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    std::cout << "Could not locate TRKR_HITTRUTHASSOC node, quit! " << std::endl;
    exit(1);
  }

  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, m_GeoNodeName);
  if (!geo)
  {
    std::cout << "Could not locate geometry node " << m_GeoNodeName << std::endl;
    exit(1);
  }
  // loop over all of the layers in the hit container
  // we need the geometry object for this layer
  if (Verbosity() > 2) std::cout << " PHG4InttHitReco: Loop over hits" << std::endl;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  for (PHG4HitContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    const int sphxlayer = hiter->second->get_detid();
    CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(geo->GetLayerGeom(sphxlayer));

    // checking ADC timing integration window cut
    // uses default values for now
    // these should depend on layer radius
    if (hiter->second->get_t(0) > m_Tmax) continue;
    if (hiter->second->get_t(1) < m_Tmin) continue;

    float time = (hiter->second->get_t(0) + hiter->second->get_t(1)) / 2.0;

    // I made this (small) diffusion up for now, we will get actual values for the Intt later
    double diffusion_width = 5.0e-04;  // diffusion radius 5 microns, in cm

    const int ladder_z_index = hiter->second->get_ladder_z_index();
    const int ladder_phi_index = hiter->second->get_ladder_phi_index();

    // What we have is a hit in the sensor. We have not yet assigned the strip(s) that were hit, we do that here
    //========================================================================

// initialize them. In case find_strip_index_values does not set them we can pick this up
    int strip_y_index_in = -99999;
    int strip_z_index_in = -99999;
    int strip_y_index_out = -99999;
    int strip_z_index_out = -99999;

    layergeom->find_strip_index_values(ladder_z_index, hiter->second->get_local_y(0), hiter->second->get_local_z(0), strip_y_index_in, strip_z_index_in);
    layergeom->find_strip_index_values(ladder_z_index, hiter->second->get_local_y(1), hiter->second->get_local_z(1), strip_y_index_out, strip_z_index_out);
    if (strip_y_index_in ==  -99999 ||
        strip_z_index_in == -99999 ||
        strip_y_index_out == -99999 ||
	strip_z_index_out == -99999)
    {
      std::cout << "setting of strip indices failed" << std::endl;
      std::cout << "strip_y_index_in: " << strip_y_index_in << std::endl;
      std::cout << "strip_z_index_in: " << strip_z_index_in << std::endl;
      std::cout << "strip_y_index_out: " << strip_y_index_out << std::endl;
      std::cout << "strip_z_index_out: " << strip_y_index_out << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
    if (Verbosity() > 5)
    {
      // check to see if we get back the positions from these strip index values
      double check_location[3] = {-1, -1, -1};
      layergeom->find_strip_center_localcoords(ladder_z_index, strip_y_index_in, strip_z_index_in, check_location);
      std::cout << " G4 entry location = " << hiter->second->get_local_x(0) << "  " << hiter->second->get_local_y(0) << "  " << hiter->second->get_local_z(0) << std::endl;
      std::cout << " Check entry location = " << check_location[0] << "  " << check_location[1] << "  " << check_location[2] << std::endl;
      layergeom->find_strip_center_localcoords(ladder_z_index, strip_y_index_out, strip_z_index_out, check_location);
      std::cout << " G4 exit location = " << hiter->second->get_local_x(1) << " " << hiter->second->get_local_y(1) << "  " << hiter->second->get_local_z(1) << std::endl;
      std::cout << " Check exit location = " << check_location[0] << "  " << check_location[1] << "  " << check_location[2] << std::endl;
    }

    // Now we find how many strips were crossed by this track, and divide the energy between them
    int minstrip_z = strip_z_index_in;
    int maxstrip_z = strip_z_index_out;
    if (minstrip_z > maxstrip_z) std::swap(minstrip_z, maxstrip_z);

    int minstrip_y = strip_y_index_in;
    int maxstrip_y = strip_y_index_out;
    if (minstrip_y > maxstrip_y) std::swap(minstrip_y, maxstrip_y);

    // Use an algorithm similar to the one for the MVTX pixels, since it facilitates adding charge diffusion
    // for now we assume small charge diffusion
    std::vector<int> vybin;
    std::vector<int> vzbin;
    //std::vector<double> vlen;
    std::vector<std::pair<double, double> > venergy;

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
    double stripenergy[13][13] = {};  // init to 0
    double stripeion[13][13] = {};    // init to 0

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
        std::cout << " segment " << i
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
                  << " segvec.Y " << gsl_vector_get(m_SegmentVec, 1) << std::endl
                  << " diffusion_radius " << diffusion_radius
                  << std::endl;

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
          double striparea_frac = PHG4Utils::circle_rectangle_intersection(y1, z1, y2, z2, gsl_vector_get(m_SegmentVec, 1), gsl_vector_get(m_SegmentVec, 2), diffusion_radius) / (M_PI * (diffusion_radius * diffusion_radius));
          // assume that the energy is deposited uniformly along the tracklet length, so that this segment gets the fraction 1/nsegments of the energy
          stripenergy[iy - minstrip_y][iz - minstrip_z] += striparea_frac * hiter->second->get_edep() / (float) nsegments;
          if (hiter->second->has_property(PHG4Hit::prop_eion))
          {
            stripeion[iy - minstrip_y][iz - minstrip_z] += striparea_frac * hiter->second->get_eion() / (float) nsegments;
          }
          if (Verbosity() > 5)
          {
            std::cout << "    strip y index " << iy << " strip z index  " << iz
                      << " strip area fraction of circle " << striparea_frac << " accumulated pixel energy " << stripenergy[iy - minstrip_y][iz - minstrip_z]
                      << std::endl;
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
          std::pair<double, double> tmppair = std::make_pair(stripenergy[iy - minstrip_y][iz - minstrip_z], stripeion[iy - minstrip_y][iz - minstrip_z]);
          venergy.push_back(tmppair);
          if (Verbosity() > 1)
          {
            std::cout << " Added ybin " << iy << " zbin " << iz << " to vectors with energy " << stripenergy[iy - minstrip_y][iz - minstrip_z] << std::endl;
          }
        }
      }
    }

    //===================================
    // End of charge sharing implementation
    //===================================

    for (unsigned int i1 = 0; i1 < vybin.size(); i1++)  // loop over all fired cells
    {
      // We add the Intt TrkrHitsets directly to the node using hitsetcontainer

      // Get the hit crossing
      int crossing = (int) (round( time / m_crossingPeriod) );
      // crossing has to fit into 5 bits
      if(crossing < -512) crossing = -512;
      if(crossing > 511) crossing = 511;
      // We need to create the TrkrHitSet if not already made - each TrkrHitSet should correspond to a sensor for the Intt ?
      // The hitset key includes the layer, the ladder_z_index (sensors numbered 0-3) and  ladder_phi_index (azimuthal location of ladder) for this hit
      TrkrDefs::hitsetkey hitsetkey = InttDefs::genHitSetKey(sphxlayer, ladder_z_index, ladder_phi_index, crossing);

      // Use existing hitset or add new one if needed
      TrkrHitSetContainer::Iterator hitsetit = hitsetcontainer->findOrAddHitSet(hitsetkey);

      // generate the key for this hit
      TrkrDefs::hitkey hitkey = InttDefs::genHitKey(vzbin[i1], vybin[i1]);
      // See if this hit already exists
      TrkrHit *hit = hitsetit->second->getHit(hitkey);
      if (!hit)
      {
        // Otherwise, create a new one
	hit = new TrkrHitv2();
        hitsetit->second->addHitSpecificKey(hitkey, hit);
      }

      // Either way, add the energy to it
      if (Verbosity() > 2)
      {
        std::cout << "add energy " << venergy[i1].first << " to intthit " << std::endl;
      }
      hit->addEnergy(venergy[i1].first * TrkrDefs::InttEnergyScaleup);

      // Add this hit to the association map
      hittruthassoc->addAssoc(hitsetkey, hitkey, hiter->first);

      if (Verbosity() > 2)
      {
        std::cout << "PHG4InttHitReco: added hit wirh hitsetkey " << hitsetkey << " hitkey " << hitkey << " g4hitkey " << hiter->first << " energy " << hit->getEnergy() << std::endl;
      }
    }
    
  }  // end loop over g4hits
  
  // print the list of entries in the association table
  if (Verbosity() > 0)
  {
    std::cout << "From PHG4InttHitReco: " << std::endl;
    hitsetcontainer->identify();
    hittruthassoc->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4InttHitReco::SetDefaultParameters()
{
  // if we ever need separate timing windows, don't patch around here!
  // use PHParameterContainerInterface which
  // provides for multiple layers/detector types
  set_default_double_param("tmax", 7020.0);   // max upper time window for extended readout
  set_default_double_param("tmin", -20.0);  // min lower time window for extended readout
  set_default_double_param("beam_crossing_period", 106.0);   

  return;
}
