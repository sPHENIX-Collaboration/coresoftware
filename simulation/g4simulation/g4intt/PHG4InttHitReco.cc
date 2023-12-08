#include "PHG4InttHitReco.h"

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeom.h>  // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainerv1.h>
#include <trackbase/ClusHitsVerbosev1.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>  // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <phparameter/PHParameterInterface.h>  // for PHParameterInterface

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Utils.h>
#include <g4main/PHG4TruthInfoContainer.h>

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


// update to make sure to clusterize clusters in loopers 

PHG4InttHitReco::PHG4InttHitReco(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_Detector("INTT")
  , m_Tmin(NAN)
  , m_Tmax(NAN)
  , m_crossingPeriod(NAN)
  , m_truth_hits { new TrkrHitSetContainerv1 }
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
  delete m_truth_hits;
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

  // get the nodes for the truth clustering
  m_truthtracks = findNode::getClass<TrkrTruthTrackContainer>(topNode, "TRKR_TRUTHTRACKCONTAINER");
  if (!m_truthtracks)
  {
    PHNodeIterator dstiter(dstNode);
    m_truthtracks = new TrkrTruthTrackContainerv1();
    auto newNode = new PHIODataNode<PHObject>(m_truthtracks, "TRKR_TRUTHTRACKCONTAINER", "PHObject");
    dstNode->addNode(newNode);
  }

  m_truthclusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!m_truthclusters)
  {
    m_truthclusters = new TrkrClusterContainerv4;
    auto newNode = new PHIODataNode<PHObject>(m_truthclusters, "TRKR_TRUTHCLUSTERCONTAINER", "PHObject");
    dstNode->addNode(newNode);
  }

  m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthinfo)
  {
    std::cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << std::endl;
    assert(m_truthinfo);
  }

  // get cluster hits verbose (the hit and energy) information
  if (record_ClusHitsVerbose) {
    // get the node
    mClusHitsVerbose = findNode::getClass<ClusHitsVerbosev1>(topNode, "Trkr_TruthClusHitsVerbose");
    if (!mClusHitsVerbose)
    {
      PHNodeIterator dstiter(dstNode);
      auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode("TRKR");
        dstNode->addNode(DetNode);
      }
      mClusHitsVerbose = new ClusHitsVerbosev1();
      auto newNode = new PHIODataNode<PHObject>(mClusHitsVerbose, "Trkr_TruthClusHitsVerbose", "PHObject");
      DetNode->addNode(newNode);
    }
  }

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

    truthcheck_g4hit(hiter->second, topNode);

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

      double hit_energy = venergy[i1].first * TrkrDefs::InttEnergyScaleup;
      hit->addEnergy(hit_energy);

      addtruthhitset(hitsetkey, hitkey, hit_energy);

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

  if (m_is_emb) {
    cluster_truthhits(topNode); // the last track was truth -- make it's clusters
    prior_g4hit = nullptr;
  }
  
  end_event_truthcluster( topNode );
  return Fun4AllReturnCodes::EVENT_OK;
} // end process_event

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

void PHG4InttHitReco::truthcheck_g4hit(PHG4Hit* g4hit, PHCompositeNode* topNode) {
  if (g4hit==nullptr) return;
  int new_trkid = g4hit->get_trkid();

  bool is_new_track = (new_trkid != m_trkid);
  if (Verbosity()>5) std::cout << PHWHERE << std::endl << " -> Checking status of PHG4Hit. Track id("<<new_trkid<<")" << std::endl;
  if (!is_new_track) {
    // check to see if it is an embedded track that meets the looper condition:
    if (m_is_emb) {
      if (prior_g4hit!=nullptr
          && (    std::abs(prior_g4hit->get_x(0)-g4hit->get_x(0)) > max_g4hitstep
               || std::abs(prior_g4hit->get_y(0)-g4hit->get_y(0)) > max_g4hitstep
             )
          )
      {
          // this is a looper track -- cluster hits up to this point already
          cluster_truthhits(topNode);
      }
      prior_g4hit = g4hit;
    }
    return;
  }
  // <- STATUS: this is a new track
  if (Verbosity()>2) std::cout << PHWHERE << std::endl << " -> Found new embedded track with id: " << new_trkid << std::endl;
  if (m_is_emb) {
    //cluster the old track
    cluster_truthhits(topNode); // cluster m_truth_hits and add m_current_track
    m_current_track = nullptr;
    prior_g4hit = nullptr;
  }
  m_trkid = new_trkid;
  m_is_emb = m_truthinfo->isEmbeded(m_trkid);
  if (m_is_emb) {
    m_current_track = m_truthtracks->getTruthTrack(m_trkid, m_truthinfo);
    prior_g4hit = g4hit;
  }
}

void PHG4InttHitReco::end_event_truthcluster ( PHCompositeNode* topNode ) {
  if (m_is_emb) {
    cluster_truthhits(topNode); // cluster m_truth_hits and add m_current_track
    m_current_track = nullptr;
    m_trkid = -1;
    m_is_emb = false;
  }
  m_hitsetkey_cnt.clear();
}

void PHG4InttHitReco::addtruthhitset(
    TrkrDefs::hitsetkey hitsetkey, 
    TrkrDefs::hitkey hitkey, 
    float neffelectrons) 
{
  if (!m_is_emb) return;
  TrkrHitSetContainer::Iterator hitsetit = m_truth_hits->findOrAddHitSet(hitsetkey);
  // See if this hit already exists
  TrkrHit *hit = nullptr;
  hit = hitsetit->second->getHit(hitkey);
  if (!hit)
  {
    // create a new one
    hit = new TrkrHitv2();
    hitsetit->second->addHitSpecificKey(hitkey, hit);
  }
  // Either way, add the energy to it  -- adc values will be added at digitization
  hit->addEnergy(neffelectrons);
}

void PHG4InttHitReco::cluster_truthhits(PHCompositeNode* topNode) {
  // -----------------------------------------------
  // Digitize, adapted from g4intt/PHG4InttDigitizer
  // -----------------------------------------------
  //
  // Note: not using digitization, because as currently implemented, the SvtxTrack clusters
  // don't use the adc weighting from the digitization code anyway.
  //
  // don't use the dead map for truth tracks
  /* TrkrHitSetContainer::ConstRange hitset_range = m_truth_hits->getHitSets(TrkrDefs::TrkrId::inttId); */
  /* for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first; */
  /*      hitset_iter != hitset_range.second; */
  /*      ++hitset_iter) */
  /* { */
  /*   // we have an itrator to one TrkrHitSet for the intt from the trkrHitSetContainer */
  /*   // get the hitset key so we can find the layer */
  /*   TrkrDefs::hitsetkey hitsetkey = hitset_iter->first; */
  /*   const int layer = TrkrDefs::getLayer(hitsetkey); */
  /*   const int ladder_phi = InttDefs::getLadderPhiId(hitsetkey); */
  /*   const int ladder_z = InttDefs::getLadderZId(hitsetkey); */

  /*   if (Verbosity() > 1) */
  /*   { */
  /*     std::cout << "PHG4InttDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << std::endl; */
  /*   } */
  /*   // get all of the hits from this hitset */
  /*   TrkrHitSet *hitset = hitset_iter->second; */
  /*   TrkrHitSet::ConstRange hit_range = hitset->getHits(); */
  /*   /1* std::set<TrkrDefs::hitkey> dead_hits;  // hits on dead channel *1/ // no dead channels implemented */
  /*   for (TrkrHitSet::ConstIterator hit_iter = hit_range.first; */
  /*        hit_iter != hit_range.second; */
  /*        ++hit_iter) */
  /*   { */
  /*     // ++m_nCells; // not really used by PHG4InttDigitizer */

  /*     TrkrHit *hit = hit_iter->second; */
  /*     TrkrDefs::hitkey hitkey = hit_iter->first; */
  /*     int strip_col = InttDefs::getCol(hitkey);  // strip z index */
  /*     int strip_row = InttDefs::getRow(hitkey);  // strip phi index */

  /*     // FIXME need energy scales here */
  /*     if (_energy_scale.count(layer) > 1) */
  /*     { */
  /*       assert(!"Error: _energy_scale has two or more keys."); */
  /*     } */
  /*     const float mip_e = _energy_scale[layer]; */

  /*     std::vector<std::pair<double, double> > vadcrange = _max_fphx_adc[layer]; */

  /*     int adc = 0; */
  /*     for (unsigned int irange = 0; irange < vadcrange.size(); ++irange) */
  /*     { */
  /*       if (hit->getEnergy() / TrkrDefs::InttEnergyScaleup >= vadcrange[irange].first * */ 
  /*           (double) mip_e && hit->getEnergy() / TrkrDefs::InttEnergyScaleup */ 
  /*           < vadcrange[irange].second * (double) mip_e) */
  /*       { */
  /*         adc = (unsigned short) irange; */
  /*       } */
  /*     } */
  /*     hit->setAdc(adc); */

  /*     if (Verbosity() > 2) */
  /*     { */
  /*       std::cout << "PHG4InttDigitizer: found hit with layer " << layer << " ladder_z " << ladder_z << " ladder_phi " << ladder_phi */
  /*                 << " strip_col " << strip_col << " strip_row " << strip_row << " adc " << hit->getAdc() << std::endl; */
  /*     } */
  /*   }  // end loop over hits in this hitset */

    /* // remove hits on dead channel in TRKR_HITSET and TRKR_HITTRUTHASSOC */
    /* for (const auto &key : dead_hits) */
    /* { */
    /*   if (Verbosity() > 2) */
    /*   { */
    /*     std::cout << " PHG4InttDigitizer: remove hit with key: " << key << std::endl; */
    /*   } */
    /*   hitset->removeHit(key); */
    /* } */
  /* }  // end loop over hitsets */

  // -----------------------------------------------
  // Cluster, adapted from intt/InttClusterizer
  // -----------------------------------------------
  if (Verbosity() > 1) std::cout << "Clustering truth clusters" << std::endl;

  //-----------
  // Clustering
  //-----------
  // get the geometry node
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container) return;


  // loop over the InttHitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange =
      m_truth_hits->getHitSets(TrkrDefs::TrkrId::inttId); // from TruthClusterizerBase

  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    // Each hitset contains only hits that are clusterizable - i.e. belong to a single sensor
    TrkrHitSet*         hitset    = hitsetitr->second;
    TrkrDefs::hitsetkey hitsetkey = hitset->getHitSetKey();

    // cluster this hitset; all pixels in it are, by definition, part of the same clusters

    if ( Verbosity() > 1 ) std::cout << "InttClusterizer found hitsetkey " << hitsetitr->first << std::endl;
    if ( Verbosity() > 2 ) hitset->identify();

    // we have a single hitset, get the info that identifies the sensor

    if (Verbosity() > 2)
      std::cout << "Filling cluster with hitsetkey " << ((int)hitsetkey) << std::endl;

    // get the bunch crossing number from the hitsetkey
    /* short int crossing = InttDefs::getTimeBucketId(hitset->getHitSetKey()); */

    // determine the size of the cluster in phi and z, useful for track fitting the cluster
    std::set<int> phibins;
    std::set<int> zbins;

    // determine the cluster position...
    double       xlocalsum   = 0.0;
    double       ylocalsum   = 0.0;
    double       zlocalsum   = 0.0;
    unsigned int clus_adc = 0.0;
    unsigned     nhits       = 0;

    // aggregate the adc values 
    double sum_adc {0};
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for ( auto ihit = hitrangei.first; ihit != hitrangei.second; ++ihit) {
      sum_adc += ihit->second->getAdc();
    }

    // tune this energy threshold in the same maner of the MVTX, namely to get the same kind of pixel sizes
    // as the SvtxTrack clusters
    /* const double threshold = sum_adc * m_truth_pixelthreshold; */
    const double threshold = sum_adc * m_pixel_thresholdrat; //FIXME -- tune this as needed
    std::map<int,unsigned int> m_iphi, m_it, m_iphiCut, m_itCut; // FIXME

    int layer = TrkrDefs::getLayer ( hitsetkey );
    CylinderGeomIntt* geom = dynamic_cast<CylinderGeomIntt*>(geom_container->GetLayerGeom(layer));

    int ladder_z_index = InttDefs::getLadderZId ( hitsetkey );

    for ( auto ihit = hitrangei.first; ihit != hitrangei.second; ++ihit)
    {
      int col = InttDefs::getCol ( ihit->first );
      int row = InttDefs::getRow ( ihit->first );
      auto adc = ihit->second->getAdc();

      if (mClusHitsVerbose) {
        std::map<int,unsigned int>& m_phi = (adc<threshold) ? m_iphiCut : m_iphi;
        std::map<int,unsigned int>& m_z   = (adc<threshold) ? m_itCut   : m_it;

        auto pnew = m_phi.try_emplace(row,adc);
        if (!pnew.second) pnew.first->second += adc;

        pnew = m_z.try_emplace(col,adc);
        if (!pnew.second) pnew.first->second += adc;
      }
      if (adc<threshold) continue;

      clus_adc += adc;
      zbins   .insert(col);
      phibins .insert(row);

      // now get the positions from the geometry
      double local_hit_location[3] = {0., 0., 0.};

      geom->find_strip_center_localcoords(ladder_z_index, row, col, local_hit_location);

      xlocalsum += local_hit_location[0];
      ylocalsum += local_hit_location[1];
      zlocalsum += local_hit_location[2];

      ++nhits;

      if (Verbosity() > 6)
      {
        std::cout << "  From  geometry object: hit x " << local_hit_location[0] 
          << " hit y " << local_hit_location[1] << " hit z " << local_hit_location[2] << std::endl;
        std::cout << "     nhits " << nhits << " clusx  = " << xlocalsum / nhits << " clusy " 
          << ylocalsum / nhits << " clusz " << zlocalsum / nhits << std::endl;
      }
      // NOTE:
      /* if (_make_e_weights[layer]) */ // these values are all false by default
      /* if ( false ) // the current implementation of the code does not weight by adc values */
      /*              // therefore the default here is to use use adc to cut the outliers and nothing else */
      /* { */
      /*   xlocalsum += local_hit_location[0] * (double) hit_adc; */
      /*   ylocalsum += local_hit_location[1] * (double) hit_adc; */
      /*   zlocalsum += local_hit_location[2] * (double) hit_adc; */
      /* } */
      /* else */
      /* { */
      /* } */
      /* if(hit_adc > clus_maxadc) clus_maxadc = hit_adc; */ //FIXME: do we want this value to be set?
      /* clus_energy += hit_adc; */
    }
    if (mClusHitsVerbose) {
      if (Verbosity()>10) {
        for (auto& hit : m_iphi) {
          std::cout << " m_phi(" << hit.first <<" : " << hit.second<<") " << std::endl;
        }
      }
      for (auto& hit : m_iphi)    mClusHitsVerbose->addPhiHit    (hit.first, (float)hit.second);
      for (auto& hit : m_it)      mClusHitsVerbose->addZHit      (hit.first, (float)hit.second);
      for (auto& hit : m_iphiCut) mClusHitsVerbose->addPhiCutHit (hit.first, (float)hit.second);
      for (auto& hit : m_itCut)   mClusHitsVerbose->addZCutHit   (hit.first, (float)hit.second);
    }


    // add this cluster-hit association to the association map of (clusterkey,hitkey)
    if (Verbosity() > 2) std::cout << "  nhits = " << nhits << std::endl;

    /* static const float invsqrt12 = 1./sqrt(12); */
    // scale factors (phi direction)
    /*
       they corresponds to clusters of size 1 and 2 in phi
       other clusters, which are very few and pathological, get a scale factor of 1
       These scale factors are applied to produce cluster pulls with width unity
       */

    /* float phierror = pitch * invsqrt12; */

    /* static constexpr std::array<double, 3> scalefactors_phi = {{ 0.85, 0.4, 0.33 }}; */
    /* if( phibins.size() == 1 && layer < 5) phierror*=scalefactors_phi[0]; */
    /* else if( phibins.size() == 2 && layer < 5) phierror*=scalefactors_phi[1]; */
    /* else if( phibins.size() == 2 && layer > 4) phierror*=scalefactors_phi[2]; */ 
    /* // z error. All clusters have a z-size of 1. */
    /* const float zerror = length * invsqrt12; */
    if (nhits == 0) continue;

    double cluslocaly = ylocalsum / nhits;
    double cluslocalz = zlocalsum / nhits;

    //if (_make_e_weights[layer]) // FIXME: this is always false for now
    /* { */
    /*   cluslocaly = ylocalsum / (double) clus_adc; */
    /*   cluslocalz = zlocalsum / (double) clus_adc; */
    /* } */
    /* else */
    /* { */
    /* } */
    if ( m_cluster_version==4 ){
      auto clus = std::make_unique<TrkrClusterv4>();
      clus->setAdc(clus_adc);
      clus->setPhiSize(phibins.size());
      clus->setZSize(1);

      if(Verbosity() > 10) clus->identify();

      clus->setLocalX(cluslocaly);
      clus->setLocalY(cluslocalz);
      // silicon has a 1-1 map between hitsetkey and surfaces. So set to 0
      clus->setSubSurfKey(0);

      m_hitsetkey_cnt.try_emplace(hitsetkey,0);
      unsigned int& cnt = m_hitsetkey_cnt[hitsetkey];
      TrkrDefs::cluskey ckey = TrkrDefs::genClusKey(hitsetkey, cnt);
      m_truthclusters->addClusterSpecifyKey(ckey, clus.release());
      m_current_track->addCluster(ckey);
      if (mClusHitsVerbose) {
        mClusHitsVerbose->push_hits(ckey);
        if (Verbosity()>10) std::cout << " ClusHitsVerbose.size (in INTT): " 
          << mClusHitsVerbose->getMap().size() << std::endl;
      }
      ++cnt;
    }  // end loop over hitsets
  }

  m_truth_hits->Reset(); 
  prior_g4hit = nullptr;
  return;
}
