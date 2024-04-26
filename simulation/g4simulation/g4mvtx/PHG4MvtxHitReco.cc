// this is the new trackbase version

#include "PHG4MvtxHitReco.h"

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxHitPruner.h>

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>  // // make iwyu happy
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>  // make iwyu happy
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssoc.h>  // make iwyu happy
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <g4tracking/TrkrTruthTrackContainerv1.h>
#include <trackbase/ClusHitsVerbosev1.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv4.h>

#include <g4detectors/PHG4CylinderGeom.h>  // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>      // for PHObject
#include <phool/PHRandomSeed.h>  // for PHRandomSeed
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <G4SystemOfUnits.hh>  // for microsecond

#include <TVector3.h>  // for TVector3, ope...

#include <cassert>  // for assert
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>  // for allocator_tra...
#include <set>     // for vector
#include <vector>  // for vector

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>  // for accessing the MVTX hot pixel file from the CDB
#include <ffarawobjects/MvtxRawEvtHeader.h>

// New headers I added

PHG4MvtxHitReco::PHG4MvtxHitReco(const std::string& name, const std::string& detector)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_detector(detector)
  , m_tmin(-5000.)
  , m_tmax(5000.)
  , m_strobe_width(5.)
  , m_strobe_separation(0.)
  , m_truth_hits{new TrkrHitSetContainerv1}
{
  if (Verbosity())
  {
    std::cout << "Creating PHG4MvtxHitReco for detector = " << detector << std::endl;
  }
  // initialize rng
  const uint seed = PHRandomSeed();
  m_rng.reset(gsl_rng_alloc(gsl_rng_mt19937));
  gsl_rng_set(m_rng.get(), seed);

  InitializeParameters();
}

int PHG4MvtxHitReco::InitRun(PHCompositeNode* topNode)
{
  UpdateParametersWithMacro();

  m_tmin = get_double_param("mvtx_tmin");
  m_tmax = get_double_param("mvtx_tmax");
  m_strobe_width = get_double_param("mvtx_strobe_width");
  m_strobe_separation = get_double_param("mvtx_strobe_separation");
  m_in_sphenix_srdo = (bool) get_int_param("mvtx_in_sphenix_srdo");

  m_extended_readout_time = m_tmax - m_strobe_width;

  // printout
  std::cout
      << "PHG4MvtxHitReco::InitRun\n"
      << " m_tmin: " << m_tmin << "ns, m_tmax: " << m_tmax << "ns\n"
      << " m_strobe_width: " << m_strobe_width << "\n"
      << " m_strobe_separation: " << m_strobe_separation << "\n"
      << " m_extended_readout_time: " << m_extended_readout_time << "\n"
      << " m_in_sphenix_srdo: " << (m_in_sphenix_srdo ? "true" : "false") << "\n"
      << std::endl;

  //! get DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  assert(dstNode);

  //! get detector run node
  auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  assert(runNode);

  PHNodeIterator runiter(runNode);
  auto runDetNode = dynamic_cast<PHCompositeNode*>(runiter.findFirst("PHCompositeNode", m_detector));
  if (!runDetNode)
  {
    runDetNode = new PHCompositeNode(m_detector);
    runNode->addNode(runDetNode);
  }
  std::string paramNodeName = "G4CELLPARAM_" + m_detector;
  SaveToNodeTree(runDetNode, paramNodeName);

  // create hitset container if needed
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    hitsetcontainer = new TrkrHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    trkrnode->addNode(newNode);
  }

  // create hit truth association if needed
  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    hittruthassoc = new TrkrHitTruthAssocv1;
    auto newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
    trkrnode->addNode(newNode);
  }

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

  if (record_ClusHitsVerbose)
  {
    // get the node
    mClusHitsVerbose = findNode::getClass<ClusHitsVerbosev1>(topNode, "Trkr_TruthClusHitsVerbose");
    if (!mClusHitsVerbose)
    {
      PHNodeIterator dstiter(dstNode);
      auto DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
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

  std::string database = CDBInterface::instance()->getUrl(
      "MVTX_DeadPixelMap");  // This is specifically for MVTX Dead Pixels
  CDBTTree* cdbttree = new CDBTTree(database);

  int NPixel = -1;
  NPixel = cdbttree->GetSingleIntValue("TotalDeadPixels");

  for (int i = 0; i < NPixel; i++)
  {
    int Layer = cdbttree->GetIntValue(i, "layer");
    int Stave = cdbttree->GetIntValue(i, "stave");
    int Chip = cdbttree->GetIntValue(i, "chip");
    int Col = cdbttree->GetIntValue(i, "col");
    int Row = cdbttree->GetIntValue(i, "row");
    if (Verbosity() > 0)
    {
      std::cout << "Layer: " << Layer << ", Stave: " << Stave << ", Chip: " << Chip << ", Row: " << Row << ", Col: " << Col << std::endl;
    }

    TrkrDefs::hitsetkey DeadPixelHitKey = MvtxDefs::genHitSetKey(Layer, Stave, Chip, 0);
    TrkrDefs::hitkey DeadHitKey = MvtxDefs::genHitKey(Col, Row);
    m_deadPixelMap.push_back({std::make_pair(DeadPixelHitKey, DeadHitKey)});
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4MvtxHitReco::process_event(PHCompositeNode* topNode)
{
  ActsGeometry* tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tgeometry)
  {
    std::cout << "Could not locate acts geometry" << std::endl;
    exit(1);
  }

  // G4Hits
  const std::string g4hitnodename = "G4HIT_" + m_detector;
  auto g4hitContainer = findNode::getClass<PHG4HitContainer>(topNode, g4hitnodename);
  assert(g4hitContainer);

  // geometry
  const std::string geonodename = "CYLINDERGEOM_" + m_detector;
  auto geoNode = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename);
  assert(geoNode);

  // Get the TrkrHitSetContainer node
  auto trkrHitSetContainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrHitSetContainer);

  // Get the TrkrHitTruthAssoc node
  auto hitTruthAssoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  assert(hitTruthAssoc);

  // Generate strobe zero relative to trigger time
  double strobe_zero_tm_start = generate_strobe_zero_tm_start();

  // assumes we want the range of accepted times to be from 0 to m_extended_readout_time
  std::pair<double, double> alpide_pulse = generate_alpide_pulse(0.0);  // this currently just returns fixed values
  double clearance = 200.0;                                             // 0.2 microsecond for luck
  m_tmax = m_extended_readout_time + alpide_pulse.first + clearance;
  m_tmin = alpide_pulse.second - clearance;

  // The above limits will select g4hit times of 0 up to m_extended_readout_time (only) with extensions by clearance
  // But we really want to select all g4hit times that will be strobed, so replace clearance with something derived from
  // the strobe start time in future

  if (Verbosity() > 0)
  {
    std::cout << " m_strobe_width " << m_strobe_width << " m_strobe_separation " << m_strobe_separation << " strobe_zero_tm_start " << strobe_zero_tm_start << " m_extended_readout_time " << m_extended_readout_time << std::endl;
  }

  // loop over all of the layers in the g4hit container
  auto layer_range = g4hitContainer->getLayers();
  for (auto layer_it = layer_range.first; layer_it != layer_range.second; ++layer_it)
  {
    // get layer
    const auto layer = *layer_it;
    assert(layer < 3);

    // we need the geometry object for this layer
    auto layergeom = dynamic_cast<CylinderGeom_Mvtx*>(geoNode->GetLayerGeom(layer));
    assert(layergeom);

    // loop over the hits in this layer
    const PHG4HitContainer::ConstRange g4hit_range = g4hitContainer->getHits(layer);

    // Get some layer parameters for later use
    double xpixw = layergeom->get_pixel_x();
    double xpixw_half = xpixw / 2.0;
    double zpixw = layergeom->get_pixel_z();
    double zpixw_half = zpixw / 2.0;
    int maxNX = layergeom->get_NX();
    int maxNZ = layergeom->get_NZ();

    // Now loop over all g4 hits for this layer
    for (auto g4hit_it = g4hit_range.first; g4hit_it != g4hit_range.second; ++g4hit_it)
    {
      // get hit
      auto g4hit = g4hit_it->second;

      truthcheck_g4hit(g4hit, topNode);

      // cout << "From PHG4MvtxHitReco: Call hit print method: " << std::endl;
      if (Verbosity() > 4)
      {
        g4hit->print();
      }

      unsigned int n_replica = 1;

      // Function returns ns
      // std::pair<double, double> alpide_pulse = generate_alpide_pulse(g4hit->get_edep());

      double lead_edge = g4hit->get_t(0) * ns + alpide_pulse.first;
      double fall_edge = g4hit->get_t(1) * ns + alpide_pulse.second;

      if (Verbosity() > 0)
      {
        std::cout << " MvtxHitReco: t0 " << g4hit->get_t(0) << " t1 " << g4hit->get_t(1) << " lead_edge " << lead_edge
                  << " fall_edge " << fall_edge << " tmin " << m_tmin << " tmax " << m_tmax << std::endl;
      }

      // check that the signal occurred witin the time window 0 to extended_readout_time, discard if not
      if (lead_edge > m_tmax or fall_edge < m_tmin)
      {
        continue;
      }

      double t0_strobe_frame = get_strobe_frame(lead_edge, strobe_zero_tm_start);
      double t1_strobe_frame = get_strobe_frame(fall_edge, strobe_zero_tm_start);
      n_replica += t1_strobe_frame - t0_strobe_frame;

      if (Verbosity() > 1)
      {
        std::cout
            << "MVTX is in strobed timing mode\n"
            << "layer " << layer << " t0(ns) " << g4hit->get_t(0) << " t1(ns) " << g4hit->get_t(1) << "\n"
            << "strobe_zero_tm_start(us): " << strobe_zero_tm_start / microsecond << "\n"
            << "strobe width(us): " << m_strobe_width / microsecond << "\n"
            << "strobe separation(us): " << m_strobe_separation / microsecond << "\n"
            << "alpide pulse start(us): " << alpide_pulse.first / microsecond << "\n"
            << "alpide pulse end(us): " << alpide_pulse.second / microsecond << "\n"
            << "tm_zero_strobe_frame: " << t0_strobe_frame << "\n"
            << "tm_one_strobe_frame: " << t1_strobe_frame << "\n"
            << "number of hit replica: " << n_replica << "\n"
            << std::endl;
      }

      assert(n_replica >= 1);

      // get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
      int stave_number = g4hit->get_property_int(PHG4Hit::prop_stave_index);
      int chip_number = g4hit->get_property_int(PHG4Hit::prop_chip_index);

      TVector3 local_in(g4hit->get_local_x(0), g4hit->get_local_y(0), g4hit->get_local_z(0));
      TVector3 local_out(g4hit->get_local_x(1), g4hit->get_local_y(1), g4hit->get_local_z(1));
      TVector3 midpoint((local_in.X() + local_out.X()) / 2.0, (local_in.Y() + local_out.Y()) / 2.0, (local_in.Z() + local_out.Z()) / 2.0);

      if (Verbosity() > 2)
      {
        std::cout
            << " world entry point position: "
            << g4hit->get_x(0) << " " << g4hit->get_y(0) << " " << g4hit->get_z(0) << "\n"
            << " world exit  point position: "
            << g4hit->get_x(1) << " " << g4hit->get_y(1) << " " << g4hit->get_z(1) << "\n"
            << " local coords of entry point from G4 "
            << g4hit->get_local_x(0) << " " << g4hit->get_local_y(0) << " " << g4hit->get_local_z(0)
            << std::endl;

        TVector3 world_in(g4hit->get_x(0), g4hit->get_y(0), g4hit->get_z(0));
        auto hskey = MvtxDefs::genHitSetKey(layer, stave_number, chip_number, 0);
        auto surf = tgeometry->maps().getSiliconSurface(hskey);

        TVector3 local_in_check = layergeom->get_local_from_world_coords(surf, tgeometry, world_in);
        std::cout
            << " local coords of entry point from geom (a check) "
            << local_in_check.X() << " " << local_in_check.Y() << " " << local_in_check.Z() << "\n"
            << " local coords of exit point from G4 "
            << g4hit->get_local_x(1) << " " << g4hit->get_local_y(1) << " " << g4hit->get_local_z(1)
            << "\n"
            << " local coords of exit point from geom (a check) "
            << local_out.X() << " " << local_out.Y() << " " << local_out.Z()
            << std::endl;
      }
      /*
      if (Verbosity() > 4)
      {
        // As a check, get the positions of the hit pixels in world coordinates from the geo object
              auto hskey = MvtxDefs::genHitSetKey(*layer,stave_number,chip_number,strobe);
              auto surf = tgeometry->maps().getSiliconSurface(hskey);

        TVector3 location_in = layergeom->get_world_from_local_coords(surf,tgeometry, local_in);
        TVector3 location_out = layergeom->get_world_from_local_coords(surf,tgeometry, local_out);

        std::cout
          << std::endl
          << "      PHG4MvtxHitReco:  Found world entry location from geometry for  "
          << " stave number " << stave_number
          << " half stave number " << half_stave_number
          << " module number " << module_number
          << " chip number " << chip_number
          << std::endl
          << " x = " << location_in.X()
          << " y = " << location_in.Y()
          << " z = " << location_in.Z()
          << " radius " << sqrt(pow(location_in.X(), 2) + pow(location_in.Y(), 2))
          << " angle " << atan(location_in.Y() / location_in.X())
          << std::endl;
          << "     PHG4MvtxHitReco: The world entry location from G4 was "
          << " x = " << g4hit->get_x(0)
          << " y = " << g4hit->get_y(0)
          << " z = " << g4hit->get_z(0)
          << " radius " << sqrt(pow(g4hit->get_x(0), 2) + pow(g4hit->get_y(0), 2))
          << " angle " << atan(g4hit->get_y(0) / g4hit->get_x(0))
          << std::endl;
          << " difference in x = " << g4hit->get_x(0) - location_in.X()
          << " difference in y = " << g4hit->get_y(0) - location_in.Y()
          << " difference in z = " << g4hit->get_z(0) - location_in.Z()
          << " difference in radius = "
          << sqrt(pow(g4hit->get_x(0), 2) + pow(g4hit->get_y(0), 2)) - sqrt(pow(location_in.X(), 2) + pow(location_in.Y(), 2))
          << " in angle = " << atan(g4hit->get_y(0) / g4hit->get_x(0)) - atan(location_in.Y() / location_in.X())
          << std::endl;
      }
  */
      // Get the pixel number of the entry location
      int pixel_number_in = layergeom->get_pixel_from_local_coords(local_in);
      // Get the pixel number of the exit location
      int pixel_number_out = layergeom->get_pixel_from_local_coords(local_out);

      if (pixel_number_in < 0 || pixel_number_out < 0)
      {
        std::cout
            << "Oops!  got negative pixel number in layer " << layergeom->get_layer()
            << " pixel_number_in " << pixel_number_in
            << " pixel_number_out " << pixel_number_out
            << " local_in = " << local_in.X() << " " << local_in.Y() << " " << local_in.Z()
            << " local_out = " << local_out.X() << " " << local_out.Y() << " " << local_out.Z()
            << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      if (Verbosity() > 3)
      {
        std::cout
            << "entry pixel number " << pixel_number_in
            << " exit pixel number " << pixel_number_out
            << std::endl;
      }

      std::vector<int> vpixel;
      std::vector<int> vxbin;
      std::vector<int> vzbin;
      std::vector<std::pair<double, double> > venergy;
      // double trklen = 0.0;

      //===================================================
      // OK, now we have found which sensor the hit is in, extracted the hit
      // position in local sensor coordinates,  and found the pixel numbers of the
      // entry point and exit point

      //====================================================
      // Beginning of charge sharing implementation
      //    Find tracklet line inside sensor
      //    Divide tracklet line into n segments (vary n until answer stabilizes)
      //    Find centroid of each segment
      //    Diffuse charge at each centroid
      //    Apportion charge between neighboring pixels
      //    Add the pixel energy contributions from different track segments together
      //====================================================

      TVector3 pathvec = local_in - local_out;

      // See figure 7.3 of the thesis by  Lucasz Maczewski (arXiv:10053.3710) for diffusion simulations in a MAPS epitaxial layer
      // The diffusion widths below were inspired by those plots, corresponding to where the probability drops off to 1/3 of the peak value
      // However note that we make the simplifying assumption that the probability distribution is flat within this diffusion width,
      // while in the simulation it is not
      // double diffusion_width_max = 35.0e-04;   // maximum diffusion radius 35 microns, in cm
      // double diffusion_width_min = 12.0e-04;   // minimum diffusion radius 12 microns, in cm
      double diffusion_width_max = 25.0e-04;  // maximum diffusion radius 35 microns, in cm
      double diffusion_width_min = 8.0e-04;   // minimum diffusion radius 12 microns, in cm

      double ydrift_max = pathvec.Y();
      int nsegments = 4;

      // we want to make a list of all pixels possibly affected by this hit
      // we take the entry and exit locations in local coordinates, and build
      // a rectangular array of pixels that encompasses both, with "nadd" pixels added all around

      int xbin_in = layergeom->get_pixel_X_from_pixel_number(pixel_number_in);
      int zbin_in = layergeom->get_pixel_Z_from_pixel_number(pixel_number_in);
      int xbin_out = layergeom->get_pixel_X_from_pixel_number(pixel_number_out);
      int zbin_out = layergeom->get_pixel_Z_from_pixel_number(pixel_number_out);

      int xbin_max, xbin_min;
      int nadd = 2;
      if (xbin_in > xbin_out)
      {
        xbin_max = xbin_in + nadd;
        xbin_min = xbin_out - nadd;
      }
      else
      {
        xbin_max = xbin_out + nadd;
        xbin_min = xbin_in - nadd;
      }

      int zbin_max, zbin_min;
      if (zbin_in > zbin_out)
      {
        zbin_max = zbin_in + nadd;
        zbin_min = zbin_out - nadd;
      }
      else
      {
        zbin_max = zbin_out + nadd;
        zbin_min = zbin_in - nadd;
      }

      // need to check that values of xbin and zbin are within the valid range
      // YCM: Fix pixel range: Xbin (row) 0 to 511, Zbin (col) 0 to 1023
      if (xbin_min < 0)
      {
        xbin_min = 0;
      }
      if (zbin_min < 0)
      {
        zbin_min = 0;
      }
      if (xbin_max >= maxNX)
      {
        xbin_max = maxNX - 1;
      }
      if (zbin_max >= maxNZ)
      {
        xbin_max = maxNZ - 1;
      }

      if (Verbosity() > 1)
      {
        std::cout << " xbin_in " << xbin_in << " xbin_out " << xbin_out << " xbin_min " << xbin_min << " xbin_max " << xbin_max << std::endl;
        std::cout << " zbin_in " << zbin_in << " zbin_out " << zbin_out << " zbin_min " << zbin_min << " zbin_max " << zbin_max << std::endl;
      }

      // skip this hit if it involves an unreasonable  number of pixels
      // this skips it if either the xbin or ybin range traversed is greater than 8 (for 8 adding two pixels at each end makes the range 12)
      if (xbin_max - xbin_min > 12 || zbin_max - zbin_min > 12)
      {
        continue;
      }

      // this hit is skipped earlier if this dimensioning would be exceeded
      double pixenergy[12][12] = {};  // init to 0
      double pixeion[12][12] = {};    // init to 0

      // Loop over track segments and diffuse charge at each segment location, collect energy in pixels
      for (int i = 0; i < nsegments; i++)
      {
        // Find the tracklet segment location
        // If there are n segments of equal length, we want 2*n intervals
        // The 1st segment is centered at interval 1, the 2nd at interval 3, the nth at interval 2n -1
        double interval = 2 * (double) i + 1;
        double frac = interval / (double) (2 * nsegments);
        TVector3 segvec(pathvec.X() * frac, pathvec.Y() * frac, pathvec.Z() * frac);
        segvec = segvec + local_out;

        //  Find the distance to the back of the sensor from the segment location
        // That projection changes only the value of y
        double ydrift = segvec.Y() - local_out.Y();

        // Caculate the charge diffusion over this drift distance
        // increases from diffusion width_min to diffusion_width_max
        double ydiffusion_radius = diffusion_width_min + (ydrift / ydrift_max) * (diffusion_width_max - diffusion_width_min);

        if (Verbosity() > 5)
        {
          std::cout
              << " segment " << i
              << " interval " << interval
              << " frac " << frac
              << " local_in.X " << local_in.X()
              << " local_in.Z " << local_in.Z()
              << " local_in.Y " << local_in.Y()
              << " pathvec.X " << pathvec.X()
              << " pathvec.Z " << pathvec.Z()
              << " pathvec.Y " << pathvec.Y()
              << " segvec.X " << segvec.X()
              << " segvec.Z " << segvec.Z()
              << " segvec.Y " << segvec.Y()
              << " ydrift " << ydrift
              << " ydrift_max " << ydrift_max
              << " ydiffusion_radius " << ydiffusion_radius
              << std::endl;
        }
        // Now find the area of overlap of the diffusion circle with each pixel and apportion the energy
        for (int ix = xbin_min; ix <= xbin_max; ix++)
        {
          for (int iz = zbin_min; iz <= zbin_max; iz++)
          {
            // Find the pixel corners for this pixel number
            int pixnum = layergeom->get_pixel_number_from_xbin_zbin(ix, iz);

            if (pixnum < 0)
            {
              std::cout
                  << " pixnum < 0 , pixnum = " << pixnum << "\n"
                  << " ix " << ix << " iz " << iz << "\n"
                  << " xbin_min " << xbin_min << " zbin_min " << zbin_min << "\n"
                  << " xbin_max " << xbin_max << " zbin_max " << zbin_max << "\n"
                  << " maxNX " << maxNX << " maxNZ " << maxNZ
                  << std::endl;
            }

            TVector3 tmp = layergeom->get_local_coords_from_pixel(pixnum);
            // note that (x1,z1) is the top left corner, (x2,z2) is the bottom right corner of the pixel - circle_rectangle_intersection expects this ordering
            double x1 = tmp.X() - xpixw_half;
            double z1 = tmp.Z() + zpixw_half;
            double x2 = tmp.X() + xpixw_half;
            double z2 = tmp.Z() - zpixw_half;

            // here segvec.X and segvec.Z are the center of the circle, and diffusion_radius is the circle radius
            // circle_rectangle_intersection returns the overlap area of the circle and the pixel. It is very fast if there is no overlap.
            double pixarea_frac = PHG4Utils::circle_rectangle_intersection(x1, z1, x2, z2, segvec.X(), segvec.Z(), ydiffusion_radius) / (M_PI * pow(ydiffusion_radius, 2));
            // assume that the energy is deposited uniformly along the tracklet length, so that this segment gets the fraction 1/nsegments of the energy
            pixenergy[ix - xbin_min][iz - zbin_min] += pixarea_frac * g4hit->get_edep() / (float) nsegments;
            if (g4hit->has_property(PHG4Hit::prop_eion))
            {
              pixeion[ix - xbin_min][iz - zbin_min] += pixarea_frac * g4hit->get_eion() / (float) nsegments;
            }
            if (Verbosity() > 5)
            {
              std::cout
                  << "    pixnum " << pixnum << " xbin " << ix << " zbin " << iz
                  << " pixel_area fraction of circle " << pixarea_frac << " accumulated pixel energy " << pixenergy[ix - xbin_min][iz - zbin_min]
                  << std::endl;
            }
          }
        }
      }  // end loop over segments

      // now we have the energy deposited in each pixel, summed over all tracklet segments. We make a vector of all pixels with non-zero energy deposited
      for (int ix = xbin_min; ix <= xbin_max; ix++)
      {
        for (int iz = zbin_min; iz <= zbin_max; iz++)
        {
          if (pixenergy[ix - xbin_min][iz - zbin_min] > 0.0)
          {
            int pixnum = layergeom->get_pixel_number_from_xbin_zbin(ix, iz);
            vpixel.push_back(pixnum);
            vxbin.push_back(ix);
            vzbin.push_back(iz);
            std::pair<double, double> tmppair = std::make_pair(pixenergy[ix - xbin_min][iz - zbin_min], pixeion[ix - xbin_min][iz - zbin_min]);
            venergy.push_back(tmppair);
            if (Verbosity() > 1)
            {
              std::cout
                  << " Added pixel number " << pixnum << " xbin " << ix
                  << " zbin " << iz << " to vectors with energy " << pixenergy[ix - xbin_min][iz - zbin_min]
                  << std::endl;
            }
          }
        }
      }

      //===================================
      // End of charge sharing implementation
      //===================================

      // loop over all fired cells for this g4hit and add them to the TrkrHitSet

      for (unsigned int i1 = 0; i1 < vpixel.size(); i1++)  // loop over all fired cells
      {
        // This is the new storage object version
        //====================================
        for (unsigned int i_rep = 0; i_rep < n_replica; i_rep++)
        {
          int strobe = t0_strobe_frame + i_rep;
          // to fit in a 5 bit field in the hitsetkey [-16,15]
          if (strobe < -16)
          {
            strobe = -16;
          }
          if (strobe >= 16)
          {
            strobe = 15;
          }

          // We need to create the TrkrHitSet if not already made - each TrkrHitSet should correspond to a chip for the Mvtx
          TrkrDefs::hitsetkey hitsetkey = MvtxDefs::genHitSetKey(layer, stave_number, chip_number, strobe);
          // Use existing hitset or add new one if needed
          TrkrHitSetContainer::Iterator hitsetit = trkrHitSetContainer->findOrAddHitSet(hitsetkey);

          // generate the key for this hit
          TrkrDefs::hitkey hitkey = MvtxDefs::genHitKey(vzbin[i1], vxbin[i1]);
          // See if this hit already exists
          TrkrHit* hit = nullptr;
          hit = hitsetit->second->getHit(hitkey);

          if (hit)
          {
            std::cout << PHWHERE << "::" << __func__
                      << " - duplicated hit, hitsetkey: " << hitsetkey
                      << " hitkey: " << hitkey << std::endl;
            continue;
          }
          const TrkrDefs::hitsetkey hitsetkeymask = MvtxDefs::genHitSetKey(layer, stave_number, chip_number, 0);

          if (std::find(m_deadPixelMap.begin(), m_deadPixelMap.end(), std::make_pair(hitsetkeymask, hitkey)) == m_deadPixelMap.end())
          {
            // create hit and insert in hitset

            if (Verbosity() > 0)
            {
              std::cout << "Layer: " << layer << ", Stave: " << (uint16_t) MvtxDefs::getStaveId(hitsetkey) << ", Chip: " << (uint16_t) MvtxDefs::getChipId(hitsetkey) << ", Row: " << (uint16_t) MvtxDefs::getRow(hitkey) << ", Col: " << (uint16_t) MvtxDefs::getCol(hitkey) << ", Strobe: " << (int) MvtxDefs::getStrobeId(hitsetkey) << ", added hit " << hitkey << " to hitset " << hitsetkey << ", hit is not associated with a dead pixel, it proceeds to create a new hit." << std::endl;
            }

            hit = new TrkrHitv2();
            hitsetit->second->addHitSpecificKey(hitkey, hit);
          }
          else
          {
            if (Verbosity() > 0)
            {
              std::cout << "Layer: " << layer << ", Stave: " << (uint16_t) MvtxDefs::getStaveId(hitsetkey) << ", Chip: " << (uint16_t) MvtxDefs::getChipId(hitsetkey) << ", Row: " << (uint16_t) MvtxDefs::getRow(hitkey) << ", Col: " << (uint16_t) MvtxDefs::getCol(hitkey) << ", Strobe: " << (int) MvtxDefs::getStrobeId(hitsetkey) << ", added hit " << hitkey << " to hitset " << hitsetkey << ", hit is not associated with a dead pixel, it proceeds to continue." << std::endl;
            }
            continue;
          }
          // Either way, add the energy to it
          double hitenergy = venergy[i1].first * TrkrDefs::MvtxEnergyScaleup;
          hit->addEnergy(hitenergy);

          addtruthhitset(hitsetkey, hitkey, hitenergy);

          if (Verbosity() > 0)
          {
            std::cout << "Layer: " << layer << ", Stave: " << (uint16_t) MvtxDefs::getStaveId(hitsetkey) << ", Chip: " << (uint16_t) MvtxDefs::getChipId(hitsetkey) << ", Row: " << (uint16_t) MvtxDefs::getRow(hitkey) << ", Col: " << (uint16_t) MvtxDefs::getCol(hitkey) << ", Strobe: " << (int) MvtxDefs::getStrobeId(hitsetkey) << ", added hit " << hitkey << " to hitset " << hitsetkey << " with energy " << hit->getEnergy() / TrkrDefs::MvtxEnergyScaleup << std::endl;
          }

          // now we update the TrkrHitTruthAssoc map - the map contains <hitsetkey, std::pair <hitkey, g4hitkey> >
          // There is only one TrkrHit per pixel, but there may be multiple g4hits
          // How do we know how much energy from PHG4Hit went into TrkrHit? We don't, have to sort it out in evaluator to save memory

          // we set the strobe ID to zero in the hitsetkey
          // we use the findOrAdd method to keep from adding identical entries
          TrkrDefs::hitsetkey bare_hitsetkey = zero_strobe_bits(hitsetkey);
          hitTruthAssoc->findOrAddAssoc(bare_hitsetkey, hitkey, g4hit_it->first);
        }
      }  // end loop over hit cells
    }    // end loop over g4hits for this layer

  }  // end loop over layers

  // print the list of entries in the association table
  if (Verbosity() > 0)
  {
    std::cout << "From PHG4MvtxHitReco: " << std::endl;
    hitTruthAssoc->identify();
  }

  // spit out the truth clusters

  end_event_truthcluster(topNode);

  if (Verbosity() > 3)
  {
    int nclusprint = -1;
    std::cout << PHWHERE << ": content of clusters " << std::endl;
    auto& tmap = m_truthtracks->getMap();
    std::cout << " Number of tracks: " << tmap.size() << std::endl;
    for (auto& _pair : tmap)
    {
      auto& track = _pair.second;

      printf("id(%2i) phi:eta:pt(", (int) track->getTrackid());
      std::cout << "phi:eta:pt(";
      printf("%5.2f:%5.2f:%5.2f", track->getPhi(), track->getPseudoRapidity(), track->getPt());
      /* Form("%5.2:%5.2:%5.2", track->getPhi(), track->getPseudoRapidity(), track->getPt()) */
      //<<track->getPhi()<<":"<<track->getPseudoRapidity()<<":"<<track->getPt()
      std::cout << ") nclusters(" << track->getClusters().size() << ") ";
      int nclus = 0;
      for (auto cluskey : track->getClusters())
      {
        std::cout << " "
                  << ((int) TrkrDefs::getHitSetKeyFromClusKey(cluskey)) << ":index(" << ((int) TrkrDefs::getClusIndex(cluskey)) << ")" << std::endl;
        ++nclus;
        if (nclusprint > 0 && nclus >= nclusprint)
        {
          std::cout << " ... ";
          break;
        }
      }
    }
    std::cout << PHWHERE << " ----- end of clusters " << std::endl;
  }

  if (m_is_emb)
  {
    cluster_truthhits(topNode);  // the last track was truth -- cluster it
    prior_g4hit = nullptr;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::pair<double, double> PHG4MvtxHitReco::generate_alpide_pulse(const double energy_deposited)
{
  // We need to translate energy deposited to num/ electrons released
  if (Verbosity() > 2)
  {
    std::cout << "energy_deposited: " << energy_deposited << std::endl;
  }
  // int silicon_band_gap = 1.12; //Band gap energy in eV
  // int Q_in = rand() % 5000 + 50;
  // int clipping_point = 110;
  // double ToT_start = Q_in < 200 ? 395.85*exp(-0.5*pow((Q_in+851.43)/286.91, 2)) : 0.5;
  // double ToT_end = Q_in < clipping_point ? 5.90*exp(-0.5*pow((Q_in-99.86)/54.80, 2)) : 5.8 - 6.4e-4 * Q_in;

  // return make_pair(ToT_start*1e3, ToT_end*1e3);
  //  Using constant alpide pulse length
  return std::make_pair<double, double>(1.5 * microsecond, 5.9 * microsecond);
}

double PHG4MvtxHitReco::generate_strobe_zero_tm_start()
{
  return -1. * gsl_rng_uniform_pos(m_rng.get()) * (m_strobe_separation + m_strobe_width);
}

int PHG4MvtxHitReco::get_strobe_frame(double alpide_time, double strobe_zero_tm_start)
{
  int strobe_frame = int((alpide_time - strobe_zero_tm_start) / (m_strobe_width + m_strobe_separation));
  strobe_frame += (alpide_time < strobe_zero_tm_start) ? -1 : 0;
  return strobe_frame;
}

void PHG4MvtxHitReco::set_timing_window(const int detid, const double tmin, const double tmax)
{
  if (false)
  {
    std::cout
        << "PHG4MvtxHitReco: Set Mvtx timing window parameters from macro for layer = "
        << detid << " to tmin = " << tmin << " tmax = " << tmax
        << std::endl;
  }

  return;
}

void PHG4MvtxHitReco::SetDefaultParameters()
{
  // cout << "PHG4MvtxHitReco: Setting Mvtx timing window defaults to tmin = -5000 and  tmax = 5000 ns" << std::endl;
  set_default_double_param("mvtx_tmin", -5000);
  set_default_double_param("mvtx_tmax", 5000);
  set_default_double_param("mvtx_strobe_width", 5 * microsecond);
  set_default_double_param("mvtx_strobe_separation", 0.);
  set_default_int_param("mvtx_in_sphenix_srdo", (int) false);
  return;
}

TrkrDefs::hitsetkey PHG4MvtxHitReco::zero_strobe_bits(TrkrDefs::hitsetkey hitsetkey)
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  unsigned int stave = MvtxDefs::getStaveId(hitsetkey);
  unsigned int chip = MvtxDefs::getChipId(hitsetkey);
  TrkrDefs::hitsetkey bare_hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip, 0);

  return bare_hitsetkey;
}

PHG4MvtxHitReco::~PHG4MvtxHitReco()
{
  delete m_truth_hits;
}

void PHG4MvtxHitReco::truthcheck_g4hit(PHG4Hit* g4hit, PHCompositeNode* topNode)
{
  int new_trkid = (g4hit == nullptr) ? -1 : g4hit->get_trkid();
  bool is_new_track = (new_trkid != m_trkid);
  if (Verbosity() > 5)
  {
    std::cout << PHWHERE << std::endl
              << " -> Checking status of PHG4Hit. Track id(" << new_trkid << ")" << std::endl;
  }
  if (!is_new_track)
  {
    /* FIXME */
    /* std::cout << " FIXME PEAR Checking g4hit : " << g4hit->get_x(0) << " " */
    /*   << g4hit->get_y(0) << " " << g4hit->get_z(0) */
    /*   << "  trackid("<<g4hit->get_trkid() << ") " << std::endl; */
    if (m_is_emb)
    {
      /* std::cout << " FIXME Checking MVTX g4hit : " << g4hit->get_x(0) << " " */
      /*   << g4hit->get_y(0) << " " << g4hit->get_z(0) */
      /*   << "  trackid("<<g4hit->get_trkid() << ") " << std::endl; */
      if (prior_g4hit != nullptr && (std::abs(prior_g4hit->get_x(0) - g4hit->get_x(0)) > max_g4hitstep || std::abs(prior_g4hit->get_y(0) - g4hit->get_y(0)) > max_g4hitstep))
      {
        // this is a looper track -- cluster hits up to this point already
        cluster_truthhits(topNode);
      }
      prior_g4hit = g4hit;
    }
    return;
  }
  // <- STATUS: this is a new track
  if (Verbosity() > 2)
  {
    std::cout << PHWHERE << std::endl
              << " -> Found new embedded track with id: " << new_trkid << std::endl;
  }
  if (m_is_emb)
  {
    cluster_truthhits(topNode);  // cluster m_truth_hits and add m_current_track
    m_current_track = nullptr;
    prior_g4hit = nullptr;
  }
  m_trkid = new_trkid;
  m_is_emb = m_truthinfo->isEmbeded(m_trkid);
  if (m_is_emb)
  {
    m_current_track = m_truthtracks->getTruthTrack(m_trkid, m_truthinfo);
    prior_g4hit = g4hit;
  }
}

void PHG4MvtxHitReco::end_event_truthcluster(PHCompositeNode* topNode)
{
  if (m_is_emb)
  {
    cluster_truthhits(topNode);  // cluster m_truth_hits and add m_current_track
    m_current_track = nullptr;
    m_trkid = -1;
    m_is_emb = false;
  }
  m_hitsetkey_cnt.clear();
}

void PHG4MvtxHitReco::addtruthhitset(
    TrkrDefs::hitsetkey hitsetkey,
    TrkrDefs::hitkey hitkey,
    float neffelectrons)
{
  if (!m_is_emb)
  {
    return;
  }
  TrkrHitSetContainer::Iterator hitsetit = m_truth_hits->findOrAddHitSet(hitsetkey);
  // See if this hit already exists
  TrkrHit* hit = nullptr;
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

void PHG4MvtxHitReco::cluster_truthhits(PHCompositeNode* topNode)
{
  // clusterize the mvtx hits in m_truth_hits, put them in m_truthclusters, and put the id's in m_current_track
  // ----------------------------------------------------------------------------------------------------
  // Digitization -- simplified from g4mvtx/PHG4MvtxDigitizer --
  // ----------------------------------------------------------------------------------------------------

  /* // We want all hitsets for the Mvtx */
  /* TrkrHitSetContainer::ConstRange hitset_range = m_truth_hits->getHitSets(TrkrDefs::TrkrId::mvtxId); */
  /* for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first; */
  /*      hitset_iter != hitset_range.second; */
  /*      ++hitset_iter) */
  /* { */
  /*   // we have an iterator to one TrkrHitSet for the mvtx from the trkrHitSetContainer */
  /*   // get the hitset key so we can find the layer */
  /*   TrkrDefs::hitsetkey hitsetkey = hitset_iter->first; */
  /*   int layer = TrkrDefs::getLayer(hitsetkey); */
  /*   // FIXME */
  /*   /1* if (Verbosity() > 1) std::cout << "PHG4MvtxDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << std::endl; *1/ */
  /*   std::cout << "PHG4MvtxDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << std::endl; */

  /*   // get all of the hits from this hitset */
  /*   TrkrHitSet *hitset = hitset_iter->second; */
  /*   TrkrHitSet::ConstRange hit_range = hitset->getHits(); */
  /*   std::set<TrkrDefs::hitkey> hits_rm; */
  /*   for (TrkrHitSet::ConstIterator hit_iter = hit_range.first; */
  /*        hit_iter != hit_range.second; */
  /*        ++hit_iter) */
  /*   { */
  /*     TrkrHit *hit = hit_iter->second; */

  /*     // Convert the signal value to an ADC value and write that to the hit */
  /*     // Unsigned int adc = hit->getEnergy() / (TrkrDefs::MvtxEnergyScaleup *_energy_scale[layer]); */
  /*     if (Verbosity() > 0) */
  /*       std::cout << "    PHG4MvtxDigitizer: found hit with key: " << hit_iter->first << " and signal " << hit->getEnergy() / TrkrDefs::MvtxEnergyScaleup << " in layer " << layer << std::endl; */
  /*     // Remove the hits with energy under threshold */
  /*     bool rm_hit = false; */
  /*     // FIXME: */
  /*     double _energy_threshold = 0.; // FIXME */
  /*     const double dummy_pixel_thickness = 0.001; */
  /*     const double mip_e = 0.003876; */
  /*     double _energy_scale =  mip_e * dummy_pixel_thickness; // FIXME: note: this doesn't actually */
  /*                                                            // matter either here or for the Svtx Tracks -- the energy is digital -- either the hit is there or it isn't */
  /*     double _max_adc = 255; */
  /*     if ((hit->getEnergy() / TrkrDefs::MvtxEnergyScaleup) < _energy_threshold) */
  /*     { */
  /*       if (Verbosity() > 0) std::cout << "         remove hit, below energy threshold of " << _energy_threshold << std::endl; */
  /*       rm_hit = true; */
  /*     } */
  /*     unsigned short adc = (unsigned short) (hit->getEnergy() / (TrkrDefs::MvtxEnergyScaleup * _energy_scale)); */
  /*     if (adc > _max_adc) adc = _max_adc; */
  /*     hit->setAdc(adc); */

  /* if (rm_hit) hits_rm.insert(hit_iter->first); */
  /*   } */

  /*   for (const auto &key : hits_rm) */
  /*   { */
  /*     if (Verbosity() > 0) std::cout << "    PHG4MvtxDigitizer: remove hit with key: " << key << std::endl; */
  /*     hitset->removeHit(key); */
  /*   } */
  /* } */

  // ----------------------------------------------------------------------------------------------------
  // Prune hits -- simplified from mvtx/MvtxHitReco
  // ----------------------------------------------------------------------------------------------------
  std::multimap<TrkrDefs::hitsetkey, TrkrDefs::hitsetkey> hitset_multimap;  // will map (bare hitset, hitset with strobe)
  std::set<TrkrDefs::hitsetkey> bare_hitset_set;                            // list of all physical sensor hitsetkeys (i.e. with strobe set to zero)

  TrkrHitSetContainer::ConstRange hitsetrange = m_truth_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);  // actually m_truth_hits can only have MVTX hits at this point...
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    auto hitsetkey = hitsetitr->first;

    // get the hitsetkey value for strobe 0
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    unsigned int stave = MvtxDefs::getStaveId(hitsetitr->first);
    unsigned int chip = MvtxDefs::getChipId(hitsetitr->first);
    auto bare_hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip, 0);

    hitset_multimap.insert(std::make_pair(bare_hitsetkey, hitsetkey));
    bare_hitset_set.insert(bare_hitsetkey);

    if (Verbosity() > 0)
    {
      std::cout << " found hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << std::endl;
    }
  }

  // Now consolidate all hits into the hitset with strobe 0, and delete the other hitsets
  //==============================================================
  for (unsigned int bare_hitsetkey : bare_hitset_set)
  {
    TrkrHitSet* bare_hitset = (m_truth_hits->findOrAddHitSet(bare_hitsetkey))->second;
    if (Verbosity() > 0)
    {
      std::cout << "         bare_hitset " << bare_hitsetkey << " initially has " << bare_hitset->size() << " hits " << std::endl;
    }

    auto bare_hitsetrange = hitset_multimap.equal_range(bare_hitsetkey);
    for (auto it = bare_hitsetrange.first; it != bare_hitsetrange.second; ++it)
    {
      auto hitsetkey = it->second;

      int strobe = MvtxDefs::getStrobeId(hitsetkey);
      if (strobe != 0)
      {
        if (Verbosity() > 0)
        {
          std::cout << "            process hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << std::endl;
        }

        // copy all hits to the hitset with strobe 0
        TrkrHitSet* hitset = m_truth_hits->findHitSet(hitsetkey);

        if (Verbosity() > 0)
        {
          std::cout << "                hitsetkey " << hitsetkey << " has strobe " << strobe << " and has " << hitset->size() << " hits,  so copy it" << std::endl;
        }

        TrkrHitSet::ConstRange hitrangei = hitset->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          auto hitkey = hitr->first;
          if (Verbosity() > 0)
          {
            std::cout << "                 found hitkey " << hitkey << std::endl;
          }
          // if it is already there, leave it alone, this is a duplicate hit
          auto tmp_hit = bare_hitset->getHit(hitkey);
          if (tmp_hit)
          {
            if (Verbosity() > 0)
            {
              std::cout << "                          hitkey " << hitkey << " is already in bare hitsest, do not copy" << std::endl;
            }
            continue;
          }

          // otherwise copy the hit over
          if (Verbosity() > 0)
          {
            std::cout << "                          copying over hitkey " << hitkey << std::endl;
          }
          auto old_hit = hitr->second;
          TrkrHit* new_hit = new TrkrHitv2();
          new_hit->setAdc(old_hit->getAdc());
          bare_hitset->addHitSpecificKey(hitkey, new_hit);
        }

        // all hits are copied over to the strobe zero hitset, remove this hitset
        m_truth_hits->removeHitSet(hitsetkey);
      }
    }
  }

  // ----------------------------------------------------------------------------------------------------
  // cluster hits -- simplified from mvtx/MvtxClusterizer
  // ----------------------------------------------------------------------------------------------------
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "WARNING: cannot find the geometry CYLINDERGEOM_MVTX" << std::endl;
    m_truth_hits->Reset();
    return;
  }

  //-----------
  // Clustering
  //-----------

  // loop over each MvtxHitSet object (chip)
  hitsetrange = m_truth_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {                                          // hitsetitr    : pair(TrkrDefs::hitsetkey, TrkrHitSet>;   TrkrHitSet : map <HitKey, TrkrHit>
    TrkrHitSet* hitset = hitsetitr->second;  // hitset : map <TrkrDefs::hitkey, TrkrHit>

    /* if (true) */
    if (Verbosity() > 0)
    {
      unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
      unsigned int stave = MvtxDefs::getStaveId(hitsetitr->first);
      unsigned int chip = MvtxDefs::getChipId(hitsetitr->first);
      unsigned int strobe = MvtxDefs::getStrobeId(hitsetitr->first);
      std::cout << "MvtxClusterizer found hitsetkey " << hitsetitr->first << " layer " << layer << " stave " << stave << " chip " << chip << " strobe " << strobe << std::endl;
    }

    if (Verbosity() > 2)
    {
      hitset->identify();
    }

    TrkrHitSet::ConstRange hitrangei = hitset->getHits();

    auto hitsetkey = hitset->getHitSetKey();
    auto ckey = TrkrDefs::genClusKey(hitsetkey, 0);  // there is only one cluster made per cluskey

    // determine the size of the cluster in phi and z
    std::set<int> phibins;
    std::set<int> zbins;

    // determine the cluster position...
    double locxsum = 0.;
    double loczsum = 0.;

    double locclusx = NAN;
    double locclusz = NAN;

    // we need the geometry object for this layer to get the global positions
    int layer = TrkrDefs::getLayer(ckey);
    auto layergeom = dynamic_cast<CylinderGeom_Mvtx*>(geom_container->GetLayerGeom(layer));
    if (!layergeom)
    {
      exit(1);
    }

    // make a tunable threshold for energy in a given hit
    //  -- percentage of someting? (total cluster energy)
    double sum_adc{0.};
    for (auto ihit = hitrangei.first; ihit != hitrangei.second; ++ihit)
    {
      sum_adc += ihit->second->getAdc();
    }
    const double threshold = sum_adc * m_pixel_thresholdrat;       // FIXME -- tune this as needed
    std::map<int, unsigned int> m_iphi, m_it, m_iphiCut, m_itCut;  // FIXME
                                                                   /* const unsigned int npixels = std::distance( hitrangei.first, hitrangei.second ); */
    // to tune this parameter: run a bunch of tracks and compare truth sizes and reco sizes,
    // should come out the same
    double npixels{0.};
    for (auto ihit = hitrangei.first; ihit != hitrangei.second; ++ihit)
    {
      const auto adc = ihit->second->getAdc();
      int col = MvtxDefs::getCol(ihit->first);
      int row = MvtxDefs::getRow(ihit->first);

      if (mClusHitsVerbose)
      {
        std::map<int, unsigned int>& m_phi = (adc < threshold) ? m_iphiCut : m_iphi;
        std::map<int, unsigned int>& m_z = (adc < threshold) ? m_itCut : m_it;

        auto pnew = m_phi.try_emplace(row, adc);
        if (!pnew.second)
        {
          pnew.first->second += adc;
        }

        pnew = m_z.try_emplace(col, adc);
        if (!pnew.second)
        {
          pnew.first->second += adc;
        }
      }
      if (adc < threshold)
      {
        continue;
      }

      // size
      npixels += 1.;
      zbins.insert(col);
      phibins.insert(row);

      // get local coordinates, in stave reference frame, for hit
      auto local_coords = layergeom->get_local_coords_from_pixel(row, col);

      /*
         manually offset position along y (thickness of the sensor),
         to account for effective hit position in the sensor, resulting from diffusion.
         Effective position corresponds to 1um above the middle of the sensor
         */
      local_coords.SetY(1e-4);

      // update cluster position
      locxsum += local_coords.X();
      loczsum += local_coords.Z();
      // add the association between this cluster key and this hitkey to the table
    }  // mapiter
    if (mClusHitsVerbose)
    {
      if (Verbosity() > 10)
      {
        for (auto& hit : m_iphi)
        {
          std::cout << " m_phi(" << hit.first << " : " << hit.second << ") " << std::endl;
        }
      }
      for (auto& hit : m_iphi)
      {
        mClusHitsVerbose->addPhiHit(hit.first, (float) hit.second);
      }
      for (auto& hit : m_it)
      {
        mClusHitsVerbose->addZHit(hit.first, (float) hit.second);
      }
      for (auto& hit : m_iphiCut)
      {
        mClusHitsVerbose->addPhiCutHit(hit.first, (float) hit.second);
      }
      for (auto& hit : m_itCut)
      {
        mClusHitsVerbose->addZCutHit(hit.first, (float) hit.second);
      }
    }

    // This is the local position
    locclusx = locxsum / npixels;
    locclusz = loczsum / npixels;

    const double pitch = layergeom->get_pixel_x();
    const double length = layergeom->get_pixel_z();
    const double phisize = phibins.size() * pitch;
    const double zsize = zbins.size() * length;

    /* if (true) { */
    if (Verbosity() > 0)
    {
      std::cout << " MvtxClusterizer: cluskey " << ckey << " layer " << layer << " rad " << layergeom->get_radius() << " phibins " << phibins.size() << " pitch " << pitch << " phisize " << phisize
                << " zbins " << zbins.size() << " length " << length << " zsize " << zsize
                << " local x " << locclusx << " local y " << locclusz
                << std::endl;
    }

    // ok force it use use cluster version v4 for now (Valgrind is not happy with application of v5)
    /* if (m_cluster_version==4){ */
    if (m_cluster_version == 4)
    {
      auto clus = std::make_unique<TrkrClusterv4>();
      clus->setAdc(npixels);
      clus->setLocalX(locclusx);
      clus->setLocalY(locclusz);

      clus->setPhiSize(phibins.size());
      clus->setZSize(zbins.size());
      // All silicon surfaces have a 1-1 map to hitsetkey.
      // So set subsurface key to 0
      clus->setSubSurfKey(0);

      if (Verbosity() > 2)
      {
        clus->identify();
      }

      // get the count of how many clusters have allready been added to this hitsetkey (possibly from other embedded tracks tracks)
      m_hitsetkey_cnt.try_emplace(hitsetkey, 0);
      unsigned int& cnt = m_hitsetkey_cnt[hitsetkey];
      ckey = TrkrDefs::genClusKey(hitsetkey, cnt);
      m_truthclusters->addClusterSpecifyKey(ckey, clus.release());
      m_current_track->addCluster(ckey);
      if (mClusHitsVerbose)
      {
        mClusHitsVerbose->push_hits(ckey);
      }
      ++cnt;
    }
    else
    {
      std::cout << PHWHERE << std::endl;
      std::cout << "Error: only cluster version 4 allowed." << std::endl;
    }  // loop over hitsets
  }
  m_truth_hits->Reset();
  prior_g4hit = nullptr;
  return;
}
