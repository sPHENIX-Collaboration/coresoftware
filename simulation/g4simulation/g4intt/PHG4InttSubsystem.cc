#include "PHG4InttSubsystem.h"
#include "PHG4InttDefs.h"
#include "PHG4InttDetector.h"
#include "PHG4InttDisplayAction.h"
#include "PHG4InttSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

#include <boost/format.hpp>

using namespace std;

//_______________________________________________________________________
PHG4InttSubsystem::PHG4InttSubsystem(const std::string &detectorname, const vpair &layerconfig)
  : PHG4DetectorGroupSubsystem(detectorname)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
  , m_LayerConfigVector(layerconfig)
  , m_DetectorType(detectorname)
{
  for (vector<pair<int, int>>::const_iterator piter = layerconfig.begin(); piter != layerconfig.end(); ++piter)
  {
    AddDetId((*piter).second);
  }

  InitializeParameters();
  // put the layer into the name so we get unique names
  // for multiple layers
  Name(detectorname);
  SuperDetector(detectorname);
}

PHG4InttSubsystem::~PHG4InttSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4InttSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4InttSubsystem::Init started" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector (detector adds its volumes to it)
  m_DisplayAction = new PHG4InttDisplayAction(Name());
  // create detector
  pair<vector<pair<int, int>>::const_iterator, vector<pair<int, int>>::const_iterator> layer_begin_end = make_pair(m_LayerConfigVector.begin(), m_LayerConfigVector.end());
  m_Detector = new PHG4InttDetector(this, topNode, GetParamsContainer(), Name(), layer_begin_end);
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->Detector(m_DetectorType);
  m_Detector->OverlapCheck(CheckOverlap());

  int active = 0;
  // initialize with support active flag (if support is active we need the absorber hit node)
  int absorberactive = GetParamsContainer()->GetParameters(PHG4InttDefs::SUPPORTPARAMS)->get_int_param("supportactive");
  int blackhole = 0;
  for (set<int>::const_iterator parcontaineriter = GetDetIds().first; parcontaineriter != GetDetIds().second; ++parcontaineriter)
  {
    if (active || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("active"))
    {
      active = 1;
    }
    if (absorberactive || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("absorberactive"))
    {
      absorberactive = 1;
    }
    if (blackhole || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("blackhole"))
    {
      blackhole = 1;
    }
  }
  if (active)
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    std::string nodename = (SuperDetector() != "NONE") ? boost::str(boost::format("G4HIT_%s") % SuperDetector()) : boost::str(boost::format("G4HIT_%s") % m_DetectorType);

    // create hit list
    PHG4HitContainer *hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
    if (!hitcontainer)
      DetNode->addNode(new PHIODataNode<PHObject>(hitcontainer = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));

    if (absorberactive)
    {
      nodename = (SuperDetector() != "NONE") ? boost::str(boost::format("G4HIT_ABSORBER_%s") % SuperDetector()) : boost::str(boost::format("G4HIT_ABSORBER_%s") % m_DetectorType);

      hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
      if (!hitcontainer)
      {
        DetNode->addNode(new PHIODataNode<PHObject>(hitcontainer = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));
      }
    }

    // create stepping action
    m_SteppingAction = new PHG4InttSteppingAction(m_Detector, GetParamsContainer(), layer_begin_end);
  }
  else
  {
    if (blackhole)
    {
      m_SteppingAction = new PHG4InttSteppingAction(m_Detector, GetParamsContainer(), layer_begin_end);
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4InttSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }

  return 0;
}

//_______________________________________________________________________
PHG4Detector *PHG4InttSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4InttSubsystem::SetDefaultParameters()
{
  // We have only two types of ladders, one with vertical strips (SEGMENTATION_Z) and one with horizontal strips (SEGMENTATION_PHI)
  // There are 4 sensors in each ladder
  //     In ladder type 0 the sensor is special and inner and outer sensors are the same.
  //     In ladder type 1 there are two different sensor types, inner and outer
  // We do not want to hard code the ladder types for the layers

  // We define default ladder types for 8 layers, but these can be changed at the macro level

  int laddertype[8] = {PHG4InttDefs::SEGMENTATION_Z,
                       PHG4InttDefs::SEGMENTATION_Z,
                       PHG4InttDefs::SEGMENTATION_PHI,
                       PHG4InttDefs::SEGMENTATION_PHI,
                       PHG4InttDefs::SEGMENTATION_PHI,
                       PHG4InttDefs::SEGMENTATION_PHI,
                       PHG4InttDefs::SEGMENTATION_PHI,
                       PHG4InttDefs::SEGMENTATION_PHI};                                    // default
  int nladder[8] = {17, 17, 15, 15, 18, 18, 21, 21};                                       // default
  double sensor_radius[8] = {6.876, 7.462, 8.987, 9.545, 10.835, 11.361, 12.676, 13.179};  // radius of center of sensor for layer default

  auto detid = GetDetIds();  // get pair of iterators to begin/end of set<int> of detids
  for (auto detiter = detid.first; detiter != detid.second; ++detiter)
  {
    set_default_int_param(*detiter, "active", 1);

    // To reconfigure the layers, all you have to do is overide the defaults for these four arrays from the tracking macro
    set_default_int_param(*detiter, "laddertype", laddertype[*detiter]);
    set_default_int_param(*detiter, "nladder", nladder[*detiter]);  // ladders per layer
    set_default_double_param(*detiter, "sensor_radius", sensor_radius[*detiter]);
    // These offsets should be kept at zero in the new design
    set_default_double_param(*detiter, "offsetphi", 0.);
    set_default_double_param(*detiter, "offsetrot", 0.);
  }

  // These are the parameters that describe the internal ladder geometry for the two ladder types
  {  // just being lazy, using namespace in this scope for less clutter
    using namespace PHG4InttDefs;
    set_default_int_param(SEGMENTATION_Z, "nstrips_phi_cell", 1);
    set_default_int_param(SEGMENTATION_Z, "nstrips_phi_sensor", 1);
    set_default_int_param(SEGMENTATION_Z, "nstrips_z_sensor_0", 128 * 5);
    set_default_int_param(SEGMENTATION_Z, "nstrips_z_sensor_1", 128 * 5);
    set_default_double_param(SEGMENTATION_Z, "fphx_x", 0.032);
    set_default_double_param(SEGMENTATION_Z, "fphx_y", 0.27);
    set_default_double_param(SEGMENTATION_Z, "fphx_z", 0.91);
    set_default_double_param(SEGMENTATION_Z, "fphx_offset_z", 0.005);
    set_default_double_param(SEGMENTATION_Z, "gap_sensor_fphx", 0.1);
    set_default_double_param(SEGMENTATION_Z, "halfladder_z", 40.00);
    set_default_double_param(SEGMENTATION_Z, "halfladder_inside_z", 23.9622);
    set_default_double_param(SEGMENTATION_Z, "hdi_copper_x", 0.0052);
    set_default_double_param(SEGMENTATION_Z, "hdi_edge_z", 0.);
    set_default_double_param(SEGMENTATION_Z, "hdi_kapton_x", 0.038);
    set_default_double_param(SEGMENTATION_Z, "hdi_y", 2.55);
    set_default_double_param(SEGMENTATION_Z, "pgs_x", 0.02);
    set_default_double_param(SEGMENTATION_Z, "sensor_edge_phi", 0.13);
    set_default_double_param(SEGMENTATION_Z, "sensor_edge_z", 0.1);
    set_default_double_param(SEGMENTATION_Z, "sensor_offset_y", 0.304);
    set_default_double_param(SEGMENTATION_Z, "stave_straight_cooler_y", 0.47);
    set_default_double_param(SEGMENTATION_Z, "stave_straight_inner_y", 0.1);
    set_default_double_param(SEGMENTATION_Z, "stave_straight_outer_y", 0.672);
    set_default_double_param(SEGMENTATION_Z, "strip_x", 0.032);
    set_default_double_param(SEGMENTATION_Z, "strip_y", 1.6);
    set_default_double_param(SEGMENTATION_Z, "strip_z_0", 0.01406);
    set_default_double_param(SEGMENTATION_Z, "strip_z_1", 0.01406);

    set_default_int_param(SEGMENTATION_Z, "carbon_stave_type", 1);
    set_default_double_param(SEGMENTATION_Z, "waterstave_straight_cooler_y", 0.47);
    set_default_double_param(SEGMENTATION_Z, "waterstave_slant_cooler_y", 1.4362);
    set_default_double_param(SEGMENTATION_Z, "waterstave_straight_outer_y", 0.672);
    set_default_double_param(SEGMENTATION_Z, "waterstave_straight_rohacell_y", 0.57181);

    set_default_int_param(SEGMENTATION_PHI, "nstrips_phi_cell", 256);
    set_default_int_param(SEGMENTATION_PHI, "nstrips_phi_sensor", 256);
    set_default_int_param(SEGMENTATION_PHI, "nstrips_z_sensor_0", 8);
    set_default_int_param(SEGMENTATION_PHI, "nstrips_z_sensor_1", 5);
    set_default_double_param(SEGMENTATION_PHI, "fphx_x", 0.032);
    set_default_double_param(SEGMENTATION_PHI, "fphx_y", 0.27);
    set_default_double_param(SEGMENTATION_PHI, "fphx_z", 0.91);
    set_default_double_param(SEGMENTATION_PHI, "fphx_offset_z", 0.005);
    set_default_double_param(SEGMENTATION_PHI, "gap_sensor_fphx", 0.1);
    set_default_double_param(SEGMENTATION_PHI, "halfladder_z", 40.00);
    set_default_double_param(SEGMENTATION_PHI, "halfladder_inside_z", 23.9622);
    set_default_double_param(SEGMENTATION_PHI, "hdi_copper_x", 0.0052);
    set_default_double_param(SEGMENTATION_PHI, "hdi_edge_z", 0.);
    set_default_double_param(SEGMENTATION_PHI, "hdi_kapton_x", 0.038);
    set_default_double_param(SEGMENTATION_PHI, "hdi_y", 3.8);
    set_default_double_param(SEGMENTATION_PHI, "pgs_x", 0.02);
    set_default_double_param(SEGMENTATION_PHI, "sensor_edge_phi", 0.13);
    set_default_double_param(SEGMENTATION_PHI, "sensor_edge_z", 0.1);
    set_default_double_param(SEGMENTATION_PHI, "sensor_offset_y", 0.);
    set_default_double_param(SEGMENTATION_PHI, "stave_straight_cooler_y", 0.47);
    set_default_double_param(SEGMENTATION_PHI, "stave_straight_inner_y", 0.344);
    set_default_double_param(SEGMENTATION_PHI, "stave_straight_outer_y", 0.522);
    set_default_double_param(SEGMENTATION_PHI, "strip_x", 0.032);
    set_default_double_param(SEGMENTATION_PHI, "strip_y", 0.0078);
    set_default_double_param(SEGMENTATION_PHI, "strip_z_0", 1.6);
    set_default_double_param(SEGMENTATION_PHI, "strip_z_1", 2.);

    set_default_int_param(SEGMENTATION_PHI, "carbon_stave_type", 1);
    set_default_double_param(SEGMENTATION_PHI, "waterstave_straight_cooler_y", 1.47684);
    set_default_double_param(SEGMENTATION_PHI, "waterstave_slant_cooler_y", 0.6345120525);
    set_default_double_param(SEGMENTATION_PHI, "waterstave_straight_outer_y", 0.33451);
    set_default_double_param(SEGMENTATION_PHI, "waterstave_straight_rohacell_y", 0.58842);

    set_default_int_param(SUPPORTPARAMS, "supportactive", 0);

    set_default_double_param(SUPPORTPARAMS, "inner_skin_inner_radius", 6.2416);
    set_default_double_param(SUPPORTPARAMS, "inner_skin_length", 50.7);
    set_default_double_param(SUPPORTPARAMS, "inner_skin_outer_radius", 6.2666);

    set_default_double_param(SUPPORTPARAMS, "outer_skin_cfcin_inner_radius", 12.0444);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_cfcin_outer_radius", 12.0694);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_cfcin_length", 50.7);

    set_default_double_param(SUPPORTPARAMS, "outer_skin_foam_inner_radius", 12.0694);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_foam_outer_radius", 12.2194);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_foam_length", 50.7);

    set_default_double_param(SUPPORTPARAMS, "outer_skin_cfcout_inner_radius", 12.2194);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_cfcout_outer_radius", 12.2444);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_cfcout_length", 50.7);

    // Endcap ring flag
    set_default_int_param(SUPPORTPARAMS, "endcap_ring_enabled", 1);
    set_default_int_param(SUPPORTPARAMS, "endcap_ring_type", 1);  // 0: Al+SS+WG, 1 : CarbonPEEK

    // Aluminum endcap ring position
    set_default_double_param(SUPPORTPARAMS, "endcap_ring_z", 24.35);

    // Aluminum endcap ring
    set_default_double_param(SUPPORTPARAMS, "endcap_Alring_inner_radius", 6.267);
    set_default_double_param(SUPPORTPARAMS, "endcap_Alring_outer_radius", 12.0444);
    set_default_double_param(SUPPORTPARAMS, "endcap_Alring_length", 0.3645);

    // Stainless steel endcap ring
    set_default_double_param(SUPPORTPARAMS, "endcap_SSring_inner_radius", 6.267);
    set_default_double_param(SUPPORTPARAMS, "endcap_SSring_outer_radius", 12.0444);
    set_default_double_param(SUPPORTPARAMS, "endcap_SSring_length", 0.0047);

    // Water Glycol endcap ring
    set_default_double_param(SUPPORTPARAMS, "endcap_WGring_inner_radius", 6.267);
    set_default_double_param(SUPPORTPARAMS, "endcap_WGring_outer_radius", 12.0444);
    set_default_double_param(SUPPORTPARAMS, "endcap_WGring_length", 0.0186);

    // CarbonPEEK endcap ring position
    set_default_double_param(SUPPORTPARAMS, "endcap_CPring_z", 24.4185);

    // CarbonPEEK endcap ring position
    set_default_double_param(SUPPORTPARAMS, "endcap_CPring_inner_radius", 6.6675);
    set_default_double_param(SUPPORTPARAMS, "endcap_CPring_outer_radius", 11.43);
    set_default_double_param(SUPPORTPARAMS, "endcap_CPring_length", 0.6370);

    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_foam_core_thickness", 0.18);
    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_inner_skin_inner_radius", 4.8);
    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_length", 42.);
    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_skin_thickness", 0.01);

    set_default_double_param(SUPPORTPARAMS, "rail_dphi", 90.);  // deg
    set_default_double_param(SUPPORTPARAMS, "rail_inner_radius", 0.45);
    set_default_double_param(SUPPORTPARAMS, "rail_length", 410);  // tpc length
    set_default_double_param(SUPPORTPARAMS, "rail_outer_radius", 0.6);
    set_default_double_param(SUPPORTPARAMS, "rail_phi_start", 45.);  // deg
    set_default_double_param(SUPPORTPARAMS, "rail_radius", 16.85);
  }

  return;
}

void PHG4InttSubsystem::Print(const string &what) const
{
  PrintDefaultParams();
  cout << endl
       << "------" << endl;
  PrintMacroParams();
  cout << endl
       << "------" << endl
       << endl;
  GetParamsContainer()->Print();
}
