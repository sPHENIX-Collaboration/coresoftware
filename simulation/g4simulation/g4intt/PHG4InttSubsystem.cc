#include "PHG4InttSubsystem.h"
#include "PHG4InttDefs.h"
#include "PHG4InttDetector.h"
#include "PHG4InttDisplayAction.h"
#include "PHG4InttSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>
#include <g4detectors/PHG4DetectorGroupSubsystem.h>  // for PHG4DetectorGrou...

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <iostream>     // for operator<<, basi...
#include <set>          // for _Rb_tree_const_i...
#include <type_traits>  // for __decay_and_stri...

class PHG4Detector;

//_______________________________________________________________________
PHG4InttSubsystem::PHG4InttSubsystem(const std::string &detectorname, const vpair &layerconfig)
  : PHG4DetectorGroupSubsystem(detectorname)
  , m_LayerConfigVector(layerconfig)
  , m_DetectorType(detectorname)
{
  for (std::vector<std::pair<int, int>>::const_iterator piter = layerconfig.begin(); piter != layerconfig.end(); ++piter)
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
  std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> layer_begin_end = std::make_pair(m_LayerConfigVector.begin(), m_LayerConfigVector.end());
  m_Detector = new PHG4InttDetector(this, topNode, GetParamsContainer(), Name(), layer_begin_end);
  m_Detector->Verbosity(Verbosity());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->Detector(m_DetectorType);
  m_Detector->OverlapCheck(CheckOverlap());

  int active = 0;
  // initialize with support active flag (if support is active we need the absorber hit node)
  int absorberactive = GetParamsContainer()->GetParameters(PHG4InttDefs::SUPPORTPARAMS)->get_int_param("supportactive");
  int blackhole = 0;
  for (std::set<int>::const_iterator parcontaineriter = GetDetIds().first; parcontaineriter != GetDetIds().second; ++parcontaineriter)
  {
    if (active || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("active"))
    {
      active = 1;
    }
    if (absorberactive || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("absorberactive") || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("supportactive"))
    {
      absorberactive = 1;
    }
    if (blackhole || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("blackhole"))
    {
      blackhole = 1;
    }
  }
  std::set<std::string> nodes;
  if (active)
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dstNode;
    if (SuperDetector() != "NONE" && !SuperDetector().empty())
    {
      PHNodeIterator iter_dst(dstNode);
      DetNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(DetNode);
      }
    }
    // create hit output nodes
    std::string detector_suffix = SuperDetector();
    if (detector_suffix == "NONE" || detector_suffix.empty())
    {
      detector_suffix = Name();
    }
    m_HitNodeName = "G4HIT_" + detector_suffix;
    nodes.insert(m_HitNodeName);
    m_AbsorberNodeName = "G4HIT_ABSORBER_" + detector_suffix;
    if (absorberactive)
    {
      nodes.insert(m_AbsorberNodeName);
    }
    for (auto nodename : nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4InttSteppingAction(m_Detector, GetParamsContainer(), layer_begin_end);
    m_SteppingAction->Verbosity(Verbosity());
    m_SteppingAction->SetHitNodeName("G4HIT", m_HitNodeName);
    m_SteppingAction->SetHitNodeName("G4HIT_ABSORBER", m_AbsorberNodeName);
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
  // We have only two types of ladders, one with vertical strips (SEGMENTATION_Z) and one with horizontal strips (SEGMENTATION_PHI)set
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
                       PHG4InttDefs::SEGMENTATION_PHI};  // default

  int nladder[8] = {17, 17, 12, 12, 16, 16, 21, 21};  // default, new 03/05/2020

  double sensor_radius[8] = {6.876, 7.462,
                             // 4 elements are those for PHG4InttDefs::SEGMENTATION_PHI, 36um subtracted to set si sensors at the place
                             // these subtractions are due to different thickness of glue for the sensors (14um) and the FPHX chips (50um)
                             7.188 - 36e-4, 7.732 - 36e-4, 9.680 - 36e-4, 10.262 - 36e-4,
                             12.676, 13.179};  // radius of center of sensor for layer default, new 30/05/2020

  double offsetphi[4] = {-0.5 * 360.0 / nladder[0+2], 0.0, -0.5 * 360.0 / nladder[2+2], 0.0 }; // the final configuration, July/09/202
 
  auto detid = GetDetIds();  // get pair of iterators to begin/end of set<int> of detids
  for (auto detiter = detid.first; detiter != detid.second; ++detiter)
  {
    set_default_int_param(*detiter, "active", 1);

    // To reconfigure the layers, all you have to do is overide the defaults for these four arrays from the tracking macro
    set_default_int_param(*detiter, "laddertype", laddertype[*detiter]);
    set_default_int_param(*detiter, "nladder", nladder[*detiter]);  // ladders per layer
    set_default_double_param(*detiter, "sensor_radius", sensor_radius[*detiter]);
    // These offsets should be kept at zero in the new design
    //  set_default_double_param(*detiter, "offsetphi", 0.);// obsolete
    set_default_double_param(*detiter, "offsetphi", offsetphi[*detiter] );
    set_default_double_param(*detiter, "offsetrot", 0.);

    // 	sitrack->set_int_param(i, "laddertype", laddertype[i]);
  }

  // These are the parameters that describe the internal ladder geometry for the two ladder types
  // SEGMENTATION_Z //////////////////////////////////////
  // int param
  set_default_int_param(PHG4InttDefs::SEGMENTATION_Z, "nstrips_phi_cell", 1);
  set_default_int_param(PHG4InttDefs::SEGMENTATION_Z, "nstrips_phi_sensor", 1);
  set_default_int_param(PHG4InttDefs::SEGMENTATION_Z, "nstrips_z_sensor_0", 128 * 5);
  set_default_int_param(PHG4InttDefs::SEGMENTATION_Z, "nstrips_z_sensor_1", 128 * 5);

  // double param
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "fphx_x", 0.032);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "fphx_y", 0.27);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "fphx_z", 0.91);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "fphx_offset_z", 0.005);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "gap_sensor_fphx", 0.1);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "halfladder_z", 40.00);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "halfladder_inside_z", 23.9622);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "hdi_copper_x", 0.0052);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "hdi_edge_z", 0.);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "hdi_kapton_x", 0.038);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "hdi_y", 2.55);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "pgs_x", 0.02);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "sensor_edge_phi", 0.13);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "sensor_edge_z", 0.1);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "sensor_offset_y", 0.304);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "strip_x", 0.032);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "strip_y", 1.6);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "strip_z_0", 0.01406);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "strip_z_1", 0.01406);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "stave_straight_cooler_x", 0.01905);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "stave_straight_cooler_y", 0.47);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "stave_slant_cooler_y", 1.4362);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "stave_straight_outer_y", 0.672);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_Z, "stave_straight_rohacell_y", 0.57181);

  // SEGMENTATION_PHI //////////////////////////////////////
  // int param
  set_default_int_param(PHG4InttDefs::SEGMENTATION_PHI, "nstrips_phi_cell", 256);
  set_default_int_param(PHG4InttDefs::SEGMENTATION_PHI, "nstrips_phi_sensor", 256);
  set_default_int_param(PHG4InttDefs::SEGMENTATION_PHI, "nstrips_z_sensor_0", 8);
  set_default_int_param(PHG4InttDefs::SEGMENTATION_PHI, "nstrips_z_sensor_1", 5);

  // double param
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "fphx_x", 0.032);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "fphx_y", 0.27);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "fphx_z", 0.91);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "fphx_offset_z", 0.005);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "gap_sensor_fphx", 0.1);

  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "si_glue_x", 0.0014);   // 14 um, don't forget to change double sensor_radius when it's changed
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "fphx_glue_x", 0.005);  // 50 um

  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "halfladder_z", 40.00);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "halfladder_inside_z", 23.9622);

  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "hdi_copper_x", 0.00376);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "hdi_edge_z", 0.);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "hdi_kapton_x", 0.038);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "hdi_y", 3.8);

  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "sensor_edge_phi", 0.13);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "sensor_edge_z", 0.1);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "sensor_offset_y", 0.);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "strip_x", 0.032);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "strip_y", 0.0078);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "strip_z_0", 1.6);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "strip_z_1", 2.);

  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "stave_straight_cooler_x", 0.03);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "stave_straight_cooler_y", 1.47684);

  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "stave_slant_cooler_y", 0.6322614829);

  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "stave_straight_outer_y", 0.33227);
  set_default_double_param(PHG4InttDefs::SEGMENTATION_PHI, "stave_straight_rohacell_y", 0.58842);

  // SUPPORTPARAMS //////////////////////////////////////
  // int param
  set_default_int_param(PHG4InttDefs::SUPPORTPARAMS, "supportactive", 0);

  // double param
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "inner_skin_inner_radius", 6.2416);
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "inner_skin_length",      50.7   );
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "inner_skin_outer_radius", 6.2666);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "inner_skin_inner_radius", 12.9667 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "inner_skin_outer_radius", 13.0175 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "inner_skin_length", 49.7);

  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_cfcin_inner_radius", 12.0444);
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_cfcin_outer_radius", 12.0694);
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_cfcin_length",       50.7   );

  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_foam_inner_radius", 12.0694);
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_foam_outer_radius", 12.2194);
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_foam_length",       50.7   );

  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_cfcout_inner_radius", 12.2194);
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_cfcout_outer_radius", 12.2444);
  // set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_cfcout_length",       50.7   );

  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_inner_radius", 23.4950 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_outer_radius", 23.5458 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "outer_skin_length", 49.7);

  // Endcap ring flag
  set_default_int_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_ring_enabled", 1);
  set_default_int_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_ring_type", 2);  // 0: Al+SS+WG, 1 : CarbonPEEK, 2(default) : new  model Jan/2021

  // Aluminum endcap ring position
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_ring_z", 24.35);

  // Aluminum endcap ring
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_Alring_inner_radius", 6.267);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_Alring_outer_radius", 12.0444);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_Alring_length", 0.3645);

  // Stainless steel endcap ring
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_SSring_inner_radius", 6.267);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_SSring_outer_radius", 12.0444);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_SSring_length", 0.0047);

  // Water Glycol endcap ring
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_WGring_inner_radius", 6.267);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_WGring_outer_radius", 12.0444);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_WGring_length", 0.0186);

  // CarbonPEEK endcap ring position
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_CPring_z", 24.4185);

  // CarbonPEEK endcap ring position
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_CPring_inner_radius", 6.6675);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_CPring_outer_radius", 11.43);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_CPring_length", 0.6370);

  ////////////////////////////////////////////////////////////////////////////////////////
  // the new endcap model
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Alring_z", 24.4185);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Alring_1_outer_radius", 11.7475);  // outer radius of the outermost part
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Cring_1_outer_radius", 11.2020);   // outer radius of the 2nd outermost part
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Alring_2_outer_radius", 9.65);     // outer radius of the 3rd outermost part, slightly shrinked from the reeeal drawing of 9.6971 cm to avoid overwlapping
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Cring_2_outer_radius", 8.7095);    // outer radius of the 4th outermost part
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Alring_3_outer_radius", 7.15);     // outer radius of the 5th outermost part, slightly shrinked from the real drawing of 7.2045 cm to avoid overlapping
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Alring_3_inner_radius", 6.5088);   // inner radius of the 5th outermost (=the outer most) part

  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Alring_length", 0.75);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "endcap_AlPEEK_Cring_length", 0.5);

  ////////////////////////////////////////////////////////////////////////////////////////
  // Survice barrel, outer
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "service_barrel_outer_inner_radius", 33.02 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "service_barrel_outer_outer_radius", 33.34 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "service_barrel_outer_length", 273.69);

  ////////////////////////////////////////////////////////////////////////////////////////
  // Support tube
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "support_tube_inner_radius", 37.47 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "support_tube_outer_radius", 38.10 / 2);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "support_tube_length", 273.69);

  // Cylinders for the bus extenders
  set_default_int_param(PHG4InttDefs::SUPPORTPARAMS, "bus_extender", 1);                           // 0: OFF, 1: ON
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "bus_extender_length", 111.0);             // in cm
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "bus_extender_ends_at", 328.5);            // z-coordinate in cm where the bus extender ends at
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "bus_extender_radius", 15.0);              // radius of the innermost layer (copper for the inner barrel)
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "bus_extender_copper_x", 48.0e-4 * 1.5);   // thickness of the copper layer of the bus extenders in cm, it's 48 um
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "bus_extender_kapton_x", 300.0e-4 * 1.5);  // thickness of the kapton layer of the bus extenders in cm, it's 300 um

  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_dphi", 90.);  // deg
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_inner_radius", 0.45);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_length", 410);  // tpc length
  //set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_length",     20   );  // tpc length
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_outer_radius", 0.6);
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_phi_start", 45.);  // deg
  //set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_radius",     16.85 );
  set_default_double_param(PHG4InttDefs::SUPPORTPARAMS, "rail_radius", (33.34 + 0.6 * 2) / 2);  // tentativevalue

  return;
}

void PHG4InttSubsystem::Print(const std::string & /*what*/) const
{
  PrintDefaultParams();
  std::cout << std::endl
            << "------" << std::endl;
  PrintMacroParams();
  std::cout << std::endl
            << "------" << std::endl
            << std::endl;
  GetParamsContainer()->Print();
}
