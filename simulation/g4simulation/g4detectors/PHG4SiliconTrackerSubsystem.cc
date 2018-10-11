#include "PHG4SiliconTrackerSubsystem.h"
#include "PHG4SiliconTrackerDefs.h"
#include "PHG4SiliconTrackerDetector.h"
#include "PHG4SiliconTrackerSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

#include <boost/format.hpp>

using namespace std;

//_______________________________________________________________________
PHG4SiliconTrackerSubsystem::PHG4SiliconTrackerSubsystem(const std::string &detectorname, const vpair &layerconfig)
  : PHG4DetectorGroupSubsystem(detectorname)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
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

//_______________________________________________________________________
int PHG4SiliconTrackerSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4SiliconTrackerSubsystem::Init started" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  pair<vector<pair<int, int>>::const_iterator, vector<pair<int, int>>::const_iterator> layer_begin_end = make_pair(m_LayerConfigVector.begin(), m_LayerConfigVector.end());
  m_Detector = new PHG4SiliconTrackerDetector(topNode, GetParamsContainer(), Name(), layer_begin_end);
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->Detector(m_DetectorType);
  m_Detector->OverlapCheck(CheckOverlap());

  int active = 0;
  // initialize with support active flag (if support is active we need the absorber hit node)
  int absorberactive = GetParamsContainer()->GetParameters(PHG4SiliconTrackerDefs::SUPPORTPARAMS)->get_int_param("supportactive");
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
    m_SteppingAction = new PHG4SiliconTrackerSteppingAction(m_Detector, GetParamsContainer(), layer_begin_end);
  }
  else
  {
    if (blackhole)
    {
      m_SteppingAction = new PHG4SiliconTrackerSteppingAction(m_Detector, GetParamsContainer(), layer_begin_end);
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4SiliconTrackerSubsystem::process_event(PHCompositeNode *topNode)
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
PHG4Detector *PHG4SiliconTrackerSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4SiliconTrackerSubsystem::SetDefaultParameters()
{
  // We have only two types of ladders, one with vertical strips (SEGMENTATION_Z) and one with horizontal strips (SEGMENTATION_PHI)
  // There are 4 sensors in each ladder
  //     In ladder type 0 the sensor is special and inner and outer sensors are the same.
  //     In ladder type 1 there are two different sensor types, inner and outer
  // We do not want to hard code the ladder types for the layers

  // We define default ladder types for 8 layers, but these can be changed at the macro level

  int laddertype[8] = {PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI};  // default
  int nladder[8] = {17,  17, 15, 15, 18, 18, 21, 21};  // default
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
    using namespace PHG4SiliconTrackerDefs;
    set_default_int_param(SEGMENTATION_Z, "nstrips_phi_cell", 1);
    set_default_int_param(SEGMENTATION_Z, "nstrips_phi_sensor", 1);
    set_default_int_param(SEGMENTATION_Z, "nstrips_z_sensor_0", 128 * 5);
    set_default_int_param(SEGMENTATION_Z, "nstrips_z_sensor_1", 128 * 5);
    set_default_double_param(SEGMENTATION_Z, "fphx_x", 0.032);
    set_default_double_param(SEGMENTATION_Z, "fphx_y", 0.27);
    set_default_double_param(SEGMENTATION_Z, "fphx_z", 0.91);
    set_default_double_param(SEGMENTATION_Z, "gap_sensor_fphx", 0.1);
    set_default_double_param(SEGMENTATION_Z, "halfladder_z", 48.0);
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

    set_default_int_param(SEGMENTATION_PHI, "nstrips_phi_cell", 256);
    set_default_int_param(SEGMENTATION_PHI, "nstrips_phi_sensor", 256);
    set_default_int_param(SEGMENTATION_PHI, "nstrips_z_sensor_0", 8);
    set_default_int_param(SEGMENTATION_PHI, "nstrips_z_sensor_1", 5);
    set_default_double_param(SEGMENTATION_PHI, "fphx_x", 0.032);
    set_default_double_param(SEGMENTATION_PHI, "fphx_y", 0.27);
    set_default_double_param(SEGMENTATION_PHI, "fphx_z", 0.91);
    set_default_double_param(SEGMENTATION_PHI, "gap_sensor_fphx", 0.1);
    set_default_double_param(SEGMENTATION_PHI, "halfladder_z", 48.0);
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

    set_default_int_param(SUPPORTPARAMS, "supportactive", 0);

    set_default_double_param(SUPPORTPARAMS, "inner_skin_inner_radius",6.385);
    set_default_double_param(SUPPORTPARAMS, "inner_skin_length",96);
    set_default_double_param(SUPPORTPARAMS, "inner_skin_outer_radius",6.4);

    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_foam_core_thickness",0.18);
    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_inner_skin_inner_radius",4.8);
    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_length",42.);
    set_default_double_param(SUPPORTPARAMS, "mvtx_shell_skin_thickness",0.01);

    set_default_double_param(SUPPORTPARAMS, "rail_dphi",60.); // deg
    set_default_double_param(SUPPORTPARAMS, "rail_inner_radius", 0.45);
    set_default_double_param(SUPPORTPARAMS, "rail_length", 410); // tpc length
    set_default_double_param(SUPPORTPARAMS, "rail_outer_radius", 0.6);
    set_default_double_param(SUPPORTPARAMS, "rail_phi_start",30.); // deg
    set_default_double_param(SUPPORTPARAMS, "rail_radius",17.5);

    set_default_double_param(SUPPORTPARAMS, "outer_skin_inner_radius",15.7);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_outer_radius",15.8);
    set_default_double_param(SUPPORTPARAMS, "outer_skin_length",96);
  }

  return;
}

void PHG4SiliconTrackerSubsystem::Print(const string &what) const
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
