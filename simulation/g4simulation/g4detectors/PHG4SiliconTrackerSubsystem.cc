#include "PHG4SiliconTrackerSubsystem.h"
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
  , detector_(0)
  , steppingAction_(nullptr)
  , layerconfig_(layerconfig)
  , detector_type(detectorname)
{
  nlayers = 0;
  for (vector<pair<int, int>>::const_iterator piter = layerconfig.begin(); piter != layerconfig.end(); ++piter)
  {
    if(verbosity > 1) cout << PHWHERE << " adding INTT layer " << (*piter).second << endl;
    nlayers++;
    AddDetId((*piter).second);
  }

  // set the default sensor radii to the current design
  double sensor_radius_inner_[4] = {6.876, 8.987, 10.835, 12.676};
  double sensor_radius_outer_[4] = {7.462, 9.545, 11.361, 13.179};
  for(int i=0;i<nlayers;i++)
    {
      sensor_radius_inner[i] = sensor_radius_inner_[i];
      sensor_radius_outer[i] = sensor_radius_outer_[i];
    }

  InitializeParameters();
  // put the layer into the name so we get unique names
  // for multiple layers
  Name(detectorname);
  SuperDetector(detectorname);
}


//_______________________________________________________________________
PHG4SiliconTrackerSubsystem::PHG4SiliconTrackerSubsystem(const double sensor_radius_inner_[], const double sensor_radius_outer_[], const std::string &detectorname, const vpair &layerconfig)
   : PHG4DetectorGroupSubsystem(detectorname)
  , detector_(0)
  , steppingAction_(nullptr)
  , layerconfig_(layerconfig)
  , detector_type(detectorname)
{
  // This constructor is an ugly way to get around a problem with the parameter class - it can be removed once that is fixed
  nlayers = 0;
  for (vector<pair<int, int>>::const_iterator piter = layerconfig.begin(); piter != layerconfig.end(); ++piter)
  {
    if(verbosity > 1) cout << PHWHERE << " adding INTT layer " << (*piter).second << endl;
    nlayers++;
    AddDetId((*piter).second);
  }

  for(int i=0;i<nlayers;i++)
    {
      sensor_radius_inner[i] = sensor_radius_inner_[i];
      sensor_radius_outer[i] = sensor_radius_outer_[i];
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
  if (verbosity > 0)
    std::cout << "PHG4SiliconTrackerSubsystem::Init started" << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4SiliconTrackerDetector(topNode, GetParamsContainer(), Name(), layerconfig_);
  detector_->SuperDetector(SuperDetector());
  detector_->Detector(detector_type);
  detector_->OverlapCheck(CheckOverlap());

  int active = 0;
  int absorberactive = 0;
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
    std::string nodename = (SuperDetector() != "NONE") ? boost::str(boost::format("G4HIT_%s") % SuperDetector()) : boost::str(boost::format("G4HIT_%s") % detector_type);

    // create hit list
    PHG4HitContainer *hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
    if (!hitcontainer)
      DetNode->addNode(new PHIODataNode<PHObject>(hitcontainer = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));

    if (absorberactive)
    {
      nodename = (SuperDetector() != "NONE") ? boost::str(boost::format("G4HIT_ABSORBER_%s") % SuperDetector()) : boost::str(boost::format("G4HIT_ABSORBER_%s") % detector_type);

      hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
      if (!hitcontainer)
      {
        DetNode->addNode(new PHIODataNode<PHObject>(hitcontainer = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));
      }
    }

    // create stepping action
    steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParamsContainer());
  }
  else
  {
    if (blackhole)
    {
      steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParamsContainer());
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4SiliconTrackerSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
    steppingAction_->SetInterfacePointers(topNode);

  return 0;
}

//_______________________________________________________________________
PHG4Detector *PHG4SiliconTrackerSubsystem::GetDetector(void) const
{
  return detector_;
}

void PHG4SiliconTrackerSubsystem::SetDefaultParameters()
{
  // We have only two types of ladders, one with vertical strips (0) and one with horizontal strips (1)
  // There are 4 sensors in each ladder
  //     In ladder type 0 the sensor is special and inner and outer sensors are the same. 
  //     In ladder type 1 there are two different sensor types, inner and outer 

  // We define default parameters for laddertype 0 or 1, all dimensions in cm
  int nstrips_phi_sensor[2] = {1, 256};  
  int nstrips_z_sensor_0[2] = {128*5, 8};  // inner sensor  
  int nstrips_z_sensor_1[2] = {128*5, 5};  // outer sensor
  double strip_y[2] = {1.6, 0.0078};
  double strip_z_0[2] = {0.01406, 1.6};
  double strip_z_1[2] = {0.01406, 2.0};
  double halfladder_z[2] = {48.0, 48.0};
  double hdi_y[2] = {2.55, 3.8};
  double stave_straight_outer_y[2] = {0.672, 0.522};
  double stave_straight_inner_y[2] = {0.1, 0.344};  // the first value is a dummy, not used, to avoid issues with making a G4Logical;Volume with y = 0
  double stave_straight_cooler_y[2] = {0.47, 0.47};
  double sensor_offset_y[2] = {0.304, 0.0};

  // We do not want to hard code the ladder types for the layers
  // We define default ladder types for 4 layers, but these can be changed at the macro level
  int laddertype[4] = {0, 1, 1, 1};
  int nladder[4] = {34, 30, 36, 42};
  // sensor radius_inner and sensor_radius_outer are set in the constructor for now, to avoid a problem with the parameter class

  for(int i=0;i<nlayers;i++)
    {
      // To reconfigure the layers, all you have to do is overide the defaults for these four arrays from the tracking macro
      set_default_int_param(i, "laddertype", laddertype[i]);
      set_default_int_param(i, "nladder", nladder[i]);  // ladders per layer
      set_default_double_param(i, "sensor_radius_inner", sensor_radius_inner[i]*cm);
      set_default_double_param(i, "sensor_radius_outer", sensor_radius_outer[i]*cm);
      //cout << " PHG4SiliconTrackerSubsystem setting default parameters to: " << endl;
      //cout << "  layer " << i << " laddertype " << laddertype[i] << " nladder " << nladder[i] 
      //   << " sensor_radius_inner " << sensor_radius_inner[i] << " sensor_radius_outer " << sensor_radius_outer[i] << endl;
      // These should be kept at zero in the new design
      set_default_double_param(i, "offsetphi", 0.);
      set_default_double_param(i, "offsetrot", 0.);
    }

  // Set the parameters for the two laddertypes. All values in cm!
  for (int i = 0; i < 2; i++)
  {
    set_default_int_param(i, "nstrips_z_sensor_0", nstrips_z_sensor_0[i]);
    set_default_int_param(i, "nstrips_z_sensor_1", nstrips_z_sensor_1[i]);
    set_default_int_param(i, "nstrips_phi_sensor", nstrips_phi_sensor[i]);
    set_default_int_param(i, "nstrips_phi_cell", nstrips_phi_sensor[i]);

    set_default_double_param(i, "stave_straight_inner_y", stave_straight_inner_y[i]*cm);
    set_default_double_param(i, "stave_straight_outer_y", stave_straight_outer_y[i]*cm);
    set_default_double_param(i, "stave_straight_cooler_y", stave_straight_cooler_y[i]*cm);

    set_default_double_param(i, "strip_y", strip_y[i]*cm);
    set_default_double_param(i, "strip_z_0", strip_z_0[i]*cm);
    set_default_double_param(i, "strip_z_1", strip_z_1[i]*cm);
    set_default_double_param(i, "strip_x", 0.032*cm);  // 320 microns deep
    set_default_double_param(i, "sensor_edge_phi", 0.13*cm);
    set_default_double_param(i, "sensor_edge_z", 0.10*cm);
    set_default_double_param(i, "sensor_offset_y", sensor_offset_y[i]*cm);
    set_default_double_param(i, "hdi_kapton_x", 0.038*cm);
    set_default_double_param(i, "hdi_copper_x", 0.0052*cm);  // effective width of all copper in ground layers and signal layers
    set_default_double_param(i, "hdi_y", hdi_y[i]*cm);
    set_default_double_param(i, "hdi_edge_z", 0.0*cm);
    set_default_double_param(i, "fphx_x", 0.032*cm); 
    set_default_double_param(i, "fphx_y", 0.27*cm);
    set_default_double_param(i, "fphx_z", 0.91*cm);
    set_default_double_param(i, "gap_sensor_fphx", 0.1*cm);
    set_default_double_param(i, "pgs_x", 0.02*cm);  // 0.2 mm
    set_default_double_param(i, "halfladder_z", halfladder_z[i]*cm);
  }

  //std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator> begin_end = GetDetIds();
  //for (set<int>::const_iterator it = begin_end.first; it != begin_end.second; ++it)
  for (int i=0; i < nlayers; ++i)
  {
    set_default_int_param(i, "active", 1);
  }

  return;
}

void PHG4SiliconTrackerSubsystem::Print(const string &what) const
{
  PrintDefaultParams();
  PrintMacroParams();
  GetParamsContainer()->Print();
}
