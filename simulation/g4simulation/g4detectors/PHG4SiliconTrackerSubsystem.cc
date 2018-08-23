#include "PHG4SiliconTrackerSubsystem.h"
#include "PHG4SiliconTrackerDetector.h"
#include "PHG4SiliconTrackerSteppingAction.h"
#include "PHG4SiliconTrackerDefs.h"

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
    cout << PHWHERE << " adding INTT layer " << (*piter).second << endl;
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
  {
    std::cout << "PHG4SiliconTrackerSubsystem::Init started" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4SiliconTrackerDetector(topNode, GetParamsContainer(), Name(), GetDetIds());
  detector_->SuperDetector(SuperDetector());
  detector_->Detector(detector_type);
  detector_->OverlapCheck(CheckOverlap());

  int active = 0;
  int absorberactive = 0;
  int blackhole = 0;
  for (set<int>::const_iterator parcontaineriter = GetDetIds().first; parcontaineriter != GetDetIds().second; ++parcontaineriter)
  {
    const PHParameters *par = GetParamsContainer()->GetParameters(*parcontaineriter);
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
    steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParamsContainer(), GetDetIds());
  }
  else
  {
    if (blackhole)
    {
      steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParamsContainer(),GetDetIds());
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


  // We do not want to hard code the ladder types for the layers
  // We define default ladder types for 4 layers, but these can be changed at the macro level
  int laddertype[4] = {PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI};
  int nladder[4] = {34, 30, 36, 42};
  // sensor radius_inner and sensor_radius_outer are set in the constructor for now, to avoid a problem with the parameter class
  cout << "layers: " << nlayers << endl;
  auto detid = GetDetIds(); // get pair of iterators to begin/end of set<int> of detids
  for (auto detiter = detid.first; detiter != detid.second; ++detiter)
  {
      set_default_int_param(*detiter,"active",1);

      // To reconfigure the layers, all you have to do is overide the defaults for these four arrays from the tracking macro
      set_default_int_param(*detiter, "laddertype", laddertype[*detiter]);
      set_default_int_param(*detiter, "nladder", nladder[*detiter]);  // ladders per layer
      set_default_double_param(*detiter, "sensor_radius_inner", sensor_radius_inner[*detiter]*cm);
      set_default_double_param(*detiter, "sensor_radius_outer", sensor_radius_outer[*detiter]*cm);
      // These offsets should be kept at zero in the new design
    set_default_double_param(*detiter,"offsetphi",0.);
    set_default_double_param(*detiter,"offsetrot",0.);
      cout << " PHG4SiliconTrackerSubsystem setting default parameters to: " << endl;
      cout << "  layer " << *detiter << " laddertype " << laddertype[*detiter] << " nladder " << nladder[*detiter] 
	   << " sensor_radius_inner " << sensor_radius_inner[*detiter] << " sensor_radius_outer " << sensor_radius_outer[*detiter] << endl;
    }
  { // just being lazy, using namespace in this scope for less clutter
    using namespace PHG4SiliconTrackerDefs; 
    set_default_int_param(SEGMENTATION_Z,"nstrips_phi_cell",1);
    set_default_int_param(SEGMENTATION_Z,"nstrips_phi_sensor",1);
    set_default_int_param(SEGMENTATION_Z,"nstrips_z_sensor_0",128*5);
    set_default_int_param(SEGMENTATION_Z,"nstrips_z_sensor_1",128*5);
    set_default_double_param(SEGMENTATION_Z,"fphx_x",0.032*cm);
    set_default_double_param(SEGMENTATION_Z,"fphx_y",0.27*cm);
    set_default_double_param(SEGMENTATION_Z,"fphx_z",0.91*cm);
    set_default_double_param(SEGMENTATION_Z,"gap_sensor_fphx",0.1*cm);
    set_default_double_param(SEGMENTATION_Z,"halfladder_z",48.0*cm);
    set_default_double_param(SEGMENTATION_Z,"hdi_copper_x",0.0052*cm);
    set_default_double_param(SEGMENTATION_Z,"hdi_edge_z",0.*cm);
    set_default_double_param(SEGMENTATION_Z,"hdi_kapton_x",0.038*cm);
    set_default_double_param(SEGMENTATION_Z,"hdi_y",2.55*cm);
    set_default_double_param(SEGMENTATION_Z,"pgs_x",0.02*cm);
    set_default_double_param(SEGMENTATION_Z,"sensor_edge_phi",0.13*cm);
    set_default_double_param(SEGMENTATION_Z,"sensor_edge_z",0.1*cm);
    set_default_double_param(SEGMENTATION_Z,"sensor_offset_y",0.304*cm);
    set_default_double_param(SEGMENTATION_Z,"stave_straight_cooler_y",0.47*cm);
    set_default_double_param(SEGMENTATION_Z,"stave_straight_inner_y",0.1*cm);
    set_default_double_param(SEGMENTATION_Z,"stave_straight_outer_y",0.672*cm);
    set_default_double_param(SEGMENTATION_Z,"strip_x",0.032*cm);
    set_default_double_param(SEGMENTATION_Z,"strip_y",1.6*cm);
    set_default_double_param(SEGMENTATION_Z,"strip_z_0",0.01406*cm);
    set_default_double_param(SEGMENTATION_Z,"strip_z_1",0.01406*cm);

    set_default_int_param(SEGMENTATION_PHI,"nstrips_phi_cell",256);
    set_default_int_param(SEGMENTATION_PHI,"nstrips_phi_sensor",256);
    set_default_int_param(SEGMENTATION_PHI,"nstrips_z_sensor_0",8);
    set_default_int_param(SEGMENTATION_PHI,"nstrips_z_sensor_1",5);
    set_default_double_param(SEGMENTATION_PHI,"fphx_x",0.032*cm);
    set_default_double_param(SEGMENTATION_PHI,"fphx_y",0.27*cm);
    set_default_double_param(SEGMENTATION_PHI,"fphx_z",0.91*cm);
    set_default_double_param(SEGMENTATION_PHI,"gap_sensor_fphx",0.1*cm);
    set_default_double_param(SEGMENTATION_PHI,"halfladder_z",48.0*cm);
    set_default_double_param(SEGMENTATION_PHI,"hdi_copper_x",0.0052*cm);
    set_default_double_param(SEGMENTATION_PHI,"hdi_edge_z",0.*cm);
    set_default_double_param(SEGMENTATION_PHI,"hdi_kapton_x",0.038*cm);
    set_default_double_param(SEGMENTATION_PHI,"hdi_y",3.8*cm);
    set_default_double_param(SEGMENTATION_PHI,"pgs_x",0.02*cm);
    set_default_double_param(SEGMENTATION_PHI,"sensor_edge_phi",0.13*cm);
    set_default_double_param(SEGMENTATION_PHI,"sensor_edge_z",0.1*cm);
    set_default_double_param(SEGMENTATION_PHI,"sensor_offset_y",0.*cm);
    set_default_double_param(SEGMENTATION_PHI,"stave_straight_cooler_y",0.47*cm);
    set_default_double_param(SEGMENTATION_PHI,"stave_straight_inner_y",0.344*cm);
    set_default_double_param(SEGMENTATION_PHI,"stave_straight_outer_y",0.522*cm);
    set_default_double_param(SEGMENTATION_PHI,"strip_x",0.032*cm);
    set_default_double_param(SEGMENTATION_PHI,"strip_y",0.0078*cm);
    set_default_double_param(SEGMENTATION_PHI,"strip_z_0",1.6*cm);
    set_default_double_param(SEGMENTATION_PHI,"strip_z_1",2.*cm);
  }


  return;
}

void PHG4SiliconTrackerSubsystem::Print(const string &what) const
{
  PrintDefaultParams();
  cout << endl << "------" << endl;
  PrintMacroParams();
  cout << endl << "------" << endl << endl;
  GetParamsContainer()->Print();
}
