#include "PHG4SiliconTrackerSubsystem.h"
#include "PHG4Parameters.h"
#include "PHG4ParametersContainer.h"
#include "PHG4SiliconTrackerDetector.h"
#include "PHG4SiliconTrackerSteppingAction.h"

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
  int nladder[4] = {20, 26, 32, 38};
  int nstrips_z_sensor_0[4] = {5, 8, 8, 8};
  int nstrips_z_sensor_1[4] = {5, 5, 5, 5};
  double halfladder_z[4] = {22.0, 26.8, 26.8, 26.8};
  double hdi_y[4] = {3.8, 4.3, 4.3, 4.3};
  double offsetrot[4] = {14.0, 14.0, 12.0, 11.5};
  double radius[4] = {6.0, 8.0, 10.0, 12.0};
  double strip_y[4] = {0.0078, 0.0086, 0.0086, 0.0086};
  double strip_z_0[4] = {1.8, 1.6, 1.6, 1.6};
  double strip_z_1[4] = {1.8, 2.0, 2.0, 2.0};

  // all values in cm!
  for (int i = 0; i < nlayers; i++)
  {
    set_default_int_param(i, "nstrips_phi_cell", 128);
    set_default_double_param(i, "fphx_x", 0.032);
    set_default_double_param(i, "fphx_y", 0.27);
    set_default_double_param(i, "fphx_z", 0.9);
    set_default_double_param(i, "gap_sensor_fphx", 0.1);
    set_default_double_param(i, "hdi_x", 0.038626);
    set_default_double_param(i, "hdi_edge_z", 0.01);
    set_default_double_param(i, "offsetphi", 0.);
    set_default_double_param(i, "pgs_x", 0.021);
    set_default_double_param(i, "sensor_edge_phi", 0.1305);
    set_default_double_param(i, "sensor_edge_z", 0.098);
    set_default_double_param(i, "stave_x", 0.023);
    set_default_double_param(i, "strip_x", 0.02);

    set_default_int_param(i, "nladder", nladder[i]);
    set_default_int_param(i, "nstrips_z_sensor_0", nstrips_z_sensor_0[i]);
    set_default_int_param(i, "nstrips_z_sensor_1", nstrips_z_sensor_1[i]);
    set_default_double_param(i, "halfladder_z", halfladder_z[i]);
    set_default_double_param(i, "hdi_y", hdi_y[i]);
    set_default_double_param(i, "offsetrot", offsetrot[i]);
    set_default_double_param(i, "radius", radius[i]);
    set_default_double_param(i, "strip_y", strip_y[i]);
    set_default_double_param(i, "strip_z_0", strip_z_0[i]);
    set_default_double_param(i, "strip_z_1", strip_z_1[i]);
  }

  /*
  set_default_int_param(0, "nladder", 22);
  set_default_int_param(1, "nladder", 26);
  set_default_int_param(2, "nladder", 32);
  set_default_int_param(3, "nladder", 38);

  set_default_int_param(0, "nstrips_z_sensor_0", 5);
  set_default_int_param(1, "nstrips_z_sensor_0", 8);
  set_default_int_param(2, "nstrips_z_sensor_0", 8);
  set_default_int_param(3, "nstrips_z_sensor_0", 8);

  set_default_int_param(0, "nstrips_z_sensor_1", 5);
  set_default_int_param(1, "nstrips_z_sensor_1", 5);
  set_default_int_param(2, "nstrips_z_sensor_1", 5);
  set_default_int_param(3, "nstrips_z_sensor_1", 5);

  set_default_double_param(0, "halfladder_z", 22.);
  set_default_double_param(1, "halfladder_z", 26.8);
  set_default_double_param(2, "halfladder_z", 26.8);
  set_default_double_param(3, "halfladder_z", 26.8);

  set_default_double_param(0, "hdi_y", 3.8);
  set_default_double_param(1, "hdi_y", 4.3);
  set_default_double_param(2, "hdi_y", 4.3);
  set_default_double_param(3, "hdi_y", 4.3);

  set_default_double_param(0, "offsetrot", 14.0);
  set_default_double_param(1, "offsetrot", 14.0);
  set_default_double_param(2, "offsetrot", 12.0);
  set_default_double_param(3, "offsetrot", 11.5);

  set_default_double_param(0, "radius", 6.);
  set_default_double_param(1, "radius", 8.);
  set_default_double_param(2, "radius", 10.);
  set_default_double_param(3, "radius", 12.);

  set_default_double_param(0, "strip_y", 0.0078);
  set_default_double_param(1, "strip_y", 0.0086);
  set_default_double_param(2, "strip_y", 0.0086);
  set_default_double_param(3, "strip_y", 0.0086);

  set_default_double_param(0, "strip_z_0", 1.8);
  set_default_double_param(1, "strip_z_0", 1.6);
  set_default_double_param(2, "strip_z_0", 1.6);
  set_default_double_param(3, "strip_z_0", 1.6);

  set_default_double_param(0, "strip_z_1", 1.8);
  set_default_double_param(1, "strip_z_1", 2.0);
  set_default_double_param(2, "strip_z_1", 2.0);
  set_default_double_param(3, "strip_z_1", 2.0);
  */

  std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator> begin_end = GetDetIds();
  for (set<int>::const_iterator it = begin_end.first; it != begin_end.second; ++it)
  {
    set_default_int_param(*it, "active", 1);
  }

  return;
}

void PHG4SiliconTrackerSubsystem::Print(const string &what) const
{
  PrintDefaultParams();
  PrintMacroParams();
  GetParamsContainer()->Print();
}
