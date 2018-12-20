#include "PHG4BlockSubsystem.h"
#include "PHG4BlockDetector.h"
#include "PHG4BlockGeomContainer.h"
#include "PHG4BlockGeomv1.h"
#include "PHG4BlockSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Utils.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4BlockSubsystem::PHG4BlockSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
  , _detector(nullptr)
  , _steppingAction(nullptr)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4BlockSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  _detector = new PHG4BlockDetector(topNode, GetParams(), Name(), GetLayer());
  _detector->SuperDetector(SuperDetector());
  _detector->OverlapCheck(CheckOverlap());
  if (GetParams()->get_int_param("active"))
  {
    ostringstream nodename;
    ostringstream geonode;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
      geonode << "BLOCKGEOM_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
      geonode << "BLOCKGEOM_" << Name();
    }

    // create hit list
    PHG4HitContainer *block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!block_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }

    block_hits->AddLayer(GetLayer());
    PHG4BlockGeomContainer *geocont = findNode::getClass<PHG4BlockGeomContainer>(topNode,
                                                                                 geonode.str());
    if (!geocont)
    {
      geocont = new PHG4BlockGeomContainer();
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geocont, geonode.str(), "PHObject");
      runNode->addNode(newNode);
    }

    PHG4BlockGeom *geom = new PHG4BlockGeomv1(GetLayer(),
                                              GetParams()->get_double_param("size_x"),
                                              GetParams()->get_double_param("size_y"),
                                              GetParams()->get_double_param("size_z"),
                                              GetParams()->get_double_param("place_x"),
                                              GetParams()->get_double_param("place_y"),
                                              GetParams()->get_double_param("place_z"),
                                              GetParams()->get_double_param("rot_z"));
    geocont->AddLayerGeom(GetLayer(), geom);

    _steppingAction = new PHG4BlockSteppingAction(_detector, GetParams());
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    _steppingAction = new PHG4BlockSteppingAction(_detector, GetParams());
  }

  return 0;
}

//_______________________________________________________________________
int PHG4BlockSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (_steppingAction)
  {
    _steppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector *
PHG4BlockSubsystem::GetDetector(void) const
{
  return _detector;
}

void PHG4BlockSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("steplimits", NAN);
  set_default_double_param("size_x", 10.);
  set_default_double_param("size_y", 10.);
  set_default_double_param("size_z", 10.);

  set_default_int_param("use_g4steps", 0);

  set_default_string_param("material", "G4_Galactic");
}
