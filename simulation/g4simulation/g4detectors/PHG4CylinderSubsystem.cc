#include "PHG4CylinderSubsystem.h"
#include "PHG4CylinderDetector.h"
#include "PHG4CylinderDisplayAction.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv1.h"
#include "PHG4CylinderSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <cmath>     // for NAN
#include <iostream>  // for operator<<, basic_ostream, endl
#include <sstream>

class PHG4CylinderGeom;
class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4CylinderSubsystem::PHG4CylinderSubsystem(const std::string &na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
{
  m_ColorArray.fill(NAN);
  InitializeParameters();
}

//_______________________________________________________________________
PHG4CylinderSubsystem::~PHG4CylinderSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4CylinderSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  // create hit list only for active layers
  double detlength = GetParams()->get_double_param("length");
  if (!isfinite(detlength) && GetParams()->get_int_param("lengthviarapidity"))
  {
    GetParams()->set_double_param("length", PHG4Utils::GetLengthForRapidityCoverage(GetParams()->get_double_param("radius") + GetParams()->get_double_param("thickness")) * 2);
    detlength = GetParams()->get_double_param("length");
  }
  else
  {
    GetParams()->set_int_param("lengthviarapidity", 0);
  }
  // use world material if material was not set so far
  if (GetParams()->get_string_param("material") == "WorldMaterial")
  {
    recoConsts *rc = recoConsts::instance();
    GetParams()->set_string_param("material", rc->get_StringFlag("WorldMaterial"));
  }
  // create display settings before detector
  PHG4CylinderDisplayAction *disp_action = new PHG4CylinderDisplayAction(Name(), GetParams());
  if (isfinite(m_ColorArray[0]) &&
      isfinite(m_ColorArray[1]) &&
      isfinite(m_ColorArray[2]) &&
      isfinite(m_ColorArray[3]))
  {
    disp_action->SetColor(m_ColorArray[0], m_ColorArray[1], m_ColorArray[2], m_ColorArray[3]);
  }
  m_DisplayAction = disp_action;

  // create detector
  m_Detector = new PHG4CylinderDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));

    string nodename;
    string geonode;
    if (SuperDetector() != "NONE")
    {
      // create super detector subnodes
      PHNodeIterator iter_dst(dstNode);
      PHCompositeNode *superSubNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(superSubNode);
      }
      dstNode = superSubNode;
      PHNodeIterator iter_run(runNode);
      superSubNode = dynamic_cast<PHCompositeNode *>(iter_run.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        runNode->addNode(superSubNode);
      }
      runNode = superSubNode;

      nodename = "G4HIT_" + SuperDetector();
      geonode = "CYLINDERGEOM_" + SuperDetector();
    }

    else
    {
      nodename = "G4HIT_" + Name();
      geonode = "CYLINDERGEOM_" + Name();
    }
    PHG4HitContainer *cylinder_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
    if (!cylinder_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(cylinder_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
    }
    cylinder_hits->AddLayer(GetLayer());
    PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonode);
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode, "PHObject");
      runNode->addNode(newNode);
    }
    PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv1(GetParams()->get_double_param("radius"), GetParams()->get_double_param("place_z") - detlength / 2., GetParams()->get_double_param("place_z") + detlength / 2., GetParams()->get_double_param("thickness"));
    geo->AddLayerGeom(GetLayer(), mygeom);
    auto *tmp = new PHG4CylinderSteppingAction(this, m_Detector, GetParams());
    tmp->HitNodeName(nodename);
    m_SteppingAction = tmp;
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new PHG4CylinderSteppingAction(this, m_Detector, GetParams());
  }
  if (m_SteppingAction)
  {
    (dynamic_cast<PHG4CylinderSteppingAction *>(m_SteppingAction))->SaveAllHits(m_SaveAllHitsFlag);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4CylinderSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4CylinderSubsystem::SetDefaultParameters()
{
  set_default_double_param("length", NAN);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("radius", NAN);
  set_default_double_param("steplimits", NAN);
  set_default_double_param("thickness", NAN);
  set_default_double_param("tmin", NAN);
  set_default_double_param("tmax", NAN);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("start_phi_rad", 0.);
  set_default_double_param("delta_phi_rad", M_PI * 2);
  set_default_int_param("lengthviarapidity", 1);
  set_default_int_param("lightyield", 0);
  set_default_int_param("use_g4steps", 0);

  // place holder, will be replaced by world material if not set by other means (macro)
  set_default_string_param("material", "WorldMaterial");
}

PHG4Detector *
PHG4CylinderSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4CylinderSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  if (!BeginRunExecuted())
  {
    cout << "Need to execute BeginRun() before parameter printout is meaningful" << endl;
    cout << "To do so either run one or more events or on the command line execute: " << endl;
    cout << "Fun4AllServer *se = Fun4AllServer::instance();" << endl;
    cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << endl;
    cout << "g4->InitRun(se->topNode());" << endl;
    cout << "PHG4CylinderSubsystem *cyl = (PHG4CylinderSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << endl;
    cout << "cyl->Print()" << endl;
    return;
  }
  GetParams()->Print();
  if (m_SteppingAction)
  {
    m_SteppingAction->Print(what);
  }
  return;
}
