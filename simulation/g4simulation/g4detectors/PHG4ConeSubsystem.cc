#include "PHG4ConeSubsystem.h"

#include "PHG4ConeDetector.h"
#include "PHG4ConeDisplayAction.h"
#include "PHG4ConeSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <cmath>  // for tan, atan, exp, M_PI
#include <sstream>

class PHG4Detector;

//_______________________________________________________________________
PHG4ConeSubsystem::PHG4ConeSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  m_ColorArray.fill(NAN);
  InitializeParameters();
}

//_______________________________________________________________________
PHG4ConeSubsystem::~PHG4ConeSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4ConeSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  // create display settings before detector
  PHG4ConeDisplayAction *disp_action = new PHG4ConeDisplayAction(Name(), GetParams());
  if (std::isfinite(m_ColorArray[0]) &&
      std::isfinite(m_ColorArray[1]) &&
      std::isfinite(m_ColorArray[2]) &&
      std::isfinite(m_ColorArray[3]))
  {
    disp_action->SetColor(m_ColorArray[0], m_ColorArray[1], m_ColorArray[2], m_ColorArray[3]);
  }
  m_DisplayAction = disp_action;

  // create detector
  m_Detector = new PHG4ConeDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));

    std::string nodename;
    if (SuperDetector() != "NONE")
    {
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
      nodename = "G4HIT_" + SuperDetector();
    }
    else
    {
      nodename = "G4HIT_" + Name();
    }
    // create hit list
    PHG4HitContainer *cone_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
    if (!cone_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(cone_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
    }
    // create stepping action
    m_SteppingAction = new PHG4ConeSteppingAction(m_Detector);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4ConeSubsystem::process_event(PHCompositeNode *topNode)
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
PHG4Detector *PHG4ConeSubsystem::GetDetector() const
{
  return m_Detector;
}

//_______________________________________________________________________
void PHG4ConeSubsystem::Set_eta_range(const double etaMin, const double etaMax)
{
  double thetaMin = 2 * atan(exp(-etaMax));
  double thetaMax = 2 * atan(exp(-etaMin));

  double z1 = get_double_param("place_z") - get_double_param("length");
  double z2 = get_double_param("place_z") + get_double_param("length");

  set_double_param("rmin1", z1 * tan(thetaMin));
  set_double_param("rmax1", z1 * tan(thetaMax));

  set_double_param("rmin2", z2 * tan(thetaMin));
  set_double_param("rmax2", z2 * tan(thetaMax));
}

void PHG4ConeSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);

  set_default_double_param("length", NAN);
  set_default_double_param("rmin1", NAN);
  set_default_double_param("rmax1", NAN);
  set_default_double_param("rmin2", NAN);
  set_default_double_param("rmax2", NAN);
  set_default_double_param("sphi", 0.);
  set_default_double_param("dphi", 360.);  // degrees
  set_default_double_param("rot_x", 0);
  set_default_double_param("rot_y", 0);
  set_default_double_param("rot_z", 0);

  set_default_string_param("material", "WorldMaterial");
}

void PHG4ConeSubsystem::SetR1(const double min, const double max)
{
  set_double_param("rmin1", min);
  set_double_param("rmax1", max);
}

void PHG4ConeSubsystem::SetR2(const double min, const double max)
{
  set_double_param("rmin2", min);
  set_double_param("rmax2", max);
}

void PHG4ConeSubsystem::SetZlength(const double a)
{
  set_double_param("length", a);
}

void PHG4ConeSubsystem::SetPhi(const double a, const double b)
{
  set_double_param("sphi", a);
  set_double_param("dphi", b);
}

void PHG4ConeSubsystem::SetPlaceZ(const double dbl)
{
  set_double_param("place_z", dbl);
}

void PHG4ConeSubsystem::SetPlace(const double place_x, const double place_y, const double place_z)
{
  set_double_param("place_x", place_x);
  set_double_param("place_y", place_y);
  set_double_param("place_z", place_z);
}

void PHG4ConeSubsystem::SetZRot(const double dbl)
{
  set_double_param("rot_z", dbl);
}

void PHG4ConeSubsystem::SetMaterial(const std::string &mat)
{
  set_string_param("material", mat);
}
