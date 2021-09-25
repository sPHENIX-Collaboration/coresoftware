/* vim: set sw=2 ft=cpp: */

#include "PHG4EPDSubsystem.h"

#include "PHG4EPDDetector.h"
#include "PHG4EPDDisplayAction.h"
#include "PHG4EPDSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4detectors/PHG4DetectorGroupSubsystem.h>  // for PHG4DetectorGroupSubsystem

#include <g4main/PHG4HitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

PHG4EPDSubsystem::PHG4EPDSubsystem(std::string const& name)
  : PHG4DetectorSubsystem(name)
{
  InitializeParameters();
}

int32_t PHG4EPDSubsystem::InitRunSubsystem(PHCompositeNode* node)
{
//  PHParametersContainer* params = GetParamsContainer();
  std::string const& name = Name();

  m_DisplayAction = new PHG4EPDDisplayAction(Name());

m_Detector = new PHG4EPDDetector(this, node, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  m_SteppingAction = new PHG4EPDSteppingAction(m_Detector, GetParams());

  if (!GetParams()->get_int_param("active"))
    return 0;

  std::string const& superdet = SuperDetector();
  std::string label = "G4HIT_" + ((superdet != "NONE") ? superdet : name);

  PHCompositeNode* dst = dynamic_cast<PHCompositeNode*>(
      PHNodeIterator(node).findFirst("PHCompositeNode", "DST"));

  PHCompositeNode* det = dynamic_cast<PHCompositeNode*>(
      PHNodeIterator(dst).findFirst("PHCompositeNode", superdet));

  if (det == nullptr)
  {
    det = new PHCompositeNode(superdet);

    dst->addNode(det);
  }

  if (findNode::getClass<PHG4HitContainer>(node, label.data()) == nullptr)
    det->addNode(new PHIODataNode<PHObject>(new PHG4HitContainer(label),
                                            label.data(), "PHObject"));

  return 0;
}

int32_t PHG4EPDSubsystem::process_event(PHCompositeNode* node)
{
  if (m_SteppingAction != nullptr)
    m_SteppingAction->SetInterfacePointers(node);

  return 0;
}

PHG4Detector* PHG4EPDSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4EPDSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_z", 300.);
}
