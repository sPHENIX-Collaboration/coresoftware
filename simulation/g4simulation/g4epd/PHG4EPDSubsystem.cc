/* vim: set sw=2 ft=cpp: */

#include "PHG4EPDSubsystem.h"

#include "PHG4EPDetector.h"
#include "PHG4EPSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4HitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

PHG4EPDSubsystem::PHG4EPDSubsystem(std::string const& name)
  : PHG4DetectorGroupSubsystem(name)
  , m_detector(nullptr)
  , m_stepaction(nullptr)
{
  InitializeParameters();
}

int32_t PHG4EPDSubsystem::InitRunSubsystem(PHCompositeNode* node) {
  PHParametersContainer* params = GetParamsContainer();
  std::string const& name = Name();

  m_detector = new PHG4EPDetector(this, node, params, name);
  m_detector->SuperDetector(SuperDetector());
  m_detector->OverlapCheck(CheckOverlap());

  m_stepaction = new PHG4EPSteppingAction(m_detector, params);
  m_stepaction->Init();

  if (!params->GetParameters(-1)->get_int_param("active"))
    return 0;

  std::string const& superdet = SuperDetector();
  std::string label = "G4HIT_" + ((superdet != "NONE") ? superdet : name);

  PHCompositeNode* dst = dynamic_cast<PHCompositeNode*>(
    PHNodeIterator(node).findFirst("PHCompositeNode", "DST"));

  PHCompositeNode* det = dynamic_cast<PHCompositeNode*>(
    PHNodeIterator(dst).findFirst("PHCompositeNode", superdet));

  if (det == nullptr) {
    det = new PHCompositeNode(superdet);

    dst->addNode(det);
  }

  if (findNode::getClass<PHG4HitContainer>(node, label.data()) == nullptr)
    det->addNode(new PHIODataNode<PHObject>(new PHG4HitContainer(label),
      label.data(), "PHObject"));

  return 0;
}

int32_t PHG4EPDSubsystem::process_event(PHCompositeNode* node) {
  if (m_stepaction != nullptr)
    m_stepaction->SetInterfacePointers(node);

  return 0;
}

PHG4Detector *PHG4EPDSubsystem::GetDetector() const {
  return m_detector;
}

PHG4SteppingAction *PHG4EPDSubsystem::GetSteppingAction() const {
  return m_stepaction;
}

void PHG4EPDSubsystem::SetDefaultParameters() {
  set_default_int_param(-1, "active", 1);
}
