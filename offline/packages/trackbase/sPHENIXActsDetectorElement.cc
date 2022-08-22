#include "sPHENIXActsDetectorElement.h"


sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{
  auto map = ctxt.getContainer();
  Acts::Transform3& transform = map->find(identifier());
  return transform;
}
