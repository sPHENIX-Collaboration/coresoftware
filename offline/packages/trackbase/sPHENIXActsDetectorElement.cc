#include "sPHENIXActsDetectorElement.h"


sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{
  auto map = ctxt.getContainer();
  auto id = identifier;
  Acts::Transform3 transform = map->find(id);
  return transform;
}
