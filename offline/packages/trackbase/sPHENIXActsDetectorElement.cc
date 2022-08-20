#include "sPHENIXActsDetectorElement.h"

sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;


const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{

  // using Transform3 = Eigen::Transform<ActsScalar, 3, Eigen::Affine>;

  auto map = ctxt.get<std::map<const Acts::GeometryIdentifier, Acts::Transform3>>();
  Acts::Transform3& transform = map.find(identifier())->second;
  return transform;
}

