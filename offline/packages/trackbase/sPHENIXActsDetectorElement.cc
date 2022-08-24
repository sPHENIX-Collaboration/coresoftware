#include "sPHENIXActsDetectorElement.h"

sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{

  // if(m_alignmentParameters.empty())
  //   {
  return TGeoDetectorElement::transform(ctxt);
  //   }

  // else
    //{
      // auto id  = surface().geometryId();
      // auto map = ctxt.get<std::map<const Acts::GeometryIdentifier, Acts::Transform3>>();

      // Acts::Transform3& transform = map.find(id)->second;

      // return transform;
      // }

}
