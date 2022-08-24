#include "sPHENIXActsDetectorElement.h"

#include <ActsGeometry.h>

sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

bool sPHENIXActsDetectorElement::use_alignment = false;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{
  // Typedef for Acts::Transform3 is: using Transform3 = Eigen::Transform<ActsScalar, 3, Eigen::Affine>;

  if(use_alignment)
    {
      std::cout << "map is filled " << std::endl;
     const std::map<Acts::GeometryIdentifier, Acts::Transform3>& map = ctxt.get<std::map<Acts::GeometryIdentifier, Acts::Transform3>>();
     //auto map = ctxt.get<std::map<Acts::GeometryIdentifier, Acts::Transform3>&>();

      std::cout << "      map size " << map.size() << std::endl; 
      Acts::GeometryIdentifier id = surface().geometryId();
      std::cout << " Get transform for identifier " << id << std::endl;
      const Acts::Transform3& transform = map.find(id)->second;
      std::cout << "          transform is: " << transform.matrix() << std::endl;

      return transform;
    }
  else
    {
      const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
      return transform;
    }

}

