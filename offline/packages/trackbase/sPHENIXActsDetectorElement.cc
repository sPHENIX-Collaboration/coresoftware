#include "sPHENIXActsDetectorElement.h"
#include "alignmentTransformationContainer.h"

#include <ActsGeometry.h>

sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

bool sPHENIXActsDetectorElement::use_alignment = false;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{

  if(alignmentTransformationContainer::use_alignment)
    {
      Acts::GeometryIdentifier id = surface().geometryId();
      //std::cout << " sPHENIXActsDetectorElement: Get alignment transform for identifier: " << id << std::endl;

      const std::map<Acts::GeometryIdentifier, Acts::Transform3>& map =  ctxt.get<std::map<Acts::GeometryIdentifier, Acts::Transform3>>();

      auto it = map.find(id);
      if(it!=map.end())
	{
	  const Acts::Transform3& transform = it->second;
	  //std::cout << "        got alignment transform: " << std::endl << transform.matrix() << std::endl;
	  return transform;
	}
      else
	{
	  std::cout << " Alignment transform not found, for identifier " << id << " use construction transform " << std::endl;
	  const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
	  //std::cout << "           construction transform: " << std::endl << transform.matrix() << std::endl;      
	  return transform;
	}
    }
  else
    {
      // return the construction transform
      const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
      return transform;
    }

}

