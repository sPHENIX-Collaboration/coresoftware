#include "ActsGeometry.h"
#include "alignmentTransformationContainer.h"
#include "sPHENIXActsDetectorElement.h"

#include <phool/phool.h>

sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{

  if(alignmentTransformationContainer::use_alignment)
    {
      Acts::GeometryIdentifier id = surface().geometryId();

      unsigned int volume = id.volume(); 
      unsigned int layer = id.layer(); 
      unsigned int sphlayer = base_layer_map.find(volume)->second + layer / 2 -1;
      unsigned int sensor = id.sensitive() - 1;  // Acts sensor ID starts at 1
 
      const alignmentTransformationContainer* transformContainer 
	= ctxt.get<alignmentTransformationContainer*>();
    
      const auto& transformVec = transformContainer->getMap();
      auto& layerVec = transformVec[sphlayer];    // get the vector of transforms for this layer
   
      if(layerVec.size() > sensor)
	{
	  //std::cout << "sPHENIXActsDetectorElement: return transform:  volume " << volume <<" Acts  layer " << layer << " sensor " << sensor
	  //	  	    << " sphenix layer " << sphlayer << " layerVec size " << layerVec.size() << std::endl;
	
	  return layerVec[sensor];
	}
      
      // if we are still here, it was not found
      std::cout << PHWHERE << " Alignment transform not found, for identifier " << id << " continuing on with ideal geometry is not ideal so we exit" << std::endl;
   
      exit(1);
   
    }
 else
    {
      // return the construction transform
      const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
      return transform;
    }

}


