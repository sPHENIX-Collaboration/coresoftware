#include "ActsGeometry.h"
#include "alignmentTransformationContainer.h"
#include "sPHENIXActsDetectorElement.h"


sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{

  if(alignmentTransformationContainer::use_alignment)
    {
      Acts::GeometryIdentifier id = surface().geometryId();

      unsigned int layer = id.layer(); 
      unsigned int volume = id.volume(); 
      unsigned int sphlayer = base_layer_map.find(volume)->second + layer / 2 -1;

      const std::map<unsigned int, std::map<Acts::GeometryIdentifier, Acts::Transform3>>& transformMap 
	= ctxt.get<std::map<unsigned int, std::map<Acts::GeometryIdentifier, Acts::Transform3>>&>();
            
      auto it = transformMap.find(sphlayer);    // get an iterator to the map for this layer
      if(it != transformMap.end())
	{
	  //std::cout << "sPHENIXActsDetectorElement: return transform for volume " << volume <<"  layer " << layer 
	  //	    << " sphenix layer " << sphlayer << " layer map size " << it->second.size() << std::endl;
	
	  auto lyr_mapit = it->second.find(id);     //iterator to the map entry for this surface
	  if(lyr_mapit != it->second.end())
	    {
	      return lyr_mapit->second;
	    }
	}
      
      // if we are still here, it was not found
      std::cout << " Alignment transform not found, for identifier " << id << " use construction transform " << std::endl;
      const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
      //std::cout << "           construction transform: " << std::endl << transform.matrix() << std::endl;      
      return transform;
    }
 else
    {
      // return the construction transform
      const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
      return transform;
    }

}


