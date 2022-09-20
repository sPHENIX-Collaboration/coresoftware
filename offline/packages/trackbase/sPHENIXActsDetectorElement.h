#ifndef TRACKBASE_SPHENIXACTSDETECTORELEMENT_H
#define TRACKBASE_SPHENIXACTSDETECTORELEMENT_H

#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>
#include <Acts/Plugins/Identification/Identifier.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>

/**
 * This class implements an sphenix detector element to build
 * the Acts geometry with. Relevant info on how to handle alignment
 * information from the GeometryContext is logically processed here
 */

class ActsGeometry;

class sPHENIXActsDetectorElement : public Acts::TGeoDetectorElement
{
 public:

  sPHENIXActsDetectorElement() = delete;
  
  sPHENIXActsDetectorElement(const Identifier& identifier, 
			     const TGeoNode& tGeoNode,
			     const TGeoMatrix& tGeoMatrix = TGeoIdentity(),
			     const std::string& axes = "XYZ", 
			     double scalor = 10.,
			     std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr) : 
    Acts::TGeoDetectorElement(identifier, tGeoNode, tGeoMatrix, 
			      axes, scalor, material) {}
    
  sPHENIXActsDetectorElement(const Identifier& identifier, 
			     const TGeoNode& tGeoNode,
			      Acts::Transform3& tgTransform,
			     std::shared_ptr<const Acts::PlanarBounds> tgBounds,
			     double tgThickness = 0.) : 
     Acts::TGeoDetectorElement(identifier, tGeoNode, tgTransform, 
			       tgBounds, tgThickness) {}
  
     
  sPHENIXActsDetectorElement(const Identifier& identifier, 
			     const TGeoNode& tGeoNode,
			     Acts::Transform3& tgTransform,
			     std::shared_ptr<const Acts::DiscBounds> tgBounds,
			     double tgThickness = 0.) :
    Acts::TGeoDetectorElement(identifier, tGeoNode, tgTransform, tgBounds, 
			       tgThickness) {}
     
  ~sPHENIXActsDetectorElement() override;

  const Acts::Transform3& transform(const Acts::GeometryContext& ctxt) const override;


private:

  std::map<unsigned int, unsigned int> base_layer_map = { {10, 0}, {12,3}, {14,7}, {16,55} };

};

std::shared_ptr<sPHENIXActsDetectorElement> sPHENIXElementFactory(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const TGeoMatrix& tGeoMatrix, const std::string& axes, double scalor,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  return std::make_shared<sPHENIXActsDetectorElement>(identifier, tGeoNode, tGeoMatrix,
                                               axes, scalor, material);
}

#endif
