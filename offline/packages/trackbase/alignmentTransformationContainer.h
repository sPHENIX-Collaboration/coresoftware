#ifndef TRACKBASE_ALIGNMENTTRANSFORMATIONCONTAINER_H
#define TRACKBASE_ALIGNMENTTRANSFORMATIONCONTAINER_H

/**
 * @file trackbase/alignmentTranformationContainer
 * @author R. Boucher
 * @date August 2022
 * @brief Storage class for alignment transformation to node tree
 */

#include "ActsGeometry.h"
#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for pair

/**
 * @brief Container class for Alignment transformations
 *
 * Association object holding transformations associated with given tracker hitset
 */
class alignmentTransformationContainer : public Acts::GeometryContext
{
 public:
  alignmentTransformationContainer();

  virtual ~alignmentTransformationContainer() {}

  void Reset();
  void identify(std::ostream& os = std::cout);
  void addTransform(Acts::GeometryIdentifier, const Acts::Transform3&);
  Acts::Transform3& getTransform(Acts::GeometryIdentifier id);
  void replaceTransform(const Acts::GeometryIdentifier id, Acts::Transform3 transform);
  const std::vector<std::vector<Acts::Transform3>>& getMap() const;
  void setMisalignmentFactor(uint8_t layer, double factor);
  const double& getMisalignmentFactor(uint8_t layer) const { return m_misalignmentFactor.find(layer)->second; }
  static bool use_alignment;

 private:
  unsigned int getsphlayer(Acts::GeometryIdentifier);

  std::map<unsigned int, unsigned int> base_layer_map = {{10, 0}, {12, 3}, {14, 7}, {16, 55}};

  std::vector<std::vector<Acts::Transform3>> transformVec;

  /// Map of TrkrDefs::Layer to misalignment factor
  std::map<uint8_t, double> m_misalignmentFactor;

  ClassDef(alignmentTransformationContainer, 1);
};

#endif  // TRACKBASE_ALIGNMENTTRANSFORMATIONCONTAINER_H
