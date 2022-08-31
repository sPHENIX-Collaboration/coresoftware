#ifndef TRACKBASE_ALIGNMENTTRANSFORMATION_H 
#define TRACKBASE_ALIGNMENTTRANSFORMATION_H
#include <map>

#include "TrkrDefs.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "alignmentTransformationContainer.h"

class PHCompositeNode;

class ActsGeometry;

class AlignmentTransformation {

 public:

  AlignmentTransformation() = default;

  ~AlignmentTransformation() {} 

  void createMap(PHCompositeNode* topNode);
  void createAlignmentTransformContainer(PHCompositeNode* topNode);

 private:

  std::string alignmentParamsFile = "";

  bool localVerbosity = false;

  Acts::Transform3 makeTransform(Surface surf, Eigen::Vector3d millepedeTranslation, Eigen::Vector3d sensorAngles);

  Acts::Transform3 makeAffineMatrix(Eigen::Matrix3d rotationMatrix, Eigen::Vector3d translationVector);

  Eigen::Matrix3d rotateToGlobal(Surface surf);

  alignmentTransformationContainer* transformMap;
  ActsGeometry* m_tGeometry;
  
  int getNodes(PHCompositeNode* topNode);

};

#endif
