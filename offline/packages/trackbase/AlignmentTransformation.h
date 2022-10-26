#ifndef TRACKBASE_ALIGNMENTTRANSFORMATION_H 
#define TRACKBASE_ALIGNMENTTRANSFORMATION_H

#include "alignmentTransformationContainer.h"
#include "TrkrDefs.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <map>


class PHCompositeNode;

class ActsGeometry;

class AlignmentTransformation {

 public:

  AlignmentTransformation() = default;

  ~AlignmentTransformation() {} 

  void createMap(PHCompositeNode* topNode);
  void createAlignmentTransformContainer(PHCompositeNode* topNode);

 private:

  std::string alignmentParamsFile = "./localAlignmentParamsFile.txt";

  int localVerbosity = 0;

  Acts::Transform3 makeTransform(Surface surf, Eigen::Vector3d millepedeTranslation, Eigen::Vector3d sensorAngles);

  Acts::Transform3 makeAffineMatrix(Eigen::Matrix3d rotationMatrix, Eigen::Vector3d translationVector);

  Eigen::Matrix3d rotateToGlobal(Surface surf);

  alignmentTransformationContainer* transformMap = NULL;
  ActsGeometry* m_tGeometry = NULL;
  
  int getNodes(PHCompositeNode* topNode);

};

#endif
