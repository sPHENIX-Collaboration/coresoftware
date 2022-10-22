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

  std::string alignmentParamsFile = "";

  int localVerbosity = 0;

  bool _test_translation = true;
  Acts::Vector3  _silicon_test_translation[7] = {
    {0.0, 0.0, 0.0},
    {0.1, 0.0, 0.0},
    {-0.05, 0.0, 0.0},
    {0.05, 0.0, 0.0},
    {-0.1, 0.0, 0.0},
    {0.05, 0.0, 0.0},
    {-0.05, 0.0, 0.0}
  }; // mm

  Acts::Transform3 makeTransform(Surface surf, Eigen::Vector3d millepedeTranslation, Eigen::Vector3d sensorAngles);

  Acts::Transform3 makeAffineMatrix(Eigen::Matrix3d rotationMatrix, Eigen::Vector3d translationVector);

  Eigen::Matrix3d rotateToGlobal(Surface surf);

  alignmentTransformationContainer* transformMap = NULL;
  ActsGeometry* m_tGeometry = NULL;
  
  int getNodes(PHCompositeNode* topNode);

};

#endif
