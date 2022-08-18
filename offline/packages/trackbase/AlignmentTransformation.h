#ifndef TRACKBASE_ALIGNMENTTRANSFORMATION_H 
#define TRACKBASE_ALIGNMENTTRANSFORMATION_H
#include <map>
#include <trackbase/TrkrDefs.h>
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

 private:

  int localVerbosity = false;

  Eigen::Matrix4d makeTransform(TrkrDefs::hitsetkey hitsetkey, Eigen::Vector3d millepedeTranslation, Eigen::Vector3d sensorAngles, unsigned int subsurfkey);

  Eigen::Matrix4d makeAffineMatrix(Eigen::Matrix3d rotationMatrix, Eigen::Vector3d translationVector);

  Eigen::Matrix3d rotateToGlobal(TrkrDefs::hitsetkey hitsetkey, unsigned int subsurfkey);

  //std::map<const TrkrDefs::hitsetkey, Eigen::Matrix4d>* transformMap;
  alignmentTransformationContainer* transformMap;
  ActsGeometry* m_tGeometry;
  
  int createNodes(PHCompositeNode* topNode);

};

#endif
