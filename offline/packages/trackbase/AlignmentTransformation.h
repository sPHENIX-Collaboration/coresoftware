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

  Eigen::Matrix4f makeTransform(TrkrDefs::hitsetkey hitsetkey, Eigen::Vector3f staveTrans,Eigen::Vector3f staveAngles, Eigen::Vector3f sensorTrans,Eigen::Vector3f sensorAngles);

  Eigen::Matrix4f mvtxTransform(TrkrDefs::hitsetkey hitsetkey, Eigen::Vector3f staveTrans, Eigen::Vector3f staveAngles, Eigen::Vector3f sensorTrans,Eigen::Vector3f sensorAngles);

  Eigen::Matrix4f inttTransform(TrkrDefs::hitsetkey hitsetkey, Eigen::Vector3f staveTrans, Eigen::Vector3f staveAngles, Eigen::Vector3f sensorTrans,Eigen::Vector3f sensorAngles);

  Eigen::Matrix4f makeAffineMatrix(Eigen::Vector3f rotationAnglesXYZ, Eigen::Vector3f translationVector);

  Eigen::Matrix3d rotateToGlobal(TrkrDefs::hitsetkey hitsetkey);

  //std::map<const TrkrDefs::hitsetkey, Eigen::Matrix4f>* transformMap;
  alignmentTransformationContainer* transformMap;
  ActsGeometry* m_tGeometry;
  
  int createNodes(PHCompositeNode* topNode);

};

#endif
