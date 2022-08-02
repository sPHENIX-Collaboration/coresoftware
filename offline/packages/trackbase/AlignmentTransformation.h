#ifndef TRACKBASE_ALIGNMENTTRANSFORMATION_H 
#define TRACKBASE_ALIGNMENTTRANSFORMATION_H
#include <map>
#include <trackbase/TrkrDefs.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>



class AlignmentTransformation {

 public:

  AlignmentTransformation() = default;
  ~AlignmentTransformation() {} 

  void createMap();

 private:

  Eigen::Matrix4f makeTransform(int layer, int stave, int chip, Eigen::Vector3f staveTrans,Eigen::Vector3f staveAngles,Eigen::Vector3f sensorTrans,Eigen::Vector3f sensorAngles);

  Eigen::Matrix4f makeAffineMatrix(Eigen::Vector3f rotationAnglesXYZ, Eigen::Vector3f translationVector);

  std::map<TrkrDefs::hitsetkey, Eigen::Matrix4f> transformMap;


};

#endif
