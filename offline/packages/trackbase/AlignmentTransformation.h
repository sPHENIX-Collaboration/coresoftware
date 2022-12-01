#ifndef TRACKBASE_ALIGNMENTTRANSFORMATION_H 
#define TRACKBASE_ALIGNMENTTRANSFORMATION_H

#include "alignmentTransformationContainer.h"
#include "TrkrDefs.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <map>
#include <random>



class PHCompositeNode;

class ActsGeometry;

class AlignmentTransformation {

 public:

  AlignmentTransformation() = default;

  ~AlignmentTransformation(){} 

  void createMap(PHCompositeNode* topNode);
  void createAlignmentTransformContainer(PHCompositeNode* topNode);

  void generateRandomPerturbations(Eigen::Vector3d angleDev, Eigen::Vector3d transformDev);

  bool perturbMVTX = false;
  bool perturbINTT = false;
  bool perturbTPC  = false;
  bool perturbMM   = false;

  Eigen::Vector3d perturbationAngles      = Eigen::Vector3d(0.0,0.0,0.0);
  Eigen::Vector3d perturbationTranslation = Eigen::Vector3d(0.0,0.0,0.0);

 void setMVTXParams(double mvtxDevs[6])
  {
    mvtxAngleDev(0) = mvtxDevs[0];
    mvtxAngleDev(1) = mvtxDevs[1];
    mvtxAngleDev(2) = mvtxDevs[2];
    mvtxTransDev(0) = mvtxDevs[3];
    mvtxTransDev(1) = mvtxDevs[4];
    mvtxTransDev(2) = mvtxDevs[5];

    perturbMVTX = true;

    if(localVerbosity)
      {
	std::cout << "perturbMVTX: "<<perturbMVTX <<" MVTX Angle Std Dev: " << mvtxAngleDev <<"MVTX Trans Std Dev:"<< mvtxTransDev<< std::endl;
      }
  }

 void setINTTParams(double inttDevs[6])
  {
    inttAngleDev(0) = inttDevs[0];
    inttAngleDev(1) = inttDevs[1];
    inttAngleDev(2) = inttDevs[2];
    inttTransDev(0) = inttDevs[3];
    inttTransDev(1) = inttDevs[4];
    inttTransDev(2) = inttDevs[5];

    perturbINTT = true;

    if(localVerbosity)
      {
	std::cout << "perturbINTT: "<<perturbINTT <<" INTT Angle Std Dev: " << inttAngleDev <<"INTT Trans Std Dev:"<< inttTransDev<< std::endl;
      }
  } 
void setTPCParams(double tpcDevs[6])
  {
    tpcAngleDev(0) = tpcDevs[0];
    tpcAngleDev(1) = tpcDevs[1];
    tpcAngleDev(2) = tpcDevs[2];
    tpcTransDev(0) = tpcDevs[3];
    tpcTransDev(1) = tpcDevs[4];
    tpcTransDev(2) = tpcDevs[5];

    perturbTPC = true;

    if(localVerbosity)
      {
	std::cout << "perturbTPC: "<<perturbTPC <<" TPC Angle Std Dev: " << tpcAngleDev <<"TPC Trans Std Dev:"<< tpcTransDev<< std::endl;
      }
  }
 void setMMParams(double mmDevs[6])
  {
    mmAngleDev(0) = mmDevs[0];
    mmAngleDev(1) = mmDevs[1];
    mmAngleDev(2) = mmDevs[2];
    mmTransDev(0) = mmDevs[3];
    mmTransDev(1) = mmDevs[4];
    mmTransDev(2) = mmDevs[5];

    perturbMM = true;

    if(localVerbosity)
      {
	std::cout << "perturbMM: "<<perturbMM <<" MM Angle Std Dev: " << mmAngleDev <<"MM Trans Std Dev:"<< mmTransDev<< std::endl;
      }
  }

 void misalignmentFactor(TrkrDefs::TrkrId id, const double factor);

 private:

  Eigen::Vector3d mvtxAngleDev;
  Eigen::Vector3d mvtxTransDev;
  Eigen::Vector3d inttAngleDev;
  Eigen::Vector3d inttTransDev;  
  Eigen::Vector3d tpcAngleDev;
  Eigen::Vector3d tpcTransDev;  
  Eigen::Vector3d mmAngleDev;
  Eigen::Vector3d mmTransDev;

  std::default_random_engine generator;

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
