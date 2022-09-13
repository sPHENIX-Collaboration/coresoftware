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

  void setMVTXParams(double mvtxAlphaDev, double mvtxBetaDev, double mvtxGammaDev, double mvtxDxDev, double mvtxDyDev, double mvtxDzDev)
  {
    mvtxAngleDev(0) = mvtxAlphaDev;
    mvtxAngleDev(1) = mvtxBetaDev;
    mvtxAngleDev(2) = mvtxGammaDev;
    mvtxTransDev(0) = mvtxDxDev;
    mvtxTransDev(1) = mvtxDyDev;
    mvtxTransDev(2) = mvtxDzDev;

    perturbMVTX = true;

    if(localVerbosity)
      {
	std::cout << "perturbMVTX: "<<perturbMVTX << " alpha in: " << mvtxAlphaDev <<" mvtxAngleDevs: "<< mvtxAngleDev<< std::endl;
	std::cout << mvtxAngleDev << " " << mvtxTransDev << " " << perturbMVTX << std::endl;
      }
  }

  void setINTTParams(double inttAlphaDev, double inttBetaDev, double inttGammaDev, double inttDxDev, double inttDyDev, double inttDzDev)
  {
    inttAngleDev(0) = inttAlphaDev;
    inttAngleDev(1) = inttBetaDev;
    inttAngleDev(2) = inttGammaDev;
    inttTransDev(0) = inttDxDev;
    inttTransDev(1) = inttDyDev;
    inttTransDev(2) = inttDzDev;

    perturbINTT = true;

    if(localVerbosity)
      {
	std::cout << "perturbINTT: "<<perturbINTT << " " << inttAlphaDev << std::endl;
	std::cout << inttAngleDev << " " << inttTransDev << " " << perturbINTT << std::endl;
      }
  }

  void setTPCParams(double tpcAlphaDev, double tpcBetaDev, double tpcGammaDev, double tpcDxDev, double tpcDyDev, double tpcDzDev)
  {
    tpcAngleDev(0) = tpcAlphaDev;
    tpcAngleDev(1) = tpcBetaDev;
    tpcAngleDev(2) = tpcGammaDev;
    tpcTransDev(0) = tpcDxDev;
    tpcTransDev(1) = tpcDyDev;
    tpcTransDev(2) = tpcDzDev;

    perturbTPC = true;

    if(localVerbosity)
      {
	std::cout << "perturbTPC: "<<perturbTPC << " " << tpcAlphaDev << std::endl;
	std::cout << "tpcAngleDevs " << tpcAngleDev << " tpctransDev "<< tpcTransDev <<std::endl;
	std::cout << tpcAngleDev << " " << tpcTransDev << " " << perturbTPC << std::endl;
      }
  }

  void setMMParams(double mmAlphaDev, double mmBetaDev, double mmGammaDev, double mmDxDev, double mmDyDev, double mmDzDev)
  {
    mmAngleDev(0) = mmAlphaDev;
    mmAngleDev(1) = mmBetaDev;
    mmAngleDev(2) = mmGammaDev;
    mmTransDev(0) = mmDxDev;
    mmTransDev(1) = mmDyDev;
    mmTransDev(2) = mmDzDev;

    perturbMM = true;

    if(localVerbosity)
      {
	std::cout << "perturbMM: "<<perturbMM << " " << mmAlphaDev << std::endl;
	std::cout << mmAngleDev << " " << mmTransDev << " " << perturbMM << std::endl;
      }
  }

 private:

  std::string alignmentParamsFile = "";
  
  Eigen::Vector3d mvtxAngleDev;
  Eigen::Vector3d mvtxTransDev;
  Eigen::Vector3d inttAngleDev;
  Eigen::Vector3d inttTransDev;  
  Eigen::Vector3d tpcAngleDev;
  Eigen::Vector3d tpcTransDev;  
  Eigen::Vector3d mmAngleDev;
  Eigen::Vector3d mmTransDev;

  bool localVerbosity = false;

  Acts::Transform3 makeTransform(Surface surf, Eigen::Vector3d millepedeTranslation, Eigen::Vector3d sensorAngles);

  Acts::Transform3 makeAffineMatrix(Eigen::Matrix3d rotationMatrix, Eigen::Vector3d translationVector);

  Eigen::Matrix3d rotateToGlobal(Surface surf);

  alignmentTransformationContainer* transformMap = NULL;
  ActsGeometry* m_tGeometry = NULL;
  
  int getNodes(PHCompositeNode* topNode);

};

#endif
