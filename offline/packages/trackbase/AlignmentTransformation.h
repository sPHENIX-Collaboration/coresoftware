#ifndef TRACKBASE_ALIGNMENTTRANSFORMATION_H
#define TRACKBASE_ALIGNMENTTRANSFORMATION_H

#include "TrkrDefs.h"
#include "alignmentTransformationContainer.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <map>
#include <random>

class PHCompositeNode;

class ActsGeometry;

/**
 * Build or populate the alignment transformation map from the node tree.
 * @param topNode Top-level node containing geometry and alignment information.
 */
 
/**
 * Create and attach an alignmentTransformationContainer to the node tree and initialize it.
 * @param topNode Top-level node where the container will be created or retrieved.
 */

/**
 * Generate and store random perturbations for rotations and translations.
 * The generated samples are stored in the object's internal perturbation vectors.
 * @param angleDev Standard deviations for rotation about X, Y, Z (radians).
 * @param transformDev Standard deviations for translation along X, Y, Z (same units as detector coordinates).
 */

/**
 * Set MVTX subsystem perturbation standard deviations and enable MVTX perturbations.
 * The array layout is [angleX, angleY, angleZ, transX, transY, transZ].
 * @param mvtxDevs Six-element array with angle (first 3) then translation (last 3) standard deviations.
 */

/**
 * Set INTT subsystem perturbation standard deviations and enable INTT perturbations.
 * The array layout is [angleX, angleY, angleZ, transX, transY, transZ].
 * @param inttDevs Six-element array with angle (first 3) then translation (last 3) standard deviations.
 */

/**
 * Set TPC subsystem perturbation standard deviations and enable TPC perturbations.
 * The array layout is [angleX, angleY, angleZ, transX, transY, transZ].
 * @param tpcDevs Six-element array with angle (first 3) then translation (last 3) standard deviations.
 */

/**
 * Set MM subsystem perturbation standard deviations and enable MM perturbations.
 * The array layout is [angleX, angleY, angleZ, transX, transY, transZ].
 * @param mmDevs Six-element array with angle (first 3) then translation (last 3) standard deviations.
 */

/**
 * Enable local verbose output.
 */

/**
 * Set the misalignment factor for a given layer.
 * @param layer Layer identifier.
 * @param factor Multiplicative misalignment factor to apply for the specified layer.
 */

/**
 * Retrieve the misalignment factor for the specified layer from the active transform container.
 * @param layer Layer identifier.
 * @returns The misalignment factor for the given layer.
 */

/**
 * Toggle use of INTT survey geometry.
 * @param sur If true, survey geometry for INTT will be used.
 */

/**
 * Toggle use of the new silicon rotation application order.
 * @param flag If true, the new silicon rotation order will be used.
 */

/**
 * Toggle always-on module tilt behavior.
 * @param flag If true, module tilt will always be applied.
 */
class AlignmentTransformation
{
 public:
  AlignmentTransformation() = default;

  ~AlignmentTransformation() {}

  void createMap(PHCompositeNode* topNode);
  void createAlignmentTransformContainer(PHCompositeNode* topNode);
  void generateRandomPerturbations(Eigen::Vector3d angleDev, Eigen::Vector3d transformDev);

  bool perturbMVTX = false;
  bool perturbINTT = false;
  bool perturbTPC = false;
  bool perturbMM = false;

  Eigen::Vector3d perturbationAngles = Eigen::Vector3d(0.0, 0.0, 0.0);
  Eigen::Vector3d perturbationAnglesGlobal = Eigen::Vector3d(0.0, 0.0, 0.0);
  Eigen::Vector3d perturbationTranslation = Eigen::Vector3d(0.0, 0.0, 0.0);

  void setMVTXParams(double mvtxDevs[6])
  {
    mvtxAngleDev(0) = mvtxDevs[0];
    mvtxAngleDev(1) = mvtxDevs[1];
    mvtxAngleDev(2) = mvtxDevs[2];
    mvtxTransDev(0) = mvtxDevs[3];
    mvtxTransDev(1) = mvtxDevs[4];
    mvtxTransDev(2) = mvtxDevs[5];

    perturbMVTX = true;

    if (localVerbosity)
    {
      std::cout << "perturbMVTX: " << perturbMVTX << " MVTX Angle Std Dev: " << mvtxAngleDev << "MVTX Trans Std Dev:" << mvtxTransDev << std::endl;
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

    if (localVerbosity)
    {
      std::cout << "perturbINTT: " << perturbINTT << " INTT Angle Std Dev: " << inttAngleDev << "INTT Trans Std Dev:" << inttTransDev << std::endl;
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

    if (localVerbosity)
    {
      std::cout << "perturbTPC: " << perturbTPC << " TPC Angle Std Dev: " << tpcAngleDev << "TPC Trans Std Dev:" << tpcTransDev << std::endl;
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

    if (localVerbosity)
    {
      std::cout << "perturbMM: " << perturbMM << " MM Angle Std Dev: " << mmAngleDev << "MM Trans Std Dev:" << mmTransDev << std::endl;
    }
  }

  void verbosity() { localVerbosity = 1; }
  void misalignmentFactor(uint8_t layer, const double factor);
  double misalignmentFactor(uint8_t layer)
  {
    return transformMap->getMisalignmentFactor(layer);
  }
  void useInttSurveyGeometry(bool sur) { use_intt_survey_geometry = sur; }
  void setUseNewSiliconRotationOrder(bool flag) { use_new_silicon_rotation_order = flag; }
  void setUseModuleTiltAlways(bool flag) { use_module_tilt_always = flag; }

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

  bool use_new_silicon_rotation_order = false;
  bool use_module_tilt_always = false;
  bool use_module_tilt = false;   // starts at false in all cases
  
  bool use_intt_survey_geometry = false;
  
  Acts::Transform3 newMakeTransform(const Surface& surf, Eigen::Vector3d& millepedeTranslation, Eigen::Vector3d& sensorAngles, Eigen::Vector3d& localFrameTranslation, Eigen::Vector3d& sensorAnglesGlobal, unsigned int trkrid, bool survey);

  Eigen::Vector3d getTpcLocalFrameTranslation(float moduleRadius, float layerRadius, Eigen::Vector3d& localRotation) const; 
  void extractModuleCenterPositions();
  double extractModuleCenter(TrkrDefs::hitsetkey hitsetkey, double sectorphi);  

  alignmentTransformationContainer* transformMap = NULL;
  alignmentTransformationContainer* transformMapTransient = NULL;
  ActsGeometry* m_tGeometry = NULL;

  int getNodes(PHCompositeNode* topNode);

  // These should be checked and updated (ADF 12/2/2025)
  float TpcModuleRadii[2][12][3] = {}; // module radial center in local coords
  unsigned int innerLayer[3] = {};
  double sectorPhi[2][12] = {};
};

#endif