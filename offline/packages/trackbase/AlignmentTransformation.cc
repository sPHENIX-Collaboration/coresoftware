#include "AlignmentTransformation.h"

#include "ActsGeometry.h"
#include "TpcDefs.h"
#include "TrkrDefs.h"

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>

#include <cmath>
#include <fstream>
#include <sstream>

void AlignmentTransformation::createMap(PHCompositeNode* topNode)
{
  localVerbosity = 0;
  // The default is to use translation parameters that are in global coordinates
  std::cout << "AlignmentTransformation: use INTT survey geometry = " << use_intt_survey_geometry << std::endl;
  std::cout << "AlignmentTransformation: localVerbosity = " << localVerbosity << std::endl;

  getNodes(topNode);

  // Use construction transforms as a reference for making the map
  if (alignmentTransformationContainer::use_alignment)
  {
    alignmentTransformationContainer::use_alignment = false;
  }

  // Define Parsing Variables
  TrkrDefs::hitsetkey hitsetkey = 0;
  //  float alpha = 0.0, beta = 0.0, gamma = 0.0, dr = 0.0, rdphi = 0.0, dz = 0.0;
  float alpha = 0.0, beta = 0.0, gamma = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;

  // load alignment constants file
  std::ifstream datafile;
  datafile.open(alignmentParamsFile);  //  looks for default file name on disk
  if (datafile.is_open())
  {
    std::cout << "AlignmentTransformation: Reading alignment parameters from disk file: "
              << alignmentParamsFile << " localVerbosity = " << localVerbosity << std::endl;
  }
  else
  {
    datafile.clear();
    // load alignment constants file from database
    alignmentParamsFile = CDBInterface::instance()->getUrl("TRACKINGALIGNMENT");
    std::cout << "AlignmentTransformation: Reading alignment parameters from database file: " << alignmentParamsFile << std::endl;
    datafile.open(alignmentParamsFile);
  }

  // check to see how many parameters per line in the file
  // If it is old, there will only be six. In that case, set the global rotation pars to zero, and issue a warning.


  ActsSurfaceMaps surfMaps = m_tGeometry->maps();
  Surface surf;

  int fileLines = 1824;
  for (int i = 0; i < fileLines; i++)
  {
    // guard against reading in bad parameter files
    std::string str; 
    std::getline(datafile, str);
    std::stringstream ss(str);
    std::string dummy;
    int count = 0;
    while(ss >> dummy)
      {
	count ++;
      }
    if(count > 6)
      {
	std::stringstream str6(str);
	str6 >>  hitsetkey >> alpha >> beta >> gamma >> dx >> dy >> dz;
      }
    else
      {
	std::cout << "WARNING: fields = " << count << " you are reading a bad alignment parameters file, quit!" << std::endl;
	  exit(1);
      }

    if(localVerbosity > 0)
      {
	std::cout  <<  hitsetkey << "  " << alpha  << "  " << beta  << "  " << gamma  << "  " << dx  << "  " << dy << "  "  << dz 
		   << std::endl;  
      }

    // Perturbation translations and angles for stave and sensor
    Eigen::Vector3d sensorAngles(alpha, beta, gamma);
    Eigen::Vector3d millepedeTranslation(dx, dy, dz);

    unsigned int trkrId = TrkrDefs::getTrkrId(hitsetkey);  // specify between detectors

    perturbationAngles = Eigen::Vector3d(0.0, 0.0, 0.0);
    perturbationAnglesGlobal = Eigen::Vector3d(0.0, 0.0, 0.0);
    perturbationTranslation = Eigen::Vector3d(0.0, 0.0, 0.0);

    if (trkrId == TrkrDefs::mvtxId)
    {
      if (perturbMVTX)
      {
        generateRandomPerturbations(mvtxAngleDev, mvtxTransDev);
        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }

      surf = surfMaps.getSiliconSurface(hitsetkey);

      Acts::Transform3 transform;
      transform = newMakeTransform(surf, millepedeTranslation, sensorAngles, false);

      Acts::GeometryIdentifier id = surf->geometryId();

      if (localVerbosity)
      {
        std::cout << " Add transform for MVTX with surface GeometryIdentifier " << id << " trkrid " << trkrId << std::endl;
        std::cout << " final mvtx transform:" << std::endl
                  << transform.matrix() << std::endl;
      }
      transformMap->addTransform(id, transform);
      transformMapTransient->addTransform(id, transform);
    }

    else if (trkrId == TrkrDefs::inttId)
    {
      if (perturbINTT)
      {
        generateRandomPerturbations(inttAngleDev, inttTransDev);
        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }

      surf = surfMaps.getSiliconSurface(hitsetkey);

      Acts::Transform3 transform;
      transform = newMakeTransform(surf, millepedeTranslation, sensorAngles, use_intt_survey_geometry);
      Acts::GeometryIdentifier id = surf->geometryId();

      if (localVerbosity)
      {
        std::cout << " Add transform for INTT with surface GeometryIdentifier " << id << " trkrid " << trkrId << std::endl;
      }

      transformMap->addTransform(id, transform);
      transformMapTransient->addTransform(id, transform);
    }

    else if (trkrId == TrkrDefs::tpcId)
    {
      if (perturbTPC)
      {
        generateRandomPerturbations(tpcAngleDev, tpcTransDev);
        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }
      unsigned int sector = TpcDefs::getSectorId(hitsetkey);
      unsigned int side = TpcDefs::getSide(hitsetkey);
      int subsurfkey_min = (1 - side) * 144 + (144 - sector * 12) - 12 - 6;
      int subsurfkey_max = subsurfkey_min + 12;
      // std::cout << " sector " << sector << " side " << side << " subsurfkey_min " << subsurfkey_min << " subsurfkey_max " << subsurfkey_max << std::endl;

      for (int subsurfkey = subsurfkey_min; subsurfkey < subsurfkey_max; subsurfkey++)
      {
        int sskey = subsurfkey;
        if (sskey < 0)
        {
          sskey += 288;
        }

        surf = surfMaps.getTpcSurface(hitsetkey, (unsigned int) sskey);

        Acts::Transform3 transform;
        transform = newMakeTransform(surf, millepedeTranslation, sensorAngles, false);
        Acts::GeometryIdentifier id = surf->geometryId();

        if (localVerbosity)
        {
          unsigned int layer = TrkrDefs::getLayer(hitsetkey);
          std::cout << " Add transform for TPC with surface GeometryIdentifier " << id << std::endl
                    << " trkrid " << trkrId << " hitsetkey " << hitsetkey << " layer " << layer << " sector " << sector << " side " << side
                    << " subsurfkey " << subsurfkey << std::endl;
          Acts::Vector3 center = surf->center(m_tGeometry->geometry().getGeoContext()) * 0.1;  // convert to cm
          std::cout << "Ideal surface center: " << std::endl
                    << center << std::endl;
          std::cout << "transform matrix: " << std::endl
                    << transform.matrix() << std::endl;
        }
        transformMap->addTransform(id, transform);
        transformMapTransient->addTransform(id, transform);
      }
    }
    else if (trkrId == TrkrDefs::micromegasId)
    {
      if (perturbMM)
      {
        generateRandomPerturbations(mmAngleDev, mmTransDev);

        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }
      surf = surfMaps.getMMSurface(hitsetkey);

      Acts::Transform3 transform;
      transform = newMakeTransform(surf, millepedeTranslation, sensorAngles,  false);
      Acts::GeometryIdentifier id = surf->geometryId();

      if (localVerbosity)
      {
        std::cout << " Add transform for Micromegas with surface GeometryIdentifier " << id << " trkrid " << trkrId << std::endl;
      }

      transformMap->addTransform(id, transform);
      transformMapTransient->addTransform(id, transform);
    }

    else
    {
      std::cout << "Error: Invalid Hitsetkey" << std::endl;
    }
  }

  // copy map into geoContext
  m_tGeometry->geometry().geoContext = transformMap;

  // map is created, now we can use the transforms
  alignmentTransformationContainer::use_alignment = true;
}

// currently used as the transform maker
Acts::Transform3 AlignmentTransformation::newMakeTransform(const Surface& surf, Eigen::Vector3d& millepedeTranslation, Eigen::Vector3d& sensorAngles, bool survey)
{
  // define null matrices
  Eigen::Vector3d nullTranslation(0, 0, 0);
  Eigen::AngleAxisd a(0, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd b(0, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd g(0, Eigen::Vector3d::UnitZ());
  Eigen::Quaternion<double> qnull = g * b * a;
  Eigen::Matrix3d nullRotation = qnull.matrix();

  // get the acts transform components
  // Note that Acts transforms local coordinates of (x,z,y) to global (x,y,z)
  Acts::Transform3 actsTransform = surf->transform(m_tGeometry->geometry().getGeoContext());
  Eigen::Matrix3d actsRotationPart = actsTransform.rotation();
  Eigen::Vector3d actsTranslationPart = actsTransform.translation();

  // Create  alignment local coordinates rotation matrix
  Eigen::AngleAxisd alpha(sensorAngles(0), Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd beta(sensorAngles(1), Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd gamma(sensorAngles(2), Eigen::Vector3d::UnitZ());
  Eigen::Quaternion<double> q = gamma * beta * alpha;
  Eigen::Matrix3d millepedeRotation = q.matrix();

  // and make affine matrices from each

  Acts::Transform3 mpRotationAffine;
  mpRotationAffine.linear() = millepedeRotation;
  mpRotationAffine.translation() = nullTranslation;

  Acts::Transform3 mpTranslationAffine;
  mpTranslationAffine.linear() = nullRotation;
  mpTranslationAffine.translation() = millepedeTranslation;
  
  Acts::Transform3 actsRotationAffine;
  actsRotationAffine.linear() = actsRotationPart;
  actsRotationAffine.translation() = nullTranslation;
  Acts::Transform3 actsTranslationAffine;
  actsTranslationAffine.linear() = nullRotation;
  actsTranslationAffine.translation() = actsTranslationPart;

  // Put them together into a combined transform
  Acts::Transform3 transform;
  //! If we read the survey parameters directly, that is the full transform
  if (survey)
  {
    //! The millepede affines will just be what was read in, which was the
    //! survey information. This should (in principle) be equivalent to
    //! the ideal position + any misalignment
    transform = mpTranslationAffine * mpRotationAffine;
  }
  else   
  {
    transform = mpTranslationAffine * actsTranslationAffine * mpRotationAffine * actsRotationAffine;
  }

  if (localVerbosity)
  {
    Acts::Transform3 actstransform = actsTranslationAffine * actsRotationAffine;

    std::cout << "newMakeTransform" << std::endl;
    std::cout << "Input translation: " << std::endl << millepedeTranslation << std::endl;
    std::cout << "Input sensorAngles: " << std::endl << sensorAngles << std::endl;
    std::cout << "mpRotationAffine: " << std::endl
              << mpRotationAffine.matrix() << std::endl;
      std::cout << "millepederotation * acts " << std::endl
              << millepedeRotation * actsRotationPart << std::endl;
    std::cout << "actsRotationAffine: " << std::endl
              << actsRotationAffine.matrix() << std::endl;
    std::cout << "actsTranslationAffine: " << std::endl
              << actsTranslationAffine.matrix() << std::endl;
    std::cout << "full acts transform " << std::endl
              << actstransform.matrix() << std::endl;
    std::cout << "mpTranslationAffine: " << std::endl
	      << mpTranslationAffine.matrix() << std::endl;
    std::cout << "Overall transform: " << std::endl
              << transform.matrix() << std::endl;
    std::cout << "overall * idealinv " << std::endl
              << (transform * actstransform.inverse()).matrix() << std::endl;
    std::cout << "overall - ideal " << std::endl;
    for (int test = 0; test < transform.matrix().rows(); test++)
    {
      for (int test2 = 0; test2 < transform.matrix().cols(); test2++)
      {
        std::cout << transform(test, test2) - actstransform(test, test2) << ", ";
      }
      std::cout << std::endl;
    }
  }

  return transform;
}

int AlignmentTransformation::getNodes(PHCompositeNode* topNode)
{
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return 0;
}

void AlignmentTransformation::misalignmentFactor(uint8_t layer, const double factor)
{
  transformMap->setMisalignmentFactor(layer, factor);
}
void AlignmentTransformation::createAlignmentTransformContainer(PHCompositeNode* topNode)
{
  //​ Get a pointer to the top of the node tree
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in AlignmentTransformation::createNodes");
  }

  transformMap = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainer");
  if (!transformMap)
  {
    transformMap = new alignmentTransformationContainer;
    auto node = new PHDataNode<alignmentTransformationContainer>(transformMap, "alignmentTransformationContainer");
    dstNode->addNode(node);
  }

  transformMapTransient = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainerTransient");
  if (!transformMapTransient)
  {
    transformMapTransient = new alignmentTransformationContainer;
    auto node = new PHDataNode<alignmentTransformationContainer>(transformMapTransient, "alignmentTransformationContainerTransient");
    dstNode->addNode(node);
  }
}

void AlignmentTransformation::generateRandomPerturbations(Eigen::Vector3d angleDev, Eigen::Vector3d transformDev)
{
  /*Creates random perturbations for the correctional parameters with a given standard deviation and mean of zero*/

  std::cout << "Generating Random Perturbations..." << std::endl;

  if (angleDev(0) != 0)
  {
    std::normal_distribution<double> distribution(0, angleDev(0));
    perturbationAngles(0) = distribution(generator);
  }
  if (angleDev(1) != 0)
  {
    std::normal_distribution<double> distribution(0, angleDev(1));
    perturbationAngles(1) = distribution(generator);
  }
  if (angleDev(2) != 0)
  {
    std::normal_distribution<double> distribution(0, angleDev(2));
    perturbationAngles(2) = distribution(generator);
  }
  if (transformDev(0) != 0)
  {
    std::normal_distribution<double> distribution(0, transformDev(0));
    perturbationTranslation(0) = distribution(generator);
  }
  if (transformDev(1) != 0)
  {
    std::normal_distribution<double> distribution(0, transformDev(1));
    perturbationTranslation(1) = distribution(generator);
  }
  if (transformDev(2) != 0)
  {
    std::normal_distribution<double> distribution(0, transformDev(2));
    perturbationTranslation(2) = distribution(generator);
  }
  if (localVerbosity)
  {
    std::cout << "randomperturbationAngles" << perturbationAngles << " randomperturbationTrans:" << perturbationTranslation << std::endl;
  }
}
