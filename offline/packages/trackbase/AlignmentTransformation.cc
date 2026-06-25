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

#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4detectors/PHG4TpcGeom.h>

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

  innerLayer[0] = 7;
  innerLayer[1] = 23;
  innerLayer[2] = 39;

  double m_moduleStepPhi = 2.0 * M_PI / 12.0;
  double m_modulePhiStart = -M_PI;
  for (int iside : {0, 1})
  {
    for (int isector = 0; isector < 12; ++isector)
    {
      sectorPhi[iside][isector] = m_modulePhiStart + m_moduleStepPhi * (double) isector;
    }
  }

  extractModuleCenterPositions();  // needed for TPC module transforms

  // Define Parsing Variables
  TrkrDefs::hitsetkey hitsetkey = 0;
  float alpha = 0.0;
  float beta = 0.0;
  float gamma = 0.0;
  float dx = 0.0;
  float dy = 0.0;
  float dz = 0.0;
  float dgrx = 0.0;
  float dgry = 0.0;
  float dgrz = 0.0;

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

  ActsSurfaceMaps surfMaps = m_tGeometry->maps();
  Surface surf;

  int linecount = 0;
  std::string str;
  while (std::getline(datafile, str))
  {
    // trim leading space characters
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch)
                                        { return !std::isspace(ch); }));

    // skip empty lines, or commented lines
    if (str.empty())
    {
      continue;
    }
    if (str.substr(0, 2) == "//")
    {
      continue;
    }
    if (str.substr(0, 1) == "#")
    {
      continue;
    }

    // try read
    std::stringstream ss(str);

    // check to see how many parameters per line in the file
    // If it is old, there may be only six. In that case, set the global rotation pars to zero, print a message.
    std::string dummy;
    int count = 0;
    while (ss >> dummy)
    {
      count++;
    }
    if (count < 9)
    {
      std::stringstream str6(str);
      str6 >> hitsetkey >> alpha >> beta >> gamma >> dx >> dy >> dz;
      if (str6.rdstate() & std::ios::failbit)
      {
        std::cout << "AlignmentTransformation::createMap - invalid line: " << str << " -------- Exiting" << std::endl;
        exit(1);
      }
      dgrx = 0;
      dgry = 0;
      dgrz = 0;

      if (linecount == 1 && localVerbosity > 0)
      {
        std::cout << PHWHERE << "The  alignment parameters file has only 6 parameters" << std::endl
                  << "     --- setting global rotation parameters to zero!" << std::endl;
      }
    }
    else
    {
      std::stringstream str9(str);
      str9 >> hitsetkey >> alpha >> beta >> gamma >> dx >> dy >> dz >> dgrx >> dgry >> dgrz;
      if (str9.rdstate() & std::ios::failbit)
      {
        std::cout << "AlignmentTransformation::createMap - invalid line: " << str << " -------- Exiting" << std::endl;
        exit(1);
      }
    }

    linecount++;

    if(localVerbosity > 0)
      {
	std::cout  <<  hitsetkey << "  " << alpha  << "  " << beta  << "  " << gamma  << "  " << dx  << "  " << dy << "  " << dz
		   << "  " << dgrx << "  " << dgry << "  " << dgrz << std::endl;
      }
    
    // Perturbation translations and angles for stave and sensor
    Eigen::Vector3d sensorAngles(alpha, beta, gamma);
    Eigen::Vector3d millepedeTranslation(dx, dy, dz);
    Eigen::Vector3d sensorAnglesGlobal(dgrx, dgry, dgrz);

    perturbationAngles = Eigen::Vector3d(0.0, 0.0, 0.0);
    perturbationAnglesGlobal = Eigen::Vector3d(0.0, 0.0, 0.0);
    perturbationTranslation = Eigen::Vector3d(0.0, 0.0, 0.0);

    unsigned int trkrId = TrkrDefs::getTrkrId(hitsetkey);  // specify between detectors

    switch (trkrId)
    {
    case TrkrDefs::mvtxId:
    {
      if (perturbMVTX)
      {
        generateRandomPerturbations(mvtxAngleDev, mvtxTransDev);
        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }

      surf = surfMaps.getSiliconSurface(hitsetkey);

      Eigen::Vector3d localFrameTranslation(0.0, 0.0, 0.0);

      Acts::Transform3 transform;
      transform = newMakeTransform(surf, millepedeTranslation, sensorAngles, localFrameTranslation, sensorAnglesGlobal, trkrId, false);

      Acts::GeometryIdentifier id = surf->geometryId();

      if (localVerbosity > 1)
      {
        std::cout << " Add transform for MVTX with surface GeometryIdentifier " << id << " trkrid " << trkrId << std::endl;
        std::cout << " final mvtx transform:" << std::endl
                  << transform.matrix() << std::endl;
      }
      transformMap->addTransform(id, transform);
      transformMapTransient->addTransform(id, transform);

      break;
    }

    case TrkrDefs::inttId:
    {
      if (perturbINTT)
      {
        generateRandomPerturbations(inttAngleDev, inttTransDev);
        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }

      surf = surfMaps.getSiliconSurface(hitsetkey);

      Eigen::Vector3d localFrameTranslation(0.0, 0.0, 0.0);

      Acts::Transform3 transform;
      transform = newMakeTransform(surf, millepedeTranslation, sensorAngles, localFrameTranslation, sensorAnglesGlobal, trkrId, use_intt_survey_geometry);
      Acts::GeometryIdentifier id = surf->geometryId();

      if (localVerbosity > 1)
      {
        std::cout << " Add transform for INTT with surface GeometryIdentifier " << id << " trkrid " << trkrId << std::endl;
      }

      transformMap->addTransform(id, transform);
      transformMapTransient->addTransform(id, transform);
      break;
    }

    case TrkrDefs::tpcId:
    {
      if (perturbTPC)
      {
        generateRandomPerturbations(tpcAngleDev, tpcTransDev);
        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }

      unsigned int nlayers = 1;
      unsigned int test_layer = TrkrDefs::getLayer(hitsetkey);
      unsigned int layer_begin = test_layer;
      if (test_layer < 3)
      {
        // This is a TPC module hitsetkey ("test_layer" will be 0, 1, 2)
        nlayers = 16;
        layer_begin = innerLayer[test_layer];
      }

      unsigned int side = TpcDefs::getSide(hitsetkey);
      unsigned int sector = TpcDefs::getSectorId(hitsetkey);
      //      std::cout << "New module hitsetkey " << hitsetkey << "test_layer " << test_layer <<  " side " << side << " sector " << sector << " nlayers " << nlayers << " layer_begin " << layer_begin << std::endl;

      // loop over layers in module
      for (unsigned int this_layer = layer_begin; this_layer < layer_begin + nlayers; ++this_layer)
      {
        TrkrDefs::hitsetkey this_hitsetkey = TpcDefs::genHitSetKey(this_layer, sector, side);

	//   std::cout << " *** module hitsetkey " << hitsetkey << " this_hitsetkey " << this_hitsetkey << " this layer " << this_layer << " side " << side << " sector " << sector << std::endl;

	// Each TPC hitsetkey has 12 fake surfaces associated with it
	// We want to make a transform for every fake surface in this hitsetkey
	// Loop over the sector phi angles for the fake surfaces and get each surface
	auto layergeom = m_tpccellgeo->GetLayerCellGeom((int) this_layer);
        auto sec_min_phi = layergeom->get_sector_min_phi();
	auto min_phi = sec_min_phi[side][sector];
        auto sec_max_phi = layergeom->get_sector_max_phi();
	auto max_phi = sec_max_phi[side][sector];
	double dphi  = (max_phi - min_phi)/12.0;
	for(int is = 0; is < 12; ++is)
	  {
	    double phis = min_phi + is*dphi + dphi/2.0;
	    double radius = layergeom->get_radius();
	    double zcenter = 51.0;
	    if(side == 0)
	      {
		zcenter *= -1;
	      }
	    Acts::Vector3 env_pos(radius*std::cos(phis), radius * std::sin(phis), zcenter);
	    Acts::Vector3 world_pos = m_tGeometry->transformTpcEnvelopeToWorld(env_pos);	    
	    unsigned short  sskey = 999;
	    Surface this_surf = m_tGeometry->get_tpc_surface_from_coords(this_hitsetkey, world_pos, sskey);
	    if(sskey == 999 || !this_surf)
	      {
		std::cout << PHWHERE << "Failed to get surface for layer " << this_layer << " side " << side << " sector " << sector << "  quit!" << std::endl;
		exit(1);
	      }
	    /*
	    std::cout << " layer " << this_layer << " radius " << radius << " phis " << phis << " min_phi " << min_phi << " max_phi " << max_phi
		      << " side " << side << " sector " << sector << " world " << world_pos.x() << "  " << world_pos.y() << "  " << world_pos.z()
		      <<"  world_radius " << sqrt(world_pos.x() * world_pos.x() + world_pos.y() * world_pos.y())
		      << " sskey " << sskey << std::endl;
	    */
	    
	    Eigen::Vector3d localFrameTranslation(0, 0, 0);
	    use_module_tilt = false;
	    if (test_layer < 4 || use_module_tilt_always)
	      {
		// get the local frame translation that puts the local surface center at the tilted position after the local rotations are applied
		unsigned int this_region = (this_layer - 7) / 16;                                           // 0-2
		Eigen::Vector3d this_center = this_surf->center(m_tGeometry->geometry().getGeoContext()) * 0.1;  // mm to cm
		//this_center includes the TPC tilt used in PHG4TpcDetector construction, transform to tpc envelope coords
		Acts::Vector3 this_center_envelope = m_tGeometry->transformTpcWorldToEnvelope(this_center);
		double this_radius = std::sqrt(this_center_envelope[0] * this_center_envelope[0] + this_center_envelope[1] * this_center_envelope[1]);
		float moduleRadius = TpcModuleRadii[side][sector][this_region];                                     // radius of the center of the module in cm
		localFrameTranslation = getTpcLocalFrameTranslation(moduleRadius, this_radius, sensorAngles) * 10;  // cm to mm
		
		// set this flag for later use 
		use_module_tilt = true;
	      }

	    Acts::Transform3 transform;
	    transform = newMakeTransform(this_surf, millepedeTranslation, sensorAngles, localFrameTranslation, sensorAnglesGlobal, trkrId, false);
	    Acts::GeometryIdentifier id = this_surf->geometryId();
	    
	    if (localVerbosity)
	      {
		std::cout << "    Add transform for TPC with surface GeometryIdentifier " << id 
			  << " trkrid " << trkrId << " hitsetkey " << this_hitsetkey << " layer " << this_layer << " sector " << sector
			  << " side " << side << std::endl;
		if(localVerbosity > 1)
		  {
		    Acts::Vector3 center = this_surf->center(m_tGeometry->geometry().getGeoContext()) * 0.1;  // convert to cm
		    std::cout << "Ideal surface center: " << std::endl
			      << center << std::endl;
		    std::cout << "transform matrix: " << std::endl
			      << transform.matrix() << std::endl;
		  }
	      }
	    transformMap->addTransform(id, transform);
	    transformMapTransient->addTransform(id, transform);
	  }
      }
      
      break;
    }
    
    case TrkrDefs::micromegasId:
    {
      if (perturbMM)
      {
        generateRandomPerturbations(mmAngleDev, mmTransDev);

        sensorAngles = sensorAngles + perturbationAngles;
        millepedeTranslation = millepedeTranslation + perturbationTranslation;
      }
      surf = surfMaps.getMMSurface(hitsetkey);

      Eigen::Vector3d localFrameTranslation(0.0, 0.0, 0.0);

      Acts::Transform3 transform;
      transform = newMakeTransform(surf, millepedeTranslation, sensorAngles, localFrameTranslation, sensorAnglesGlobal, trkrId, false);
      Acts::GeometryIdentifier id = surf->geometryId();

      if (localVerbosity > 1)
      {
        std::cout << " Add transform for Micromegas with surface GeometryIdentifier " << id << " trkrid " << trkrId << std::endl;
      }

      transformMap->addTransform(id, transform);
      transformMapTransient->addTransform(id, transform);
      break;
    }

    default:
    {
      std::cout << "AlignmentTransformation::createMap - Invalid Hitsetkey: " << hitsetkey << std::endl;
      break;
    }
    }
  }

  // copy map into geoContext
  Acts::GeometryContext gctx{transformMap};
  m_tGeometry->geometry().geoContext = gctx;

  std::cout << " AlignmentTransformation processed " << linecount << " input lines " << std::endl;

  // map is created, now we can use the transforms
  alignmentTransformationContainer::use_alignment = true;
}

// currently used as the transform maker
Acts::Transform3 AlignmentTransformation::newMakeTransform(const Surface& surf, Eigen::Vector3d& millepedeTranslation, Eigen::Vector3d& sensorAngles, Eigen::Vector3d& localFrameTranslation, Eigen::Vector3d& sensorAnglesGlobal, unsigned int trkrid, bool survey)
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
  auto actsTransform = surf->localToGlobalTransform(m_tGeometry->geometry().getGeoContext());
  Eigen::Matrix3d actsRotationPart = actsTransform.rotation();
  Eigen::Vector3d actsTranslationPart = actsTransform.translation();

  // Create  alignment local coordinates rotation matrix
  // the measurement is in the local XY plane, which becomes the global xz plane
  Eigen::AngleAxisd alpha(sensorAngles(0), Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd beta(sensorAngles(1), Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd gamma(sensorAngles(2), Eigen::Vector3d::UnitZ());
  Eigen::Quaternion<double> q = gamma * beta * alpha;
  Eigen::Matrix3d millepedeRotationLocal = q.matrix();

  // Create alignment global coordinates rotation matrix
  Eigen::AngleAxisd grx(sensorAnglesGlobal(0), Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd gry(sensorAnglesGlobal(1), Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd grz(sensorAnglesGlobal(2), Eigen::Vector3d::UnitZ());
  Eigen::Quaternion<double> gqr = grz * gry * grx;
  Eigen::Matrix3d millepedeRotationGlobal = gqr.matrix();

  // and make affine matrices from each

  // careful! what axes are what in local coordinates?
  Acts::Transform3 mpLocalRotationAffine;
  mpLocalRotationAffine.linear() = millepedeRotationLocal;
  mpLocalRotationAffine.translation() = nullTranslation;

  Acts::Transform3 mpLocalTranslationAffine;
  mpLocalTranslationAffine.linear() = nullRotation;
  mpLocalTranslationAffine.translation() = localFrameTranslation;

  Acts::Transform3 mpGlobalRotationAffine;
  mpGlobalRotationAffine.linear() = millepedeRotationGlobal;
  mpGlobalRotationAffine.translation() = nullTranslation;

  Acts::Transform3 mpGlobalTranslationAffine;
  mpGlobalTranslationAffine.linear() = nullRotation;
  mpGlobalTranslationAffine.translation() = millepedeTranslation;

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
      transform = mpGlobalTranslationAffine * mpGlobalRotationAffine * mpLocalRotationAffine;
    }
  else
    {
      // not survey. this is the normal usage

      if (trkrid == TrkrDefs::tpcId)
	{
	  if(use_module_tilt)
	    {
	      // use module tilt transforms with local rotation followed by local translation 
	      transform = mpGlobalTranslationAffine * mpGlobalRotationAffine * actsTranslationAffine * actsRotationAffine * mpLocalTranslationAffine * mpLocalRotationAffine;
	    }
	  else
	    {
	      // backward compatibility for old alignment params sets
	      transform = mpGlobalTranslationAffine * mpGlobalRotationAffine * actsTranslationAffine * mpLocalRotationAffine * actsRotationAffine;
	    }
	}
      else
	{
	  // silicon and TPOT	  
	  if(use_new_silicon_rotation_order)
	    {
	      // use new transform order for silicon as well as TPC
	      transform = mpGlobalTranslationAffine * mpGlobalRotationAffine * actsTranslationAffine * actsRotationAffine * mpLocalTranslationAffine * mpLocalRotationAffine;
	    }
	  else
	    {
	      // needed for backward compatibility to existing local rotation parmeter sets in silicon
	      transform = mpGlobalTranslationAffine * mpGlobalRotationAffine * actsTranslationAffine * mpLocalRotationAffine * actsRotationAffine;
	    }
	}
    }
  
  if (localVerbosity > 1)
    {
      Acts::Transform3 actstransform = actsTranslationAffine * actsRotationAffine;
      
      std::cout << "newMakeTransform" << std::endl;
      std::cout << "Input sensorAngles: " << std::endl
		<< sensorAngles << std::endl;
      std::cout << "Input sensorAnglesGlobal: " << std::endl
		<< sensorAnglesGlobal << std::endl;
      std::cout << "Input translation: " << std::endl
		<< millepedeTranslation << std::endl;
      std::cout << "mpLocalRotationAffine: " << std::endl
		<< mpLocalRotationAffine.matrix() << std::endl;
      std::cout << "mpLocalTranslationAffine: " << std::endl
		<< mpLocalTranslationAffine.matrix() << std::endl;
      std::cout << "actsRotationAffine: " << std::endl
		<< actsRotationAffine.matrix() << std::endl;
      std::cout << "actsTranslationAffine: " << std::endl
		<< actsTranslationAffine.matrix() << std::endl;
      std::cout << "mpRotationGlobalAffine: " << std::endl
		<< mpGlobalRotationAffine.matrix() << std::endl;
      std::cout << "mpTranslationGlobalAffine: " << std::endl
		<< mpGlobalTranslationAffine.matrix() << std::endl;
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

Eigen::Vector3d AlignmentTransformation::getTpcLocalFrameTranslation(float moduleRadius, float layerRadius, Eigen::Vector3d& localRotation) const
{
  // everything in cm here

  float Rdiff = layerRadius - moduleRadius;

  // alpha local translation around X axis
  float alpha = localRotation(0);
  float dx = 0.0;
  float dy = -Rdiff * std::sin(alpha);
  float dz = -Rdiff * (1 - std::cos(alpha));

  // beta local translation for rotation around Y axis
  float beta = localRotation(1);
  dx += Rdiff * std::sin(beta);
  dy += 0.0;
  dz += -Rdiff * (1 - std::cos(beta));

  // gamma local translation around Z axis
  float gamma = localRotation(2);
  dx += -Rdiff * std::sin(gamma);
  dy += -Rdiff * (1 - std::cos(gamma));
  dz += 0.0;

  if (localVerbosity > 1)
  {
    std::cout << " alpha, beta, gamma " << alpha << "  " << beta << "  " << gamma << " radius " << moduleRadius << " Rdiff " << Rdiff
              << " dx, dy dz " << dx << "  " << dy << "  " << dz << std::endl;
  }

  Eigen::Vector3d localTranslation(dx, dy, dz);

  return localTranslation;
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

  m_tpccellgeo = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!m_tpccellgeo)
  {
    std::cout << PHWHERE << " unable to find DST node TPCGEOMCONTAINER" << std::endl;
    exit(1);
  }
  
  return 0;
}

void AlignmentTransformation::misalignmentFactor(uint8_t layer, const double factor)
{
  transformMap->setMisalignmentFactor(layer, factor);
  transformMapTransient->setMisalignmentFactor(layer, factor);
}
void AlignmentTransformation::createAlignmentTransformContainer(PHCompositeNode* topNode)
{
  // ​ Get a pointer to the top of the node tree
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
    auto* node = new PHDataNode<alignmentTransformationContainer>(transformMap, "alignmentTransformationContainer");
    dstNode->addNode(node);
  }

  transformMapTransient = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainerTransient");
  if (!transformMapTransient)
  {
    transformMapTransient = new alignmentTransformationContainer;
    auto* node = new PHDataNode<alignmentTransformationContainer>(transformMapTransient, "alignmentTransformationContainerTransient");
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
  if (localVerbosity > 1)
  {
    std::cout << "randomperturbationAngles" << perturbationAngles << " randomperturbationTrans:" << perturbationTranslation << std::endl;
  }
}

void AlignmentTransformation::extractModuleCenterPositions()
{
  if (localVerbosity > 1)
  {
    std::cout << "Extracting TPC module center radii:" << std::endl;
  }

  for (int iside = 0; iside < 2; ++iside)
  {
    for (int iregion = 0; iregion < 3; ++iregion)
    {
      unsigned int lin = innerLayer[iregion];
      unsigned int lout = lin + 15;

      for (int isector = 0; isector < 12; ++isector)
      {
        double sectorphi = sectorPhi[iside][isector];

        TrkrDefs::hitsetkey hitsetkey_in = TpcDefs::genHitSetKey(lin, isector, iside);
	double surf_rad_in = extractModuleCenter(hitsetkey_in, sectorphi);

        TrkrDefs::hitsetkey hitsetkey_out = TpcDefs::genHitSetKey(lout, isector, iside);
	double surf_rad_out = extractModuleCenter(hitsetkey_out, sectorphi);
	
        double mod_radius = (surf_rad_in + surf_rad_out) / 2.0;
        TpcModuleRadii[iside][isector][iregion] = mod_radius;

	if (localVerbosity > 1)
	  {
          std::cout << "  hitsetkey_in " << hitsetkey_in << " lin " << lin << " sector " << isector << " side " << iside << " region " << iregion << std::endl;
	  std::cout << "  hitsetkey_out " << hitsetkey_out << " lout " << lout << " sector " << isector << " side " << iside  << " region " << iregion << std::endl;
	  std::cout << "          module radius " << mod_radius << std::endl;
	  }

      }
    }
  }
}

double AlignmentTransformation::extractModuleCenter(TrkrDefs::hitsetkey hitsetkey, double sectorphi)
{
  // We want the module center position from the ideal geometry in the tpc envelope frame

  // the radius and z are not used, only the phi value
  double x = std::cos(sectorphi + 0.01) * 10.0;
  double y = std::sin(sectorphi + 0.01) * 10.0;
  double z = 0.0;

  Acts::Vector3 world(x, y, z);
  TrkrDefs::subsurfkey subsurfkey = 0;

  //   std::cout << "extractModuleCenter: sectorphi " << sectorphi << " world " << world(0) << "  " << world(1) << "  " << world(2) << std::endl;
    
  // Note: the "world" position here is in pre-tilt tpc envelope coordinates, not global coordinates
  // But, get_tpc_surface_from_coords() expects a global position as input, so we convert to world coordinates
  Acts::Vector3 world_envelope = m_tGeometry->transformTpcEnvelopeToWorld(world);
  Surface surface = m_tGeometry->get_tpc_surface_from_coords(hitsetkey, world_envelope, subsurfkey);
  if (!surface)
  {
    std::cout << PHWHERE << "Failed to find surface, quit " << std::endl;
    exit(1);
  }

  Eigen::Vector3d surf_center = surface->center(m_tGeometry->geometry().getGeoContext());
  surf_center /= 10.0;  // convert from mm to cm
  // convert to tpc envelope coords
  Acts::Vector3 surf_center_envelope = m_tGeometry->transformTpcWorldToEnvelope(surf_center);
  double surf_radius = std::sqrt(surf_center_envelope[0] * surf_center_envelope[0] + surf_center_envelope[1] * surf_center_envelope[1]);

  return surf_radius;
}
