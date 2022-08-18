#include "AlignmentTransformation.h"
#include <g4mvtx/PHG4MvtxDefs.h>
// #include <mvtx/CylinderGeom_Mvtx.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>

#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
//#include <fun4all/Fun4AllReturnCodes.h>
//#include <g4intt/G4_Intt.C>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <mvtx/SegmentationAlpide.h>  // for Alpide constants


#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <trackbase/ActsGeometry.h>



void AlignmentTransformation::createMap(PHCompositeNode* topNode)
{ 
  std::cout << "Entering createMap..." << std::endl;

  createNodes(topNode);

  // Define Parsing Variables
  TrkrDefs::hitsetkey hitsetkey = 0;
  float alpha = 0.0, beta = 0.0, gamma = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;
  
  // load alignment constants file
  std::ifstream datafile("data.txt");

  // to loop through all lines in file ---- make sure this is correct number for loop 
  //int fileLines = std::count(std::istreambuf_iterator<char>(datafile), std::istreambuf_iterator<char>(), '\n');  
  // int fileLines = 432; // mvtx lines
  int fileLines = 656;
  for (int i=0; i<fileLines; i++)
     {
      datafile >> hitsetkey >> alpha >> beta >> gamma >> dx >> dy >> dz;
      std::cout << i << " " <<hitsetkey << " " <<alpha<< " " <<beta<< " " <<gamma<< " " <<dx<< " " <<dy<< " " <<dz << std::endl;

      // Perturbation translations and angles for stave and sensor
      Eigen::Vector3f staveTranslation (0.0,0.0,0.0);
      Eigen::Vector3f staveAngles (0.0,0.0,0.0);
      Eigen::Vector3f sensorTranslation (0.0,0.0,0.0);
      Eigen::Vector3f sensorAngles (0.0,0.0,0);  

      Eigen::Matrix4f transform = makeTransform(hitsetkey, staveTranslation, staveAngles, sensorTranslation, sensorAngles);

      Eigen::Matrix3d globalRotation = rotateToGlobal(hitsetkey);
      std::cout << " global Rotation: " << globalRotation<< std::endl;
      //Eigen::Vector4f test (0.5,0.0,0.6,1);
      //Eigen::Vector4f firstPoint = transform*test;
      //std::cout << firstPoint << std::endl;

      //Eigen::Vector4f finalPoint = transform.inverse()*firstPoint;
      //std::cout << finalPoint << std::endl;
	   
      transformMap->addTransform(hitsetkey,transform);
	  
      //std::cout << "hitsetkey: " << hitsetkey << std::endl; 
      std::cout << "transform: "<< transform << std::endl;
           
     }

   transformMap->identify();
} 


Eigen::Matrix3d AlignmentTransformation::rotateToGlobal(TrkrDefs::hitsetkey hitsetkey)
{  

//unsigned int trkrId = TrkrDefs::getTrkrId(hitsetkey); // specify between detectors
  
  auto surfMaps = m_tGeometry->maps();
  Surface surf = surfMaps.getSiliconSurface(hitsetkey);;

  // if(trkrId == TrkrDefs::mvtxId or trkrId == TrkrDefs::inttId)
  //   {
  //     surf = surfMaps.getSiliconSurface(hitsetkey);
  //   }
  //else if(trkrId == TrkrDefs::tpcId)
  //   {
  //     surf = surfMaps.getTpcSurface(hitsetkey);
  //   }
  // else
  //   {
  //     std::cout << "Invalid hitsetkey";
  //     exit(1);
  //   }


  //surface *surf = surfMaps->getSiliconSurface(hitsetkey);

  Eigen::Vector3d ylocal (0,1,0);
  Eigen::Vector3d sensorNormal = surf->normal(m_tGeometry->geometry().geoContext);
  //Eigen::Vector3d sensorCenter = surf->center(m_tGeometry->geoContext());

  sensorNormal = sensorNormal/sensorNormal.norm(); // make unit vector 
  double sinTheta     = ylocal.dot(sensorNormal);
  double cosTheta     = ylocal.cross(sensorNormal).norm();

  Eigen::Vector3d vectorRejection = (sensorNormal - (ylocal.dot(sensorNormal))*sensorNormal)/(sensorNormal - (ylocal.dot(sensorNormal))*sensorNormal).norm();
  Eigen::Vector3d perpVector      =  sensorNormal.cross(ylocal);


  //Initialize and fill matrices
  Eigen::Matrix3d fInverse;
  fInverse.col(0) = ylocal;  // also try F.row(0) = A, then use inverse
  fInverse.col(1) = vectorRejection;
  fInverse.col(2) = perpVector;
  
  Eigen::Matrix3d G;
  G(0,0) =  cosTheta;
  G(0,1) = -sinTheta;
  G(0,2) =  0;
  G(1,0) =  sinTheta;
  G(1,1) =  cosTheta;
  G(1,2) =  0;
  G(2,1) =  0;
  G(2,2) =  0;
  G(2,2) =  1;

  Eigen::Matrix3d globalRotation = fInverse*G*(fInverse.inverse());
 
  //ylocal is A, sensor Normal is b
  return globalRotation;
}



Eigen::Matrix4f AlignmentTransformation::makeAffineMatrix(Eigen::Vector3f rotationAnglesXYZ, Eigen::Vector3f translationVector)
{
  // Creates 4x4 affine matrix given rotation angles about each axis and translationVector 
    
  Eigen::Transform<float, 3, Eigen::Affine> affineMatrix;
  affineMatrix = Eigen::Translation<float, 3>(translationVector);
  affineMatrix.rotate(Eigen::AngleAxis<float>(rotationAnglesXYZ(0), Eigen::Vector3f::UnitX()));
  affineMatrix.rotate(Eigen::AngleAxis<float>(rotationAnglesXYZ(1), Eigen::Vector3f::UnitY()));
  affineMatrix.rotate(Eigen::AngleAxis<float>(rotationAnglesXYZ(2), Eigen::Vector3f::UnitZ()));
  return affineMatrix.matrix();
}


Eigen::Matrix4f AlignmentTransformation::mvtxTransform(TrkrDefs::hitsetkey hitsetkey, Eigen::Vector3f staveTrans, Eigen::Vector3f staveAngles, Eigen::Vector3f sensorTrans, Eigen::Vector3f sensorAngles)
{
  std::cout << "Creating MVTX Transform..." << std::endl;

  double loc_sensor_in_chip_data[3]           = {0.058128, -0.0005, 0.0};   // mvtx_stave_v1.gdml
  double inner_loc_chip_in_module_data[9][3]  = {{0.0275, -0.02075, -12.060},{0.0275, -0.02075, -9.0450},{0.0275, -0.02075, -6.0300},{0.0275, -0.02075, -3.0150}, 
						  {0.0275,-0.02075,0.0},{0.0275,-0.02075,3.0150},{0.0275,-0.02075,6.0300},{0.0275,-0.02075, 9.0450},{0.0275, -0.02075, 12.060}};
  double inner_loc_halfstave_in_stave_data[3] = {-0.0275, 0.01825, 0.0};    //only one module in stave  
  unsigned int layer                          = TrkrDefs::getLayer(hitsetkey);
  unsigned int stave                          = MvtxDefs::getStaveId(hitsetkey);
  unsigned int chip                           = MvtxDefs::getChipId(hitsetkey);


  Eigen::Vector3f moduleInStave (inner_loc_halfstave_in_stave_data[0],inner_loc_halfstave_in_stave_data[1],inner_loc_halfstave_in_stave_data[2]);
  Eigen::Vector3f sensorInChip  (loc_sensor_in_chip_data[0], loc_sensor_in_chip_data[1], loc_sensor_in_chip_data[2]);
  Eigen::Vector3f chipInModule  (inner_loc_chip_in_module_data[chip][0],inner_loc_chip_in_module_data[chip][1],inner_loc_chip_in_module_data[chip][2]); 
  Eigen::Vector3f implicitTrans = moduleInStave + sensorInChip + chipInModule;  // sensor translation implicit to mvtx
    

  int staveNum   = PHG4MvtxDefs::mvtxdat[layer][5];       // Number of staves per layer
  float rMin     = PHG4MvtxDefs::mvtxdat[layer][0]*0.1;   // convert to cm
  float rMid     = PHG4MvtxDefs::mvtxdat[layer][1]*0.1;   // convert to cm
  float rMax     = PHG4MvtxDefs::mvtxdat[layer][2]*0.1;   // convert to cm
  float phi_0    = PHG4MvtxDefs::mvtxdat[layer][4];
  float phi_tilt = std::asin((rMax*rMax - rMin*rMin) / (2*rMid*SegmentationAlpide::SensorSizeRows));
  float phi      = 2*M_PI * stave/staveNum + phi_0;
  float xRmid    = rMid * cos(phi);  
  float yRmid    = rMid * sin(phi);


  // Create Eigen Matrices for sensor to stave and stave to global transformation
  Eigen::Vector3f sensorTranslation = sensorTrans + implicitTrans;
  Eigen::Vector3f sensorRotation    = sensorAngles; 
  Eigen::Matrix4f sensorTransform   = AlignmentTransformation::makeAffineMatrix(sensorRotation, sensorTranslation);
  Eigen::Vector3f staveTranslation  (staveTrans(0) + xRmid, staveTrans(1) + yRmid, staveTrans(2));
  Eigen::Vector3f staveRotation     (staveAngles(0), staveAngles(1), staveAngles(2) + phi + M_PI/2 + phi_tilt);
  Eigen::Matrix4f staveTransform    = AlignmentTransformation::makeAffineMatrix(staveRotation, staveTranslation);
  Eigen::Matrix4f combinedTransform = staveTransform*sensorTransform;


  //if(Verbosity()>-1)
  //{
      // float testradius = sqrt(xRmid*xRmid + yRmid*yRmid);

      // std::cout <<"hitsetkey: " <<hitsetkey<< " phi: "<<phi<< " phi_0: "<< phi_0 << " stave: " << stave << " stavenum: " << staveNum << " rmid: "<< rMid << " sensorTranslation: "<<sensorTranslation<<" xrmid: " <<  xRmid << " yrmid: "<< yRmid <<" testRadius: "<< testradius<<" layer: "<< layer << std::endl;

      // std::cout << " sensorTransform: " << sensorTransform << std::endl;
      // std::cout << "implicit translation: " << implicitTrans << std::endl;
      // std::cout << " staveTransform: " << staveTransform << std::endl; 
      // std::cout << "phiTilt: " << phi_tilt << " sensW: " << SegmentationAlpide::SensorSizeRows <<std::endl;

  // }

  return combinedTransform;
}


Eigen::Matrix4f AlignmentTransformation::inttTransform(TrkrDefs::hitsetkey hitsetkey, Eigen::Vector3f staveTrans, Eigen::Vector3f staveAngles, Eigen::Vector3f sensorTrans, Eigen::Vector3f sensorAngles)
{
  std::cout << "Creating INTT Transform..." << std::endl;

  double cm                = 1;  // check: Must be modified to correct unit transformaiton
  unsigned int layer       = TrkrDefs::getLayer(hitsetkey);
  unsigned int ladderzId   = InttDefs::getLadderZId(hitsetkey);
  unsigned int ladderphiId = InttDefs::getLadderPhiId(hitsetkey);
  std::cout <<"layer: "<<layer<< " ladderz: "<< ladderzId << " ladderphi: " << ladderphiId << " " << std::endl;
  int nladder[4]          = {12, 12, 16, 16};
  double sensor_radius[4] = {7.188 - 36e-4, 7.732 - 36e-4, 9.680 - 36e-4, 10.262 - 36e-4};
  double offsetphi[4]     = {0.0, 0.5 * 2*M_PI / nladder[1], 0.0, 0.5 * 2*M_PI / nladder[3]};
  int inttlayer           = layer - 3;
  double dphi             = 2*M_PI/nladder[inttlayer];
  const double phi        = offsetphi[inttlayer] + dphi * ladderphiId;  // if offsetphi is zero we start at zero, icopy goes up to number of ladders
  std::cout << "phi: " << phi <<std::endl;


  // All variables needed for Tv_Si_x
  double fphx_x                 = 0.032;
  double fphx_glue_x            = 0.005;
  double hdi_kapton_x           = 0.038;
  double hdi_copper_x           = 0.00376;
  const double stave_thickness  = 0.03;                               // stave thickness
  double si_glue_x              = 0.0014;
  const double strip_x          = 0.032;
  double cooler_gap_x           = 0.3 * cm;                           // id of cooling tube in cm
  double cooler_wall            = stave_thickness;                    // outer wall thickness of cooling tube
  double cooler_x               = cooler_gap_x + 2.0 * cooler_wall;   // thickness of the formed sheet, the flat sheet, and the gap b/w the sheets
  double stave_x                = cooler_x;
  const double ladder_x         = stave_x + hdi_kapton_x + hdi_copper_x + fphx_glue_x + fphx_x;
  const double TVstave_x        = ladder_x / 2. - stave_x / 2.;
  const double TVhdi_kapton_x   = TVstave_x - stave_x / 2. - hdi_kapton_x / 2.;
  const double TVhdi_copper_x   = TVhdi_kapton_x - hdi_kapton_x / 2. - hdi_copper_x / 2.;
  const double siactive_x       = strip_x;
  const double TVsi_glue_x      = TVhdi_copper_x - hdi_copper_x / 2. - si_glue_x / 2.;
  const double TVSi_x           = TVsi_glue_x - si_glue_x / 2. - siactive_x / 2.;
  double sensor_offset_x_ladder = 0.0 - TVSi_x;                       // ladder center is at x = 0.0 by construction. Sensor is at lower x, so TVSi_x is negative

  // PROBLEMS
  // Also there are most likely more translations or rotations that are needed
  // to correctly transform double check 
  
  int itype;
  double strip_z;
  int nstrips_z_sensor;
  if(ladderzId %2 == 0)
    {
      itype = 0;
      strip_z = 1.6;
      nstrips_z_sensor = 8;
    }
  else
    {
      itype = 1;
      strip_z = 2.0;
      nstrips_z_sensor = 5;
    }


  double sensor_edge_z = 0.1;
  double hdi_edge_z    = 0;
  double siactive_z    = strip_z * nstrips_z_sensor;
  double sifull_z      = siactive_z + 2.0 * sensor_edge_z * cm;
  double hdi_z         = sifull_z + hdi_edge_z * cm;
  double m_PosZ        = (itype == 0) ? hdi_z / 2. : hdi_z + hdi_z / 2.;  // location of center of ladder in Z

  if(ladderzId == 0) // There are only even ladderzid in ouput make sure this makes sense 
    {
      m_PosZ = m_PosZ*-1;
    }
  else if(ladderzId == 2)
    {
      m_PosZ = m_PosZ*1;
    }


  double radius      = sensor_radius[inttlayer];
  radius            += sensor_offset_x_ladder;
  const double xRmid = radius * cos(phi);
  const double yRmid = radius * sin(phi);


  // incoroporate below translation from line 1024 phg4 intt
  //const double fRotate = phi + offsetrot;  // rotate in its own frame to make sensor perp to radial vector (p), then additionally rotate to account for ladder phi


  Eigen::Vector3f implicitTrans (sensor_offset_x_ladder,0,m_PosZ); // sensor translation implicit to intt
  Eigen::Vector3f sensorTranslation = sensorTrans + implicitTrans;
  Eigen::Matrix4f sensorTransform   = AlignmentTransformation::makeAffineMatrix(sensorAngles, sensorTranslation);
  Eigen::Vector3f staveTranslation (staveTrans(0) + xRmid, staveTrans(1) + yRmid, staveTrans(2));
  Eigen::Vector3f staveRot         (staveAngles(0), staveAngles(1), staveAngles(2) + phi + M_PI/2);
  Eigen::Matrix4f staveTransform    = AlignmentTransformation::makeAffineMatrix(staveRot, staveTranslation);
  Eigen::Matrix4f combinedTransform = staveTransform*sensorTransform;


  std::cout << "m_PosZ:" << m_PosZ << " implicit trans: " << implicitTrans<<std::endl;

  
  return combinedTransform;
}


Eigen::Matrix4f AlignmentTransformation::makeTransform(TrkrDefs::hitsetkey hitsetkey, Eigen::Vector3f staveTrans, Eigen::Vector3f staveAngles, Eigen::Vector3f sensorTrans, Eigen::Vector3f sensorAngles)
{
  unsigned int trkrId = TrkrDefs::getTrkrId(hitsetkey); // specify between detectors

  if(trkrId == TrkrDefs::mvtxId)
    {
      // Eigen::Vector3d nullVector (0,0,0);
      //make affine matrix with millepede sensor rotations ad null transaltion z placement 
      //Eigen::Matrix4d millepedeRotation = AlignmentTransformation::makeAffineMatrix(sensorAngles,nullVector);
      // make affine matrix with global rotations and millepede perturbation + center vector perturbation
      
      
      Eigen::Matrix4f transformation = AlignmentTransformation::mvtxTransform(hitsetkey,staveTrans,staveAngles,sensorTrans,sensorAngles);
      return transformation;
    }
  else if(trkrId == TrkrDefs::inttId)
    {
      Eigen::Matrix4f transformation = AlignmentTransformation::inttTransform(hitsetkey,staveTrans,staveAngles,sensorTrans,sensorAngles);
      return transformation;
    }
  else
    {
      std::cout << "Invalid hitsetkey";
      exit(1);
    }
}


int AlignmentTransformation::createNodes(PHCompositeNode* topNode)
{
  //​ Get a pointer to the top of the node tree
  PHNodeIterator iter(topNode);
 
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in AlignmentTransformation::createNodes");
  }

  transformMap = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainer");
  if(!transformMap)
    {
      transformMap = new alignmentTransformationContainer;
      auto node    = new PHDataNode<alignmentTransformationContainer>(transformMap, "alignmentTransformationContainer");
      dstNode->addNode(node);
    }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }

 return 0; 
}



  // return combinedTransform;
//   // Number of sensors layers and staves will most likely change 
//   //look at jon4 placement code
//   //where is sensor in ladder/stave
//   //take into account sensor type
//   //there are four layers
//   //two types of each sensor inner and outer
//   //Transform strip to stave coordiantes
//   // then stave to layer (right radiua and angle)

//   //ladder placement is rotation angle then radisu or translation  
