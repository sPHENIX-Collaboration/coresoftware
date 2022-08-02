#include "AlignmentTransformation.h"
#include <g4mvtx/PHG4MvtxDefs.h>
// #include <mvtx/CylinderGeom_Mvtx.h>
#include <trackbase/MvtxDefs.h>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>



void AlignmentTransformation::createMap()
{ 
  std::cout << "Entering createMap..." << std::endl;

  // Define Parsing Variables
  TrkrDefs::hitsetkey hitsetkey = 0;
  float alpha = 0.0, beta = 0.0, gamma = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;
  
  // load alignment constants file
  std::ifstream datafile("data.txt");

  // to loop through all lines in file ---- make sure this is correct number for loop 
  //int fileLines = std::count(std::istreambuf_iterator<char>(datafile), std::istreambuf_iterator<char>(), '\n');  
  int fileLines = 432; 

  for (int i=0; i<fileLines; i++)
     {
      datafile >> hitsetkey >> alpha >> beta >> gamma >> dx >> dy >> dz;
      std::cout << i << " " <<hitsetkey << " " <<alpha<< " " <<beta<< " " <<gamma<< " " <<dx<< " " <<dy<< " " <<dz << std::endl;

      unsigned int trkrId = TrkrDefs::getTrkrId(hitsetkey); // specify between detectors
      if (trkrId == TrkrDefs::mvtxId)
	{
	  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
	  unsigned int stave = MvtxDefs::getStaveId(hitsetkey);
	  unsigned int chip  = MvtxDefs::getChipId(hitsetkey);
	  	  
	  std::cout << "layer: " << layer << " stave: " << stave << " chip: " << chip << std::endl;

	  // Perturbation translations and angles for stave and sensor
	  Eigen::Vector3f staveTranslation (0.0,0.0,0.0);
	  Eigen::Vector3f staveAngles (0.0,0.0,0.0);
	  Eigen::Vector3f sensorTranslation (0.0,0.0,0.0);
	  Eigen::Vector3f sensorAngles (0.0,0.0,0.0);  

	  // now use these transform method for mvtx
	  Eigen::Matrix4f transform = makeTransform(layer, stave, chip, staveTranslation, staveAngles, sensorTranslation, sensorAngles);
	  
	  Eigen::Vector4f test (0.5,0.0,0.6,1);
	  Eigen::Vector4f firstPoint = transform*test;
	  std::cout << firstPoint << std::endl;

	  Eigen::Vector4f finalPoint = transform.inverse()*firstPoint;
	  std::cout << finalPoint << std::endl;
	   
	  transformMap.insert(std::make_pair(hitsetkey,transform));
	  
	  std::cout << "hitsetkey: " << hitsetkey << std::endl; 
	  std::cout << "transform: "<< transform << std::endl;

	}
    }
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




// ask about making methods for mvtx, intt, and tpc
Eigen::Matrix4f AlignmentTransformation::makeTransform(int layer, int stave, int chip, Eigen::Vector3f staveTrans, Eigen::Vector3f staveAngles, Eigen::Vector3f sensorTrans, Eigen::Vector3f sensorAngles, bool mvtx = true)
{
  if(mvtx==true)
    {
      double loc_sensor_in_chip_data[3]            = {0.058128, -0.0005, 0.0};   // mvtx_stave_v1.gdml
      double inner_loc_chip_in_module_data[9][3]   = {{0.0275, -0.02075, -12.060},{0.0275, -0.02075, -9.0450},{0.0275, -0.02075, -6.0300},{0.0275, -0.02075, -3.0150}, 
						   {0.0275,-0.02075,0.0},{0.0275,-0.02075,3.0150},{0.0275,-0.02075,6.0300},{0.0275,-0.02075, 9.0450},{0.0275, -0.02075, 12.060}};
      double inner_loc_halfstave_in_stave_data[3]  = {-0.0275, 0.01825, 0.0};  //only one module in stave
  
      Eigen::Vector3f moduleInStave (inner_loc_halfstave_in_stave_data[0],inner_loc_halfstave_in_stave_data[1],inner_loc_halfstave_in_stave_data[2]);
      Eigen::Vector3f sensorInChip (loc_sensor_in_chip_data[0], loc_sensor_in_chip_data[1], loc_sensor_in_chip_data[2]);
      Eigen::Vector3f chipInModule (inner_loc_chip_in_module_data[chip][0],inner_loc_chip_in_module_data[chip][1],inner_loc_chip_in_module_data[chip][2]);

      Eigen::Vector3f implicitTrans = moduleInStave + sensorInChip + chipInModule;  // sensor translation implicit to mvtx
  
      int staveNum = PHG4MvtxDefs::mvtxdat[layer][5];           // Number of staves per layer
      float rMid   = PHG4MvtxDefs::mvtxdat[layer][1]*0.1;       // convert to cm
      float phi_0  = PHG4MvtxDefs::mvtxdat[layer][4];
      float phi    = 2*M_PI * stave/staveNum + phi_0;
      float xRmid  = rMid * cos(phi);  
      float yRmid  = rMid * sin(phi);
    }
  else if(mvtx==false)
    {

      // sensor center
      const double phi = offsetphi + dphi * icopy;  // if offsetphi is zero we start at zero

      double radius;
      // Make each layer at a single radius - i.e. what was formerly a sub-layer is now considered a layer
      //m_SensorRadius[inttlayer] = params1->get_double_param("sensor_radius") * cm;

      radius = m_SensorRadius[inttlayer];
      radius += sensor_offset_x_ladder;

      double p = 0.0;
      if (laddertype == PHG4InttDefs::SEGMENTATION_Z)
        {
          // The Z sensitive ladders have the sensors offset in y relative to the ladder center
          // We have to slightly rotate the ladder in its own frame to make the radial vector to the sensor center normal to the sensor face
          p = atan(sensor_offset_y / radius);
          // then we adjust the distance to the center of the ladder to put the sensor at the requested distance from the center of the barrel
          radius /= cos(p);
        }

      //    The sensors have no tilt in the new design
      //    The type 1 ladders have the sensor at the center of the ladder in phi, so that is easy
      //    The type 0 ladders are more complicated because the sensor center is perpendicular to the radial vector and the sensor is not at the ladder center
      //         We made the stave box symmetric in y around the sensor center to simplify things


      double sensor_offset_x_ladder = 0.0 - TVSi_x;  // ladder center is at x = 0.0 by construction. Sensor is at lower x, so TVSi_x is negative
      

      //still confused by icopy and p
      // these describe the center of the ladder volume, placing it so that the center of the sensor is at phi = dphi * icopy, and at the correct radius
      const double posx = radius * cos(phi - p);
      const double posy = radius * sin(phi - p);

      const double xRmid = radius * cos(phi - p);
      const double yRmid = radius * sin(phi - p);


      // const double fRotate = p + (phi - p) + offsetrot; // rotate in its own frame to make sensor perp to radial vector (p), then additionally rotate to account for ladder phi
      
      if (itype != 0)
        {
          // We have added the outer sensor above, now we add the HDI extension tab to the end of the outer sensor HDI
          const double posz_ext = (hdi_z_arr[inttlayer][0] + hdi_z) + hdiext_z / 2.;
	  const double m_posZ_new  = m_posZ + posz_ext
        }


      Eigen::Vector3f implicitTrans; //sensor translation implicit to intt 
    }
  else
    {
      std::cout << "Boolean value mvtx must be either true or false";
      Eigen::Matrix4f failedRun; 
      return failedRun;
    }


  Eigen::Vector3f staveTranslation (staveTrans(0) + xRmid, staveTrans(1) + yRmid, staveTrans(2));
  Eigen::Vector3f staveRot (staveAngles(0), staveAngles(1), staveAngles(2) + phi);
  Eigen::Matrix4f staveTransform = AlignmentTransformation::makeAffineMatrix(staveRot, staveTranslation);

  Eigen::Vector3f sensorTranslation = sensorTrans + implicitTrans;
  Eigen::Matrix4f sensorTransform   = AlignmentTransformation::makeAffineMatrix(sensorAngles, sensorTranslation);
  Eigen::Matrix4f combinedTransform = staveTransform*sensorTransform;

  return combinedTransform;
}


//   // Number of sensors layers and staves will most likely change 

//   //look at jon4 placement code

//   //where is sensor in ladder/stave
//   //take into account sensor type

//   //there are four layers
//   //two types of each sensor inner and outer

//   //Transform strip to stave coordiantes
//   // then stave to layer (right radiua and angle)


//   //ladder placement is rotation angle then radisu or translation
  
 
