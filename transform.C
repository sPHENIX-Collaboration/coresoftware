//#include <cmath>
//#include <iostream>

#include <Math/Vector3D.h>
#include <Math/GenVector/DisplacementVector3D.h>
#include <Math/VectorUtil.h>
#include <Math/RotationZYX.h>
#include <Math/Transform3D.h>

//#include <string>
//#include <bits/stdc++.h> 
#include <Math/Point3Dfwd.h>
//#include <Math/PositionVector3D.h>
//#include <Math/Point3D.h>



void transform()
{
  //sensor in chip +chipin module +module in half stave +half stave in stave = sensor in stave
  double loc_sensor_in_chip_data[3]          = {0.058128, -0.0005, 0.0};  // mvtx_stave_v1.gdml
  double inner_loc_chip_in_module_data[9][3] = {{0.0275, -0.02075, -12.060},{0.0275, -0.02075, -9.0450},{0.0275, -0.02075, -6.0300},{0.0275, -0.02075, -3.0150}, 
                                             {0.0275, -0.02075, 0.0},{0.0275, -0.02075, 3.0150},{0.0275, -0.02075, 6.0300},{0.0275, -0.02075, 9.0450},{0.0275, -0.02075, 12.060}};
  
  // double inner_loc_module_in_halfstave_data[3] = {0.0, 0.0, 0.0};       //unnecessary
  double inner_loc_halfstave_in_stave_data[3]  = {-0.0275, 0.01825, 0.0};  //only one module in stave
  
  ROOT::Math::XYZVector moduleInStave (inner_loc_halfstave_in_stave_data[0],inner_loc_halfstave_in_stave_data[1],inner_loc_halfstave_in_stave_data[2]);
  ROOT::Math::XYZVector sensorInChip (loc_sensor_in_chip_data[0], loc_sensor_in_chip_data[1], loc_sensor_in_chip_data[2]);

  double sensorAngles[9][3]       = {{30.0,60.0,90.0},{30.0,60.0,90.0},{30.0,60.0,90.0},{30.0,60.0,90.0},{30.0,60.0,90.0}
			             ,{30.0,60.0,90.0},{30.0,60.0,90.0},{30.0,60.0,90.0},{30.0,60.0,90.0}}; // angles from Millepede for sensor alignment
  double sensorTranslations[9][3] = {{0.1,0.1,0.1},{0.1,0.1,0.1},{0.1,0.1,0.1},{0.1,0.1,0.1},{0.1,0.1,0.1},{0.1,0.1,0.1},{0.1,0.1,0.1},{0.1,0.1,0.1},{0.1,0.1,0.1}};
  double staveTranslation[3][3]   = {{0.1,0.1,0.1},{0.2,0.2,0.2},{0.3,0.3,0.3}};
  double staveRotation[3][3]      = {{30.0,60.0,90.0},{30.0,60.0,90.0},{30.0,60.0,90.0}};

  for(int j=0;j<=2;j++) // loop over staves
  {

    ROOT::Math::XYZVector staveTrans (staveTranslation[j][0],staveTranslation[j][1],staveTranslation[j][2]);
    ROOT::Math::RotationZYX staveRot (staveRotation[j][0]*M_PI/180.0,staveRotation[j][1]*M_PI/180.0,staveRotation[j][2]*M_PI/180.0);
    ROOT::Math::Transform3D staveTransform (staveRot,staveTrans);
    
    std::cout << staveTransform << std::endl;

    for(int i=0;i<=8;i++)
      {
	ROOT::Math::XYZVector chipInModule (inner_loc_chip_in_module_data[i][0],inner_loc_chip_in_module_data[i][1],inner_loc_chip_in_module_data[i][2]);
	//std::cout << chipInModule << std::endl;

	ROOT::Math::XYZVector sensorInStaveVec = sensorInChip + chipInModule + moduleInStave;
	ROOT::Math::XYZPoint sensorInStave (sensorInStaveVec.X(),sensorInStaveVec.Y(),sensorInStaveVec.Z());
    
	ROOT::Math::XYZVector translation (sensorTranslations[i][0], sensorTranslations[i][1], sensorTranslations[i][2]);
	ROOT::Math::RotationZYX rotation (sensorAngles[i][0]*M_PI/180.0, sensorAngles[i][1]*M_PI/180.0, sensorAngles[i][2]*M_PI/180.0);

	ROOT::Math::Transform3D sensorTransform (rotation,translation);
	//std::cout << transform << std::endl;
	ROOT::Math::XYZPoint transformedSensorInStave = sensorTransform*sensorInStave;
	//std::cout << transformedSensorInStave << std::endl;
	
	ROOT::Math::XYZPoint transformedStave = staveTransform*transformedSensorInStave;
	// std::cout << transformedStave << std::endl;


	staveTransform.Invert();
	sensorTransform.Invert();

	ROOT::Math::XYZPoint test = staveTransform*transformedStave;
	ROOT::Math::XYZPoint test2 = sensorTransform*test;
	std::cout << test2 << std::endl; 

              

      }
  }




}
