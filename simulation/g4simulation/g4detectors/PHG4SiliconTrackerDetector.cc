#include "PHG4SiliconTrackerDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv4.h"

#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4ReflectionFactory.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <cmath>
#include <sstream>

using namespace std;

static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes

PHG4SiliconTrackerDetector::PHG4SiliconTrackerDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr  ):
  PHG4Detector(Node, dnam),
  envelope_inner_radius(26.0*mm),
  envelope_outer_radius(880*mm),
  envelope_z(2300*mm + no_overlap),
  place_in_x(0 * cm),
  place_in_y(0 * cm),
  place_in_z(0 * cm),
  x_rot(0),
  y_rot(0),
  z_rot(0),
  active(0),
  absorberactive(0),
  layer(lyr),
  layer_nominal_radius(800*mm),
  sensor_z(98.0*mm / 2.0),
  sensor_x(0.32*mm / 2.0),
  strip_x(0.320*mm / 2.0),
  strip_y(0.060*mm / 2.00),
  strip_z(8.0*mm / 2.0),
  roc_y(12.0*mm / 2.0),
  roc_cu_x(0.072*mm / 2.0),
  roc_ins_x(0.256*mm / 2.0),
  tab_z(20.0*mm / 2.0),
  tab_y(15.0*mm / 2.0),
  tab_cu_x(0.072*mm / 2.0),  
  tab_ins_x(0.256*mm / 2.0),  
  fpga_z(10.0*mm / 2.0),
  fpga_y(10.0*mm / 2.0),
  fpga_x(0.3*mm / 2.0),
  stave_y(8.0*mm / 2.0),
  stave_inner_y(6.0*mm / 2.0),
  stave_x(4.0*mm / 2.0),
  stave_inner_x(3.2*mm / 2.0),
  chip_z(6.4*mm / 2.0),
  chip_y(9.1*mm / 2.0),
  chip_x(0.3*mm / 2.0),  
  radius_stagger(15.0*mm),
  N_staggers(2),
  maxrap(1.1),
  overlap_fraction(0.15),
  strip_tilt(0),
  sensor_edge_phi(2.92 * mm),
  sensor_edge_z(1.0 * mm),
  N_strips_per_column(1536),
  option_double_layer(false),
  add_lower_roc(true),
  N_sensors_in_layer(0)
{

  //================================
  // The geometry of the components
  //
  // Sensors:
  // depth is 320 microns
  // material is Si (G4_Si)
  // outer envelope is 98 mm x 98 mm
  // inactive edges are: (98-92.16) / 2 = 2.92 mm in phi, (98 - 96) / 2 = 1 mm in z
  // The smallest block has active area 7.68 mm in phi x 8 mm in z  (128 strps x 60 microns x 8 mm long)
  //
  // Large (s2 layer)  - 12 blocks in phi, 12 blocks in z
  //    active: 92.16 mm in phi x 96 mm in z
  //    outer dimensions:  98 mm in phi x 98 mm in z
  //
  // Medium (s1 layer)  - 6 blocks in phi, 12 blocks in z
  //     active: 46.08 mm in phi x 96 mm in z
  //     outer: 48 mm x 98 mm
  //
  // Small (s0 layer): 2 blocks in phi, 12 blocks in z
  //     active: 15.35 mm x 96 mm
  //     outer:  16 mm x 98 mm
  //  
  // ROC card:
  //     12 mm in phi x 96 mm in z
  //     72 microns Cu
  //     256 microns insulator (34 cm L_R)
  //
  // FPGA tab:
  //     15 mm in phi x 20 mm in z
  //     72 microns Cu
  //     256 microns insulator (34 cm L_R)
  //
  // Chip:
  //     9.1 mm in phi x 6.4 mm in z x 300 microns
  //     material is Si (G4_Si)
  //     
  // FPGA:
  //     dimensions - placeholder for now
  //     10 mm in phi x 10 mm in z
  //     material is Si (G4_Si)
  //
  // Readout bus:
  //     dimensions?
  //     materail?
  //     not implemented, Yasuyuki says too little mass to worry about
  // 
  // Stave:
  //     external: 4 mm in phi x 7 mm in z
  //     internal;  3.2 mm in phi x 6 mm in z
  //     stave material: CFC
  //     coolant: Novec
  //============================

  //===============================
  // Novec:
  //    C_4F_9OCH_3
  //    liquid density1400 kg/m^3 at 25 C, 1475 kg/m^3 at 0 C
  //    
  //===============================

  ostringstream name;
  name.str("");
  name << "layer_" << layer;
  layer_string = name.str().c_str();
  if (verbosity > 0) cout << layer_string.c_str() << endl;
}

PHG4SiliconTrackerDetector::~PHG4SiliconTrackerDetector()
{}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4SiliconTrackerDetector::IsInSiliconTracker(G4VPhysicalVolume * volume) const
{
  // Is this volume one of the sensor strips?

  if (volume->GetName().find("sensor_strip") != string::npos)
    {
      // Check to see if this strip is in the layer belonging to this instance of SiliconTrackerDetector
      // layer_string contains "layer_n", where n = the layer number for this instance
      if (volume->GetName().find(layer_string.c_str()) != string::npos)
      return 1;
    }
  
  return 0;
}

void
PHG4SiliconTrackerDetector::Construct( G4LogicalVolume* logicWorld )
{
  // This is called from PHG4PhenixDetector::Construct()

  if(verbosity>0)
    cout << endl << "PHG4SiliconTrackerDetector::Construct called for layer " << layer << endl;

  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  ConstructSiliconTracker(logicWorld);

  // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int
PHG4SiliconTrackerDetector::ConstructSiliconTracker(G4LogicalVolume* trackerenvelope)
{
  if(verbosity>0)
    {
      cout << " PHG4SiliconTrackerDetector::ConstructSiliconTracker:" <<endl;
      cout << " Layer " << layer << endl;
      cout << " add_lower_roc " << add_lower_roc << endl;
      cout << " option_double_layer " << option_double_layer << endl;
      cout << " layer_nominal_radius " << layer_nominal_radius << endl;
      cout << " radius_stagger  " << radius_stagger << endl;
      cout << " N_strips_per_column " << N_strips_per_column << endl;
      cout << " strip_tilt " << strip_tilt << endl;
      cout << endl;
    }


  // Construct the outer layer of the tracker 
  G4double clearance = 0.0001*mm;

  //===================================
  // Make all of the pieces of the ladder section here
  // These are all defined so that they are in the correct relative 
  // location to each other when the sensor is at 0,0,0 
  //===================================

  // Define the sensor envelope box. It should contain the strips, with clearance all around
  // NOTE: box dimensions are twice the argument values!
  // we calculate the sensor y dimension from the strop parameters and the edge dimensions
  // recall that the sensor_y and strip_y are half the length of the side, and there are 2 edges 
  // we want clearance between all strips
  // There are 12 columns of sensors per sensor envelope, 128, or 256, or 512, etc.  strips per column
  strip_z_spacing = strip_z * 2.0  + clearance;  // distance between strip centers (box dimenstions are 1/2 of side)
  strip_y_spacing = strip_y * 2.0  + clearance;
  N_strip_columns = (int) ( (sensor_z)* 2.0 / strip_z_spacing);
  // calculate the y extent of the sensor including clearances between all components, then halve it to get the box y
  sensor_y = (sensor_edge_phi + clearance + (double) N_strips_per_column * strip_y_spacing + sensor_edge_phi) / 2.0;
  G4double sensor_envelope_x = sensor_x + clearance;
  G4double sensor_envelope_y = sensor_y + clearance;
  G4double sensor_envelope_z = (clearance + sensor_edge_z + clearance + (double) N_strip_columns * strip_z_spacing + sensor_edge_z + clearance)  / 2.0;

  if(verbosity>0)  
    {
      cout << "sensor_envelope_y is calculated to be " << sensor_envelope_y 
	   << " from N_strips_per_column " << N_strips_per_column 
	   << " and strip y " << strip_y
	   << " and sensor_edge_phi " << sensor_edge_phi
	   << endl;
      cout   << " sensor_envelope_z is calculated to be " << sensor_envelope_z 
	     << " from N_strip_columns " << N_strip_columns 
	     << " and strip_z " << strip_z
	     << " and sensor_edge_z " << sensor_edge_z
	     << endl;
    }

  G4VSolid *sensor_envelope_box = new G4Box("sensor_envelope_box", sensor_envelope_x, sensor_envelope_y, sensor_envelope_z);

  // Define a strip box
  // strips are 60 microns pitch x 8 mm x 320 microns thick
  // but we make the box a little longer in z so that tilted strips can be truncated to the right length
  G4VSolid *strip_box = new G4Box("strip_box", strip_x, strip_y, strip_z * 1.10);

  // Define the box for the ROC  - this consists of two parts
  // The first part holds the chips (12 (128 channel) readout chips per roc)
  // Each part has layers of Cu and layers of insulator
  // Define a box for each
  roc_z = sensor_z;
  G4VSolid *roc_cu_box = new G4Box("roc_cu_box", roc_cu_x, roc_y, roc_z);
  G4VSolid *roc_ins_box = new G4Box("roc_ins_box", roc_ins_x, roc_y, roc_z);

  // The second part of the ROC is a tab that holds the FPGA
  G4VSolid *tab_cu_box = new G4Box("tab_cu_box", tab_cu_x, tab_y, tab_z);
  G4VSolid *tab_ins_box = new G4Box("tab_ins_box", tab_ins_x, tab_y, tab_z);

  // Define the FPGA that sits on the tab
  G4VSolid *fpga_box = new G4Box("fpga_box", fpga_x, fpga_y, fpga_z);

  // The bus carries power and signals
  // dimensions? material?

  // The stave is a rectangular cross section tube attached to the back of the main part of the ROC
  stave_z = sensor_z;
  G4VSolid *stave_box = new G4Box("stave_box", stave_x, stave_y, stave_z);
  G4VSolid *stave_inner_box = new G4Box("stave_inner_box", stave_inner_x, stave_inner_y, stave_z);

  // The readout chip is a little different
  // There are 12 per ROC
  int N_chips_per_roc = 12;
  G4VSolid *chip_box = new G4Box("chip_box", chip_x, chip_y, chip_z);

  // We have all of the basic volumes defined. Now we want to place them in the barrel
  //============================================================

 ostringstream name;

//=======================================
  // First we create the sensor envelope and populate it
  // with strips once, and the dead edge areas. Then we will 
  // place copies as needed to populate the entire layer.
  //======================================= 

  // Create the sensor envelope volume 
  name.str("");
  name << "sensor_envelope_layer_" << layer;
  G4LogicalVolume *sensor_envelope_volume = new G4LogicalVolume(sensor_envelope_box,G4Material::GetMaterial("G4_AIR"),name.str().c_str(), 0, 0, 0);
  G4VisAttributes* sensorVisAtt = new G4VisAttributes();
  sensorVisAtt->SetVisibility(false);
  sensorVisAtt->SetForceSolid(false);
  sensorVisAtt->SetColour(G4Colour::White());
  sensor_envelope_volume->SetVisAttributes(sensorVisAtt);

  // create the helper assembly volume to contain all of the strips that will be placed in the sensor
  G4AssemblyVolume *sensor_assembly = new G4AssemblyVolume();

  // Populate the (so far unplaced) sensor assembly with placed strips
  //-----------------------------------------------------------------------------------
  // Here is where we need to implement stereo strips
  // Stereo strips have to be rotated around the x axis by the tilt angle
  // Then we have to cut off the part of the strip that is outside the sensor envelope
  // This must be done at the solid stage
  // We use a G4IntersectionSolid, which keeps the part that is inside both solids

  // Define the box that defines the parts of the tilted strips to be kept in a column
  // this will be used to discard the parts that do not belong in the sensor
  G4double column_box_x = sensor_x + clearance;
  // This should be only the sensitive area of the sensor, including clearance between strips
  //G4double column_box_y = (double) N_strips_per_column * (strip_y + clearance) - clearance;
  G4double column_box_y = ( (double) N_strips_per_column * strip_y_spacing -clearance) / 2.0;
  G4double column_box_z = strip_z;
  G4VSolid *column_box = new G4Box("column_box", column_box_x, column_box_y, column_box_z);

  for(int icolumn = 0; icolumn<N_strip_columns; icolumn++)
    {
      // We want to position the strip in the same way that it is done in the geometry object
      // we refer everything to the center of the sensor envelope
      G4double strip_z_offset = 0.0;
      if(N_strip_columns%2)
	strip_z_offset  =  ((double) (icolumn - N_strip_columns/2) ) * strip_z_spacing;
      else
	strip_z_offset  =  ((double) (icolumn - N_strip_columns/2) + 0.5) * strip_z_spacing; 

      for(int istrip = 0; istrip<N_strips_per_column; istrip++)
	{
	  // get the offset of the center of the strip relative to the center of the sensor envelope
	  // remember that the box dimensions are half the side length
	  // move the strip to the end of the sensor envelope, then move it back by edge z and strip z to get it inside the sensor, then add icolumn * strip length in z (i.e. 2.0*strip_z)
	  // for each new column. Add clearance between all components.
	  //G4double strip_z_offset = -sensor_envelope_z + clearance + sensor_edge_z  + clearance + (double) icolumn * strip_z_spacing  + strip_z;
	  // Note that this centers the strips on the middle of the sensor in y
	  //G4double strip_y_offset =  - sensor_envelope_y +  clearance + sensor_edge_phi + clearance + (double) istrip * strip_y_spacing + strip_y;

	  // We want to position the strip in the same way that it is done in the geometry object
	  // we refer everything to the center of the sensor
	  double strip_y_offset = 0.0;
	  if(N_strips_per_column%2)
	    strip_y_offset = (double) (istrip - N_strips_per_column/2) * strip_y_spacing;
	  else
	    strip_y_offset = ((double) (istrip - N_strips_per_column/2) + 0.5) * strip_y_spacing;

	  // assume all strips are tilted, make the cuts here that define the part of the strip that is kept
	  // define the transformation that tilts the strip (Ra) and puts it at the correct y position in the frame of the column box (Ta)
	  // Apart from the z dimension, this is the same transformation as we will use later to place the truncated strip in the sensor envelope
	  G4RotationMatrix Ra1;
	  Ra1.rotateX(strip_tilt);
	  G4ThreeVector Ta1;
	  Ta1.setX(0.0);
	  Ta1.setY(strip_y_offset);
	  Ta1.setZ(0.0);
	  G4Transform3D Tr1(Ra1,Ta1);
	  // now make the cut that removes any part of the strip that is outside the column box 
	  G4VSolid *tilted_strip_box =  new G4IntersectionSolid("tilted_strip_box",column_box,strip_box, Tr1);

	  // Note that the origin and coordinates of the combined solid are the same as those of the column box
	  // Therefore the tilted strips are already placed correctly in y

	  // name the tilted strip and make its logical volume
	  name.str("");
	  name << "sensor_strip_layer_" << layer << "_address_" << icolumn << "_" << istrip;
	  G4LogicalVolume *strip_volume = new G4LogicalVolume(tilted_strip_box,G4Material::GetMaterial("G4_Si"),name.str().c_str(), 0, 0, 0);
	  G4VisAttributes* stripVisAtt = new G4VisAttributes();
	  stripVisAtt->SetVisibility(true);
	  stripVisAtt->SetForceSolid(true);
	  if(istrip % 2)
	    stripVisAtt->SetColour(G4Colour::Green());
	  else
	    stripVisAtt->SetColour(G4Colour::Red());
	  strip_volume->SetVisAttributes(stripVisAtt);

	  // Add the strip to the sensor assembly in the correct y position relative to the sensor cente
	  G4RotationMatrix Ra;
	  G4ThreeVector Ta;
	  Ta.setX(0.0); 
	  Ta.setY(0.0) ; 
	  Ta.setZ(strip_z_offset);
	  G4Transform3D Tr(Ra,Ta);
	  sensor_assembly->AddPlacedVolume(strip_volume, Tr);
	}
    }

  // Now add the inactive edges of the sensor
  //-----------------------------------------------------
  // The phi edges  
  G4VSolid *sensor_phi_edge_box = new G4Box("sensor_phi_edge_box", sensor_x, sensor_edge_phi/2.0, sensor_z);
  G4LogicalVolume *sensor_phi_edge_volume = new G4LogicalVolume(sensor_phi_edge_box,G4Material::GetMaterial("G4_Si"),"sensor_phi_edge_volume", 0, 0, 0);
  G4VisAttributes* edge_phi_VisAtt = new G4VisAttributes();
  edge_phi_VisAtt->SetVisibility(true);
  edge_phi_VisAtt->SetForceSolid(true);
  edge_phi_VisAtt->SetColour(G4Colour::Yellow());
  sensor_phi_edge_volume->SetVisAttributes(edge_phi_VisAtt);

  G4RotationMatrix Ra2;
  G4ThreeVector Ta2;
  Ta2.setX(0.0); 
  Ta2.setY(sensor_envelope_y - clearance - sensor_edge_phi/2.0) ; 
  Ta2.setZ(0.0);
  G4Transform3D Tr2(Ra2,Ta2);
  sensor_assembly->AddPlacedVolume(sensor_phi_edge_volume, Tr2);

  G4RotationMatrix Ra3;
  G4ThreeVector Ta3;
  Ta3.setX(0.0); 
  Ta3.setY(-sensor_envelope_y + clearance + sensor_edge_phi/2.0) ; 
  Ta3.setZ(0.0);
  G4Transform3D Tr3(Ra3,Ta3);
  sensor_assembly->AddPlacedVolume(sensor_phi_edge_volume, Tr3);

  // The z edges
  G4VSolid *sensor_z_edge_box = new G4Box("sensor_z_edge_box", sensor_x, sensor_y - sensor_edge_phi, sensor_edge_z/2.0);
  G4LogicalVolume *sensor_z_edge_volume = new G4LogicalVolume(sensor_z_edge_box,G4Material::GetMaterial("G4_Si"),"sensor_z_edge_volume", 0, 0, 0);
  G4VisAttributes* edge_z_VisAtt = new G4VisAttributes();
  edge_z_VisAtt->SetVisibility(true);
  edge_z_VisAtt->SetForceSolid(true);
  edge_z_VisAtt->SetColour(G4Colour::Yellow());
  sensor_z_edge_volume->SetVisAttributes(edge_z_VisAtt);
  //

  G4RotationMatrix Ra4;
  G4ThreeVector Ta4;
  Ta4.setX(0.0); 
  Ta4.setY(0.0) ; 
  Ta4.setZ(sensor_envelope_z - clearance - sensor_edge_z/2.0);
  G4Transform3D Tr4(Ra4,Ta4);
  sensor_assembly->AddPlacedVolume(sensor_z_edge_volume, Tr4);
  //
  G4RotationMatrix Ra5;
  G4ThreeVector Ta5;
  Ta5.setX(0.0); 
  Ta5.setY(0.0);
  Ta5.setZ(-sensor_envelope_z + clearance + sensor_edge_z/2.0) ; 
  G4Transform3D Tr5(Ra5,Ta5);
  sensor_assembly->AddPlacedVolume(sensor_z_edge_volume, Tr5);

  // now imprint the sensor assembly in the (still unplaced) sensor_envelope
  G4RotationMatrix Ra6;
  Ra6.rotateZ(0.0);
  G4ThreeVector Ta6;
  Ta6.setX(0.0); 
  Ta6.setY(0.0) ; 
  Ta6.setZ(0.0);
  G4Transform3D Tr6(Ra6,Ta6);
  sensor_assembly->MakeImprint(sensor_envelope_volume, Tr6,0,overlapcheck);
  // This takes a very long time to do an overlapcheck for the sensor, use this to avoid it
  //sensor_assembly->MakeImprint(sensor_envelope_volume, *Tr,0,false);

  //==================================
  // done with filling the sensor envelope volume
  //==================================

  //======================================================
  // Now we create a volume to hold the upper ROCS, chips, fpga's and stave
  //======================================================

  // Create a box that will hold the upper or lower ROCs and staves
  // It has to be able to accomodate the stave plus ROCs + chips + fpga on the pper or lower side
  // these are all half the side length
  G4double roc_envelope_z = sensor_envelope_z;   // already includes clearance
  G4double roc_envelope_y = clearance/2.0 + roc_y + clearance/2.0  + tab_y + clearance/2.0;
  G4double roc_envelope_x = clearance/2.0 + chip_x + clearance/2.0 + roc_cu_x + clearance/2.0 + roc_ins_x + clearance/2.0 + stave_x + clearance/2.0;
  if(option_double_layer)
    {
      // this is a second layer, so it  will not have its own stave, we want only the roc and the chip
      roc_envelope_x -= stave_x + clearance/2.0;
    }
  G4VSolid *roc_envelope_box = new G4Box("roc_envelope_box",roc_envelope_x, roc_envelope_y, roc_envelope_z);
  //
  // Create the upper ROC envelope volume 
  name.str("");
  name  << "upper_ROC_envelope_volume_layer_" <<  layer;
  if (verbosity>0) cout << " upper ROC envelope name = " << name.str().c_str() << endl;
  G4LogicalVolume *upper_ROC_envelope_volume = new G4LogicalVolume(roc_envelope_box,G4Material::GetMaterial("G4_AIR"),name.str().c_str(), 0, 0, 0);
 G4VisAttributes* roc_envelope_VisAtt = new G4VisAttributes();
  roc_envelope_VisAtt->SetVisibility(false);
  roc_envelope_VisAtt->SetForceSolid(false);
  roc_envelope_VisAtt->SetForceWireframe(false);
  roc_envelope_VisAtt->SetColour(G4Colour::White());
  upper_ROC_envelope_volume->SetVisAttributes(roc_envelope_VisAtt);

  // create the assembly helper for filling the upper ROC envelope
  G4AssemblyVolume *upper_roc_assembly = new G4AssemblyVolume();

  // Create the actual upper ROC card logical volumes and place them in the assembly helper
  // There are two volumes - the Cu and the insulator
  name.str("");
  name << "upper_roc_cu_logical_volume";
  G4LogicalVolume *single_upper_roc_cu_volume = new G4LogicalVolume(roc_cu_box,G4Material::GetMaterial("G4_Cu"),name.str().c_str(), 0, 0, 0);
  G4VisAttributes* roc_cu_VisAtt = new G4VisAttributes();
  roc_cu_VisAtt->SetVisibility(true);
  roc_cu_VisAtt->SetForceSolid(true);
  roc_cu_VisAtt->SetColour(G4Colour::Yellow());
  single_upper_roc_cu_volume->SetVisAttributes(roc_cu_VisAtt);
  //
  name.str("");
  name << "upper_roc_ins_logical_volume";
  G4LogicalVolume *single_upper_roc_ins_volume = new G4LogicalVolume(roc_ins_box,G4Material::GetMaterial("G4_Si"),name.str().c_str(), 0, 0, 0);
  G4VisAttributes* roc_ins_VisAtt = new G4VisAttributes();
  roc_ins_VisAtt->SetVisibility(true);
  roc_ins_VisAtt->SetForceSolid(true);
  roc_ins_VisAtt->SetColour(G4Colour::Yellow());
  single_upper_roc_ins_volume->SetVisAttributes(roc_ins_VisAtt);
  //
  // The upper roc goes at the bottom of the upper roc envelope in y
  // and in x between the chip and stave
  //
  // first the Cu
  G4RotationMatrix Ra7; 
  G4ThreeVector Ta7;
  G4double xloc = -roc_envelope_x + clearance  + chip_x * 2.0 + clearance + roc_cu_x;
  Ta7.setX(xloc);
  Ta7.setZ(0.);
  Ta7.setY( -roc_envelope_y  + clearance + roc_y);
  G4Transform3D Tr7(Ra7,Ta7);
  upper_roc_assembly->AddPlacedVolume(single_upper_roc_cu_volume,Tr7);
  //
  // now the insulator
  G4RotationMatrix Ra8; 
  G4ThreeVector Ta8;
  xloc = -roc_envelope_x + clearance  + chip_x * 2.0 + clearance + roc_cu_x * 2.0 + clearance + roc_ins_x;
  Ta8.setX(xloc);
  Ta8.setZ(0.);
  Ta8.setY( -roc_envelope_y  + clearance + roc_y);
  G4Transform3D Tr8(Ra8,Ta8);
  upper_roc_assembly->AddPlacedVolume(single_upper_roc_ins_volume,Tr8);

  // The upper tab volume goes at the top of the upper roc envelope
  // this also has a Cu and an insulator box
  //
  // first the Cu
  name.str("");
  name << "upper_tab_cu_logical_volume";
  G4LogicalVolume *single_upper_tab_cu_volume = new G4LogicalVolume(tab_cu_box,G4Material::GetMaterial("G4_Si"),name.str().c_str(), 0, 0, 0);
  G4VisAttributes* tab_cu_VisAtt = new G4VisAttributes();
  tab_cu_VisAtt->SetVisibility(true);
  tab_cu_VisAtt->SetForceSolid(true);
  tab_cu_VisAtt->SetColour(G4Colour::Yellow());
  single_upper_tab_cu_volume->SetVisAttributes(tab_cu_VisAtt);
  G4RotationMatrix Ra9; 
  G4ThreeVector Ta9;
  xloc = -roc_envelope_x + clearance + fpga_x * 2.0 + clearance + tab_cu_x;
  Ta9.setX(xloc); 
  Ta9.setY( +roc_envelope_y - clearance - tab_y);
  Ta9.setZ( 0.);
  G4Transform3D Tr9(Ra9,Ta9);
  upper_roc_assembly->AddPlacedVolume(single_upper_tab_cu_volume,Tr9);
  //
  // Then the insulator
  name.str("");
  name << "upper_tab_ins_logical_volume";
  G4LogicalVolume *single_upper_tab_ins_volume = new G4LogicalVolume(tab_ins_box,G4Material::GetMaterial("G4_Si"),name.str().c_str(), 0, 0, 0);
  G4VisAttributes* tab_ins_VisAtt = new G4VisAttributes();
  tab_ins_VisAtt->SetVisibility(true);
  tab_ins_VisAtt->SetForceSolid(true);
  tab_ins_VisAtt->SetColour(G4Colour::Yellow());
  single_upper_tab_ins_volume->SetVisAttributes(tab_ins_VisAtt);
  G4RotationMatrix Ra10; 
  G4ThreeVector Ta10;
  xloc = -roc_envelope_x + clearance + fpga_x * 2.0 + clearance + tab_cu_x * 2.0 + clearance + tab_ins_x;
  Ta10.setX(xloc); 
  Ta10.setY( +roc_envelope_y - clearance - tab_y);
  Ta10.setZ( 0.);
  G4Transform3D Tr10(Ra10,Ta10);
  upper_roc_assembly->AddPlacedVolume(single_upper_tab_ins_volume,Tr10);
  
  // The upper fpga goes in the middle of the upper tab in y
  // and in front of the tab in x
  name.str("");
  name << "upper_fpga_logical_volume";
  G4LogicalVolume *single_upper_fpga_volume = new G4LogicalVolume(fpga_box,G4Material::GetMaterial("G4_Si"),name.str().c_str(), 0, 0, 0);
  G4VisAttributes* fpgaVisAtt = new G4VisAttributes();
  fpgaVisAtt->SetVisibility(true);
  fpgaVisAtt->SetForceSolid(true);
  fpgaVisAtt->SetColour(G4Colour::Red());
  single_upper_fpga_volume->SetVisAttributes(fpgaVisAtt);
  G4RotationMatrix Ra11; 
  G4ThreeVector Ta11;
  xloc = -roc_envelope_x  + clearance + fpga_x; 
  Ta11.setX(xloc); 
  Ta11.setY(roc_envelope_y - clearance - tab_y);
  Ta11.setZ(0.);
  G4Transform3D Tr11(Ra11,Ta11);
  upper_roc_assembly->AddPlacedVolume(single_upper_fpga_volume,Tr11);

  // Now we add the stave, unless this is a double layer - in which case it does not have its own stave
  // give the stave its own logical volume, which will hold both the tube and the coolant volumes
  G4VSolid *stave_volume_box = new G4Box("stave_volume_box",stave_x, stave_y, stave_z);
  name.str("");
  name << "stave_volume_layer_" << layer;
  G4LogicalVolume *stave_volume = new G4LogicalVolume(stave_volume_box,G4Material::GetMaterial("G4_AIR"),name.str().c_str(), 0, 0, 0);
  //
  // The stave is a hollow rectangular tube, filled with Novec
  // We create the stave tube as a solid here using a subtraction boolean between the inner and outer stave boxes:
  // The subtraction solid retains the coordinate system of the first solid - i.e. the outer stave box
  //G4RotationMatrix Ra;
  G4RotationMatrix Ra12(0.0, 0.0, 0.0);
  G4ThreeVector Ta12(0.0, 0.0, 0.0);
  G4VSolid *stave_tube =  new G4SubtractionSolid("stave_tube",stave_box,stave_inner_box, &Ra12, Ta12);
  G4LogicalVolume *stave_tube_volume = new G4LogicalVolume(stave_tube,G4Material::GetMaterial("G4_C"),"stave_tube_volume", 0, 0, 0);
  G4VisAttributes* staveVisAtt = new G4VisAttributes();
  staveVisAtt->SetVisibility(true);
  staveVisAtt->SetForceSolid(true);
  staveVisAtt->SetColour(G4Colour::Yellow());
  stave_tube_volume->SetVisAttributes(staveVisAtt);
  //
  G4Transform3D Tr12(Ra12,Ta12);
  if( !option_double_layer )
    new G4PVPlacement(Tr12, stave_tube_volume, "stave_tube_volume", stave_volume, false, 0, overlapcheck);
  //
  // Now add the Novec as a box inside the stave tube
  G4VSolid *stave_coolant_box = new G4Box("stave_coolant_box", stave_inner_x - clearance, stave_inner_y - clearance, sensor_z);
  //
  // define Novec material:   C_{4} F_{9} O C H_{3} = 5 C atoms, 9 F atoms, 1 O atom, 3 H atoms
  G4Element* elH  = new G4Element("Hydrogen","H", 1.0, 1.01*g/mole);
  G4Element* elC  = new G4Element("Carbon","C", 6.0, 12.0*g/mole);
  G4Element* elO  = new G4Element("Oxygen","O", 8.0, 16.0*g/mole);
  G4Element* elF  = new G4Element("Fluorine","F", 9.0, 19.0*g/mole);
  G4Material* Novec = new G4Material("Novec", 1.40*g/cm3, 4);   // Novec is 1400 kg/m^3
  Novec->AddElement(elC , 5);
  Novec->AddElement(elF , 9);
  Novec->AddElement(elO, 1);
  Novec->AddElement(elH, 3);
  //  
  G4LogicalVolume *stave_coolant_volume = new G4LogicalVolume(stave_coolant_box,Novec,"stave_coolant_volume", 0, 0, 0);
  G4VisAttributes* coolantVisAtt = new G4VisAttributes();
  coolantVisAtt->SetVisibility(true);
  coolantVisAtt->SetForceSolid(true);
  coolantVisAtt->SetColour(G4Colour::White());
  stave_coolant_volume->SetVisAttributes(coolantVisAtt);
  // place it at the same location as the stave tube - it sits inside it
  G4RotationMatrix Ra13; 
  G4ThreeVector Ta13;
  G4Transform3D Tr13(Ra13,Ta13);
  if( !option_double_layer )
    new G4PVPlacement(Tr13, stave_coolant_volume, "stave_coolant_volume", stave_volume, false, 0, overlapcheck);
  
  // now place the stave volume inside the ladder_segment_volume, if this is a norm,al layer
  G4RotationMatrix Ra14(0.0, 0.0, 0.0);
  G4ThreeVector Ta14;
  G4double stx =   -roc_envelope_x + clearance + chip_x * 2.0 + clearance + roc_cu_x * 2.0 + clearance + roc_ins_x * 2.0 + clearance + stave_x; 
  G4double sty =  -roc_envelope_y  + clearance + stave_y;
  Ta14.setX(stx);
  Ta14.setY(sty);
  Ta14.setZ(0.);
  G4Transform3D Tr14(Ra14,Ta14);
  if( !option_double_layer )
    upper_roc_assembly->AddPlacedVolume(stave_volume,Tr14);
  //
  // Now loop over the 12 chips per ROC
  for(int i=0;i<N_chips_per_roc;i++)
    {
      G4LogicalVolume *chip_volume_up = new G4LogicalVolume(chip_box,G4Material::GetMaterial("G4_Si"),"chip_volume_up", 0,0,0);
      G4VisAttributes *chipVisAtt = new G4VisAttributes();
      chipVisAtt->SetVisibility(true);
      chipVisAtt->SetForceSolid(true);
      chipVisAtt->SetColour(G4Colour::Red());
      chip_volume_up->SetVisAttributes(chipVisAtt); 
      // Position the chip on the front face of the upper ROC card, just inside the envelope
      // and in the middle of the roc
      G4double chip_x_loc = -roc_envelope_x + clearance + chip_x;
      G4double chip_y_loc = -roc_envelope_y  + roc_y + clearance;
      G4double chip_spacing = (roc_z * 2.0) / (double) N_chips_per_roc;
      // move all the way to negative z on the roc, add the first chip at 1/2 chip spacing, add the rest at chip spacing + clearance from that
      G4double chip_z_loc = -roc_z + chip_spacing / 2.0  + (double) i * chip_spacing + (double) i * clearance;
      G4RotationMatrix Ra15;
      G4ThreeVector Ta15;
      Ta15.setX(chip_x_loc); 
      Ta15.setY(chip_y_loc);
      Ta15.setZ(chip_z_loc);
      G4Transform3D Tr15(Ra15,Ta15);
      upper_roc_assembly->AddPlacedVolume( chip_volume_up,Tr15 );
    }

  // now imprint the upper roc assembly in the (still unplaced) upper roc envelope
  G4RotationMatrix Ra16;
  // If this is a double layer, we rotate the upper ROC through 180 degrees
  // around the y axis so that the chips face outward in x
  if( option_double_layer)
    Ra16.rotateY(180.0*deg);
  G4ThreeVector Ta16;
  Ta16.setX(0.0); 
  Ta16.setY(0.0) ; 
  Ta16.setZ(0.0);
  G4Transform3D Tr16(Ra16,Ta16);
  upper_roc_assembly->MakeImprint(upper_ROC_envelope_volume, Tr16, 0, overlapcheck);
    
  //======================================================
  // The lower ROCS, chips, fpga's and stave are made by rotating the upper 
  // roc assembly around the x axis. This is done during placement of the 
  // ROC envelope in the lower ROC location 
  //======================================================

  //=====================================================
  //  We now have three logical volumes that together make up a ladder section:
  //     upper_ROC_envelope_volume
  //     lower_ROC_envelope_volume (will be made by rotation and translation of upperROC_envelope_volume)
  //     sensor_envelope_volume
  // We can combine them into a single volume that has all three assembled 
  // inside it - call it ladder_segment_volume
  //=====================================================

  // Create a volume to hold the ladder segment
  // the roc envelopes are bigger in x than the sensor, so the x dimension has to accomodate them
  // if this is a double layer, roc_envelope volume has already had the stave width removed
  G4double ladder_segment_x = roc_envelope_x +  clearance;
  G4double ladder_segment_y = sensor_envelope_y + clearance + 2.0 * roc_envelope_y  + 2.0 * clearance;
  // if we do not add the lower roc, we reduce the size of the ladder segment accordingly
  if( !add_lower_roc )
    ladder_segment_y = sensor_envelope_y + clearance + 1.0 * roc_envelope_y  + 1.0 * clearance;
  G4double ladder_segment_z = roc_envelope_z + clearance;
  G4Box *ladder_segment_box = new G4Box("ladder_segment_box", ladder_segment_x, ladder_segment_y, ladder_segment_z);
  name.str("");
  name << "ladder_segment_volume_layer_" << layer;
  G4LogicalVolume *ladder_segment_volume = new G4LogicalVolume(ladder_segment_box, G4Material::GetMaterial("G4_AIR"), name.str().c_str(), 0,0,0);
  G4VisAttributes* ladder_segment_VisAtt = new G4VisAttributes();
  ladder_segment_VisAtt->SetVisibility(false);
  ladder_segment_VisAtt->SetForceSolid(false);
  ladder_segment_VisAtt->SetColour(G4Colour::Red());
  ladder_segment_volume->SetVisAttributes(ladder_segment_VisAtt);

  // create the assembly helper for filling the ladder segment volume
  G4AssemblyVolume *ladder_segment_assembly = new G4AssemblyVolume();

  // Place the sensor envelope in the ladder segment envelope
  //--------------------------------------------------------------------------
  // It is centered in z and in y (if there are two ROCs)
  // For a normal layer, in x, move it all the way to the front of the ladder segment envelope, then
  // move it out by the chip width and half the sensor envelope width
  sensor_x_offset = -ladder_segment_x + clearance + chip_x * 2.0 + clearance + sensor_envelope_x; 
  if(option_double_layer)
    {
      // if this is a double layer insread, the roc is in front of the chip and there is no stave
      // so the sensor should be all the way to the front of the ladder_segment_envelope
      sensor_x_offset -= chip_x * 2.0 + clearance; 
    }
  G4RotationMatrix Ra17;
  Ra17.rotateZ(0.0);
  G4ThreeVector Ta17;
  Ta17.setX(sensor_x_offset);
  // if upper and lower ROCs are in place, the sensor is centered in y, otherwise it is not!
  sensor_y_offset = 0.0;
  if( !add_lower_roc)
    sensor_y_offset = -ladder_segment_y + clearance + sensor_envelope_y;
  Ta17.setY(sensor_y_offset); 
  Ta17.setZ(0.0);
  G4Transform3D Tr17(Ra17,Ta17);
  ladder_segment_assembly->AddPlacedVolume(sensor_envelope_volume, Tr17);

  // Place the upper roc envelope in the ladder segment envelope
  //--------------------------------------------------------------------------
  // Remember that the upper roc envelope already has the chips and fpga orientation right for normmal or double layers
  // It is centered in z, in y, move it all the way to the top
  // In x, move it all the way to the front of the ladder segment envelope, then
  // move it out by half the upper ROC envelope width
  G4RotationMatrix Ra18;
  G4ThreeVector Ta18;
  Ta18.setX(-ladder_segment_x + clearance + roc_envelope_x); 
  Ta18.setY(+ladder_segment_y - clearance - roc_envelope_y) ; 
  Ta18.setZ(0.0);
  G4Transform3D Tr18(Ra18,Ta18);
  ladder_segment_assembly->AddPlacedVolume(upper_ROC_envelope_volume, Tr18);

  if( add_lower_roc )
    {
      // Place the lower roc envelope in the ladder segment envelope
      //--------------------------------------------------------------------------
      // Rotate it 180 degress about the x axis to get the orientation right (i.e. chips facing large x)
      // This works for a normal layer or a double layer
      G4RotationMatrix Ra;
      Ra.rotateX(180.0 * deg);
      G4ThreeVector Ta;
      Ta.setX(-ladder_segment_x + clearance + roc_envelope_x); 
      Ta.setY(-ladder_segment_y + clearance + roc_envelope_y) ; 
      Ta.setZ(0.0);
      G4Transform3D Tr(Ra,Ta);
      ladder_segment_assembly->AddPlacedVolume(upper_ROC_envelope_volume, Tr);
    }

  // Now instantiate the assembly in the new ladder_segment_volume
  G4RotationMatrix Ra19;
  Ra19.rotateZ(0.0);
  G4ThreeVector Ta19;
  Ta19.setX(0.0); 
  Ta19.setY(0.0) ; 
  Ta19.setZ(0.0);
  G4Transform3D Tr19(Ra19,Ta19);
  if (verbosity>0) cout << "Place ladder_segment_assemly (upper_ROC_enevelope_vlume + sensor_envelope_volume) in the ladder_segment_volume" << endl;
  ladder_segment_assembly->MakeImprint(ladder_segment_volume, Tr19, 0, overlapcheck);
    
  //=========================================
  // Now we populate the whole layer with the ladder segments
  // that we created above
  //==========================================

  // Get the number of sensors in z in this barrel
  // We want to place enough sensors in a line in z to cover the rapidity range
  // maxrap is specified in the constructor - it will normally be 1.1
  segment_z_spacing = ladder_segment_z * 2.0 + clearance;
  G4double theta_max = 2.0 * atan(exp(-maxrap));
  G4double layer_length = 2.0 * layer_nominal_radius * tan(M_PI/2.0 - theta_max);
  // do not allow the inner layers to be less than 40 cm in z
  if(layer_length < 40.0 * cm)
    layer_length = 40.0 * cm; 
  layer_NZ = (int) ( layer_length / segment_z_spacing );

  // get the interval in phi between sensors in the barrel
  // we want to consider only the sensor active phi here, so we subtract the inactive edge
  G4double ladder_segment_y_offset = (sensor_y -sensor_edge_phi) * 2.0 - (sensor_y - sensor_edge_phi)  * 2.0 * overlap_fraction;
  layer_NPHI = (int) (2.0 * M_PI * layer_nominal_radius / ladder_segment_y_offset);
  segment_phi_spacing = 2.0 * M_PI / (double) layer_NPHI;

  if(verbosity>0)
    cout << " sensor y " << sensor_y
	 << " ladder_segment_y_offset " << ladder_segment_y_offset
	 << " ladder_segment_y " << ladder_segment_y
	 << " layer_nominal_radius " << layer_nominal_radius 
	 << " layer_NPHI " << layer_NPHI
	 << " segment_phi_spacing " << segment_phi_spacing
	 << endl;

  if(verbosity > 0)
		   cout << "Layer " << layer << ":  enter loop over  " << layer_NZ << " z positions " 
		   << " and "  << layer_NPHI << " phi positions " 
		   << endl;
  
  // if this is a double layer, the nominal radius is the same as the first layer
  // we now increase it by the stave thickness so that this layer is on the other side of the common stave
  // we did not do this earlier because we want the phi spacing to match the first layer
  if(option_double_layer)
    layer_nominal_radius += stave_x * 2.0+clearance;
  
  for(int iz = 0;iz<layer_NZ;iz++)
    {
      // step through all of the z locations in turn
      G4double z_spacing = ladder_segment_z * 2.0 + clearance;
      G4double z_location =  (double) (iz-layer_NZ/2) * z_spacing;
      
      // At each z location. step through all phi values
      int istagger = 0; 
     for (int iphi=0; iphi<layer_NPHI; iphi++)
	{
	  // Place the ladder segment envelopes at the correct z and phi
	  G4double phi_rotation = (double) iphi * segment_phi_spacing;

	  // this determines the stggered layer radius
	  istagger = iphi % N_staggers;

	  // We need to stagger the radii at alternate phi values by radius_stagger, since the ladders overlap in phi
	  // The number of staggers is an input number, since it has to be the same for both parts of a double layer!
	  G4double R_layer;
	  R_layer = layer_nominal_radius + (double) istagger * radius_stagger;

	  if(verbosity > 1)
	    cout << " iphi " << iphi
		 << " istagger " << istagger
		 << " N_staggers " << N_staggers
		 << " R_layer " << R_layer
		 << endl;
	  
	  // Check that there is still enough space for this ladder segment to not touch the first one at phi = 0
	  // this requires there to be two phi half-widths between the center of this ladder and the one at phi = 0
	  // A ladder segment has a phi *half width* at this radius of: 
	  G4double ladder_segment_phi = ladder_segment_y / R_layer;
	  if(2.0 * M_PI - phi_rotation <= 2.0 * ladder_segment_phi && istagger == 0)
	    {
	      if (verbosity > 0) cout << " dropping ladder segment for iphi = " << iphi << " to avoid overlap " << endl;
	      break;
	    }
	
	  // there is no version of AddPlacedVolume  that does the translation followed by the rotation
	  // The ladder segment is first rotated in phi
	  G4RotationMatrix Ra;
	  Ra.rotateZ(phi_rotation);
	  // Then translated as folows
	  G4ThreeVector Ta;
	  Ta.setX(R_layer * cos(phi_rotation)); 
	  Ta.setY(R_layer * sin(phi_rotation)) ; 
	  Ta.setZ( z_location );
	  G4Transform3D Tr(Ra,Ta);

	  G4int copynum = 1000 * iz + iphi;
	  // name the ladder_segment
	  name.str("");
	  name << "layer_" << layer << "_ladder_segment_" << iz << "_" << iphi;
	  new G4PVPlacement(Tr, ladder_segment_volume, name.str().c_str(), trackerenvelope, false, copynum, overlapcheck);

	  // This ladder segment was actually added, so increment the sensor count for this layer 
	  N_sensors_in_layer++;
	} 
    }
  if(verbosity > 0)
    cout << "This layer has a total of " << N_sensors_in_layer << " sensors" << endl;

  return 0;
}

int
PHG4SiliconTrackerDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
{
  static int i = 0;
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume,G4Material::GetMaterial("G4_POLYSTYRENE"),"DISPLAYLOGICAL", 0, 0, 0);
  G4VisAttributes* visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  switch(i)
    {
    case 0:
      visattchk->SetColour(G4Colour::Red());
      i++;
      break;
    case 1:
      visattchk->SetColour(G4Colour::Magenta());
      i++;
      break;
    case 2:
      visattchk->SetColour(G4Colour::Yellow());
      i++;
      break;
    case 3:
      visattchk->SetColour(G4Colour::Blue());
      i++;
      break;
    case 4:
      visattchk->SetColour(G4Colour::Cyan());
      i++;
      break;
    default:
      visattchk->SetColour(G4Colour::Green());
      i=0;
      break;
    }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm,G4ThreeVector(0,0,0),checksolid,"DISPLAYVOL",logvol, 0, false, overlapcheck);
  return 0;
}

void
PHG4SiliconTrackerDetector::AddGeometryNode()
{
  if (active)
    {
      ostringstream geonode;
      if (superdetector != "NONE")
	{
           geonode << "CYLINDERGEOM_" << superdetector;
        }
      else
        {
          geonode << "CYLINDERGEOM_" << detector_type << "_" << layer;
        }
      PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonode.str().c_str());
      if (!geo)
        {
          geo = new PHG4CylinderGeomContainer();
          PHNodeIterator iter( topNode );
          PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.str().c_str(), "PHObject");
          runNode->addNode(newNode);
        }
      // here in the detector class we have internal units, convert to cm
      // before putting into the geom object
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv4(N_sensors_in_layer, layer_NZ, N_strips_per_column, N_strip_columns, N_staggers,
							layer_nominal_radius/cm, radius_stagger/cm, segment_z_spacing/cm, segment_phi_spacing/rad,
							sensor_x_offset/cm, sensor_y_offset/cm, strip_z_spacing/cm, strip_y_spacing/cm, strip_x/cm, strip_tilt/rad);

      geo->AddLayerGeom(layer, mygeom);
      geo->identify();
    }
}
