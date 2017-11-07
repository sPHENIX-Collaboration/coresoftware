#include "PHG4MapsDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_MAPS.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4GDMLParser.hh>
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
#include <memory>

using namespace std;

//static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes

PHG4MapsDetector::PHG4MapsDetector( PHCompositeNode *Node,  PHG4Parameters *parameters, const std::string &dnam ):
  PHG4Detector(Node, dnam),
  //envelope_inner_radius(26.0*mm),
  //envelope_outer_radius(880*mm),
  //envelope_z(2300*mm + no_overlap),
  place_in_x(0 * cm),
  place_in_y(0 * cm),
  place_in_z(0 * cm),
  x_rot(0),
  y_rot(0),
  z_rot(0),
  active(parameters->get_int_param("active")),
  absorberactive(parameters->get_int_param("absorberactive")),
  layer(parameters->get_int_param("layer")),
  blackhole(parameters->get_int_param("blackhole")),
  stave_type(parameters->get_int_param("stave_type")),
  layer_nominal_radius(parameters->get_double_param("layer_nominal_radius")),
  N_staves(parameters->get_int_param("N_staves")),
  phistep(NAN),
  phitilt(parameters->get_double_param("phitilt")),
  pixel_x(parameters->get_double_param("pixel_x")),
  pixel_z(parameters->get_double_param("pixel_z")),
  pixel_thickness(parameters->get_double_param("pixel_thickness")),
  stave_geometry_file(parameters->get_string_param("stave_geometry_file"))
{
//  verbosity = 2;

  if(verbosity > 0)
    cout << "PHG4MapsDetector constructor called" << endl;

  if(verbosity > 10)
    cout << " cm " << cm << " mm " << mm << endl;

  ostringstream name;
  name.str("");
  name << "ITSUSensor" << layer;
  layer_string = name.str().c_str();

  if (verbosity > 0) 
    cout << "PHG4MapsDetector constructor: making layer with sensor string: " << layer_string.c_str() << "using stave_type " << stave_type << endl;
}

PHG4MapsDetector::~PHG4MapsDetector()
{}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4MapsDetector::IsInMaps(G4VPhysicalVolume * volume) const
{
  // Is this volume one of the sensors?

  if (volume->GetName().find("ITSUSensor") != string::npos)
    {
      // Check to see if this strip is in the layer belonging to this instance of MapsDetector
      // layer_string contains "layer_n", where n = the layer number for this instance
      if (volume->GetName().find(layer_string.c_str()) != string::npos)
      return 1;
    }
  
  return 0;
}

void
PHG4MapsDetector::Construct( G4LogicalVolume* logicWorld )
{
  // This is called from PHG4PhenixDetector::Construct()

  if(verbosity>0)
    cout << endl << "PHG4MapsDetector::Construct called for layer " << layer << endl;

  if(layer_nominal_radius < 0 || stave_type < 0)
    {
      cout << PHWHERE << "layer radius or stave type undefined, quit!" << endl;
      exit(1);
    }

  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  // this reads in the ITS stave geometry from a file and constructs the layer from it 
  ConstructMaps(logicWorld);

    // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int
PHG4MapsDetector::ConstructMaps(G4LogicalVolume* trackerenvelope)
{
  if(verbosity>0)
    {
      cout << " PHG4MapsDetector::ConstructMaps:" <<endl;
      cout << " Constructing Layer " << layer << endl;
      cout << endl;
    }


  //===================================
  // Import the stave physical volume here
  //===================================

  // import the staves from the gemetry file
  std::unique_ptr<G4GDMLReadStructure>  reader (new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.Read(stave_geometry_file, false);
  
  // figure out which assembly we want
  char assemblyname[500];
  sprintf(assemblyname, "ITSUStave%i",layer);

  if (Verbosity())
    cout << "Geting the stave assembly named " << assemblyname << endl;
  G4AssemblyVolume* av_ITSUStave = reader->GetAssembly(assemblyname);

//  if (reader) delete reader;

  //=========================================
  // Now we populate the whole layer with the staves
  //==========================================

  // Step through all phi values

  // The phistep is determined by the arc length spacing needed bewtween staves, which depends on the stave type
  // The numbers below are inferred from the arc length spacing in ITS.gdml, and reduced a small amount (< 1%) to make sure the default radius cases work
  // These arcstep numbers matter, that is why they are hardcoded here - leave them alone!
  double arcstep;
  if(stave_type == 0)
    {
      arcstep = 12.25;  // inner barrel, mm
    }
  else if(stave_type == 1)
    {
      arcstep = 56.8;  // middle barrel, mm
    }
  else
    {
      // stave_type is 2
      arcstep = 54.2;  // outer barrel, mm
    }

  if(N_staves < 0)
    {
      // The user did not specify how many staves to use for this layer, so we calculate the best value
      // We choose a phistep that fits N_staves into this radius, but with an arclength separation of AT LEAST arcstep
      // ideally, the radius would be chosen so that numstaves = N_staves exactly, to get the closest spacing of staves possible without overlaps
      double numstaves = 2.0 * M_PI * layer_nominal_radius / arcstep;  // this is just to print out
      N_staves = int(2.0 * M_PI * layer_nominal_radius / arcstep);  // this is the number of staves used  

      // Also suggest the ideal radius for this layer
      if (verbosity > 0)
	{
	  cout << " Calculated N_staves for layer " << layer 
	       << " layer_nominal_radius " << layer_nominal_radius
	       << " ITS arcstep " << arcstep
	       << " circumference divided by arcstep  " << numstaves 
	       << " N_staves " << N_staves
	       << endl;
	  cout << "A radius for this layer of " << (double) N_staves * arcstep / (2.0 * M_PI)  + 0.01 << " or " 
	       <<  (double) (N_staves+1) * arcstep / (2.0 * M_PI) + 0.01 << " would produce  perfect stave spacing" << endl;
	}
    }
  phistep = 2.0 * M_PI / (double) N_staves;  // this produces even stave spacing
  double z_location = 0.0;

  // The tilt for stave type 0 (usually layers 0-2) is phitilt = 0.304;   // radians, equivalent to 17.4 degrees, now input from parameters
  if(stave_type != 0)
    phitilt = 0.0;

  if (Verbosity())
    {
      cout << " layer " << layer 
	   << " layer_nominal_radius " << layer_nominal_radius
	   << " N_staves " << N_staves
	   << " phistep " << phistep
	   << " phtilt " << phitilt
	   << endl;
    }

  // The stave starts out at (0,0,0) and oriented so that the sensors face upward in y
  // So we need to rotate the sensor 90 degrees before placing it using phi_offset
  // note that the gdml file uses a negative phi_offset - different coord system, apparently - the following works
  double phi_offset =  M_PI /2.0;

  for (int iphi=0; iphi<N_staves; iphi++)
    {
      // Place the ladder segment envelopes at the correct z and phi 
      // This is the azimuthal angle at which we place the stave
      G4double phi_rotation = (double) iphi * phistep;

      G4RotationMatrix Ra;
      G4ThreeVector Ta;
      
      if(verbosity >0)
	cout << "phi_offset = " << phi_offset << " iphi " << iphi << " phi_rotation = " << phi_rotation << " phitilt " << phitilt << endl;

      // It  is first rotated in phi by the azimuthal angle phi_rotation, plus the 90 degrees needed to point the face of the sensor  at the origin,  plus the tilt (if a tilt is appropriate)
      
      Ra.rotateZ(phi_rotation + phi_offset + phitilt); // note - if this is layer 0-2, phitilt is the additional tilt for clearance. Otherwise it is zero
      // Then translated as follows
      
      Ta.setX(layer_nominal_radius * cos(phi_rotation)); 
      Ta.setY(layer_nominal_radius * sin(phi_rotation)) ; 
      Ta.setZ( z_location );
      
      if(verbosity > 0)
	cout << " iphi " << iphi << " phi_rotation " << phi_rotation 
	     << " x " << layer_nominal_radius * cos(phi_rotation)
	     << " y " <<  layer_nominal_radius * sin(phi_rotation)
	     << " z " << z_location 
	     << endl;      	  
      
      
      G4Transform3D Tr(Ra,Ta);
      
      av_ITSUStave->MakeImprint(trackerenvelope, Tr, 0, overlapcheck);
    } 
  
  if(verbosity > 0)
    cout << "This layer has a total of " << N_staves << " staves" << endl;

  SetDisplayProperty(av_ITSUStave);
  
  return 0;
}

void PHG4MapsDetector::SetDisplayProperty( G4AssemblyVolume* av)
{
  
  //  cout <<"SetDisplayProperty - G4AssemblyVolume w/ TotalImprintedVolumes "<<av->TotalImprintedVolumes()
  //   <<"/"<<av->GetImprintsCount()<<endl;
  
  std::vector<G4VPhysicalVolume*>::iterator it = av->GetVolumesIterator();
  
  int nDaughters = av->TotalImprintedVolumes();
  for(int i = 0; i < nDaughters; ++i, ++it)
    {
      //  cout <<"SetDisplayProperty - AV["<<i<<"] = "<<(*it)->GetName()<<endl;
      G4VPhysicalVolume* pv =(*it);
      
      G4LogicalVolume* worldLogical = pv->GetLogicalVolume();
      SetDisplayProperty(worldLogical);
    }
}

void PHG4MapsDetector::SetDisplayProperty( G4LogicalVolume* lv)
{
  string material_name(
  lv->GetMaterial()->GetName());

  if (Verbosity() >= 5)
    cout << "SetDisplayProperty - LV " << lv->GetName() << " built with "
        << material_name << endl;

  G4VisAttributes* matVis= new G4VisAttributes();
  if (material_name.find("SI") != std::string::npos)
    {
      PHG4Utils::SetColour(matVis, "G4_Si");
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
      if (Verbosity() >= 5)
      cout << "SetDisplayProperty - LV " << lv->GetName() << " display with G4_Si" << endl;
    }
  else if (material_name.find("KAPTON") != std::string::npos)
    {
      PHG4Utils::SetColour(matVis, "G4_KAPTON");
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
      if (Verbosity() >= 5)
      cout << "SetDisplayProperty - LV " << lv->GetName() << " display with G4_KAPTON" << endl;
    }
  else if (material_name.find("ALUMINUM") != std::string::npos)
    {
      PHG4Utils::SetColour(matVis, "G4_Al");
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
      if (Verbosity() >= 5)
      cout << "SetDisplayProperty - LV " << lv->GetName() << " display with G4_Al" << endl;
    }
  else if (material_name.find("Carbon") != std::string::npos)
    {
      matVis->SetColour(0.5,0.5,0.5,.25);
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
      if (Verbosity() >= 5)
      cout << "SetDisplayProperty - LV " << lv->GetName() << " display with Gray" << endl;
    }
  else if (material_name.find("M60J3K") != std::string::npos)
    {
      matVis->SetColour(0.25,0.25,0.25,.25);
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
      if (Verbosity() >= 5)
      cout << "SetDisplayProperty - LV " << lv->GetName() << " display with Gray" << endl;
    }
  else if (material_name.find("WATER") != std::string::npos)
    {
      matVis->SetColour(0.0,0.5,0.0,.25);
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
      if (Verbosity() >= 5)
      cout << "SetDisplayProperty - LV " << lv->GetName() << " display with WATER" << endl;
    }
  else
    {
      matVis->SetColour(.2,.2,.7,.25);
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
    }
  lv->SetVisAttributes(matVis);
  
  int nDaughters = lv->GetNoDaughters();
  for(int i = 0; i < nDaughters; ++i)
    {
      G4VPhysicalVolume* pv = lv->GetDaughter(i);
      
      // cout <<"SetDisplayProperty - PV["<<i<<"] = "<<pv->GetName()<<endl;
      
      
      G4LogicalVolume* worldLogical = pv->GetLogicalVolume();
      SetDisplayProperty(worldLogical);
    }
}

void
PHG4MapsDetector::AddGeometryNode()
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
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeom_MAPS(layer, stave_type, N_staves, layer_nominal_radius/cm, phistep/rad, phitilt/rad, pixel_x, pixel_z, pixel_thickness);
      geo->AddLayerGeom(layer, mygeom);
      if (Verbosity())
        geo->identify();
    }
}
