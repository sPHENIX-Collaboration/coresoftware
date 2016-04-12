#include "PHG4MapsDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_MAPS.h"

#include "TMath.h"

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

using namespace std;

//static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes

PHG4MapsDetector::PHG4MapsDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr, const int in_stave_type  ):
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
  active(0),
  absorberactive(0),
  layer(lyr),
  stave_type(in_stave_type),
  layer_nominal_radius(-1)
{
  verbosity = 2;

  if(verbosity > 0)
    cout << "PHG4MapsDetector constructor called" << endl;

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
  G4GDMLReadStructure * reader = new G4GDMLReadStructure();
  G4GDMLParser gdmlParser(reader);
  gdmlParser.Read("./ITS.gdml");
  
  // figure out which assembly we want
  char assemblyname[500];
  sprintf(assemblyname, "ITSUStave%i",layer);

  cout << "Geting the stave assembly named " << assemblyname << endl;
  G4AssemblyVolume* av_ITSUStave = reader->GetAssembly(assemblyname);

  /*
  G4AssemblyVolume* av_ITSUStave0 = reader->GetAssembly("ITSUStave0");
  G4AssemblyVolume* av_ITSUStave1 = reader->GetAssembly("ITSUStave1");
  G4AssemblyVolume* av_ITSUStave2 = reader->GetAssembly("ITSUStave2");
  G4AssemblyVolume* av_ITSUStave3 = reader->GetAssembly("ITSUStave3");
  G4AssemblyVolume* av_ITSUStave4 = reader->GetAssembly("ITSUStave4");
  G4AssemblyVolume* av_ITSUStave5 = reader->GetAssembly("ITSUStave5");
  G4AssemblyVolume* av_ITSUStave6 = reader->GetAssembly("ITSUStave6");
  */

  //=========================================
  // Now we populate the whole layer with the staves
  //==========================================

  /*      
  // In the ITS we should have these numbers for how staves are built (from ITS.gdml file)
   lyr rad   L   staves     modules                                                                                                     chips/module
   0  23    290    12    1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   1  31    290    16   1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   2  39    290    20   1 (x=0, y=0, z=0)                                                                                     9 (x=0, y=-0.00875, z=- 12.04, -9,03, -6.02, -3.01, 0, 3.01, 6.02, 9.03, 12.04) 
   3  194  900    24   4 (x=0, y=-0.06075, z= -31.605, -10.535, 10.535, 31.605)                    14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03)      
   4  247  900    30   4 (x=0, y=-0.06075, z= -31.605, -10.535, 10.535, 31.605)                    14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03)      
   5  253 1500   42   7 (x=0, y=-0.06075, z = -63.21, -42.14, -21.07, 0.0, 21.07, 42.14, 63.21)   14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03) 
   6  405 1500  48  7 (x=0, y=-0.06075, z = -63.21, -42.14, -21.07, 0.0, 21.07, 42.14, 63.21)   14  (x = -0.755 or +0.755, y= -0.00825, z = -9.03, -9.03, -6.02, -6.02, -3.01, -3.01, 0, 0, 3.01, 3.01, 6.02, 6.02, 9.03, 9.03) 

   sensor is in chip at (x=0, y=-0.0016, z=0)
   3-6 half-staves are in 3-6 staves at:  (x = -1.29, y = +2.067, z = 0)  or (x = +1.29 cm, y = 2.243, z = 0)

   // layers 0,1,2 have one stave with 1 module and 7 chips in that module
   // where layers 3 and 4  have two half-staves with 4 modules and 7 chips/module
   // layers 5 and 6 have two half staves with 7 modules and 7 chips/module 
   // module coordinates are in the stave (0-2) or in the half staves (3-6) 
   // chip coordinates are in the module
   // 
  */

  // Step through all phi values
  double arcstep;
  if(stave_type == 0)
    {
      arcstep = 12.0;  // inner barrel, mm
    }
  else if(stave_type == 1)
    {
      arcstep = 50.5;  // middle barrel, mm
    }
  else
    {
      // stave_type is 2
      arcstep = 52.0;  // outer barrel, mm
    }

  N_staves = int(2.0 * M_PI * layer_nominal_radius / arcstep);
  phistep = 2.0 * M_PI / (double) N_staves;
  double z_location = 0.0;
  
  // this is the tilt for stave type 0 (usually layers 0-2)
  phitilt = 0.304;   // radians, equivalent to 17.4 degrees
  if(stave_type != 0)
    phitilt = 0.0;

  if(verbosity > 0)
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
  double phi_offset =  M_PI /2.0;

  for (int iphi=0; iphi<N_staves; iphi++)
  //for (int iphi=0; iphi<1; iphi++)
    {
      // Place the ladder segment envelopes at the correct z and phi 
      // This is the azimuthal angle at which we place the stave
      G4double phi_rotation = (double) iphi * phistep;

      G4RotationMatrix Ra;
      G4ThreeVector Ta;
      double notest = true;

      if(notest)
	{
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
	}
      else
	{
	  // testing!!!
	  
	  double phitest = -88.9770648826 * M_PI / 180.0;
	  double xtest = 2.24456394319 * cm;
	  double ytest = 0.661462549917 * cm;
	  double ztest = 0.0 * cm;
	  
	  cout << "phitest = " << phitest << " xtest " << xtest << " ytest " << ytest << " ztest " << ztest << endl;
	  
	  // It  is first rotated in phi by the azimuthal angle phi_rotation, plus the 90 degrees needed to point the face of the sensor  at the origin,  plus the tilt (if a tilt is appropriate)

	  Ra.rotateZ(phitest); 
	  // Then translated as follows
	  Ta.setX(xtest); 
	  Ta.setY(ytest) ; 
	  Ta.setZ( ztest);
	}
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
  
  //cout <<"SetDisplayProperty - LV "<<lv->GetName()<<endl;
  
  
  G4VisAttributes* matVis = new G4VisAttributes();
  matVis->SetColour(.2,.2,.7,.25);
  matVis->SetVisibility(true);
  matVis->SetForceSolid(true);
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
  // How to handle geometry for MAPS pixels:
  // 1) Find the center of the sensor that contains the hit.
  //     
  //
  //
  //
  //
  //

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
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeom_MAPS(layer, stave_type, N_staves, layer_nominal_radius/cm, phistep/rad, phitilt/rad);

      geo->AddLayerGeom(layer, mygeom);
      geo->identify();
    }
}
