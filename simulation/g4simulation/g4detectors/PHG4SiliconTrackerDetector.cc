#include "PHG4SiliconTrackerDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Siladders.h"
#include "PHG4SiliconTrackerParameterisation.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4PVParameterised.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <TMath.h>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

using namespace std;

PHG4SiliconTrackerDetector::PHG4SiliconTrackerDetector(PHCompositeNode *Node, PHParametersContainer *parameters, const std::string &dnam, const vpair &layerconfig)
  : PHG4Detector(Node, dnam)
  , paramscontainer(parameters)
{
  layerconfig_ = layerconfig;

  layermin_ = layerconfig_.front().first;
  layermax_ = layerconfig_.back().first;
  nlayer_ = layerconfig_.size();

  for (unsigned int ilayer = 0; ilayer < nlayer_; ++ilayer)
  {
    const int inttlayer = layerconfig_[ilayer].second;
    if (inttlayer < 0 || inttlayer >= 4)
    {
      assert(!"PHG4SiliconTrackerDetector: check INTT ladder layer.");
    }
  }
  PHParametersContainer::ConstRange begin_end = paramscontainer->GetAllParameters();
  for (PHParametersContainer::ConstIterator iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    PHParameters *par = iter->second;
    IsActive[iter->first] = par->get_int_param("active");
    IsAbsorberActive[iter->first] = par->get_int_param("absorberactive");
  }
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4SiliconTrackerDetector::IsInSiliconTracker(G4VPhysicalVolume *volume) const
{
  // Is this volume one of the sensor strips?
  // just checking if the pointer to the logical volume is in the set
  // of our active/active absorber ones makes sure we are in an active volume
  // name parsing is a bad idea since this is called for all steps
  // and we would have to trust that people give different names
  // to their volumes
  G4LogicalVolume *logvol = volume->GetLogicalVolume();
  if (!absorberlogvols.empty() && absorberlogvols.find(logvol) != absorberlogvols.end())
  {
    return -1;
  }
  if (activelogvols.find(logvol) != activelogvols.end())
  {
    return 1;
  }

  return 0;
}

void PHG4SiliconTrackerDetector::Construct(G4LogicalVolume *logicWorld)
{
  if (Verbosity() > 0)
    std::cout << "PHG4SiliconTrackerDetector::Construct called for layers " << layermin_ << " to " << layermax_ << std::endl;

  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  ConstructSiliconTracker(logicWorld);

  // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int PHG4SiliconTrackerDetector::ConstructSiliconTracker(G4LogicalVolume *trackerenvelope)
{
  // We have an arbitray number of layers (nlayer_)
  // We have 2 types of ladders (vertical strips and horizontal strips)
  // We have 2 types of sensors (inner and outer)

  double hdi_z_[nlayer_][2];
  // we loop over layers. All layers have only one laddertype
  for (unsigned int ilayer = 0; ilayer < nlayer_; ++ilayer)
    {
      const int sphxlayer = layerconfig_[ilayer].first;
      const int inttlayer = layerconfig_[ilayer].second;

      // get the parameters for this layer
      const PHParameters *params1 = paramscontainer->GetParameters(inttlayer);
      const int laddertype = params1->get_int_param("laddertype");
      double layer_radius_inner = params1->get_double_param("ladder_radius_inner") * cm;
      double layer_radius_outer = params1->get_double_param("ladder_radius_outer") * cm;

      // Look up all remaining parameters by the laddertype for this layer
      const PHParameters *params = paramscontainer->GetParameters(laddertype);
      const G4double strip_x = params->get_double_param("strip_x") * cm;
      const G4double strip_y = params->get_double_param("strip_y") * cm;
      const int nstrips_phi_sensor = params->get_int_param("nstrips_phi_sensor");
      const int nladders_layer = params->get_int_param("nladders");
      const G4double offsetphi = params->get_double_param("offsetphi") * deg;
      const G4double offsetrot = params->get_double_param("offsetrot") * deg;
      const G4double hdi_y = params->get_double_param("hdi_y") * cm;
      double hdi_x = params->get_double_param("hdi_x") * cm;
      double fphx_x = params->get_double_param("fphx_x") * cm;
      double fphx_y = params->get_double_param("fphx_y") * cm;
      double fphx_z = params->get_double_param("fphx_z") * cm;
      double pgs_x = params->get_double_param("pgs_x") * cm;
      double halfladder_z = params->get_double_param("halfladder_z") * cm;
      double stave_straight_outer_y = params->get_double_param("stave_straight_outer_y") * cm;
      double stave_straight_inner_y = params->get_double_param("stave_straight_inner_y") * cm;
      double stave_straight_cooler_y = params->get_double_param("stave_straight_cooler_y") * cm;

      // We loop over inner, then outer, sensors, where  itype specifies the inner or outer sensor type
      // The rest of this loop will construct and put in place a section of a ladder corresponding to the 
      // Z range of this sensor only
      for (int itype = 0; itype < 2; ++itype)
	{
	  if (!(itype >= 0 && itype <= 1))
	    {
	      assert(!"Error: check ladder type.");
	    }
	  G4double strip_z;
	  int nstrips_z_sensor;
	  switch (itype)
	    {
	    case 0:
	      strip_z = params->get_double_param("strip_z_0") * cm;
	      nstrips_z_sensor = params->get_int_param("nstrips_z_sensor_0");
	      break;
	    case 1:
	      strip_z = params->get_double_param("strip_z_1") * cm;
	      nstrips_z_sensor = params->get_int_param("nstrips_z_sensor_1");
	      break;
	    default:
	      cout << "invalid itype " << itype << endl;
	      exit(1);
	    }

	  // ----- Step 1 ---------------------------------------------------------------------------------------------
	  // We make the volumes for Si-sensor, FPHX, HDI, PGS sheet, and stave components
	  // We add them to the ladder later
	  //============================================================
      
	  // Create a volume for the Si-strip. Use half widths to make the box
	  G4VSolid *strip_box = new G4Box(boost::str(boost::format("strip_box_%d_%d") % sphxlayer %itype).c_str(), 
					  strip_x / 2., strip_y / 2. - strip_y / 20000., strip_z / 2. - strip_z / 2. / 10000.);
	  G4LogicalVolume *strip_volume = new G4LogicalVolume(strip_box, G4Material::GetMaterial("G4_Si"), 
							      boost::str(boost::format("strip_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsActive.find(inttlayer))->second > 0)
	    {
	      activelogvols.insert(strip_volume);
	    }
	  G4VisAttributes *strip_vis = new G4VisAttributes();
	  strip_vis->SetVisibility(false);
	  strip_vis->SetForceSolid(false);
	  strip_vis->SetColour(G4Colour::White());
	  strip_volume->SetVisAttributes(strip_vis);
      
	  // Create Si-sensor active volume
	  const double siactive_x = strip_x;
	  const double siactive_y = strip_y * nstrips_phi_sensor; 
	  const double siactive_z = strip_z * nstrips_z_sensor;  
      
	  G4VSolid *siactive_box = new G4Box(boost::str(boost::format("siactive_box_%d_%d") % sphxlayer % itype).c_str(), siactive_x, siactive_y, siactive_z);
	  G4LogicalVolume *siactive_volume = new G4LogicalVolume(siactive_box, G4Material::GetMaterial("G4_Si"), 
								 boost::str(boost::format("siactive_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
      
	  // Now make a G4PVParameterised containing all of the strips in a sensor 
	  // this works for ladder 0 because there is only one strip type - all cells are identical
	  G4VPVParameterisation *stripparam = new PHG4SiliconTrackerStripParameterisation(nstrips_phi_sensor, nstrips_z_sensor, strip_y, strip_z);
	  new G4PVParameterised(boost::str(boost::format("siactive_%d_%d") % sphxlayer % itype).c_str(), 
				strip_volume, siactive_volume, kZAxis, nstrips_phi_sensor * nstrips_z_sensor, stripparam, false);  // overlap check too long.
      
      
	  // Si-sensor full (active+inactive) area
	  const double sifull_x = siactive_x;
	  const double sifull_y = siactive_y + 2.0*params->get_double_param("sensor_edge_phi") * cm; 
	  const double sifull_z = siactive_z + 2.0*params->get_double_param("sensor_edge_z") * cm;    
      	  G4VSolid *sifull_box = new G4Box(boost::str(boost::format("sifull_box_%d_%d") % sphxlayer % itype).c_str(), sifull_x / 2., sifull_y/2.0, sifull_z/2.0);
      
	  // Si-sensor inactive area
      	  G4VSolid *siinactive_box = new G4SubtractionSolid(boost::str(boost::format("siinactive_box_%d_%d") % sphxlayer % itype).c_str(), 
							    sifull_box, siactive_box, 0, G4ThreeVector(0, 0, 0));
	  G4LogicalVolume *siinactive_volume = new G4LogicalVolume(siinactive_box, G4Material::GetMaterial("G4_Si"), 
								   boost::str(boost::format("siinactive_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(siinactive_volume);
	    }
	  G4VisAttributes *siinactive_vis = new G4VisAttributes();
	  siinactive_vis->SetVisibility(true);
	  siinactive_vis->SetForceSolid(true);
	  siinactive_vis->SetColour(G4Colour::Gray());
	  siinactive_volume->SetVisAttributes(siinactive_vis);
      
	  // Make the HDI volume
	  // This makes a HDI volume that matches this sensor in Z length
	  const G4double hdi_z = sifull_z + params->get_double_param("hdi_edge_z") * cm;  
	  hdi_z_[ilayer][itype] = hdi_z;

	  G4VSolid *hdi_box = new G4Box(boost::str(boost::format("hdi_box_%d_%d") % sphxlayer % itype).c_str(), hdi_x / 2., hdi_y / 2., hdi_z/2.0);
	  G4LogicalVolume *hdi_volume = new G4LogicalVolume(hdi_box, G4Material::GetMaterial("FPC"), 
							    boost::str(boost::format("hdi_box_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(hdi_volume);
	    }

	  // This is the part of the HDI that extends beyond the sensor
	  const G4double hdiext_z = (itype == 0) ? 0.000001 : halfladder_z  - hdi_z_[ilayer][0] - hdi_z;  // need to assign nonzero value for itype=0
	  G4VSolid *hdiext_box = new G4Box(boost::str(boost::format("hdiext_box_%d_%s") % sphxlayer % itype).c_str(), hdi_x / 2., hdi_y / 2., hdiext_z / 2.0);
	  G4LogicalVolume *hdiext_volume = new G4LogicalVolume(hdiext_box, G4Material::GetMaterial("FPC"), 
							       boost::str(boost::format("hdiext_box_%d_%s") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(hdiext_volume);
	    }
	  G4VisAttributes *hdi_vis = new G4VisAttributes();
	  hdi_vis->SetVisibility(true);
	  hdi_vis->SetForceSolid(true);
	  hdi_vis->SetColour(G4Colour::Magenta());
	  hdi_volume->SetVisAttributes(hdi_vis);
	  hdiext_volume->SetVisAttributes(hdi_vis);

	  // FPHX
        
	  G4VSolid *fphx_box = new G4Box(boost::str(boost::format("fphx_box_%d_%d") % sphxlayer % itype).c_str(), fphx_x / 2., fphx_y / 2., fphx_z / 2.);
	  G4LogicalVolume *fphx_volume = new G4LogicalVolume(fphx_box, G4Material::GetMaterial("G4_Si"), 
							     boost::str(boost::format("fphx_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(fphx_volume);
	    }
	  G4VisAttributes *fphx_vis = new G4VisAttributes();
	  fphx_vis->SetVisibility(true);
	  fphx_vis->SetForceSolid(true);
	  fphx_vis->SetColour(G4Colour::Gray());
	  fphx_volume->SetVisAttributes(fphx_vis);

	  const double gap_sensor_fphx = params->get_double_param("gap_sensor_fphx") * cm;
      
	  //  FPHX Container
	  // make a container for the FPHX chips needed for this sensor, and  then place them in the container
	  const double fphxcontainer_x = fphx_x / 2.;
	  const double fphxcontainer_y = fphx_y / 2.;
	  const double fphxcontainer_z = hdi_z / 2.0;

	  G4VSolid *fphxcontainer_box = new G4Box(boost::str(boost::format("fphxcontainer_box_%d_%d") % sphxlayer % itype).c_str(), 
						  fphxcontainer_x, fphxcontainer_y, fphxcontainer_z);
	  G4LogicalVolume *fphxcontainer_volume = new G4LogicalVolume(fphxcontainer_box, G4Material::GetMaterial("G4_AIR"), 
								      boost::str(boost::format("fphxcontainer_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(fphxcontainer_volume);
	    }
	  G4VisAttributes *fphxcontainer_vis = new G4VisAttributes();
	  fphxcontainer_vis->SetVisibility(false);
	  fphxcontainer_vis->SetForceSolid(false);
	  fphxcontainer_volume->SetVisAttributes(fphxcontainer_vis);

	  // Install multiple FPHX volumes in the FPHX container volume 
	  // one FPHX chip per cell - each cell is 128 channels
	  const double offsetx = 0.;
	  const double offsety = 0.;
	  int ncopy;
	  double offsetz, cell_length_z;
	  // For layer 0, we have 5 cells per sensor, but the strips are vertical, so we have to treat it specially
	  if(nstrips_z_sensor == 1)
	    {
	      ncopy = nstrips_z_sensor / 128.0;
	    }
	  else
	    {
	      ncopy = nstrips_z_sensor;
	    }
	  cell_length_z = strip_z * nstrips_z_sensor / ncopy;
	  offsetz = (ncopy % 2 == 0) ? -2. * cell_length_z / 2. * double(ncopy / 2) + cell_length_z / 2. : -2. * cell_length_z / 2. * double(ncopy / 2);

	  G4VPVParameterisation *fphxparam = new PHG4SiliconTrackerFPHXParameterisation(offsetx, +offsety, offsetz, 2. * strip_z / 2., ncopy);
	  new G4PVParameterised(boost::str(boost::format("fphxcontainer_%d_%d") % sphxlayer % itype).c_str(), 
				fphx_volume, fphxcontainer_volume, kZAxis, ncopy, fphxparam, OverlapCheck());

	  // PGS   - this is the carbon sheet that the HDI sits on. It forms the wall of the cooling tube that cools the HDI
        
	  const double pgs_y = hdi_y;
	  const double pgs_z = hdi_z;

	  G4VSolid *pgs_box = new G4Box(boost::str(boost::format("pgs_box_%d_%d") % sphxlayer % itype).c_str(), pgs_x / 2., pgs_y / 2., pgs_z / 2.);
	  G4LogicalVolume *pgs_volume = new G4LogicalVolume(pgs_box, G4Material::GetMaterial("G4_C"), 
							    boost::str(boost::format("pgs_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(pgs_volume);
	    }
	  // The part that extends beyond this sensor, see above for hdiext
	  G4VSolid *pgsext_box = new G4Box(boost::str(boost::format("pgsext_box_%d_%s") % sphxlayer % itype).c_str(), pgs_x / 2., pgs_y / 2., hdiext_z / 2.);
	  G4LogicalVolume *pgsext_volume = new G4LogicalVolume(pgsext_box, G4Material::GetMaterial("G4_C"), 
							       boost::str(boost::format("pgsext_volume_%d_%s") % sphxlayer % itype).c_str(), 0, 0, 0);

	  G4VisAttributes *pgs_vis = new G4VisAttributes();
	  pgs_vis->SetVisibility(true);
	  pgs_vis->SetForceSolid(true);
	  pgs_vis->SetColour(G4Colour::Yellow());
	  pgs_volume->SetVisAttributes(pgs_vis);
	  pgsext_volume->SetVisAttributes(pgs_vis);

	  // Carbon stave. This is the formed sheet that sits on the PGS and completes the cooling tube
	  // Formed from straight sections and sections of a tube of radius 2.3 mm. All have wall thickness of 0.3 mm.
	  // These are different for laddertype 0 and 1, but they use some common elements.

	  // The curved section is made from a G4Cons, which is a generalized section of a cone
	  // Two curved sections combined should move the inner wall to be 2.0 mm away from the PGS, then 2 more sections bring it back
	  // For 1 of these tube sections,starting at 90 degrees, take decrease in y=1 mm at avge R. 
	  // If avge R = 2.15 mm, dtheta = invcos( (R - y)/R ) = invcos(1.15/2.15) = 53.49 deg
	  // The extent along the x axis is then R*sin(dtheta) = 1.728 mm, so two sections combined have dx = 3.456 mm length in y 
	  const double Rcmin = 0.20; //cm  inner radius of curved section, same at both ends
	  const double Rcmax = 0.23; //cm  outer radius of curved section, same at both ends
	  double Rcavge = (Rcmax+Rcmin)/2.0;
	  double dphic = TMath::ACos( (Rcavge-1.0) / Rcavge);
	  const double phic_begin = TMath::Pi() / 2.;
	  const double stave_z = pgs_z;

	  // make a single curved section and reuse it later
	  G4Cons *stave_curve_cons = new G4Cons(boost::str(boost::format("stave_curve_%d_%d") %sphxlayer % itype).c_str(), 
						Rcmin, Rcmax, Rcmin, Rcmax, stave_z / 2., phic_begin, dphic);
	  G4LogicalVolume *stave_curve_volume = new G4LogicalVolume(stave_curve_cons, G4Material::GetMaterial("G4_C"), 
								    boost::str(boost::format("stave_curve_volume_%d_%s") % sphxlayer % itype).c_str(), 0, 0, 0);

	  G4Cons *stave_curve_ext_cons = new G4Cons(boost::str(boost::format("stave_curve_%d_%d") %sphxlayer % itype).c_str(), 
						    Rcmin, Rcmax, Rcmin, Rcmax, hdiext_z / 2., phic_begin, dphic);
	  G4LogicalVolume *stave_curve_ext_volume = new G4LogicalVolume(stave_curve_ext_cons, G4Material::GetMaterial("G4_C"), 
									boost::str(boost::format("stave_curve_ext_volume_%d_%s") % sphxlayer % itype).c_str(), 0, 0, 0);

	  // Make several straight sections for use in making the stave

	  // Outer straight sections of stave
	  double stave_wall_thickness = 0.3;
	  G4VSolid *stave_straight_outer_box = new G4Box(boost::str(boost::format("stave_straight_outer_box_%d_%d") % sphxlayer % itype).c_str(), 
							 stave_wall_thickness / 2., stave_straight_outer_y / 2., stave_z / 2.);
	  G4LogicalVolume *stave_straight_outer_volume = new G4LogicalVolume(stave_straight_outer_box, G4Material::GetMaterial("G4_C"), 
									     boost::str(boost::format("stave_straight_outer_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  G4VSolid *stave_straight_outer_ext_box = new G4Box(boost::str(boost::format("stave_straight_outer_ext_box_%d_%s") % sphxlayer % itype).c_str(), 
							     stave_wall_thickness/2., stave_straight_outer_y/2., hdiext_z / 2.);
	  G4LogicalVolume *stave_straight_outer_ext_volume = new G4LogicalVolume(stave_straight_outer_ext_box, G4Material::GetMaterial("G4_C"), 
										 boost::str(boost::format("stave_straight_outer_ext_volume_%d_%s") % sphxlayer % itype).c_str(), 0, 0, 0);

	  // connects cooling tubes together, only needed for laddertype 1
	  G4VSolid *stave_straight_inner_box = new G4Box(boost::str(boost::format("stave_straight_inner_box_%d_%d") % sphxlayer % itype).c_str(), 
							 stave_wall_thickness / 2., stave_straight_inner_y / 2., stave_z / 2.);
	  G4LogicalVolume *stave_straight_inner_volume = new G4LogicalVolume(stave_straight_inner_box, G4Material::GetMaterial("G4_C"), 
									     boost::str(boost::format("stave_straight_inner_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  G4VSolid *stave_straight_inner_ext_box = new G4Box(boost::str(boost::format("stave_straight_inner_ext_box_%d_%d") % sphxlayer % itype).c_str(), 
							     stave_wall_thickness / 2., stave_straight_inner_y / 2., hdiext_z / 2.);
	  G4LogicalVolume *stave_straight_inner_ext_volume = new G4LogicalVolume(stave_straight_inner_ext_box, G4Material::GetMaterial("G4_C"), 
										 boost::str(boost::format("stave_straight_inner_ext_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);

	  //Top surface of cooler tube
	  G4VSolid *stave_straight_cooler_box = new G4Box(boost::str(boost::format("stave_straight_cooler_box_%d_%d") % sphxlayer % itype).c_str(), 
							  stave_wall_thickness / 2., stave_straight_cooler_y / 2., stave_z / 2.);
	  G4LogicalVolume *stave_straight_cooler_volume = new G4LogicalVolume(stave_straight_cooler_box, G4Material::GetMaterial("G4_C"), 
									      boost::str(boost::format("stave_straight_cooler_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  G4VSolid *stave_straight_cooler_ext_box = new G4Box(boost::str(boost::format("stave_straight_cooler_ext_box_%d_%d") % sphxlayer % itype).c_str(), 
							      stave_wall_thickness / 2., stave_straight_cooler_y / 2., hdiext_z / 2.);
	  G4LogicalVolume *stave_straight_cooler_ext_volume = new G4LogicalVolume(stave_straight_cooler_ext_box, G4Material::GetMaterial("G4_C"), 
										  boost::str(boost::format("stave_straight_cooler_ext_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);

	  G4VisAttributes *stave_curve_vis = new G4VisAttributes();
	  stave_curve_vis->SetVisibility(true);
	  stave_curve_vis->SetForceSolid(true);
	  stave_curve_vis->SetColour(G4Colour::Yellow());
	  stave_curve_volume->SetVisAttributes(stave_curve_vis);
	  stave_curve_ext_volume->SetVisAttributes(stave_curve_vis);
	  stave_straight_cooler_volume->SetVisAttributes(stave_curve_vis);
	  stave_straight_cooler_ext_volume->SetVisAttributes(stave_curve_vis);
	  stave_straight_inner_volume->SetVisAttributes(stave_curve_vis);
	  stave_straight_inner_ext_volume->SetVisAttributes(stave_curve_vis);
	  stave_straight_outer_volume->SetVisAttributes(stave_curve_vis);
	  stave_straight_outer_ext_volume->SetVisAttributes(stave_curve_vis);

	  // Now we combine the elements of a stave defined above into a stave
	  // Create a stave volume to install the stave sections into. The volume has to be big enouigh to contain the cooling tube
	  double cooler_gap_x = 0.2;  // id of cooling tube in cm
	  double cooler_wall = 0.03;   // outer wall thickness of cooling tube 
	  double cooler_x = cooler_gap_x + cooler_wall; 
	  const double stave_x = cooler_x;  // we do not include the PGS in the stave volume
	  const double stave_y = hdi_y;
	  G4VSolid *stave_box = new G4Box(boost::str(boost::format("stave_box_%d_%d") % sphxlayer % itype).c_str(), stave_x / 2., stave_y / 2., stave_z / 2.);
	  G4LogicalVolume *stave_volume = new G4LogicalVolume(stave_box, G4Material::GetMaterial("G4_Air"), 
							      boost::str(boost::format("stave_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(stave_volume);
	    }
	  G4VSolid *staveext_box = new G4Box(boost::str(boost::format("staveext_box_%d_%s") % sphxlayer % itype).c_str(), stave_x / 2., stave_y / 2., hdiext_z / 2.);
	  G4LogicalVolume *staveext_volume = new G4LogicalVolume(staveext_box, G4Material::GetMaterial("G4_Air"), 
								 boost::str(boost::format("staveext_volume_%d_%s") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(staveext_volume);
	    }

	  // Assemble the elements into the stave volume and the stave extension volume
	  if(laddertype == 0)
	    {
	      // only one cooling tube in laddertype 0
	      // Place the straight sections. We add the middle, then above x axis, then below x axis
	      double x_off_str[3] = {Rcavge, (Rcmax - Rcmin) / 2., (Rcmax - Rcmin) / 2.};
	      double y_off_str[3] = 
		{
		  // Rcavge*TMath::Sin(dphic / 2.) is half the width of a curved section in y
		  0.0,
		  stave_straight_cooler_y / 2. + Rcavge*TMath::Sin( dphic / 2. ) + 2.0 * Rcavge*TMath::Sin(dphic / 2.) + stave_straight_outer_y / 2.0,
		  -stave_straight_cooler_y / 2. - Rcavge*TMath::Sin( dphic / 2. ) - 2.0 * Rcavge*TMath::Sin(dphic / 2.) - stave_straight_outer_y / 2.0	  
		};

   	      for(int i=0;i<3;i++)
		{	  	  
		  if(i==0)
		    {
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_volume, 
					boost::str(boost::format("stave_straight_cooler_%d_%d_%d") % i % sphxlayer % itype).c_str(), stave_volume, false, 0, OverlapCheck());
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_ext_volume, 
					boost::str(boost::format("stave_straight_cooler_ext_%d_%d_%s") % i % sphxlayer % itype).c_str(), staveext_volume, false, 0, OverlapCheck());
		    }
		  else
		    {
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_volume, 
					boost::str(boost::format("stave_straight_outer_%d_%d_%d") % i % sphxlayer % itype).c_str(), stave_volume, false, 0, OverlapCheck());
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_ext_volume, 
					boost::str(boost::format("stave_straight_outer_ext_%d_%d_%s") % i % sphxlayer % itype).c_str(), staveext_volume, false, 0, OverlapCheck());
		    }
		}
	      // The cooler curved sections are made using 2 curved sections in a recurve on each side of the cooler straight section 
	      double frot[4] = {-TMath::Pi()/2., -TMath::Pi()/2.,  TMath::Pi() / 2. + dphic, TMath::Pi()/2.};
	      double x_off_cooler[4] = 
		{
		  Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge - Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge - Rcavge*TMath::Cos(dphic / 2.) 
		};
	      double y_off_cooler[4] = 
		{
		  stave_straight_cooler_y / 2. + Rcavge*TMath::Sin( dphic / 2. ),  
		  - stave_straight_cooler_y / 2. - Rcavge*TMath::Sin( dphic / 2. ), 
		  stave_straight_cooler_y / 2. + Rcavge*TMath::Sin( dphic / 2. ) + 2.0 * Rcavge*TMath::Sin(dphic / 2.),  
		  - stave_straight_cooler_y / 2. - Rcavge*TMath::Sin( dphic / 2. ) - 2.0 * Rcavge*TMath::Sin(dphic / 2.)
		};
	      
	      for(int icurve=0;icurve<4;icurve++)
		{
		  G4RotationMatrix *rot = new G4RotationMatrix();
		  rot->rotateZ(frot[icurve]);
		  new G4PVPlacement(rot, G4ThreeVector(x_off_cooler[icurve], y_off_cooler[icurve], 0.0), stave_curve_volume, 
				    boost::str(boost::format("stave_curve_%d_%d_%d") % icurve % sphxlayer % itype).c_str(), stave_volume, false, 0, OverlapCheck());
		  new G4PVPlacement(rot, G4ThreeVector(x_off_cooler[icurve], y_off_cooler[icurve], 0.0), stave_curve_ext_volume, 
				    boost::str(boost::format("stave_curve_ext_%d_%d_%s") % icurve % sphxlayer % itype).c_str(), staveext_volume, false, 0, OverlapCheck());
		}
	    }
	  else
	    {
	      // The type 1 ladder has two cooling tubes
	      // Place the straight sections, do the extension at the same time
	      // we alternate  positive and negative y values here
	      double x_off_str[5] = 
		{
		  (Rcmax-Rcmin) / 2., // against the PGS
		  (Rcmax+Rcmin) / 2.,   // top of cooler
		  (Rcmax+Rcmin) / 2.,   // top of cooler
		  (Rcmax-Rcmin) / 2., 
		  (Rcmax-Rcmin) / 2.
		};
	      double y_off_str[5] = 
		{
		  0.0,
		  stave_straight_inner_y/2. + 2. * Rcavge*TMath::Sin( dphic / 2. ) + stave_straight_cooler_y/2.,  // top of cooler
		  - stave_straight_inner_y/2. - 2. * Rcavge*TMath::Sin( dphic / 2. ) - stave_straight_cooler_y/2.,  // top of cooler
		  stave_straight_inner_y/2. + 2. * Rcavge*TMath::Sin( dphic / 2.) + stave_straight_cooler_y + 2. *  Rcavge*TMath::Sin( dphic / 2. ) + stave_straight_outer_y/2.,
		  - stave_straight_inner_y/2. - 2. * Rcavge*TMath::Sin( dphic / 2.) - stave_straight_cooler_y - 2. *  Rcavge*TMath::Sin( dphic / 2. ) - stave_straight_outer_y/2.,
		};
	      
	      for(int i=0;i<5;i++)
		{
		  if(i==0)
		    {
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_inner_volume, 
					boost::str(boost::format("stave_straight_inner_%d_%d_%d") % i % sphxlayer % itype).c_str(), stave_volume, false, 0, OverlapCheck());
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_inner_ext_volume, 
					boost::str(boost::format("stave_straight_inner_ext_%d_%d_%s") % i % sphxlayer % itype).c_str(), staveext_volume, false, 0, OverlapCheck());
		    }
		  else if(i==1 || i==2)
		    {
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_volume, 
					boost::str(boost::format("stave_straight_cooler_%d_%d_%d") % i % sphxlayer % itype).c_str(), stave_volume, false, 0, OverlapCheck());
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_ext_volume, 
					boost::str(boost::format("stave_straight_cooler_ext_%d_%d_%s") % i % sphxlayer % itype).c_str(), staveext_volume, false, 0, OverlapCheck());
		    }
		  else
		    {
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_volume, 
					boost::str(boost::format("stave_straight_outer_%d_%d_%d") % i % sphxlayer % itype).c_str(), stave_volume, false, 0, OverlapCheck());
		      new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_ext_volume, 
					boost::str(boost::format("stave_straight_outer_ext_%d_%d_%s") % i % sphxlayer % itype).c_str(), staveext_volume, false, 0, OverlapCheck());
		    }
		}
	  
	      // Place the curved sections
	      // hre we do all above the x axis, then all below the x axis, each in order of increasing y
	      double frot[8] = 
		{
		  TMath::Pi()/2., 
		  -TMath::Pi()/2.,  
		  - TMath::Pi()/2. + dphic, 
		  TMath::Pi(),
		  TMath::Pi()/2., 
		  -TMath::Pi()/2.,  
		  - TMath::Pi()/2. + dphic, 
		  TMath::Pi()
		};
	      double x_off_curve[8] = 
		{
		  Rcavge - Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge - Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge - Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge*TMath::Cos(dphic / 2.), 
		  Rcavge - Rcavge*TMath::Cos(dphic / 2.) 
		};
	      double y_off_curve[8] = 
		{
		  stave_straight_inner_y / 2. + Rcavge*TMath::Sin( dphic / 2.),  
		  stave_straight_inner_y / 2. + 2. * Rcavge*TMath::Sin( dphic / 2. ) + Rcavge*TMath::Sin(dphic / 2.),  
		  stave_straight_inner_y / 2. + 4. * Rcavge*TMath::Sin( dphic / 2. ) + stave_straight_cooler_y + Rcavge*TMath::Sin(dphic / 2.),  
		  stave_straight_inner_y / 2. + 6. * Rcavge*TMath::Sin( dphic / 2. ) + stave_straight_cooler_y + Rcavge*TMath::Sin(dphic / 2.),  
		  - stave_straight_inner_y / 2. - 6. * Rcavge*TMath::Sin( dphic / 2. ) - stave_straight_cooler_y - Rcavge*TMath::Sin(dphic / 2.),  
		  - stave_straight_inner_y / 2. - 4. * Rcavge*TMath::Sin( dphic / 2. ) - stave_straight_cooler_y - Rcavge*TMath::Sin(dphic / 2.),  
		  - stave_straight_inner_y / 2. - 2. * Rcavge*TMath::Sin( dphic / 2. ) - Rcavge*TMath::Sin(dphic / 2.),  
		  - stave_straight_inner_y / 2. - Rcavge*TMath::Sin( dphic / 2.)  
		};
	  
	      for(int i=0;i<8;i++)
		{
		  G4RotationMatrix *rot = new G4RotationMatrix();
		  rot->rotateZ(frot[i]);
		  new G4PVPlacement(rot, G4ThreeVector(x_off_curve[i], y_off_curve[i], 0.0), stave_curve_volume, boost::str(boost::format("stave_curve_%d_%d_%d") 
															  % i % sphxlayer % itype).c_str(), stave_volume, false, 0, OverlapCheck());
		  new G4PVPlacement(rot, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_curve_ext_volume, boost::str(boost::format("stave_curve_ext_%d_%d_%s") 
															 % i % sphxlayer % itype).c_str(), staveext_volume, false, 0, OverlapCheck());
		}
	    }
      

	  // ----- Step 2 ------------------------------------------------------------------------------------
	  // We place Si-sensor, FPHX, HDI, PGS sheet, and stave in the ladder  volume.
	  // ======================================================

	  // Make the ladder volume first
	  // We are still in the loop over inner or outer sensors. This is the ladder volume corresponding to this sensor. The FPHX is taller than the sensor in x.
	  const double ladder_x = stave_x  + pgs_x  + hdi_x  + fphx_x;
	  const double ladder_y = hdi_y;
	  const double ladder_z = hdi_z;
	  G4VSolid *ladder_box = new G4Box(boost::str(boost::format("ladder_box_%d_%d") % sphxlayer % itype).c_str(), ladder_x / 2., ladder_y / 2., ladder_z / 2.);
	  G4LogicalVolume *ladder_volume = new G4LogicalVolume(ladder_box, G4Material::GetMaterial("G4_AIR"), boost::str(boost::format("ladder_volume_%d_%d") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(ladder_volume);
	    }
	  G4VSolid *ladderext_box = new G4Box(boost::str(boost::format("ladderext_box_%d_%s") % sphxlayer % itype).c_str(), ladder_x / 2., ladder_y / 2., hdiext_z / 2.);
	  G4LogicalVolume *ladderext_volume = new G4LogicalVolume(ladderext_box, G4Material::GetMaterial("G4_AIR"), boost::str(boost::format("ladderext_volume_%d_%s") % sphxlayer % itype).c_str(), 0, 0, 0);
	  if ((IsAbsorberActive.find(inttlayer))->second > 0)
	    {
	      absorberlogvols.insert(ladderext_volume);
	    }
	  G4VisAttributes *ladder_vis = new G4VisAttributes();
	  ladder_vis->SetVisibility(true);
	  ladder_vis->SetForceSolid(true);
	  ladder_vis->SetColour(G4Colour::Cyan());
	  ladder_volume->SetVisAttributes(ladder_vis);
	  ladderext_volume->SetVisAttributes(ladder_vis);

	  // Now add the components of the ladder to the ladder volume
	  // The sensor is closest to the beam pipe, the stave cooler is furthest away
	  // Note that the cooler has been assembled in the stave volume with the top at larger x, so the sensor will be at smaller x
	  // That will be the configuration when the ladder is at phi = 180 degrees, the positive x direction
	  // So we start at the most positive x value and add the stave first

	  // Carbon stave        
	  const double TVstave_x = ladder_x / 2. - stave_x / 2.;
	  new G4PVPlacement(0, G4ThreeVector(TVstave_x, 0.0, 0.0), stave_volume, boost::str(boost::format("stave_%d_%d") % sphxlayer % itype).c_str(), 
			    ladder_volume, false, 0, OverlapCheck());
	  new G4PVPlacement(0, G4ThreeVector(TVstave_x, 0.0, 0.0), staveext_volume, boost::str(boost::format("staveext_%d_%s") % sphxlayer % itype).c_str(), 
			    ladderext_volume, false, 0, OverlapCheck());

	  // PGS
	  const double TVpgs_x = TVstave_x - stave_x / 2. - pgs_x / 2.;
	  new G4PVPlacement(0, G4ThreeVector(TVpgs_x, 0.0, 0.0), pgs_volume, boost::str(boost::format("pgs_%d_%d") % sphxlayer % itype).c_str(), 
			    ladder_volume, false, 0, OverlapCheck());
	  new G4PVPlacement(0, G4ThreeVector(TVpgs_x, 0.0, 0.0), pgsext_volume, boost::str(boost::format("pgsext_%d_%s") % sphxlayer % itype).c_str(), 
			    ladderext_volume, false, 0, OverlapCheck());

	  // HDI        
	  const double TVhdi_x = TVpgs_x - pgs_x / 2. - hdi_x / 2.;
	  new G4PVPlacement(0, G4ThreeVector(TVhdi_x, 0.0, 0.0), hdi_volume, boost::str(boost::format("hdi_%d_%d") % sphxlayer % itype).c_str(), ladder_volume, false, 0, OverlapCheck());
	  new G4PVPlacement(0, G4ThreeVector(TVhdi_x, 0.0, 0.0), hdiext_volume, boost::str(boost::format("hdiext_%d_%s") % sphxlayer % itype).c_str(), ladderext_volume, false, 0, OverlapCheck());

	  // Si-sensor        
	  const double TVSi_x = TVhdi_x - hdi_x / 2. - siactive_x / 2.;
	  new G4PVPlacement(0, G4ThreeVector(TVSi_x, 0.0, 0.0), siinactive_volume, boost::str(boost::format("siinactive_%d_%d") % sphxlayer % itype).c_str(), ladder_volume, false, 0, OverlapCheck());
	  new G4PVPlacement(0, G4ThreeVector(TVSi_x, 0.0, 0.0), siactive_volume, boost::str(boost::format("siactive_%d_%d") % sphxlayer % itype).c_str(), ladder_volume, false, 0, OverlapCheck());

	  // FPHX        
	  const double TVfphx_x = TVhdi_x - hdi_x / 2. - fphx_x / 2.;
	  const double TVfphx_y = sifull_y + gap_sensor_fphx + fphx_y / 2.;
	  new G4PVPlacement(0, G4ThreeVector(TVfphx_x, +TVfphx_y, 0.0), fphxcontainer_volume, boost::str(boost::format("fphxcontainerp_%d_%d") % sphxlayer % itype).c_str(), ladder_volume, false, 0, OverlapCheck());
	  new G4PVPlacement(0, G4ThreeVector(TVfphx_x, -TVfphx_y, 0.0), fphxcontainer_volume, boost::str(boost::format("fphxcontainerm_%d_%d") % sphxlayer % itype).c_str(), ladder_volume, false, 0, OverlapCheck());

	  // ----- Step 3 --------------------------------------------------------------------------------------------------------------------
	  // We install the section of ladder for this sensor at all requested phi values and at positive and negative Z
	  //========================================================================

	  // Distribute Ladders in phi
	  // We are still in the loops over layer and sensor type, we will place copies of the ladder section for this sensor
	  //  at all ladder phi values, and at both positive and negative Z values.

	  // The ladders have no tilt in the new design, they are alternately offset in radius to get overlap in phi        
	  // radius values are for the center of the sensor, we need the offset from center of ladder to center of sensor
	  double sensor_offset = 0.0 - TVSi_x; // ladder center is at 0.0 by construction

	  const double dphi = 2 * TMath::Pi() / nladders_layer;

	  // there is no single radius for a layer
	  eff_radius[ilayer] = layer_radius_inner + sensor_offset;
	  eff_radius_alternate[ilayer] = layer_radius_outer + sensor_offset;
	  posz[ilayer][itype] = (itype == 0) ? hdi_z : hdi_z_[ilayer][0] + hdi_z;
	  strip_x_offset[ilayer] = sensor_offset;

	  for (G4int icopy = 0; icopy < nladders_layer; icopy++)
	    {
	      const double phi = offsetphi + dphi * (double) icopy;  // offsetphi is zero by default, so we start at zero
	      double radius = eff_radius[ilayer];
	      if(icopy%2)
		radius = eff_radius_alternate[ilayer];  // every odd numbered copy is placed at the larger radius
	      const double posx = radius * cos(phi);
	      const double posy = radius * sin(phi);
	      const double fRotate = phi + offsetrot; // no initial rotation, since we assembled the ladder in phi = 0 orientation 

	      G4RotationMatrix *ladderrotation = new G4RotationMatrix();
	      ladderrotation->rotateZ(-fRotate);

	      // place the copy at its ladder phi value, and at positive and negative Z
	      new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, -posz[ilayer][itype]), ladder_volume, boost::str(boost::format("ladder_%d_%d_%d_%d") 
															   % sphxlayer % inttlayer % itype % icopy).c_str(), trackerenvelope, false, 0, OverlapCheck());
	      new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, +posz[ilayer][itype]), ladder_volume, boost::str(boost::format("ladder_%d_%d_%d_%d") 
															   % sphxlayer % inttlayer % itype % icopy).c_str(), trackerenvelope, false, 0, OverlapCheck());
	
	      if (itype != 0)
		{  
		  // The HDI extension tab has to be added for the outer sensor
		  const G4double posz_ext = 2. * (hdi_z_[ilayer][0] + hdi_z) + hdiext_z;
		  new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, -posz_ext), ladderext_volume, boost::str(boost::format("ladderext_%d_%d_%d_%d") 
														       % sphxlayer % inttlayer % itype % icopy).c_str(), trackerenvelope, false, 0, OverlapCheck());
		  new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, +posz_ext), ladderext_volume, boost::str(boost::format("ladderext_%d_%d_%d_%d") 
														       % sphxlayer % inttlayer % itype % icopy).c_str(), trackerenvelope, false, 0, OverlapCheck());
		  cout <<"            posz_ext " << posz_ext << endl;
		}

	      cout << "Ladder copy " << icopy << " radius " << radius << " phi " << phi << " itype " << itype << " posz " << posz[ilayer][itype] 
		   << " fRotate " << fRotate << " posx " << posx << " posy " << posy << " sensor_offset " << sensor_offset 
		   << endl;

	    } // end loop over ladder copy placement in phi and positive and negative Z
	} // end loop over inner or outer sensor
    } // end loop over layers
  
  return 0;
}

int PHG4SiliconTrackerDetector::DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm)
{
  static int i = 0;
  G4LogicalVolume *checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  switch (i)
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
    i = 0;
    break;
  }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm, G4ThreeVector(0, 0, 0), checksolid, "DISPLAYVOL", logvol, 0, false, OverlapCheck());
  return 0;
}

void PHG4SiliconTrackerDetector::AddGeometryNode()
{
  int active = 0;
  map<int, int>::const_iterator iter;
  for (iter = IsActive.begin(); iter != IsActive.end(); ++iter)
  {
    if (iter->second > 0)
    {
      active = 1;
      break;
    }
  }
  if (active)
  {
    std::string geonode = (superdetector != "NONE") ? boost::str(boost::format("CYLINDERGEOM_%s") % superdetector) : boost::str(boost::format("CYLINDERGEOM_%s") % detector_type);

    PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode.c_str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.c_str(), "PHObject");
      runNode->addNode(newNode);
    }

    for (unsigned int ilayer = 0; ilayer < nlayer_; ++ilayer)
    {
      const int sphxlayer = layerconfig_[ilayer].first;
      const int inttlayer = layerconfig_[ilayer].second;

      const PHParameters *params = paramscontainer->GetParameters(inttlayer);
      // parameters are in cm, so no conversion needed here to get to cm (*cm/cm)
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeom_Siladders(
          sphxlayer,
          params->get_double_param("strip_x"),
          params->get_double_param("strip_y"),
          params->get_double_param("strip_z_0"),
          params->get_double_param("strip_z_1"),
          params->get_int_param("nstrips_z_sensor_0"),
          params->get_int_param("nstrips_z_sensor_1"),
          params->get_int_param("nstrips_phi_sensor"),
          params->get_int_param("nladder"),
          posz[ilayer][0] / cm,
          posz[ilayer][1] / cm,
          eff_radius[ilayer] / cm,
          eff_radius_alternate[ilayer] / cm,
          strip_x_offset[ilayer] / cm,
          params->get_double_param("offsetphi") * deg / rad,
          params->get_double_param("offsetrot") * deg / rad);
      geo->AddLayerGeom(sphxlayer, mygeom);
      if (Verbosity() > 0)
        geo->identify();
    }
  }
}
