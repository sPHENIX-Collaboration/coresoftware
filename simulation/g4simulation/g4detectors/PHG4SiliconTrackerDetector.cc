#include "PHG4SiliconTrackerDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Siladders.h"
#include "PHG4SiliconTrackerParameterisation.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4PVParameterised.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>

#include <boost/format.hpp>
#include <boost/foreach.hpp>

using namespace std;

PHG4SiliconTrackerDetector::PHG4SiliconTrackerDetector(PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam, const vpair &layerconfig):
    PHG4Detector(Node, dnam),
    params(parameters),
    active(params->get_int_param("active")),
    absorberactive(params->get_int_param("absorberactive")),
    blackhole(params->get_int_param("blackhole"))
{
  layerconfig_ = layerconfig;

  layermin_ = layerconfig_.front().first;
  layermax_ = layerconfig_.back().first;
  nlayer_   = layerconfig_.size();

  for (unsigned int ilayer=0; ilayer<nlayer_; ++ilayer)
    {
      const int inttlayer = layerconfig_[ilayer].second;
      std::cout << " PHG4SiliconTrackerDetector constructor: ilayer " << ilayer << " inttlayer " << inttlayer << std::endl;
      if (inttlayer < 0 || inttlayer >= 5)
        assert(!"PHG4SiliconTrackerDetector: check INTT ladder layer.");
    }
}

PHG4SiliconTrackerDetector::~PHG4SiliconTrackerDetector()
{}

//_______________________________________________________________
//_______________________________________________________________
int PHG4SiliconTrackerDetector::IsInSiliconTracker(G4VPhysicalVolume * volume) const
  {
    // Is this volume one of the sensor strips?
    //   cout << "volume name " << volume->GetName() << endl;
     G4LogicalVolume *logvol = volume->GetLogicalVolume();
    //   cout << "logical volume: " << logvol->GetName() << endl;
    //   cout << "pointer: " << hex << logvol << dec << endl;
    if (absorberactive)
      {
	int foundit = 0;
	if (absorberlogvols.find(logvol) != absorberlogvols.end())
	  {
	    cout << "found absorber logvol" << logvol->GetName() << endl;
	    foundit = 1;
	  }
	// inactive strip strip and other parts
        if ((volume->GetName().find("siinactive") != std::string::npos) ||
	    (volume->GetName().find("hdi")        != std::string::npos) ||
	    (volume->GetName().find("fphx")       != std::string::npos) ||
	    (volume->GetName().find("pgs")        != std::string::npos) ||
	    (volume->GetName().find("stave")      != std::string::npos) ||
	    (volume->GetName().find("ladder")     != std::string::npos))
          {
            return -1;
          }
	if (foundit)
	  {
	    cout << "absorber logvol not assigned" << endl;
	  }
      }
    if (active)
      {
	int foundit = 0;
	if (activelogvols.find(logvol) != activelogvols.end())
	  {
	    foundit = 1;
	    return 1;
	    //	    cout << "found active logvol " << logvol->GetName() << endl;
	  }
	// active strip strip
	if (volume->GetName().find("siactive") != std::string::npos)
	  {
	if (activelogvols.find(logvol) == activelogvols.end())
	  {
	    cout << "volume name " << volume->GetName() << endl;
        G4LogicalVolume *logv = volume->GetLogicalVolume();
	cout << "logv: " << logv << " name " << logv->GetName() << endl;
	    BOOST_FOREACH(G4LogicalVolume *lvol, activelogvols)
	      {
		cout << "logvol " << lvol << " name: " << lvol->GetName() << endl;
	      }
	  }
	    return 2;
	  }
	if (foundit)
	  {
	    cout << "active logvol not assigned" << endl;
	  }
      }

    return 0;
  }

void PHG4SiliconTrackerDetector::Construct( G4LogicalVolume* logicWorld )
{
  if(verbosity>0)
    std::cout << "PHG4SiliconTrackerDetector::Construct called for layers " << layermin_ << " to " << layermax_ << std::endl;

  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  ConstructSiliconTracker(logicWorld);

  // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int PHG4SiliconTrackerDetector::ConstructSiliconTracker(G4LogicalVolume* trackerenvelope)
{
  double hdi_z_[nlayer_][2];

  for (unsigned int ilayer=0; ilayer<nlayer_; ++ilayer)
    for (int itype=0; itype<2; ++itype)
      {
        if (!(itype >= 0 && itype <= 1))
          assert(!"Error: check ladder type.");

        const int sphxlayer = layerconfig_[ilayer].first;
        const int inttlayer = layerconfig_[ilayer].second;

	std::cout << " PHG4SiliconTrackerDetector::ConstrctSiliconTracker:  sphxlayer " << sphxlayer << " inttlayer " << inttlayer << std::endl;

        const G4double strip_x = arr_strip_x[inttlayer];
        const G4double strip_y = arr_strip_y[inttlayer];
        const int nstrips_phi_cell = arr_nstrips_phi_cell[inttlayer];

        const int nladders_layer = arr_nladders_layer[inttlayer];
        const G4double radius    = arr_radius[inttlayer];
        const G4double offsetphi = arr_offsetphi[inttlayer];
        const G4double offsetrot = arr_offsetrot[inttlayer];
        const G4double hdi_y     = arr_hdi_y[inttlayer];

        const G4double strip_z     = (inttlayer==0) ?          arr_strip_z[0][itype] :          arr_strip_z[1][itype];
        const int nstrips_z_sensor = (inttlayer==0) ? arr_nstrips_z_sensor[0][itype] : arr_nstrips_z_sensor[1][itype];

        /*----- Step 1 -----
         * We make volume for Si-sensor, FPHX, HDI, PGS sheet, and stave.
         * Then we make ladder volume large enough to enclose the above volume.
         */

        /*
         * Si-strip
         */
        G4VSolid *strip_box = new G4Box(boost::str(boost::format("strip_box_%d_%d") %sphxlayer %itype).c_str(), strip_x, strip_y, strip_z);
        G4LogicalVolume *strip_volume = new G4LogicalVolume(strip_box, G4Material::GetMaterial("G4_Si"), boost::str(boost::format("strip_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	activelogvols.insert(strip_volume);
        G4VisAttributes *strip_vis = new G4VisAttributes();
        strip_vis->SetVisibility(false);
        strip_vis->SetForceSolid(false);
        strip_vis->SetColour(G4Colour::White());
        strip_volume->SetVisAttributes(strip_vis);

        /*
         * Si-sensor active area
         */
        const double siactive_x = strip_x; // 0.24mm/2
        const double siactive_y = strip_y * 2.*(double)nstrips_phi_cell; // (0.078mm * 2*128)/2 = 0.078mm * 128
        const double siactive_z = strip_z *    (double)nstrips_z_sensor; // (20mm * 5or8)/2 = 10mm * 5or8

        G4VSolid *siactive_box = new G4Box(boost::str(boost::format("siactive_box_%d_%d") %sphxlayer %itype).c_str(), siactive_x, siactive_y, siactive_z);
        G4LogicalVolume *siactive_volume = new G4LogicalVolume(siactive_box, G4Material::GetMaterial("G4_AIR"), boost::str(boost::format("siactive_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	//	activelogvols.insert(siactive_volume);

        G4VPVParameterisation *stripparam = new PHG4SiliconTrackerStripParameterisation(nstrips_phi_cell * 2, nstrips_z_sensor, strip_y * 2., strip_z * 2.);
        new G4PVParameterised(boost::str(boost::format("siactive_%d_%d") %sphxlayer %itype).c_str(), strip_volume, siactive_volume, kZAxis, nstrips_phi_cell * 2 * nstrips_z_sensor, stripparam, false); // overlap check too long.

        /*
         * Si-sensor full (active+inactive) area
         */
        const double sifull_x = siactive_x; // 0.24mm/2
        const double sifull_y = siactive_y + sensor_edge_phi; // (1.305mm  + 0.078mm * 2*128 + 1.305mm)/2 = 0.078mm * 128 + 1.305mm
        const double sifull_z = siactive_z + sensor_edge_z;   // (0.98mm + 20mm * 5 + 0.98mm)/2 = 10mm * 5 + 0.98mm

        G4VSolid *sifull_box = new G4Box(boost::str(boost::format("sifull_box_%d_%d") %sphxlayer %itype).c_str(), sifull_x, sifull_y, sifull_z);

        /*
         * Si-sensor inactive area
         */
        G4VSolid *siinactive_box = new G4SubtractionSolid(boost::str(boost::format("siinactive_box_%d_%d") %sphxlayer %itype).c_str(), sifull_box, siactive_box, 0, G4ThreeVector(0, 0, 0));
        G4LogicalVolume *siinactive_volume = new G4LogicalVolume(siinactive_box, G4Material::GetMaterial("G4_Si"), boost::str(boost::format("siinactive_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(siinactive_volume);

        G4VisAttributes *siinactive_vis = new G4VisAttributes();
        siinactive_vis->SetVisibility(true);
        siinactive_vis->SetForceSolid(true);
        siinactive_vis->SetColour(G4Colour::Gray());
        siinactive_volume->SetVisAttributes(siinactive_vis);

        /*
         * HDI
         */
        const G4double hdi_z = sifull_z + hdi_edge_z; // hdi_edge_z = 100micron
        hdi_z_[ilayer][itype] = hdi_z;

        G4VSolid *hdi_box = new G4Box(boost::str(boost::format("hdi_box_%d_%d") %sphxlayer %itype).c_str(), hdi_x, hdi_y, hdi_z);
        G4LogicalVolume *hdi_volume = new G4LogicalVolume(hdi_box, G4Material::GetMaterial("FPC"), boost::str(boost::format("hdi_box_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(hdi_volume);

        const G4double hdi_ext_z = (itype==0) ? 0.000001 : arr_halfladder_z[ilayer] - hdi_z_[ilayer][0] - hdi_z; // need to assign nonzero value for itype=0
        G4VSolid *hdi_ext_box = new G4Box(boost::str(boost::format("hdi_ext_box_%d_%s") %sphxlayer %itype).c_str(), hdi_x, hdi_y, hdi_ext_z);
        G4LogicalVolume *hdi_ext_volume = new G4LogicalVolume(hdi_ext_box, G4Material::GetMaterial("FPC"), boost::str(boost::format("hdi_ext_box_%d_%s") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(hdi_ext_volume);
        G4VisAttributes *hdi_vis = new G4VisAttributes();
        hdi_vis->SetVisibility(true);
        hdi_vis->SetForceSolid(true);
        hdi_vis->SetColour(G4Colour::Magenta());
        hdi_volume->SetVisAttributes(hdi_vis);

        /*
         * FPHX
         */
        G4VSolid *fphx_box = new G4Box(boost::str(boost::format("fphx_box_%d_%d") %sphxlayer %itype).c_str(), fphx_x, fphx_y, fphx_z);
        G4LogicalVolume *fphx_volume = new G4LogicalVolume(fphx_box, G4Material::GetMaterial("G4_Si"), boost::str(boost::format("fphx_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(fphx_volume);
        G4VisAttributes *fphx_vis = new G4VisAttributes();
        fphx_vis->SetVisibility(true);
        fphx_vis->SetForceSolid(true);
        fphx_vis->SetColour(G4Colour::Gray());
        fphx_volume->SetVisAttributes(fphx_vis);

        const double gap_sensor_fphx = 1.0*mm;

        /*
         * FPHX Container
         */
        const double fphxcontainer_x = fphx_x;
        const double fphxcontainer_y = fphx_y;
        const double fphxcontainer_z = hdi_z;

        G4VSolid *fphxcontainer_box = new G4Box(boost::str(boost::format("fphxcontainer_box_%d_%d") %sphxlayer %itype).c_str(), fphxcontainer_x, fphxcontainer_y, fphxcontainer_z);
        G4LogicalVolume *fphxcontainer_volume = new G4LogicalVolume(fphxcontainer_box, G4Material::GetMaterial("G4_AIR"), boost::str(boost::format("fphxcontainer_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(fphxcontainer_volume);
        G4VisAttributes *fphxcontainer_vis = new G4VisAttributes();
        fphxcontainer_vis->SetVisibility(false);
        fphxcontainer_vis->SetForceSolid(false);
        fphxcontainer_volume->SetVisAttributes(fphxcontainer_vis);

        const int ncopy = nstrips_z_sensor;
        const double offsetx = 0.;
        const double offsety = 0.;
        const double offsetz = (ncopy%2==0) ? -2.*strip_z*double(ncopy/2) + strip_z : -2.*strip_z*double(ncopy/2);

        G4VPVParameterisation *fphxparam = new PHG4SiliconTrackerFPHXParameterisation(offsetx, +offsety, offsetz, 2.*strip_z, ncopy);
        new G4PVParameterised(boost::str(boost::format("fphxcontainer_%d_%d") %sphxlayer %itype).c_str(), fphx_volume, fphxcontainer_volume, kZAxis, ncopy, fphxparam, overlapcheck);

        /*
         * PGS
         */
        const double pgs_y = sifull_y + gap_sensor_fphx + 2.*fphx_y; // between FPHX's outer edges
        const double pgs_z = hdi_z;

        G4VSolid *pgs_box = new G4Box(boost::str(boost::format("pgs_box_%d_%d") %sphxlayer %itype).c_str(), pgs_x, pgs_y, pgs_z);
        //G4LogicalVolume *pgs_volume = new G4LogicalVolume(pgs_box, Copper, boost::str(boost::format("pgs_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	G4LogicalVolume *pgs_volume = new G4LogicalVolume(pgs_box,  G4Material::GetMaterial("G4_C"), boost::str(boost::format("pgs_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(pgs_volume);
        G4VSolid *pgs_ext_box = new G4Box(boost::str(boost::format("pgs_ext_box_%d_%s") %sphxlayer %itype).c_str(), pgs_x, pgs_y, hdi_ext_z);
        //G4LogicalVolume *pgs_ext_volume = new G4LogicalVolume(pgs_ext_box, Copper, boost::str(boost::format("pgs_ext_volume_%d_%s") %sphxlayer %itype).c_str(), 0, 0, 0);
	G4LogicalVolume *pgs_ext_volume = new G4LogicalVolume(pgs_ext_box, G4Material::GetMaterial("G4_C"), boost::str(boost::format("pgs_ext_volume_%d_%s") %sphxlayer %itype).c_str(), 0, 0, 0);

        G4VisAttributes *pgs_vis = new G4VisAttributes();
        pgs_vis->SetVisibility(true);
        pgs_vis->SetForceSolid(true);
        pgs_vis->SetColour(G4Colour::Yellow());
        pgs_volume->SetVisAttributes(pgs_vis);

        /*
         * Carbon stave
         */
        const double stave_y = pgs_y;
        const double stave_z = hdi_z;

        G4VSolid *stave_box = new G4Box(boost::str(boost::format("stave_box_%d_%d") %sphxlayer %itype).c_str(), stave_x, stave_y, stave_z);
        //G4LogicalVolume *stave_volume = new G4LogicalVolume(stave_box, Copper, boost::str(boost::format("stave_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	G4LogicalVolume *stave_volume = new G4LogicalVolume(stave_box,  G4Material::GetMaterial("G4_C"), boost::str(boost::format("stave_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(stave_volume);
        G4VSolid *stave_ext_box = new G4Box(boost::str(boost::format("stave_ext_box_%d_%s") %sphxlayer %itype).c_str(), stave_x, stave_y, hdi_ext_z);
        //G4LogicalVolume *stave_ext_volume = new G4LogicalVolume(stave_ext_box,  Copper, boost::str(boost::format("stave_ext_volume_%d_%s") %sphxlayer %itype).c_str(), 0, 0, 0);
        G4LogicalVolume *stave_ext_volume = new G4LogicalVolume(stave_ext_box,  G4Material::GetMaterial("G4_C"), boost::str(boost::format("stave_ext_volume_%d_%s") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(stave_ext_volume);
        G4VisAttributes *stave_vis = new G4VisAttributes();
        stave_vis->SetVisibility(true);
        stave_vis->SetForceSolid(true);
        stave_vis->SetColour(G4Colour::Yellow());
        stave_volume->SetVisAttributes(stave_vis);

        /*
         * Ladder
         */
        const double ladder_x = stave_x + pgs_x + hdi_x + fphx_x;
        const double ladder_y = hdi_y;
        const double ladder_z = hdi_z;
        G4VSolid *ladder_box = new G4Box(boost::str(boost::format("ladder_box_%d_%d") %sphxlayer %itype).c_str(), ladder_x, ladder_y, ladder_z);
        G4LogicalVolume *ladder_volume = new G4LogicalVolume(ladder_box, G4Material::GetMaterial("G4_AIR"), boost::str(boost::format("ladder_volume_%d_%d") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(ladder_volume);
        G4VSolid *ladder_ext_box = new G4Box(boost::str(boost::format("ladder_ext_box_%d_%s") %sphxlayer %itype).c_str(), ladder_x, ladder_y, hdi_ext_z);
        G4LogicalVolume *ladder_ext_volume = new G4LogicalVolume(ladder_ext_box, G4Material::GetMaterial("G4_AIR"), boost::str(boost::format("ladder_ext_volume_%d_%s") %sphxlayer %itype).c_str(), 0, 0, 0);
	absorberlogvols.insert(ladder_ext_volume);
        G4VisAttributes *ladder_vis = new G4VisAttributes();
        ladder_vis->SetVisibility(true);
        ladder_vis->SetForceSolid(true);
        ladder_vis->SetColour(G4Colour::Cyan());
        ladder_volume->SetVisAttributes(ladder_vis);

        /*----- Step 2 -----
         * We place Si-sensor, FPHX, HDI, PGS sheet, and stave at the ladder
         volume.
         -----*/

        /*
         * Carbon stave
         */
        const double TVstave_x = -ladder_x + stave_x;
        new G4PVPlacement(0, G4ThreeVector(TVstave_x, 0.0, 0.0), stave_volume,     boost::str(boost::format("stave_%d_%d")     %sphxlayer %itype).c_str(), ladder_volume,     false, 0, overlapcheck);
        new G4PVPlacement(0, G4ThreeVector(TVstave_x, 0.0, 0.0), stave_ext_volume, boost::str(boost::format("stave_ext_%d_%s") %sphxlayer %itype).c_str(), ladder_ext_volume, false, 0, overlapcheck);

        /*
         * PGS
         */
        const double TVpgs_x = TVstave_x + stave_x + pgs_x;
        new G4PVPlacement(0, G4ThreeVector(TVpgs_x, 0.0, 0.0), pgs_volume,     boost::str(boost::format("pgs_%d_%d")     %sphxlayer %itype).c_str(), ladder_volume,     false, 0, overlapcheck);
        new G4PVPlacement(0, G4ThreeVector(TVpgs_x, 0.0, 0.0), pgs_ext_volume, boost::str(boost::format("pgs_ext_%d_%s") %sphxlayer %itype).c_str(), ladder_ext_volume, false, 0, overlapcheck);

        /*
         * HDI
         */
        const double TVhdi_x = TVpgs_x + pgs_x + hdi_x;
        new G4PVPlacement(0, G4ThreeVector(TVhdi_x, 0.0, 0.0), hdi_volume,     boost::str(boost::format("hdi_%d_%d")     %sphxlayer %itype).c_str(), ladder_volume,     false, 0, overlapcheck);
        new G4PVPlacement(0, G4ThreeVector(TVhdi_x, 0.0, 0.0), hdi_ext_volume, boost::str(boost::format("hdi_ext_%d_%s") %sphxlayer %itype).c_str(), ladder_ext_volume, false, 0, overlapcheck);

        /*
         * Si-sensor
         */
        const double TVSi_x = TVhdi_x + hdi_x + siactive_x;
        new G4PVPlacement(0, G4ThreeVector(TVSi_x, 0.0, 0.0), siinactive_volume, boost::str(boost::format("siinactive_%d_%d") %sphxlayer %itype).c_str(), ladder_volume, false, 0, overlapcheck);
        new G4PVPlacement(0, G4ThreeVector(TVSi_x, 0.0, 0.0),   siactive_volume, boost::str(boost::format("siactive_%d_%d")   %sphxlayer %itype).c_str(), ladder_volume, false, 0, overlapcheck);

        /*
         * FPHX
         */
        const double TVfphx_x = TVhdi_x + hdi_x + fphx_x;
        const double TVfphx_y = sifull_y + gap_sensor_fphx + fphx_y;
        new G4PVPlacement(0, G4ThreeVector(TVfphx_x, +TVfphx_y, 0.0), fphxcontainer_volume, boost::str(boost::format("fphxcontainerp_%d_%d") %sphxlayer %itype).c_str(), ladder_volume, false, 0, overlapcheck);

        new G4PVPlacement(0, G4ThreeVector(TVfphx_x, -TVfphx_y, 0.0), fphxcontainer_volume, boost::str(boost::format("fphxcontainerm_%d_%d") %sphxlayer %itype).c_str(), ladder_volume, false, 0, overlapcheck);

        /*----- Step 3 -----
         * We make cylinder volume in each layer and then install the silicon
         ladders
         -----*/

        /*
         * Ladder
         */
        const double dphi      = CLHEP::twopi/(double)nladders_layer;
        eff_radius[ilayer]     = radius + ladder_x - 2.*(fphx_x-strip_x);
        posz[ilayer][itype]    = (itype == 0) ? hdi_z : 2. * hdi_z_[ilayer][0] + hdi_z;
        strip_x_offset[ilayer] = ladder_x - 2.*(fphx_x-strip_x) - strip_x;

        for (G4int icopy = 0; icopy < nladders_layer; icopy++)
          {
            const double phi  = offsetphi + dphi * (double)icopy;
            const double posx = eff_radius[ilayer] * cos(phi);
            const double posy = eff_radius[ilayer] * sin(phi);
            const double fRotate = phi + offsetrot + CLHEP::pi;

            G4RotationMatrix *ladderrotation = new G4RotationMatrix();
            ladderrotation->rotateZ(-fRotate);

            int sitype[2] = {-1};
            if (itype == 0)
              {
                sitype[0] = 1;
                sitype[1] = 2;
              }
            else
              {
                sitype[0] = 0;
                sitype[1] = 3;
              }

            new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, -posz[ilayer][itype]), ladder_volume, boost::str(boost::format("ladder_%d_%d_%d_%d") %sphxlayer %inttlayer %sitype[0] % icopy).c_str(), trackerenvelope, false, 0, overlapcheck);
            new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, +posz[ilayer][itype]), ladder_volume, boost::str(boost::format("ladder_%d_%d_%d_%d") %sphxlayer %inttlayer %sitype[1] % icopy).c_str(), trackerenvelope, false, 0, overlapcheck);

            if (itype != 0)
              { // HDI tab
                const G4double posz_ext = 2.*(hdi_z_[ilayer][0] + hdi_z) + hdi_ext_z;
                new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, -posz_ext), ladder_ext_volume, boost::str(boost::format("ladder_ext_%d_%d_%d_%d") %sphxlayer %inttlayer %sitype[0] % icopy).c_str(), trackerenvelope, false, 0, overlapcheck);
                new G4PVPlacement(ladderrotation, G4ThreeVector(posx, posy, +posz_ext), ladder_ext_volume, boost::str(boost::format("ladder_ext_%d_%d_%d_%d") %sphxlayer %inttlayer %sitype[1] % icopy).c_str(), trackerenvelope, false, 0, overlapcheck);
              }

          }
      }

  return 0;
}

int PHG4SiliconTrackerDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix *rotm)
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

void PHG4SiliconTrackerDetector::AddGeometryNode()
{
  if (active)
    {
	  std::string geonode = (superdetector != "NONE") ? boost::str(boost::format("CYLINDERGEOM_%s") %superdetector) : boost::str(boost::format("CYLINDERGEOM_%s") %detector_type %nlayer_);

      PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonode.c_str());
      if (!geo)
        {
          geo = new PHG4CylinderGeomContainer();
          PHNodeIterator iter(topNode);
          PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.c_str(), "PHObject");
          runNode->addNode(newNode);
        }

      for (unsigned int ilayer=0; ilayer<nlayer_; ++ilayer)
        {
          const int sphxlayer = layerconfig_[ilayer].first;
          const int inttlayer = layerconfig_[ilayer].second;

          const double strip_z0 = (ilayer==0) ? arr_strip_z[0][0] : arr_strip_z[1][0];
          const double strip_z1 = (ilayer==0) ? arr_strip_z[0][1] : arr_strip_z[1][1];
          const int nstrips_z_sensor0 = (ilayer==0) ? arr_nstrips_z_sensor[0][0] : arr_nstrips_z_sensor[1][0];
          const int nstrips_z_sensor1 = (ilayer==0) ? arr_nstrips_z_sensor[0][1] : arr_nstrips_z_sensor[1][1];

          PHG4CylinderGeom *mygeom = new PHG4CylinderGeom_Siladders(
				       sphxlayer,
                                       arr_strip_x[inttlayer]/cm,
                                       arr_strip_y[inttlayer]/cm,
                                       strip_z0/cm,
                                       strip_z1/cm,
                                       nstrips_z_sensor0,
                                       nstrips_z_sensor1,
                                       arr_nstrips_phi_cell[inttlayer],
                                       arr_nladders_layer[inttlayer],
                                       posz[ilayer][0]/cm,
                                       posz[ilayer][1]/cm,
                                       eff_radius[ilayer]/cm,
                                       strip_x_offset[ilayer]/cm,
                                       arr_offsetphi[inttlayer]/rad,
                                       arr_offsetrot[inttlayer]/rad
                                     );
          geo->AddLayerGeom(sphxlayer, mygeom);
	  if(verbosity>0)
	    geo->identify();
        }
    }
}
