#include "PHG4INTTDetector.h"
#include "PHG4CylinderGeomINTT.h"
#include "PHG4INTTParameterisation.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVParameterised.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VisAttributes.hh>

#include <array>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <cmath>

using namespace std;

PHG4INTTDetector::PHG4INTTDetector(PHCompositeNode *Node, PHParametersContainer *parameters, const std::string &dnam, const pair<vector<pair<int, int>>::const_iterator, vector<pair<int, int>>::const_iterator> &layer_b_e)
  : PHG4Detector(Node, dnam)
  , m_ParamsContainer(parameters)
  , m_IsSupportActive(0)
  , m_LayerBeginEndIteratorPair(layer_b_e)
{
  for (auto layeriter = m_LayerBeginEndIteratorPair.first; layeriter != m_LayerBeginEndIteratorPair.second; ++layeriter)
  {
    int layer = layeriter->second;
    const PHParameters *par = m_ParamsContainer->GetParameters(layer);
    m_IsActiveMap.insert(make_pair(layer, par->get_int_param("active")));
    m_IsAbsorberActiveMap.insert(make_pair(layer, par->get_int_param("absorberactive")));
  }
  const PHParameters *par = m_ParamsContainer->GetParameters(PHG4INTTDefs::SUPPORTPARAMS);
  m_IsSupportActive = par->get_int_param("supportactive");
  fill_n(&m_PosZ[0][0], sizeof(m_PosZ) / sizeof(double), NAN);
  fill_n(m_SensorRadius, sizeof(m_SensorRadius) / sizeof(double), NAN);
  fill_n(m_StripOffsetX, sizeof(m_StripOffsetX) / sizeof(double), NAN);
}

//_______________________________________________________________
int PHG4INTTDetector::IsInINTT(G4VPhysicalVolume *volume) const
{
  // Is this volume one of the sensor strips?
  // just checking if the pointer to the logical volume is in the set
  // of our active/active absorber ones makes sure we are in an active volume
  // name parsing is a bad idea since this is called for all steps
  // and we would have to trust that people give different names
  // to their volumes
  G4LogicalVolume *logvol = volume->GetLogicalVolume();
  if (!m_PassiveVolumeTuple.empty() && m_PassiveVolumeTuple.find(logvol) != m_PassiveVolumeTuple.end())
  {
    return -1;
  }
  if (m_ActiveLogVols.find(logvol) != m_ActiveLogVols.end())
  {
    return 1;
  }

  return 0;
}

void PHG4INTTDetector::Construct(G4LogicalVolume *logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4INTTDetector::Construct called for layers " << endl;
    for (auto layeriter = m_LayerBeginEndIteratorPair.first; layeriter != m_LayerBeginEndIteratorPair.second; ++layeriter)
    {
      cout << "layer " << layeriter->second << endl;
    }
  }
  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  ConstructINTT(logicWorld);

  // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int PHG4INTTDetector::ConstructINTT(G4LogicalVolume *trackerenvelope)
{
  // We have an arbitray number of layers (nlayer_) up to 8
  // We have 2 types of ladders (vertical strips and horizontal strips)
  // We have 2 types of sensors (inner and outer)
  array<array<double, 2>, 8> hdi_z_arr;
  // we loop over layers. All layers have only one laddertype
  for (auto layeriter = m_LayerBeginEndIteratorPair.first; layeriter != m_LayerBeginEndIteratorPair.second; ++layeriter)
  {
    int inttlayer = layeriter->second;
    // get the parameters for this layer
    const PHParameters *params1 = m_ParamsContainer->GetParameters(inttlayer);
    const int laddertype = params1->get_int_param("laddertype");
    const double offsetphi = (params1->get_double_param("offsetphi") * deg) / rad;  // use rad internally
    double offsetrot = (params1->get_double_param("offsetrot") * deg) / rad;        // offsetrot is specified in deg, we convert to rad here
    m_SensorRadius[inttlayer] = params1->get_double_param("sensor_radius") * cm;
    const int nladders_layer = params1->get_int_param("nladder");

    // Look up all remaining parameters by the laddertype for this layer
    const PHParameters *params = m_ParamsContainer->GetParameters(laddertype);
    const double strip_x = params->get_double_param("strip_x") * cm;
    const double strip_y = params->get_double_param("strip_y") * cm;
    const int nstrips_phi_sensor = params->get_int_param("nstrips_phi_sensor");
    const double sensor_offset_y = params->get_double_param("sensor_offset_y") * cm;
    const double hdi_y = params->get_double_param("hdi_y") * cm;
    double hdi_kapton_x = params->get_double_param("hdi_kapton_x") * cm;
    double hdi_copper_x = params->get_double_param("hdi_copper_x") * cm;
    double fphx_x = params->get_double_param("fphx_x") * cm;
    double fphx_y = params->get_double_param("fphx_y") * cm;
    double fphx_z = params->get_double_param("fphx_z") * cm;
    double fphx_offset_z = params->get_double_param("fphx_offset_z") * cm;
    double pgs_x = params->get_double_param("pgs_x") * cm;
    double halfladder_inside_z = params->get_double_param("halfladder_inside_z") * cm;
    double stave_straight_outer_y = params->get_double_param("stave_straight_outer_y") * cm;
    double stave_straight_inner_y = params->get_double_param("stave_straight_inner_y") * cm;
    double stave_straight_cooler_y = params->get_double_param("stave_straight_cooler_y") * cm;

    if (Verbosity() > 0)
    {
      cout << "Constructing INTT layer: " << endl;
      cout << "  layer " << inttlayer << " laddertype " << laddertype << " nladders_layer " << nladders_layer
           << " sensor_radius " << m_SensorRadius[inttlayer] << " offsetphi " << offsetphi << " rad "
           << " offsetphi " << offsetphi * rad / deg << " deg "
           << endl;
    }
    // We loop over inner, then outer, sensors, where  itype specifies the inner or outer sensor
    // The rest of this loop will construct and put in place a section of a ladder corresponding to the Z range of this sensor only
    for (int itype = 0; itype < 2; ++itype)
    {
      if (!(itype >= 0 && itype <= 1))
      {
        assert(!"Error: check ladder type.");
      }
      double strip_z;
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

      // Create Si-sensor active volume
      const double siactive_x = strip_x;
      const double siactive_y = strip_y * nstrips_phi_sensor;
      const double siactive_z = strip_z * nstrips_z_sensor;
      G4VSolid *siactive_box = new G4Box((boost::format("siactive_box_%d_%d") % inttlayer % itype).str(), siactive_x / 2, siactive_y / 2., siactive_z / 2.);
      G4LogicalVolume *siactive_volume = new G4LogicalVolume(siactive_box, G4Material::GetMaterial("G4_Si"),
                                                             boost::str(boost::format("siactive_volume_%d_%d") % inttlayer % itype).c_str(), 0, 0, 0);
      if ((m_IsActiveMap.find(inttlayer))->second > 0)
      {
        m_ActiveLogVols.insert(siactive_volume);
      }
      G4VisAttributes siactive_vis;
      siactive_vis.SetVisibility(true);
      siactive_vis.SetForceSolid(true);
      siactive_vis.SetColour(G4Colour::White());
      siactive_volume->SetVisAttributes(siactive_vis);

      // We do not subdivide the sensor in G4. We will assign hits to strips in the stepping action, using the geometry object

      // Si-sensor full (active+inactive) area
      const double sifull_x = siactive_x;
      const double sifull_y = siactive_y + 2.0 * params->get_double_param("sensor_edge_phi") * cm;
      const double sifull_z = siactive_z + 2.0 * params->get_double_param("sensor_edge_z") * cm;
      G4VSolid *sifull_box = new G4Box((boost::format("sifull_box_%d_%d") % inttlayer % itype).str(), sifull_x / 2., sifull_y / 2.0, sifull_z / 2.0);

      // Si-sensor inactive area
      G4VSolid *siinactive_box = new G4SubtractionSolid((boost::format("siinactive_box_%d_%d") % inttlayer % itype).str(),
                                                        sifull_box, siactive_box, 0, G4ThreeVector(0, 0, 0));
      G4LogicalVolume *siinactive_volume = new G4LogicalVolume(siinactive_box, G4Material::GetMaterial("G4_Si"),
                                                               (boost::format("siinactive_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(siinactive_volume, make_tuple(inttlayer, PHG4INTTDefs::SI_INACTIVE)));
      }
      G4VisAttributes siinactive_vis;
      siinactive_vis.SetVisibility(true);
      siinactive_vis.SetForceSolid(true);
      siinactive_vis.SetColour(G4Colour::Red());
      siinactive_volume->SetVisAttributes(siinactive_vis);

      // Make the HDI Kapton and copper volumes

      // This makes HDI volumes that matche this sensor in Z length
      const double hdi_z = sifull_z + params->get_double_param("hdi_edge_z") * cm;
      hdi_z_arr[inttlayer][itype] = hdi_z;
      G4VSolid *hdi_kapton_box = new G4Box((boost::format("hdi_kapton_box_%d_%d") % inttlayer % itype).str(), hdi_kapton_x / 2., hdi_y / 2., hdi_z / 2.0);
      G4LogicalVolume *hdi_kapton_volume = new G4LogicalVolume(hdi_kapton_box, G4Material::GetMaterial("G4_KAPTON"),
                                                               (boost::format("hdi_kapton_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(hdi_kapton_volume, make_tuple(inttlayer, PHG4INTTDefs::HDI_KAPTON)));
      }
      G4VSolid *hdi_copper_box = new G4Box((boost::format("hdi_copper_box_%d_%d") % inttlayer % itype).str(), hdi_copper_x / 2., hdi_y / 2., hdi_z / 2.0);
      G4LogicalVolume *hdi_copper_volume = new G4LogicalVolume(hdi_copper_box, G4Material::GetMaterial("G4_Cu"),
                                                               (boost::format("hdi_copper_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(hdi_copper_volume, make_tuple(inttlayer, PHG4INTTDefs::HDI_COPPER)));
      }

      // This is the part of the HDI that extends beyond the sensor inside the endcap ring
      const double hdiext_z = (itype == 0) ? 0.000001 : halfladder_inside_z - hdi_z_arr[inttlayer][0] - hdi_z;  // need to assign nonzero value for itype=0
      G4VSolid *hdiext_kapton_box = new G4Box((boost::format("hdiext_kapton_box_%d_%s") % inttlayer % itype).str(),
                                              hdi_kapton_x / 2., hdi_y / 2., hdiext_z / 2.0);
      G4LogicalVolume *hdiext_kapton_volume = new G4LogicalVolume(hdiext_kapton_box, G4Material::GetMaterial("G4_KAPTON"),  // was "FPC"
                                                                  (boost::format("hdiext_kapton_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(hdiext_kapton_volume, make_tuple(inttlayer, PHG4INTTDefs::HDIEXT_KAPTON)));
      }
      G4VSolid *hdiext_copper_box = new G4Box((boost::format("hdiext_copper_box_%d_%s") % inttlayer % itype).str(),
                                              hdi_copper_x / 2., hdi_y / 2., hdiext_z / 2.0);
      G4LogicalVolume *hdiext_copper_volume = new G4LogicalVolume(hdiext_copper_box, G4Material::GetMaterial("G4_Cu"),
                                                                  (boost::format("hdiext_copper_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(hdiext_copper_volume, make_tuple(inttlayer, PHG4INTTDefs::HDIEXT_COPPER)));
      }
      G4VisAttributes hdi_kapton_vis;
      hdi_kapton_vis.SetVisibility(true);
      hdi_kapton_vis.SetForceSolid(true);
      hdi_kapton_vis.SetColour(G4Colour::Yellow());
      hdi_kapton_volume->SetVisAttributes(hdi_kapton_vis);
      hdiext_kapton_volume->SetVisAttributes(hdi_kapton_vis);
      G4VisAttributes hdi_copper_vis;
      hdi_copper_vis.SetVisibility(true);
      hdi_copper_vis.SetForceSolid(true);
      hdi_copper_vis.SetColour(G4Colour::White());
      hdi_copper_volume->SetVisAttributes(hdi_copper_vis);
      hdiext_copper_volume->SetVisAttributes(hdi_copper_vis);

      // FPHX
      G4VSolid *fphx_box = new G4Box((boost::format("fphx_box_%d_%d") % inttlayer % itype).str(), fphx_x / 2., fphx_y / 2., fphx_z / 2.);
      G4LogicalVolume *fphx_volume = new G4LogicalVolume(fphx_box, G4Material::GetMaterial("G4_Si"),
                                                         (boost::format("fphx_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(fphx_volume, make_tuple(inttlayer, PHG4INTTDefs::FPHX)));
      }

      G4VisAttributes fphx_vis;
      fphx_vis.SetVisibility(true);
      fphx_vis.SetForceSolid(true);
      fphx_vis.SetColour(G4Colour::Blue());
      fphx_volume->SetVisAttributes(fphx_vis);

      const double gap_sensor_fphx = params->get_double_param("gap_sensor_fphx") * cm;

      //  FPHX Container
      // make a container for the FPHX chips needed for this sensor, and  then place them in the container
      G4VSolid *fphxcontainer_box = new G4Box((boost::format("fphxcontainer_box_%d_%d") % inttlayer % itype).str(),
                                              fphx_x / 2., fphx_y / 2., hdi_z / 2.);
      G4LogicalVolume *fphxcontainer_volume = new G4LogicalVolume(fphxcontainer_box, G4Material::GetMaterial("G4_AIR"),
                                                                  (boost::format("fphxcontainer_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      // Install multiple FPHX volumes in the FPHX container volume
      // one FPHX chip per cell - each cell is 128 channels
      const double offsetx = 0.;
      const double offsety = 0.;
      int ncopy;
      double offsetz, cell_length_z;

      if (laddertype == PHG4INTTDefs::SEGMENTATION_Z)  // vertical strips
      {
        // For laddertype 0, we have 5 cells per sensor, but the strips are vertical, so we have to treat it specially
        ncopy = nstrips_z_sensor / 128.0;
      }
      else if (laddertype == PHG4INTTDefs::SEGMENTATION_PHI)
      {
        ncopy = nstrips_z_sensor;
      }
      else
      {
        cout << PHWHERE << "invalid laddertype " << laddertype << endl;
        gSystem->Exit(1);
        // this is just to make the optimizer happy which otherwise complains about possibly
        // uninitialized variables. It doesn't know gSystem->Exit(1) quits,
        // this exit here terminates the program for it
        exit(1);
      }
      cell_length_z = strip_z * nstrips_z_sensor / ncopy;
      offsetz = (ncopy % 2 == 0) ? -2. * cell_length_z / 2. * double(ncopy / 2) + cell_length_z / 2. + fphx_offset_z : -2. * cell_length_z / 2. * double(ncopy / 2) + fphx_offset_z;
      G4VPVParameterisation *fphxparam = new PHG4INTTFPHXParameterisation(offsetx, +offsety, offsetz, 2. * cell_length_z / 2., ncopy);
      new G4PVParameterised((boost::format("fphxcontainer_%d_%d") % inttlayer % itype).str(),
                            fphx_volume, fphxcontainer_volume, kZAxis, ncopy, fphxparam, OverlapCheck());

      // PGS   - this is the carbon sheet that the HDI sits on. It forms the wall of the cooling tube that cools the HDI

      const double pgs_y = hdi_y;
      const double pgs_z = hdi_z;
      G4VSolid *pgs_box = new G4Box((boost::format("pgs_box_%d_%d") % inttlayer % itype).str(), pgs_x / 2., pgs_y / 2., pgs_z / 2.);
      G4LogicalVolume *pgs_volume = new G4LogicalVolume(pgs_box, G4Material::GetMaterial("CFRP_INTT"),
                                                        (boost::format("pgs_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(pgs_volume, make_tuple(inttlayer, PHG4INTTDefs::PGS)));
      }
      // The part that extends beyond this sensor, see above for hdiext
      G4VSolid *pgsext_box = new G4Box((boost::format("pgsext_box_%d_%s") % inttlayer % itype).str(), pgs_x / 2., pgs_y / 2., hdiext_z / 2.);
      G4LogicalVolume *pgsext_volume = new G4LogicalVolume(pgsext_box, G4Material::GetMaterial("CFRP_INTT"),
                                                           (boost::format("pgsext_volume_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(pgsext_volume, make_tuple(inttlayer, PHG4INTTDefs::PGSEXT)));
      }
      G4VisAttributes pgs_vis;
      pgs_vis.SetVisibility(true);
      pgs_vis.SetForceSolid(true);
      pgs_vis.SetColour(G4Colour::Red());
      pgs_volume->SetVisAttributes(pgs_vis);
      pgsext_volume->SetVisAttributes(pgs_vis);

      // Carbon stave. This is the formed sheet that sits on the PGS and completes the cooling tube
      // Formed from straight sections and sections of a tube of radius 2.3 mm. All have wall thickness of 0.3 mm.
      // These are different for laddertype PHG4INTTDefs::SEGMENTATION_Z  and
      // PHG4INTTDefs::SEGMENTATION_PHI, but they use some common elements.

      // The curved section is made from a G4Tubs, which is a generalized section of a cylinder
      // Two curved sections combined should move the inner wall to be 2.0 mm away from the PGS, then 2 more sections bring it back
      // For 1 of these tube sections,starting at 90 degrees, take decrease in y=1 mm at avge R.
      // If avge R = 2.15 mm, dtheta = invcos( (R - y)/R ) = invcos(1.15/2.15) = 53.49 deg
      // The extent along the x axis is then R*sin(dtheta) = 1.728 mm, so two sections combined have dx = 3.456 mm length in y
      const double Rcmin = 0.20 * cm;  // 2 mm  inner radius of curved section, same at both ends
      const double Rcmax = 0.23 * cm;  //cm  outer radius of curved section, same at both ends
      double Rcavge = (Rcmax + Rcmin) / 2.0;
      double dphi_c = acos((Rcavge - Rcmin / 2.) / Rcavge);
      const double stave_z = pgs_z;

      // makecurved sections for cooler tube
      const double phic_begin[4] = {M_PI - dphi_c, -dphi_c, 0.0, M_PI};
      const double dphic[4] = {dphi_c, dphi_c, dphi_c, dphi_c};

      G4Tubs *stave_curve_cons[4];
      G4Tubs *stave_curve_ext_cons[4];
      G4LogicalVolume *stave_curve_volume[4];
      G4LogicalVolume *stave_curve_ext_volume[4];

      for (int i = 0; i < 4; i++)
      {
        stave_curve_cons[i] = new G4Tubs((boost::format("stave_curve_cons_%d_%d_%d") % inttlayer % itype % i).str(),
                                         Rcmin, Rcmax, stave_z / 2., phic_begin[i], dphic[i]);
        stave_curve_volume[i] = new G4LogicalVolume(stave_curve_cons[i], G4Material::GetMaterial("CFRP_INTT"),
                                                    (boost::format("stave_curve_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);
        if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
        {
          m_PassiveVolumeTuple.insert(make_pair(stave_curve_volume[i], make_tuple(inttlayer, PHG4INTTDefs::STAVE_CURVE)));
        }
        stave_curve_ext_cons[i] = new G4Tubs((boost::format("stave_curve_ext_cons_%d_%d_%d") % inttlayer % itype % i).str(),
                                             Rcmin, Rcmax, hdiext_z / 2., phic_begin[i], dphic[i]);
        stave_curve_ext_volume[i] = new G4LogicalVolume(stave_curve_ext_cons[i], G4Material::GetMaterial("CFRP_INTT"),
                                                        (boost::format("stave_curve_ext_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);
        if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
        {
          m_PassiveVolumeTuple.insert(make_pair(stave_curve_ext_volume[i], make_tuple(inttlayer, PHG4INTTDefs::STAVEEXT_CURVE)));
        }
        G4VisAttributes stave_curve_vis;
        stave_curve_vis.SetVisibility(true);
        stave_curve_vis.SetForceSolid(true);
        stave_curve_vis.SetColour(G4Colour::White());
        stave_curve_volume[i]->SetVisAttributes(stave_curve_vis);
        stave_curve_ext_volume[i]->SetVisAttributes(stave_curve_vis);
      }

      // we will need the length in y of the curved section as it is installed in the stave box
      double curve_length_y = Rcavge * sin(dphi_c);

      // Make several straight sections for use in making the stave

      // Outer straight sections of stave
      double stave_wall_thickness = 0.03 * cm;
      G4VSolid *stave_straight_outer_box = new G4Box((boost::format("stave_straight_outer_box_%d_%d") % inttlayer % itype).str(),
                                                     stave_wall_thickness / 2., stave_straight_outer_y / 2., stave_z / 2.);
      G4LogicalVolume *stave_straight_outer_volume = new G4LogicalVolume(stave_straight_outer_box, G4Material::GetMaterial("CFRP_INTT"),
                                                                         (boost::format("stave_straight_outer_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(stave_straight_outer_volume, make_tuple(inttlayer, PHG4INTTDefs::STAVE_STRAIGHT_OUTER)));
      }
      G4VSolid *stave_straight_outer_ext_box = new G4Box((boost::format("stave_straight_outer_ext_box_%d_%s") % inttlayer % itype).str(),
                                                         stave_wall_thickness / 2., stave_straight_outer_y / 2., hdiext_z / 2.);
      G4LogicalVolume *stave_straight_outer_ext_volume = new G4LogicalVolume(stave_straight_outer_ext_box, G4Material::GetMaterial("CFRP_INTT"),
                                                                             (boost::format("stave_straight_outer_ext_volume_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(stave_straight_outer_ext_volume, make_tuple(inttlayer, PHG4INTTDefs::STAVEEXT_STRAIGHT_OUTER)));
      }
      // connects cooling tubes together, only needed for laddertype PHG4INTTDefs::SEGMENTATION_PHI, for laddertype PHG4INTTDefs::SEGMENTATION_Z we just make a dummy
      G4VSolid *stave_straight_inner_box = new G4Box((boost::format("stave_straight_inner_box_%d_%d") % inttlayer % itype).str(),
                                                     stave_wall_thickness / 2., stave_straight_inner_y / 2., stave_z / 2.);
      G4LogicalVolume *stave_straight_inner_volume = new G4LogicalVolume(stave_straight_inner_box, G4Material::GetMaterial("CFRP_INTT"),
                                                                         (boost::format("stave_straight_inner_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(stave_straight_inner_volume, make_tuple(inttlayer, PHG4INTTDefs::STAVE_STRAIGHT_INNER)));
      }
      G4VSolid *stave_straight_inner_ext_box = new G4Box((boost::format("stave_straight_inner_ext_box_%d_%d") % inttlayer % itype).str(),
                                                         stave_wall_thickness / 2., stave_straight_inner_y / 2., hdiext_z / 2.);
      G4LogicalVolume *stave_straight_inner_ext_volume = new G4LogicalVolume(stave_straight_inner_ext_box, G4Material::GetMaterial("CFRP_INTT"),
                                                                             (boost::format("stave_straight_inner_ext_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(stave_straight_inner_ext_volume, make_tuple(inttlayer, PHG4INTTDefs::STAVEEXT_STRAIGHT_INNER)));
      }
      //Top surface of cooler tube
      G4VSolid *stave_straight_cooler_box = new G4Box((boost::format("stave_straight_cooler_box_%d_%d") % inttlayer % itype).str(),
                                                      stave_wall_thickness / 2., stave_straight_cooler_y / 2., stave_z / 2.);
      G4LogicalVolume *stave_straight_cooler_volume = new G4LogicalVolume(stave_straight_cooler_box, G4Material::GetMaterial("CFRP_INTT"),
                                                                          (boost::format("stave_straight_cooler_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(stave_straight_cooler_volume, make_tuple(inttlayer, PHG4INTTDefs::STAVE_STRAIGHT_COOLER)));
      }
      G4VSolid *stave_straight_cooler_ext_box = new G4Box((boost::format("stave_straight_cooler_ext_box_%d_%d") % inttlayer % itype).str(),
                                                          stave_wall_thickness / 2., stave_straight_cooler_y / 2., hdiext_z / 2.);
      G4LogicalVolume *stave_straight_cooler_ext_volume = new G4LogicalVolume(stave_straight_cooler_ext_box, G4Material::GetMaterial("CFRP_INTT"),
                                                                              (boost::format("stave_straight_cooler_ext_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(stave_straight_cooler_ext_volume, make_tuple(inttlayer, PHG4INTTDefs::STAVEEXT_STRAIGHT_COOLER)));
      }
      G4VisAttributes stave_vis;
      stave_vis.SetVisibility(true);
      stave_vis.SetForceSolid(true);
      stave_vis.SetColour(G4Colour::White());
      stave_straight_cooler_volume->SetVisAttributes(stave_vis);
      stave_straight_cooler_ext_volume->SetVisAttributes(stave_vis);
      if (laddertype == PHG4INTTDefs::SEGMENTATION_PHI)
      {
        stave_straight_inner_volume->SetVisAttributes(stave_vis);
      }
      if (laddertype == PHG4INTTDefs::SEGMENTATION_PHI)
      {
        stave_straight_inner_ext_volume->SetVisAttributes(stave_vis);
      }
      stave_straight_outer_volume->SetVisAttributes(stave_vis);
      stave_straight_outer_ext_volume->SetVisAttributes(stave_vis);

      // Now we combine the elements of a stave defined above into a stave
      // Create a stave volume to install the stave sections into. The volume has to be big enouigh to contain the cooling tube
      double cooler_gap_x = 0.2 * cm;  // id of cooling tube in cm
      double cooler_wall = 0.03 * cm;  // outer wall thickness of cooling tube
      double cooler_x = cooler_gap_x + cooler_wall;
      const double stave_x = cooler_x;  // we do not include the PGS in the stave volume
      const double stave_y = hdi_y;
      G4VSolid *stave_box = new G4Box((boost::format("stave_box_%d_%d") % inttlayer % itype).str(), stave_x / 2., stave_y / 2., stave_z / 2.);
      G4LogicalVolume *stave_volume = new G4LogicalVolume(stave_box, G4Material::GetMaterial("G4_AIR"),
                                                          (boost::format("stave_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      G4VSolid *staveext_box = new G4Box((boost::format("staveext_box_%d_%d") % inttlayer % itype).str(), stave_x / 2., stave_y / 2., hdiext_z / 2.);
      G4LogicalVolume *staveext_volume = new G4LogicalVolume(staveext_box, G4Material::GetMaterial("G4_AIR"),
                                                             (boost::format("staveext_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      G4VisAttributes stave_box_vis;
      stave_box_vis.SetVisibility(false);
      stave_box_vis.SetForceSolid(false);
      stave_volume->SetVisAttributes(stave_box_vis);
      staveext_volume->SetVisAttributes(stave_box_vis);

      // Assemble the elements into the stave volume and the stave extension volume
      // They are place relative to the center of the stave box. Thus the offset of the center of the segment is relative to the center of the satev box.
      // But we want the segment to be located relative to the lowest x limit of the stave box.
      if (laddertype == PHG4INTTDefs::SEGMENTATION_Z)
      {
        // only one cooling tube in laddertype 0
        // Place the straight sections. We add the middle, then above x axis, then below x axis
        double x_off_str[3] =
            {
                Rcavge - stave_x / 2.,
                (Rcmax - Rcmin) / 2. - stave_x / 2.,
                (Rcmax - Rcmin) / 2. - stave_x / 2.};
        double y_off_str[3] =
            {
                0.0,
                stave_straight_cooler_y / 2. + 2.0 * curve_length_y + stave_straight_outer_y / 2.0,
                -stave_straight_cooler_y / 2. - 2.0 * curve_length_y - stave_straight_outer_y / 2.0};

        for (int i = 0; i < 3; i++)
        {
          if (i == 0)
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_volume,
                              (boost::format("stave_straight_cooler_%d_%d_%d") % i % inttlayer % itype).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_ext_volume,
                              (boost::format("stave_straight_cooler_ext_%d_%d_%d") % i % inttlayer % itype).str(), staveext_volume, false, 0, OverlapCheck());
          }
          else
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_volume,
                              (boost::format("stave_straight_outer_%d_%d_%d") % i % inttlayer % itype).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_ext_volume,
                              (boost::format("stave_straight_outer_ext_%d_%d_%d") % i % inttlayer % itype).str(), staveext_volume, false, 0, OverlapCheck());
          }
        }
        // The cooler curved sections are made using 2 curved sections in a recurve on each side of the cooler straight section
        // The tube sections used here have the origin of their volume at their center of rotation. Rcavge
        //      Each curve section is moved to the center of the stave volume by a translation of +/- Rcavge
        //      Then it is moved to the outside or the inside of the stave volume by a translation of +/-  cooler_gap_x / 2.
        // we start at lowest y and work up in y

        double x_off_cooler[4] =
            {
                Rcavge - cooler_gap_x / 2.,
                -Rcavge + cooler_gap_x / 2.,
                -Rcavge + cooler_gap_x / 2.,
                Rcavge - cooler_gap_x / 2.};
        double y_off_cooler[4] =
            {
                -stave_straight_cooler_y / 2. - 2. * curve_length_y,
                -stave_straight_cooler_y / 2.,
                stave_straight_cooler_y / 2.,
                stave_straight_cooler_y / 2. + 2. * curve_length_y};

        for (int i = 0; i < 4; i++)
        {
          new G4PVPlacement(0, G4ThreeVector(x_off_cooler[i], y_off_cooler[i], 0.0), stave_curve_volume[i],
                            (boost::format("stave_curve_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector(x_off_cooler[i], y_off_cooler[i], 0.0), stave_curve_ext_volume[i],
                            (boost::format("stave_curve_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
        }
      }
      else if (laddertype == PHG4INTTDefs::SEGMENTATION_PHI)  // The type PHG4INTTDefs::SEGMENTATION_PHI ladder has two cooling tubes
      {
        // First place the straight sections, do the extension at the same time
        // we alternate  positive and negative y values here
        double x_off_str[5] =
            {
                (Rcmax - Rcmin) / 2. - stave_x / 2.,  // against the PGS
                (Rcmax + Rcmin) / 2. - stave_x / 2.,  // top of cooler
                (Rcmax + Rcmin) / 2. - stave_x / 2.,  // top of cooler
                (Rcmax - Rcmin) / 2. - stave_x / 2.,
                (Rcmax - Rcmin) / 2 - stave_x / 2.};
        double y_off_str[5] =
            {
                0.0,                                                                                // center section against PGS
                stave_straight_inner_y / 2. + 2. * curve_length_y + stave_straight_cooler_y / 2.,   // top of cooler
                -stave_straight_inner_y / 2. - 2. * curve_length_y - stave_straight_cooler_y / 2.,  // top of cooler
                stave_straight_inner_y / 2. + 2. * curve_length_y + stave_straight_cooler_y + 2. * curve_length_y + stave_straight_outer_y / 2.,
                -stave_straight_inner_y / 2. - 2. * curve_length_y - stave_straight_cooler_y - 2. * curve_length_y - stave_straight_outer_y / 2.,
            };

        for (int i = 0; i < 5; i++)
        {
          if (i == 0)
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_inner_volume,
                              (boost::format("stave_straight_inner_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_inner_ext_volume,
                              (boost::format("stave_straight_inner_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
          }
          else if (i == 1 || i == 2)
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_volume,
                              (boost::format("stave_straight_cooler_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_ext_volume,
                              (boost::format("stave_straight_cooler_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
          }
          else
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_volume,
                              (boost::format("stave_straight_outer_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_ext_volume,
                              (boost::format("stave_straight_outer_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
          }
        }

        // Place the curved sections
        // here we do all above the x axis, then all below the x axis, each in order of increasing y

        double x_off_curve[8] =
            {
                // below x axis, increasing in y
                Rcavge - cooler_gap_x / 2.,
                -Rcavge + cooler_gap_x / 2.,
                -Rcavge + cooler_gap_x / 2.,
                Rcavge - cooler_gap_x / 2.,
                // above x axis, increasing in y
                Rcavge - cooler_gap_x / 2.,
                -Rcavge + cooler_gap_x / 2.,
                -Rcavge + cooler_gap_x / 2.,
                Rcavge - cooler_gap_x / 2.};
        double y_off_curve[8] =
            {
                // below the x axis, increasing in y
                -stave_straight_inner_y / 2. - 2. * curve_length_y - stave_straight_cooler_y - 2. * curve_length_y,
                -stave_straight_inner_y / 2. - 2. * curve_length_y - stave_straight_cooler_y,
                -stave_straight_inner_y / 2. - 2. * curve_length_y,
                -stave_straight_inner_y / 2.,
                // above the x axis, increasing in y
                stave_straight_inner_y / 2.,
                stave_straight_inner_y / 2. + 2. * curve_length_y,
                stave_straight_inner_y / 2. + 2. * curve_length_y + stave_straight_cooler_y,
                stave_straight_inner_y / 2. + 2. * curve_length_y + stave_straight_cooler_y + 2. * curve_length_y};

        for (int i = 0; i < 8; i++)
        {
          // initially this was ivol = i; if (ivol > 3) ivol = i -4;
          // but that triggered a false positive in coverity
          // this should make it happy and reduce the noise in the coverity reports
          int ivol;
          if (i > 3)
          {
            ivol = i - 4;
          }
          else
          {
            ivol = i;
          }
          new G4PVPlacement(0, G4ThreeVector(x_off_curve[i], y_off_curve[i], 0.0), stave_curve_volume[ivol], (boost::format("stave_curve_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector(x_off_curve[i], y_off_curve[i], 0.0), stave_curve_ext_volume[ivol], (boost::format("stave_curve_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
        }
      }
      else
      {
        cout << PHWHERE << "invalid laddertype " << laddertype << endl;
        gSystem->Exit(1);
      }

      // ----- Step 2 ------------------------------------------------------------------------------------
      // We place Si-sensor, FPHX, HDI, PGS sheet, and stave in the ladder  volume.
      // ======================================================

      // Make the ladder volume first
      // We are still in the loop over inner or outer sensors. This is the ladder volume corresponding to this sensor. The FPHX is taller than the sensor in x.
      const double ladder_x = stave_x + pgs_x + hdi_kapton_x + hdi_copper_x + fphx_x;
      double ladder_y = hdi_y;

      G4VSolid *ladder_box = new G4Box((boost::format("ladder_box_%d_%d") % inttlayer % itype).str(), ladder_x / 2., ladder_y / 2., hdi_z / 2.);
      G4LogicalVolume *ladder_volume = new G4LogicalVolume(ladder_box, G4Material::GetMaterial("G4_AIR"), (boost::format("ladder_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      G4VSolid *ladderext_box = new G4Box((boost::format("ladderext_box_%d_%s") % inttlayer % itype).str(), ladder_x / 2., ladder_y / 2., hdiext_z / 2.);
      G4LogicalVolume *ladderext_volume = new G4LogicalVolume(ladderext_box, G4Material::GetMaterial("G4_AIR"), (boost::format("ladderext_%d_%d_%d") % inttlayer % inttlayer % itype).str(), 0, 0, 0);

      G4VisAttributes ladder_vis;
      ladder_vis.SetVisibility(false);
      ladder_vis.SetForceSolid(false);
      ladder_vis.SetColour(G4Colour::Cyan());
      ladder_volume->SetVisAttributes(ladder_vis);
      ladderext_volume->SetVisAttributes(ladder_vis);

      // Now add the components of the ladder to the ladder volume
      // The sensor is closest to the beam pipe, the stave cooler is furthest away
      // Note that the cooler has been assembled in the stave volume with the top at larger x, so the sensor will be at smaller x
      // That will be the configuration when the ladder is at phi = 0 degrees, the positive x direction

      // We start at the most positive x value and add the stave first

      // Carbon stave
      double TVstave_y = 0.0;
      const double TVstave_x = ladder_x / 2. - stave_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVstave_x, TVstave_y, 0.0), stave_volume, (boost::format("stave_%d_%d") % inttlayer % itype).str(),
                        ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVstave_x, TVstave_y, 0.0), staveext_volume, (boost::format("staveext_%d_%s") % inttlayer % itype).str(),
                        ladderext_volume, false, 0, OverlapCheck());

      // PGS
      const double TVpgs_x = TVstave_x - stave_x / 2. - pgs_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVpgs_x, TVstave_y, 0.0), pgs_volume, (boost::format("pgs_%d_%d") % inttlayer % itype).str(),
                        ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVpgs_x, TVstave_y, 0.0), pgsext_volume, (boost::format("pgsext_%d_%s") % inttlayer % itype).str(),
                        ladderext_volume, false, 0, OverlapCheck());

      // HDI Kapton
      const double TVhdi_kapton_x = TVpgs_x - pgs_x / 2. - hdi_kapton_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVhdi_kapton_x, TVstave_y, 0.0), hdi_kapton_volume, (boost::format("hdikapton_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVhdi_kapton_x, TVstave_y, 0.0), hdiext_kapton_volume, (boost::format("hdiextkapton_%d_%s") % inttlayer % itype).str(), ladderext_volume, false, 0, OverlapCheck());

      // HDI copper
      const double TVhdi_copper_x = TVhdi_kapton_x - hdi_kapton_x / 2. - hdi_copper_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVhdi_copper_x, TVstave_y, 0.0), hdi_copper_volume, (boost::format("hdicopper_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVhdi_copper_x, TVstave_y, 0.0), hdiext_copper_volume, (boost::format("hdiextcopper_%d_%s") % inttlayer % itype).str(), ladderext_volume, false, 0, OverlapCheck());

      // Si-sensor
      double TVSi_y = 0.0;
      // sensor is not centered in y in the ladder volume for the Z sensitive ladders
      if (laddertype == PHG4INTTDefs::SEGMENTATION_Z)
        TVSi_y = +sensor_offset_y;
      const double TVSi_x = TVhdi_copper_x - hdi_copper_x / 2. - siactive_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVSi_x, TVSi_y, 0.0), siinactive_volume,
                        (boost::format("siinactive_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVSi_x, TVSi_y, 0.0), siactive_volume,
                        (boost::format("siactive_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());

      // FPHX container
      const double TVfphx_x = TVhdi_copper_x - hdi_copper_x / 2. - fphx_x / 2.;
      double TVfphx_y = sifull_y / 2. + gap_sensor_fphx + fphx_y / 2.;
      if (laddertype == PHG4INTTDefs::SEGMENTATION_Z)
        TVfphx_y -= sensor_offset_y;
      // laddertype PHG4INTTDefs::SEGMENTATION_Z has only one FPHX, laddertype PHG4INTTDefs::SEGMENTATION_PHI has two
      if (laddertype == PHG4INTTDefs::SEGMENTATION_PHI)
      {
        new G4PVPlacement(0, G4ThreeVector(TVfphx_x, +TVfphx_y, 0.0), fphxcontainer_volume, (boost::format("fphxcontainerp_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      }
      new G4PVPlacement(0, G4ThreeVector(TVfphx_x, -TVfphx_y, 0.0), fphxcontainer_volume, (boost::format("fphxcontainerm_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());

      // ----- Step 3 --------------------------------------------------------------------------------------------------------------------
      // We install the section of ladder for this sensor at all requested phi values and at positive and negative Z
      //========================================================================

      // Distribute Ladders in phi
      // We are still in the loops over layer and sensor type, we will place copies of the ladder section for this sensor
      //  at all ladder phi values, and at both positive and negative Z values.

      // given radius values are for the center of the sensor, we need the x offset from center of ladder to center of sensor so we can place the ladder
      double sensor_offset_x_ladder = 0.0 - TVSi_x;  // ladder center is at x = 0.0 by construction. Sensor is at lower x, so TVSi_x is negative

      const double dphi = 2 * M_PI / nladders_layer;

      m_PosZ[inttlayer][itype] = (itype == 0) ? hdi_z / 2. : hdi_z_arr[inttlayer][0] + hdi_z / 2.;  // location of center of ladder in Z
      m_StripOffsetX[inttlayer] = sensor_offset_x_ladder;

      // The sensors have no tilt in the new design
      //    The type 1 ladders have the sensor at the center of the ladder in phi, so that is easy
      //    The type 0 ladders are more complicated because the sensor center is perpendicular to the radial vector and the sensor is not at the ladder center
      //         We made the stave box symmetric in y around the sensor center to simplify things

      for (int icopy = 0; icopy < nladders_layer; icopy++)
      {
        // sensor center
        const double phi = offsetphi + dphi * icopy;  // if offsetphi is zero we start at zero

        double radius;
        // Make each layer at a single radius - i.e. what was formerly a sub-layer is now considered a layer
        radius = m_SensorRadius[inttlayer];
        radius += sensor_offset_x_ladder;

        double p = 0.0;
        if (laddertype == PHG4INTTDefs::SEGMENTATION_Z)
        {
          // The Z sensitive ladders have the sensors offset in y relative to the ladder center
          // We have to slightly rotate the ladder in its own frame to make the radial vector to the sensor center normal to the sensor face
          p = atan(sensor_offset_y / radius);
          // then we adjust the distance to the center of the ladder to put the sensor at the requested distance from the center of the barrel
          radius /= cos(p);
        }

        // these describe the center of the ladder volume, placing it so that the center of the sensor is at phi = dphi * icopy, and at the correct radius
        const double posx = radius * cos(phi - p);
        const double posy = radius * sin(phi - p);
        const double fRotate = p + (phi - p) + offsetrot;  // rotate in its own frame to make sensor perp to radial vector (p), then additionally rotate to account for ladder phi
        G4RotationMatrix ladderrotation;
        ladderrotation.rotateZ(fRotate);

        // this placement version rotates the ladder in its own frame by fRotate, then translates the center to posx, posy, +/- m_PosZ
        auto pointer_negz = new G4PVPlacement(G4Transform3D(ladderrotation, G4ThreeVector(posx, posy, -m_PosZ[inttlayer][itype])), ladder_volume,
                                              (boost::format("ladder_%d_%d_%d_negz") % inttlayer % itype % icopy).str(), trackerenvelope, false, 0, OverlapCheck());
        auto pointer_posz = new G4PVPlacement(G4Transform3D(ladderrotation, G4ThreeVector(posx, posy, +m_PosZ[inttlayer][itype])), ladder_volume,
                                              (boost::format("ladder_%d_%d_%d_posz") % inttlayer % itype % icopy).str(), trackerenvelope, false, 0, OverlapCheck());
        if (m_IsActiveMap.find(inttlayer) != m_IsActiveMap.end())
        {
          m_ActiveVolumeTuple.insert(make_pair(pointer_negz, make_tuple(inttlayer, itype, icopy, -1)));
          m_ActiveVolumeTuple.insert(make_pair(pointer_posz, make_tuple(inttlayer, itype, icopy, 1)));
        }

        // The net effect of the above manipulations for the Z sensitive ladders is that the center of the sensor is at dphi * icopy and at the requested radius
        // That us all that the geometry object needs to know, so no changes to that are necessary

        if (itype != 0)
        {
          // We have added the outer sensor above, now we add the HDI extension tab to the end of the outer sensor HDI
          //const double posz_ext = (hdi_z_arr[inttlayer][0] + hdi_z) + hdiext_z / 2.;
          const double posz_ext = (hdi_z_arr[inttlayer][0] + hdi_z) + hdiext_z / 2.;

          new G4PVPlacement(G4Transform3D(ladderrotation, G4ThreeVector(posx, posy, -posz_ext)), ladderext_volume,
                            (boost::format("ladderext_%d_%d_%d_negz") % inttlayer % itype % icopy).str(), trackerenvelope, false, 0, OverlapCheck());
          new G4PVPlacement(G4Transform3D(ladderrotation, G4ThreeVector(posx, posy, +posz_ext)), ladderext_volume,
                            (boost::format("ladderext_%d_%d_%d_posz") % inttlayer % itype % icopy).str(), trackerenvelope, false, 0, OverlapCheck());
        }

        if (Verbosity() > 100)
          cout << "   Ladder copy " << icopy << " radius " << radius << " phi " << phi << " itype " << itype << " posz " << m_PosZ[inttlayer][itype]
               << " fRotate " << fRotate << " posx " << posx << " posy " << posy
               << endl;

      }  // end loop over ladder copy placement in phi and positive and negative Z
    }    // end loop over inner or outer sensor
  }      // end loop over layers

  // Finally, we add some support material for the silicon detectors

  //
  /*
    4 rails, which are 12mm OD and 9mm ID tubes at a radius of 168.5 mm.  They are spaced equidistantly in phi.
    The rails run along the entire length of the TPC and even stick out of the TPC, but I think for the moment you don't have to put the parts that stick out in the simulation.
    An inner skin with a ID at 62.416 mm and a thickness of 0.250 mm.
    An outer skin with a ID at 120.444 mm and a sandwich of 0.25 mm cfc, 1.5 mm foam and 0.25 mm cfc.
    
    All of the above are carbon fiber.
  */
  const PHParameters *supportparams = m_ParamsContainer->GetParameters(PHG4INTTDefs::SUPPORTPARAMS);

  // rails
  G4Tubs *rail_tube = new G4Tubs((boost::format("si_support_rail")).str(),
                                 supportparams->get_double_param("rail_inner_radius") * cm,
                                 supportparams->get_double_param("rail_outer_radius") * cm,
                                 supportparams->get_double_param("rail_length") * cm / 2.0,
                                 -M_PI, 2.0 * M_PI);
  G4LogicalVolume *rail_volume = new G4LogicalVolume(rail_tube, G4Material::GetMaterial("CFRP_INTT"),
                                                     "rail_volume", 0, 0, 0);
  if (m_IsSupportActive > 0)
  {
    m_PassiveVolumeTuple.insert(make_pair(rail_volume, make_tuple(PHG4INTTDefs::SUPPORT_DETID, PHG4INTTDefs::SUPPORT_RAIL)));
  }
  G4VisAttributes rail_vis;
  rail_vis.SetVisibility(true);
  rail_vis.SetForceSolid(true);
  rail_vis.SetColour(G4Colour::Cyan());
  rail_volume->SetVisAttributes(rail_vis);

  double rail_dphi = supportparams->get_double_param("rail_dphi") * deg / rad;
  double rail_phi_start = supportparams->get_double_param("rail_phi_start") * deg / rad;
  double rail_radius = supportparams->get_double_param("rail_radius") * cm;
  for (int i = 0; i < 4; i++)
  {
    double phi = rail_phi_start + i * rail_dphi;

    // place a copy at each rail phi value
    const double posx = rail_radius * cos(phi);
    const double posy = rail_radius * sin(phi);

    new G4PVPlacement(0, G4ThreeVector(posx, posy, 0.0), rail_volume,
                      (boost::format("si_support_rail_%d") % i).str(), trackerenvelope, false, 0, OverlapCheck());
  }

  // Outer skin
  G4Tubs *outer_skin_cfcin_tube = new G4Tubs("si_outer_skin_cfcin",
                                       supportparams->get_double_param("outer_skin_cfcin_inner_radius") * cm,
                                       supportparams->get_double_param("outer_skin_cfcin_outer_radius") * cm,
                                       supportparams->get_double_param("outer_skin_cfcin_length") * cm / 2.,
                                       -M_PI, 2.0 * M_PI);
  G4LogicalVolume *outer_skin_cfcin_volume = new G4LogicalVolume(outer_skin_cfcin_tube, G4Material::GetMaterial("CFRP_INTT"),
                                                           "outer_skin_cfcin_volume", 0, 0, 0);

  G4Tubs *outer_skin_foam_tube = new G4Tubs("si_outer_skin_foam",
                                       supportparams->get_double_param("outer_skin_foam_inner_radius") * cm,
                                       supportparams->get_double_param("outer_skin_foam_outer_radius") * cm,
                                       supportparams->get_double_param("outer_skin_foam_length") * cm / 2.,
                                       -M_PI, 2.0 * M_PI);
  G4LogicalVolume *outer_skin_foam_volume = new G4LogicalVolume(outer_skin_foam_tube, G4Material::GetMaterial("ROHACELL_FOAM_110"),
                                                           "outer_skin_foam_volume", 0, 0, 0);

  G4Tubs *outer_skin_cfcout_tube = new G4Tubs("si_outer_skin_cfcout",
                                       supportparams->get_double_param("outer_skin_cfcout_inner_radius") * cm,
                                       supportparams->get_double_param("outer_skin_cfcout_outer_radius") * cm,
                                       supportparams->get_double_param("outer_skin_cfcout_length") * cm / 2.,
                                       -M_PI, 2.0 * M_PI);
  G4LogicalVolume *outer_skin_cfcout_volume = new G4LogicalVolume(outer_skin_cfcout_tube, G4Material::GetMaterial("CFRP_INTT"),
                                                           "outer_skin_cfcout_volume", 0, 0, 0);
  if (m_IsSupportActive > 0)
  {
    m_PassiveVolumeTuple.insert(make_pair(outer_skin_cfcin_volume, make_tuple(PHG4INTTDefs::SUPPORT_DETID, PHG4INTTDefs::INTT_OUTER_SKIN)));
    m_PassiveVolumeTuple.insert(make_pair(outer_skin_foam_volume, make_tuple(PHG4INTTDefs::SUPPORT_DETID, PHG4INTTDefs::INTT_OUTER_SKIN)));
    m_PassiveVolumeTuple.insert(make_pair(outer_skin_cfcout_volume, make_tuple(PHG4INTTDefs::SUPPORT_DETID, PHG4INTTDefs::INTT_OUTER_SKIN)));
  }
  outer_skin_cfcin_volume->SetVisAttributes(rail_vis);
  outer_skin_foam_volume->SetVisAttributes(rail_vis);
  outer_skin_cfcout_volume->SetVisAttributes(rail_vis);
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), outer_skin_cfcin_volume,
                    "si_support_outer_skin_cfcin", trackerenvelope, false, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), outer_skin_foam_volume,
                    "si_support_outer_skin_foam", trackerenvelope, false, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), outer_skin_cfcout_volume,
                    "si_support_outer_skin_cfcout", trackerenvelope, false, 0, OverlapCheck());

  // Inner skin

  G4Tubs *inner_skin_tube = new G4Tubs("si_inner_skin",
                                       supportparams->get_double_param("inner_skin_inner_radius") * cm,
                                       supportparams->get_double_param("inner_skin_outer_radius") * cm,
                                       supportparams->get_double_param("inner_skin_length") * cm / 2.,
                                       -M_PI, 2.0 * M_PI);
  G4LogicalVolume *inner_skin_volume = new G4LogicalVolume(inner_skin_tube, G4Material::GetMaterial("CFRP_INTT"),
                                                           "inner_skin_volume", 0, 0, 0);
  if (m_IsSupportActive > 0)
  {
    m_PassiveVolumeTuple.insert(make_pair(inner_skin_volume, make_tuple(PHG4INTTDefs::SUPPORT_DETID, PHG4INTTDefs::INTT_INNER_SKIN)));
  }
  inner_skin_volume->SetVisAttributes(rail_vis);
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), inner_skin_volume,
                    "si_support_inner_skin", trackerenvelope, false, 0, OverlapCheck());

  // Endcap ring in simulations = Endcap rings + endcap staves 
  
  // Aluminum ring
  G4Tubs *endcap_Al_ring = new G4Tubs("endcap_Al_ring",
                                       supportparams->get_double_param("endcap_Alring_inner_radius") * cm,
                                       supportparams->get_double_param("endcap_Alring_outer_radius") * cm,
                                       supportparams->get_double_param("endcap_Alring_length") * cm / 2.,
                                       -M_PI, 2.0 * M_PI);

  G4LogicalVolume *endcap_Al_ring_volume = new G4LogicalVolume(endcap_Al_ring, G4Material::GetMaterial("Al6061T6"),
                                                               "endcap_Al_ring_volume", 0, 0, 0);
  
  // Stainlees steal ring
  G4Tubs *endcap_SS_ring = new G4Tubs("endcap_SS_ring",
                                      supportparams->get_double_param("endcap_SSring_inner_radius") * cm,
                                      supportparams->get_double_param("endcap_SSring_outer_radius") * cm,
                                      supportparams->get_double_param("endcap_SSring_length") * cm / 2.,
                                      -M_PI, 2.0 * M_PI);

  G4LogicalVolume *endcap_SS_ring_volume = new G4LogicalVolume(endcap_SS_ring, G4Material::GetMaterial("SS316"),
                                                               "endcap_SS_ring_volume", 0, 0, 0);
  
  // Water Glycol ring
  G4Tubs *endcap_WG_ring = new G4Tubs("endcap_WG_ring",
                                       supportparams->get_double_param("endcap_WGring_inner_radius") * cm,
                                       supportparams->get_double_param("endcap_WGring_outer_radius") * cm,
                                       supportparams->get_double_param("endcap_WGring_length") * cm / 2.,
                                       -M_PI, 2.0 * M_PI);

  G4LogicalVolume *endcap_WG_ring_volume = new G4LogicalVolume(endcap_WG_ring, G4Material::GetMaterial("WaterGlycol_INTT"),
                                                               "endcap_WG_ring_volume", 0, 0, 0);

  // Place endcap rings
  G4RotationMatrix ringrotation;
  ringrotation.rotateY( 0 * deg / rad );

  double endcap_ring_z = supportparams->get_double_param("endcap_ring_z") * cm; 
  for(int i = 0; i < 2; i++) // i=0 : positive z, i=1 negative z
  {

    endcap_ring_z = (i==0) ?  endcap_ring_z : -1.0 * endcap_ring_z;

    double width_WGring_z = supportparams->get_double_param("endcap_WGring_length") * cm;
    double width_SSring_z = supportparams->get_double_param("endcap_SSring_length") * cm;
    double width_Alring_z = supportparams->get_double_param("endcap_Alring_length") * cm;

    for(int j = 0; j < 2; j++) // j=0 : positive side z, j=1 negative side z
    {
      width_WGring_z = (j==0) ? width_WGring_z : -1.0 * width_WGring_z; 
      width_SSring_z = (j==0) ? width_SSring_z : -1.0 * width_SSring_z;
      width_Alring_z = (j==0) ? width_Alring_z : -1.0 * width_Alring_z;
    
      double cent_WGring_z = endcap_ring_z + width_WGring_z / 2.;
      double cent_SSring_z = endcap_ring_z + width_WGring_z + width_SSring_z / 2.;
      double cent_Alring_z = endcap_ring_z + width_WGring_z + width_SSring_z + width_Alring_z / 2.;


      new G4PVPlacement(G4Transform3D(ringrotation, G4ThreeVector(0, 0, cent_WGring_z)),
          endcap_WG_ring_volume,
          (boost::format("endcap_WG_ring_pv_%d_%d") %i %j).str(),
          trackerenvelope, false, 0, OverlapCheck());

      new G4PVPlacement(G4Transform3D(ringrotation, G4ThreeVector(0, 0, cent_SSring_z)),
          endcap_SS_ring_volume,
          (boost::format("endcap_SS_ring_pv_%d_%d") %i %j).str(),
          trackerenvelope, false, 0, OverlapCheck());

      new G4PVPlacement(G4Transform3D(ringrotation, G4ThreeVector(0, 0, cent_Alring_z)),
          endcap_Al_ring_volume,
          (boost::format("endcap_Al_ring_pv_%d_%d") %i %j).str(),
          trackerenvelope, false, 0, OverlapCheck());
    }
  }



  //=======================================================
  // Add an outer shell for the MVTX - move this to the MVTX detector module
  //=======================================================
  // A Rohacell foam sandwich made of 0.1 mm thick CFRP skin and 1.8 mm Rohacell 110 foam core, it has a density of 110 kg/m**3.
  double skin_thickness = supportparams->get_double_param("mvtx_shell_skin_thickness") * cm;
  double foam_core_thickness = supportparams->get_double_param("mvtx_shell_foam_core_thickness") * cm;
  double mvtx_shell_length = supportparams->get_double_param("mvtx_shell_length") * cm;

  double mvtx_shell_inner_skin_inner_radius = supportparams->get_double_param("mvtx_shell_inner_skin_inner_radius") * cm;

  double mvtx_shell_foam_core_inner_radius = mvtx_shell_inner_skin_inner_radius + skin_thickness;
  double mvtx_shell_outer_skin_inner_radius = mvtx_shell_foam_core_inner_radius + foam_core_thickness;

  G4Tubs *mvtx_shell_inner_skin_tube = new G4Tubs("mvtx_shell_inner_skin",
                                                  mvtx_shell_inner_skin_inner_radius, mvtx_shell_inner_skin_inner_radius + skin_thickness, mvtx_shell_length / 2.0, -M_PI, 2.0 * M_PI);
  G4LogicalVolume *mvtx_shell_inner_skin_volume = new G4LogicalVolume(mvtx_shell_inner_skin_tube, G4Material::GetMaterial("CFRP_INTT"),
                                                                      "mvtx_shell_inner_skin_volume", 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), mvtx_shell_inner_skin_volume,
                    "mvtx_shell_inner_skin", trackerenvelope, false, 0, OverlapCheck());
  mvtx_shell_inner_skin_volume->SetVisAttributes(rail_vis);

  G4Tubs *mvtx_shell_foam_core_tube = new G4Tubs("mvtx_shell_foam_core",
                                                 mvtx_shell_foam_core_inner_radius, mvtx_shell_foam_core_inner_radius + foam_core_thickness, mvtx_shell_length / 2.0, -M_PI, 2.0 * M_PI);
  G4LogicalVolume *mvtx_shell_foam_core_volume = new G4LogicalVolume(mvtx_shell_foam_core_tube, G4Material::GetMaterial("ROHACELL_FOAM_110"),
                                                                     "mvtx_shell_foam_core_volume", 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), mvtx_shell_foam_core_volume,
                    "mvtx_shell_foam_core", trackerenvelope, false, 0, OverlapCheck());
  mvtx_shell_foam_core_volume->SetVisAttributes(rail_vis);

  G4Tubs *mvtx_shell_outer_skin_tube = new G4Tubs("mvtx_shell_outer_skin",
                                                  mvtx_shell_outer_skin_inner_radius, mvtx_shell_outer_skin_inner_radius + skin_thickness, mvtx_shell_length / 2.0, -M_PI, 2.0 * M_PI);
  G4LogicalVolume *mvtx_shell_outer_skin_volume = new G4LogicalVolume(mvtx_shell_outer_skin_tube, G4Material::GetMaterial("CFRP_INTT"),
                                                                      "mvtx_shell_outer_skin_volume", 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), mvtx_shell_outer_skin_volume,
                    "mvtx_shell_outer_skin", trackerenvelope, false, 0, OverlapCheck());
  mvtx_shell_outer_skin_volume->SetVisAttributes(rail_vis);
  return 0;
}

void PHG4INTTDetector::AddGeometryNode()
{
  int active = 0;
  map<int, int>::const_iterator iter;
  for (iter = m_IsActiveMap.begin(); iter != m_IsActiveMap.end(); ++iter)
  {
    if (iter->second > 0)
    {
      active = 1;
      break;
    }
  }
  if (active)
  {
    std::string geonode = (m_SuperDetector != "NONE") ? (boost::format("CYLINDERGEOM_%s") % m_SuperDetector).str() : (boost::format("CYLINDERGEOM_%s") % m_DetectorType).str();

    PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode);
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode, "PHObject");
      runNode->addNode(newNode);
    }

    for (auto layeriter = m_LayerBeginEndIteratorPair.first; layeriter != m_LayerBeginEndIteratorPair.second; ++layeriter)
    {
      const int sphxlayer = layeriter->first;
      const int inttlayer = layeriter->second;
      int ilayer = inttlayer;
      const PHParameters *params_layer = m_ParamsContainer->GetParameters(inttlayer);
      const int laddertype = params_layer->get_int_param("laddertype");
      // parameters are stored in cm per our convention
      const PHParameters *params = m_ParamsContainer->GetParameters(laddertype);
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomINTT(
          sphxlayer,
          params->get_double_param("strip_x"),
          params->get_double_param("strip_y"),
          params->get_double_param("strip_z_0"),
          params->get_double_param("strip_z_1"),
          params->get_int_param("nstrips_z_sensor_0"),
          params->get_int_param("nstrips_z_sensor_1"),
          params->get_int_param("nstrips_phi_sensor"),
          params_layer->get_int_param("nladder"),
          m_PosZ[ilayer][0] / cm,  // m_PosZ uses G4 internal units, needs to be converted to cm
          m_PosZ[ilayer][1] / cm,
          m_SensorRadius[ilayer] / cm,
          0.0,
          params_layer->get_double_param("offsetphi") * deg / rad,  // expects radians
          params_layer->get_double_param("offsetrot") * deg / rad);
      geo->AddLayerGeom(sphxlayer, mygeom);
      if (Verbosity() > 0)
      {
        geo->identify();
      }
    }
  }
}

map<G4VPhysicalVolume *, std::tuple<int, int, int, int>>::const_iterator
PHG4INTTDetector::get_ActiveVolumeTuple(G4VPhysicalVolume *physvol) const
{
  auto iter = m_ActiveVolumeTuple.find(physvol);
  if (iter == m_ActiveVolumeTuple.end())
  {
    cout << PHWHERE << " Volume " << physvol->GetName() << " not in active volume tuple" << endl;
    gSystem->Exit(1);
  }
  return iter;
}

map<G4LogicalVolume *, std::tuple<int, int>>::const_iterator
PHG4INTTDetector::get_PassiveVolumeTuple(G4LogicalVolume *logvol) const
{
  auto iter = m_PassiveVolumeTuple.find(logvol);
  if (iter == m_PassiveVolumeTuple.end())
  {
    cout << PHWHERE << " Volume " << logvol->GetName() << " not in passive volume tuple" << endl;
    gSystem->Exit(1);
  }
  return iter;
}
