#include "PHG4InttDetector.h"

#include "PHG4InttDefs.h"  // for SEGMENTATION_Z
#include "PHG4InttDisplayAction.h"
#include "PHG4InttFPHXParameterisation.h"

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>      // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4GenericTrap.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVParameterised.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>        // for G4TwoVector
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/geomdefs.hh>           // for kZAxis

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/format.hpp>
#pragma GCC diagnostic pop

#include <algorithm>  // for fill_n
#include <array>
#include <cassert>  // for assert
#include <cmath>
#include <cstdlib>   // for exit, NULL
#include <iostream>  // for operator<<, basic...

class G4VPVParameterisation;
class G4VSolid;

PHG4InttDetector::PHG4InttDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParametersContainer *parameters, const std::string &dnam, const std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> &layer_b_e)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4InttDisplayAction *>(subsys->GetDisplayAction()))
  , m_ParamsContainer(parameters)
  , m_LayerBeginEndIteratorPair(layer_b_e)
{
  for (auto layeriter = m_LayerBeginEndIteratorPair.first; layeriter != m_LayerBeginEndIteratorPair.second; ++layeriter)
  {
    int layer = layeriter->second;
    const PHParameters *par = m_ParamsContainer->GetParameters(layer);
    m_IsActiveMap.insert(std::make_pair(layer, par->get_int_param("active")));
    m_IsAbsorberActiveMap.insert(std::make_pair(layer, par->get_int_param("absorberactive")));
  }
  const PHParameters *par = m_ParamsContainer->GetParameters(PHG4InttDefs::SUPPORTPARAMS);
  m_IsSupportActive = par->get_int_param("supportactive");
  m_IsEndcapActive = par->get_int_param("endcap_ring_enabled");
  std::fill_n(&m_PosZ[0][0], sizeof(m_PosZ) / sizeof(double), NAN);
  std::fill_n(m_SensorRadius, sizeof(m_SensorRadius) / sizeof(double), NAN);
  std::fill_n(m_StripOffsetX, sizeof(m_StripOffsetX) / sizeof(double), NAN);
}

//_______________________________________________________________
int PHG4InttDetector::IsInIntt(G4VPhysicalVolume *volume) const
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

void PHG4InttDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4InttDetector::Construct called for layers " << std::endl;
    for (auto layeriter = m_LayerBeginEndIteratorPair.first; layeriter != m_LayerBeginEndIteratorPair.second; ++layeriter)
    {
      std::cout << "layer " << layeriter->second << std::endl;
    }
  }
  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  ConstructIntt(logicWorld);

  // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int PHG4InttDetector::ConstructIntt(G4LogicalVolume *trackerenvelope)
{
  recoConsts *rc = recoConsts::instance();  // use for worldmaterial in a few places
  // We have an arbitray number of layers (nlayer_) up to 8
  // We have 2 types of ladders (vertical strips and horizontal strips)
  // We have 2 types of sensors (inner and outer)
  std::array<std::array<double, 2>, 8> hdi_z_arr;
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

    double si_glue_x = params->get_double_param("si_glue_x") * cm;
    double fphx_glue_x = params->get_double_param("fphx_glue_x") * cm;
    double halfladder_inside_z = params->get_double_param("halfladder_inside_z") * cm;

    if (Verbosity() > 0)
    {
      std::cout << "Constructing Intt layer: " << std::endl;
      std::cout << "  layer " << inttlayer << " laddertype " << laddertype << " nladders_layer " << nladders_layer
                << " sensor_radius " << m_SensorRadius[inttlayer] << " offsetphi " << offsetphi << " rad "
                << " offsetphi " << offsetphi * rad / deg << " deg "
                << std::endl;
    }
    // We loop over inner, then outer (wrt the beam-axis), sensors, where  itype specifies the inner or outer sensor
    // The rest of this loop will construct and put in place a section of a ladder corresponding to the Z range of this sensor only
    for (int itype = 0; itype < 2; ++itype)
    {
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
        std::cout << "invalid itype " << itype << std::endl;
        exit(1);
      }

      // ----- Step 1 ------------------------------------------------------
      // We make the volumes for Si-sensor, FPHX, HDI, and stave components
      // We add them to the ladder later
      //====================================================================

      // Create Si-sensor active volume
      const double siactive_x = strip_x;
      const double siactive_y = strip_y * nstrips_phi_sensor;
      const double siactive_z = strip_z * nstrips_z_sensor;
      G4VSolid *siactive_box = new G4Box((boost::format("siactive_box_%d_%d") % inttlayer % itype).str(), siactive_x / 2, siactive_y / 2., siactive_z / 2.);
      G4LogicalVolume *siactive_volume = new G4LogicalVolume(siactive_box, GetDetectorMaterial("G4_Si"),
                                                             boost::str(boost::format("siactive_volume_%d_%d") % inttlayer % itype).c_str(), 0, 0, 0);
      if ((m_IsActiveMap.find(inttlayer))->second > 0)
      {
        m_ActiveLogVols.insert(siactive_volume);
      }
      m_DisplayAction->AddVolume(siactive_volume, "SiActive");
      // We do not subdivide the sensor in G4. We will assign hits to strips in the stepping action, using the geometry object

      // Si-sensor full (active+inactive) area
      const double sifull_x = siactive_x;
      const double sifull_y = siactive_y + 2.0 * params->get_double_param("sensor_edge_phi") * cm;
      const double sifull_z = siactive_z + 2.0 * params->get_double_param("sensor_edge_z") * cm;
      G4VSolid *sifull_box = new G4Box((boost::format("sifull_box_%d_%d") % inttlayer % itype).str(), sifull_x / 2., sifull_y / 2.0, sifull_z / 2.0);

      // Si-sensor inactive area
      G4VSolid *siinactive_box = new G4SubtractionSolid((boost::format("siinactive_box_%d_%d") % inttlayer % itype).str(),
                                                        sifull_box, siactive_box, 0, G4ThreeVector(0, 0, 0));
      G4LogicalVolume *siinactive_volume = new G4LogicalVolume(siinactive_box, GetDetectorMaterial("G4_Si"),
                                                               (boost::format("siinactive_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(siinactive_volume, std::make_tuple(inttlayer, PHG4InttDefs::SI_INACTIVE)));
      }
      m_DisplayAction->AddVolume(siinactive_volume, "SiInActive");

      // Glue for Si-sensor full area
      G4VSolid *si_glue_box = new G4Box((boost::format("si_glue_box_%d_%d") % inttlayer % itype).str(), si_glue_x / 2., sifull_y / 2.0, sifull_z / 2.0);

      G4LogicalVolume *si_glue_volume = new G4LogicalVolume(si_glue_box, GetDetectorMaterial("SilverEpoxyGlue_INTT"),
                                                            (boost::format("si_glue_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(siinactive_volume, std::make_tuple(inttlayer, PHG4InttDefs::SI_GLUE)));
      }
      m_DisplayAction->AddVolume(si_glue_volume, "SiGlue");

      // Make the HDI Kapton and copper volumes
      // This makes HDI volumes that matche this sensor in Z length
      const double hdi_z = sifull_z + params->get_double_param("hdi_edge_z") * cm;
      hdi_z_arr[inttlayer][itype] = hdi_z;
      G4VSolid *hdi_kapton_box = new G4Box((boost::format("hdi_kapton_box_%d_%d") % inttlayer % itype).str(), hdi_kapton_x / 2., hdi_y / 2., hdi_z / 2.0);
      G4LogicalVolume *hdi_kapton_volume = new G4LogicalVolume(hdi_kapton_box, GetDetectorMaterial("G4_KAPTON"),
                                                               (boost::format("hdi_kapton_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(hdi_kapton_volume, std::make_tuple(inttlayer, PHG4InttDefs::HDI_KAPTON)));
      }
      G4VSolid *hdi_copper_box = new G4Box((boost::format("hdi_copper_box_%d_%d") % inttlayer % itype).str(), hdi_copper_x / 2., hdi_y / 2., hdi_z / 2.0);
      G4LogicalVolume *hdi_copper_volume = new G4LogicalVolume(hdi_copper_box, GetDetectorMaterial("G4_Cu"),
                                                               (boost::format("hdi_copper_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(hdi_copper_volume, std::make_tuple(inttlayer, PHG4InttDefs::HDI_COPPER)));
      }
      m_DisplayAction->AddVolume(hdi_kapton_volume, "HdiKapton");
      m_DisplayAction->AddVolume(hdi_copper_volume, "HdiCopper");

      // This is the part of the HDI that extends beyond the sensor inside the endcap ring
      const double hdiext_z = (itype == 0) ? 0.000001 : halfladder_inside_z - hdi_z_arr[inttlayer][0] - hdi_z;  // need to assign nonzero value for itype=0
      G4VSolid *hdiext_kapton_box = new G4Box((boost::format("hdiext_kapton_box_%d_%s") % inttlayer % itype).str(),
                                              hdi_kapton_x / 2., hdi_y / 2., hdiext_z / 2.0);
      G4LogicalVolume *hdiext_kapton_volume = new G4LogicalVolume(hdiext_kapton_box, GetDetectorMaterial("G4_KAPTON"),  // was "FPC"
                                                                  (boost::format("hdiext_kapton_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(hdiext_kapton_volume, std::make_tuple(inttlayer, PHG4InttDefs::HDIEXT_KAPTON)));
      }
      G4VSolid *hdiext_copper_box = new G4Box((boost::format("hdiext_copper_box_%d_%s") % inttlayer % itype).str(),
                                              hdi_copper_x / 2., hdi_y / 2., hdiext_z / 2.0);
      G4LogicalVolume *hdiext_copper_volume = new G4LogicalVolume(hdiext_copper_box, GetDetectorMaterial("G4_Cu"),
                                                                  (boost::format("hdiext_copper_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(hdiext_copper_volume, std::make_tuple(inttlayer, PHG4InttDefs::HDIEXT_COPPER)));
      }
      m_DisplayAction->AddVolume(hdiext_kapton_volume, "HdiKapton");
      m_DisplayAction->AddVolume(hdiext_copper_volume, "HdiCopper");

      // FPHX
      G4VSolid *fphx_box = new G4Box((boost::format("fphx_box_%d_%d") % inttlayer % itype).str(), fphx_x / 2., fphx_y / 2., fphx_z / 2.);
      G4LogicalVolume *fphx_volume = new G4LogicalVolume(fphx_box, GetDetectorMaterial("G4_Si"),
                                                         (boost::format("fphx_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(fphx_volume, std::make_tuple(inttlayer, PHG4InttDefs::FPHX)));
      }
      m_DisplayAction->AddVolume(fphx_volume, "FPHX");

      const double gap_sensor_fphx = params->get_double_param("gap_sensor_fphx") * cm;

      //  FPHX Container
      // make a container for the FPHX chips needed for this sensor, and  then place them in the container
      G4VSolid *fphxcontainer_box = new G4Box((boost::format("fphxcontainer_box_%d_%d") % inttlayer % itype).str(),
                                              fphx_x / 2., fphx_y / 2., hdi_z / 2.);
      G4LogicalVolume *fphxcontainer_volume = new G4LogicalVolume(fphxcontainer_box, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")),
                                                                  (boost::format("fphxcontainer_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      m_DisplayAction->AddVolume(fphxcontainer_volume, "FPHXContainer");

      // Install multiple FPHX volumes in the FPHX container volume
      // one FPHX chip per cell - each cell is 128 channels
      const double fphx_offsetx = 0.;
      const double fphx_offsety = 0.;
      int ncopy;
      double offsetz, cell_length_z;

      if (laddertype == PHG4InttDefs::SEGMENTATION_Z)  // vertical strips
      {
        // For laddertype 0, we have 5 cells per sensor, but the strips are vertical, so we have to treat it specially
        ncopy = nstrips_z_sensor / 128.0;
      }
      else if (laddertype == PHG4InttDefs::SEGMENTATION_PHI)
      {
        ncopy = nstrips_z_sensor;
      }
      else
      {
        std::cout << PHWHERE << "invalid laddertype " << laddertype << std::endl;
        gSystem->Exit(1);
        // this is just to make the optimizer happy which otherwise complains about possibly
        // uninitialized variables. It doesn't know gSystem->Exit(1) quits,
        // this exit here terminates the program for it
        exit(1);
      }
      cell_length_z = strip_z * nstrips_z_sensor / ncopy;
      offsetz = (ncopy % 2 == 0) ? -2. * cell_length_z / 2. * double(ncopy / 2) + cell_length_z / 2. + fphx_offset_z : -2. * cell_length_z / 2. * double(ncopy / 2) + fphx_offset_z;

      G4VPVParameterisation *fphxparam = new PHG4InttFPHXParameterisation(fphx_offsetx, +fphx_offsety, offsetz, 2. * cell_length_z / 2., ncopy);
      new G4PVParameterised((boost::format("fphxcontainer_%d_%d") % inttlayer % itype).str(),
                            fphx_volume, fphxcontainer_volume, kZAxis, ncopy, fphxparam, OverlapCheck());

      // Glue for FPHX, silver powder epoxy, impletemented in the same way as FPHX
      G4VSolid *fphx_glue_box = new G4Box((boost::format("fphx_glue_box_%d_%d") % inttlayer % itype).str(), fphx_glue_x / 2., fphx_y / 2., fphx_z / 2.);

      G4LogicalVolume *fphx_glue_volume = new G4LogicalVolume(fphx_glue_box, GetDetectorMaterial("SilverEpoxyGlue_INTT"),
                                                              (boost::format("fphx_glue_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(fphx_glue_volume, std::make_tuple(inttlayer, PHG4InttDefs::FPHX_GLUE)));
      }
      m_DisplayAction->AddVolume(fphx_glue_volume, "FPHXGlue");

      //  Glue of FPHX Container
      // make a container for the glue of FPHX chips, and then place them in the container
      G4VSolid *fphx_gluecontainer_box = new G4Box((boost::format("fphx_gluecontainer_box_%d_%d") % inttlayer % itype).str(),
                                                   fphx_glue_x / 2., fphx_y / 2., hdi_z / 2.);
      G4LogicalVolume *fphx_gluecontainer_volume = new G4LogicalVolume(fphx_gluecontainer_box, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")),
                                                                       (boost::format("fphx_gluecontainer_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      // Parameters for FPHX glue for G4VPVParameterisation are the same as FPGX's, so reuse them!
      G4VPVParameterisation *fphx_glueparam = new PHG4InttFPHXParameterisation(fphx_offsetx, +fphx_offsety, offsetz, 2. * cell_length_z / 2., ncopy);

      new G4PVParameterised((boost::format("glue_fphxcontainer_%d_%d") % inttlayer % itype).str(),
                            fphx_glue_volume, fphx_gluecontainer_volume, kZAxis, ncopy, fphx_glueparam, OverlapCheck());
      m_DisplayAction->AddVolume(fphx_gluecontainer_volume, "FPHXGlueContainer");

      double stave_x = 0.;
      double stave_y = 0.;
      G4LogicalVolume *stave_volume = NULL;
      G4LogicalVolume *staveext_volume = NULL;

      // Carbon stave. This consists of the formed sheet, cooling water, the water tube, glue for the tube,
      // rohacell foam to fill space around the tube, and the flat CFRP sheet, which completes the outer shell surrounds

      // Rohacel foam and cooling water pipe inside. Formed from straight sections and sections of a tube of
      // radius 3.1905 mm. All have wall thickness of 0.1905 mm.
      const double stave_thickness = params->get_double_param("stave_straight_cooler_x") * cm;  // stave thickness
      const double Rcmin = 0.30 * cm;                                                           // inner radius of curved section, same at both ends
      const double Rcmax = Rcmin + stave_thickness;                                             // outer radius of curved section, same at both ends
      double Rcavge = (Rcmax + Rcmin) / 2.0;
      double dphi_c = 23.19859051 * M_PI / 180.;  // phi of the curved section
      const double stave_z = hdi_z;

      // Make CFC structure
      //// Make curved sections
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
        stave_curve_volume[i] = new G4LogicalVolume(stave_curve_cons[i], GetDetectorMaterial("CFRP_INTT"),
                                                    (boost::format("stave_curve_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);
        if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
        {
          m_PassiveVolumeTuple.insert(std::make_pair(stave_curve_volume[i], std::make_tuple(inttlayer, PHG4InttDefs::STAVE_CURVE)));
        }
        stave_curve_ext_cons[i] = new G4Tubs((boost::format("stave_curve_ext_cons_%d_%d_%d") % inttlayer % itype % i).str(),
                                             Rcmin, Rcmax, hdiext_z / 2., phic_begin[i], dphic[i]);
        stave_curve_ext_volume[i] = new G4LogicalVolume(stave_curve_ext_cons[i], GetDetectorMaterial("CFRP_INTT"),
                                                        (boost::format("stave_curve_ext_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);
        if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
        {
          m_PassiveVolumeTuple.insert(std::make_pair(stave_curve_ext_volume[i], std::make_tuple(inttlayer, PHG4InttDefs::STAVEEXT_CURVE)));
        }
        m_DisplayAction->AddVolume(stave_curve_volume[i], "StaveCurve");
        m_DisplayAction->AddVolume(stave_curve_ext_volume[i], "StaveCurve");
      }

      // we will need the length in y of the curved section as it is installed in the stave box
      double curve_length_y = Rcavge * sin(dphi_c);

      // Make several straight sections for use in making the stave
      double stave_straight_outer_y = params->get_double_param("stave_straight_outer_y") * cm;
      double stave_straight_cooler_y = params->get_double_param("stave_straight_cooler_y") * cm;
      double rohacell_straight_y = params->get_double_param("stave_straight_rohacell_y") * cm;

      // Outer straight sections of stave
      G4VSolid *stave_straight_outer_box = new G4Box((boost::format("stave_straight_outer_box_%d_%d") % inttlayer % itype).str(),
                                                     stave_thickness / 2., stave_straight_outer_y / 2., stave_z / 2.);
      G4LogicalVolume *stave_straight_outer_volume = new G4LogicalVolume(stave_straight_outer_box, GetDetectorMaterial("CFRP_INTT"),
                                                                         (boost::format("stave_straight_outer_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_straight_outer_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVE_STRAIGHT_OUTER)));
      }
      G4VSolid *stave_straight_outer_ext_box = new G4Box((boost::format("stave_straight_outer_ext_box_%d_%s") % inttlayer % itype).str(),
                                                         stave_thickness / 2., stave_straight_outer_y / 2., hdiext_z / 2.);
      G4LogicalVolume *stave_straight_outer_ext_volume = new G4LogicalVolume(stave_straight_outer_ext_box, GetDetectorMaterial("CFRP_INTT"),
                                                                             (boost::format("stave_straight_outer_ext_volume_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_straight_outer_ext_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVEEXT_STRAIGHT_OUTER)));
      }

      //Top surface of stave
      G4VSolid *stave_straight_cooler_box = new G4Box((boost::format("stave_straight_cooler_box_%d_%d") % inttlayer % itype).str(),
                                                      stave_thickness / 2., stave_straight_cooler_y / 2., stave_z / 2.);
      G4LogicalVolume *stave_straight_cooler_volume = new G4LogicalVolume(stave_straight_cooler_box, GetDetectorMaterial("CFRP_INTT"),
                                                                          (boost::format("stave_straight_cooler_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_straight_cooler_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVE_STRAIGHT_COOLER)));
      }
      G4VSolid *stave_straight_cooler_ext_box = new G4Box((boost::format("stave_straight_cooler_ext_box_%d_%d") % inttlayer % itype).str(),
                                                          stave_thickness / 2., stave_straight_cooler_y / 2., hdiext_z / 2.);
      G4LogicalVolume *stave_straight_cooler_ext_volume = new G4LogicalVolume(stave_straight_cooler_ext_box, GetDetectorMaterial("CFRP_INTT"),
                                                                              (boost::format("stave_straight_cooler_ext_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_straight_cooler_ext_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVEEXT_STRAIGHT_COOLER)));
      }

      // Slant straight sections of stave
      double stave_slant_cooler_y = params->get_double_param("stave_slant_cooler_y") * cm;
      G4VSolid *stave_slant_cooler_box = new G4Box((boost::format("stave_slant_cooler_box_%d_%d") % inttlayer % itype).str(),
                                                   stave_thickness / 2., stave_slant_cooler_y / 2., stave_z / 2.);
      G4LogicalVolume *stave_slant_cooler_volume = new G4LogicalVolume(stave_slant_cooler_box, GetDetectorMaterial("CFRP_INTT"),
                                                                       (boost::format("stave_slant_cooler_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_slant_cooler_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVE_STRAIGHT_COOLER)));
      }
      G4VSolid *stave_slant_cooler_ext_box = new G4Box((boost::format("stave_lant_cooler_ext_box_%d_%d") % inttlayer % itype).str(),
                                                       stave_thickness / 2., stave_slant_cooler_y / 2., hdiext_z / 2.);
      G4LogicalVolume *stave_slant_cooler_ext_volume = new G4LogicalVolume(stave_slant_cooler_ext_box, GetDetectorMaterial("CFRP_INTT"),
                                                                           (boost::format("stave_slant_cooler_ext_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_slant_cooler_ext_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVEEXT_STRAIGHT_COOLER)));
      }

      // Flat CFRP sheet on the bottom of the stave structure. It was introduced instead of PGS
      G4VSolid *stave_bottom_cooler_box = new G4Box((boost::format("stave_bottom_cooler_box_%d_%d") % inttlayer % itype).str(),
                                                    stave_thickness / 2., hdi_y / 2., stave_z / 2.);

      G4LogicalVolume *stave_bottom_cooler_volume = new G4LogicalVolume(stave_bottom_cooler_box, GetDetectorMaterial("CFRP_INTT"),
                                                                        (boost::format("stave_bottom_cooler_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_bottom_cooler_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVE_BOTTOM_COOLER)));  // should be changed soon
      }

      G4VSolid *stave_bottom_cooler_ext_box = new G4Box((boost::format("stave_bottom_cooler_ext_box_%d_%s") % inttlayer % itype).str(), stave_thickness / 2., hdi_y / 2., hdiext_z / 2.);
      G4LogicalVolume *stave_bottom_cooler_ext_volume = new G4LogicalVolume(stave_bottom_cooler_ext_box, GetDetectorMaterial("CFRP_INTT"),
                                                                            (boost::format("stave_bottom_cooler_ext_volume_%d_%s") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(std::make_pair(stave_bottom_cooler_ext_volume, std::make_tuple(inttlayer, PHG4InttDefs::STAVEEXT_BOTTOM_COOLER)));
      }

      m_DisplayAction->AddVolume(stave_straight_cooler_volume, "StaveCooler");
      m_DisplayAction->AddVolume(stave_straight_cooler_ext_volume, "StaveCooler");
      m_DisplayAction->AddVolume(stave_straight_outer_volume, "StaveStraightOuter");
      m_DisplayAction->AddVolume(stave_straight_outer_ext_volume, "StaveStraightOuter");
      m_DisplayAction->AddVolume(stave_slant_cooler_volume, "StaveCooler");
      m_DisplayAction->AddVolume(stave_slant_cooler_ext_volume, "StaveCooler");
      m_DisplayAction->AddVolume(stave_bottom_cooler_volume, "StaveCooler");
      m_DisplayAction->AddVolume(stave_bottom_cooler_ext_volume, "StaveCooler");

      // cooling pipe + water inside + glue outside
      const double Rpmin = 0.10 * cm;  // inner radius of cooling pipe section, same at both ends
      const double Rpmax = 0.15 * cm;  // outer radius of cooling pipe section, same at both ends
      G4VSolid *stave_glue_box = new G4Box((boost::format("stave_glue_box_%d_%d") % inttlayer % itype).str(), 3. / 2, 3. / 2., stave_z / 2.);
      G4LogicalVolume *stave_glue_volume = new G4LogicalVolume(stave_glue_box, GetDetectorMaterial("Epoxy"),
                                                               (boost::format("stave_glue_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      G4VSolid *staveext_glue_box = new G4Box((boost::format("staveext_glue_box_%d_%d") % inttlayer % itype).str(), 3. / 2., 3. / 2., hdiext_z / 2.);
      G4LogicalVolume *staveext_glue_volume = new G4LogicalVolume(staveext_glue_box, GetDetectorMaterial("Epoxy"),
                                                                  (boost::format("staveext_glue_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      m_DisplayAction->AddVolume(stave_glue_volume, "StaveGlueBox");
      m_DisplayAction->AddVolume(staveext_glue_volume, "StaveGlueBox");

      G4VSolid *stave_pipe_cons = new G4Tubs((boost::format("stave_pipe_cons_%d_%d") % inttlayer % itype).str(),
                                             Rpmin, Rpmax, stave_z / 2., 0, 2.0 * M_PI);
      G4LogicalVolume *stave_pipe_volume = new G4LogicalVolume(stave_pipe_cons, GetDetectorMaterial("CFRP_INTT"),
                                                               (boost::format("stave_pipe_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      G4VSolid *staveext_pipe_cons = new G4Tubs((boost::format("staveext_pipe_cons_%d_%d") % inttlayer % itype).str(),
                                                Rpmin, Rpmax, hdiext_z / 2., 0, 2.0 * M_PI);
      G4LogicalVolume *staveext_pipe_volume = new G4LogicalVolume(staveext_pipe_cons, GetDetectorMaterial("CFRP_INTT"),
                                                                  (boost::format("staveext_pipe_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      m_DisplayAction->AddVolume(stave_pipe_volume, "StavePipe");
      m_DisplayAction->AddVolume(staveext_pipe_volume, "StavePipe");

      G4VSolid *stave_water_cons = new G4Tubs((boost::format("stave_water_cons_%d_%d") % inttlayer % itype).str(),
                                              0., Rpmin, stave_z / 2., 0, 2.0 * M_PI);
      G4LogicalVolume *stave_water_volume = new G4LogicalVolume(stave_water_cons, GetDetectorMaterial("G4_WATER"),
                                                                (boost::format("stave_water_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      G4VSolid *staveext_water_cons = new G4Tubs((boost::format("staveext_water_cons_%d_%d") % inttlayer % itype).str(),
                                                 0., Rpmin, hdiext_z / 2., 0, 2.0 * M_PI);
      G4LogicalVolume *staveext_water_volume = new G4LogicalVolume(staveext_water_cons, GetDetectorMaterial("G4_WATER"),
                                                                   (boost::format("staveext_water_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      m_DisplayAction->AddVolume(stave_water_volume, "StaveWater");
      m_DisplayAction->AddVolume(staveext_water_volume, "StaveWater");

      //rohacell foam
      //straight boxes
      G4VSolid *rohacell_straight_cons = new G4Box((boost::format("rohacell_straight_cons_%d_%d") % inttlayer % itype).str(), 3. / 2, rohacell_straight_y / 2., stave_z / 2.);
      G4LogicalVolume *rohacell_straight_volume = new G4LogicalVolume(rohacell_straight_cons, GetDetectorMaterial("ROHACELL_FOAM_51"),
                                                                      (boost::format("rohacell_straight_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      G4VSolid *rohacellext_straight_cons = new G4Box((boost::format("rohacellext_straight_cons_%d_%d") % inttlayer % itype).str(), 3. / 2, rohacell_straight_y / 2., hdiext_z / 2.);
      G4LogicalVolume *rohacellext_straight_volume = new G4LogicalVolume(rohacellext_straight_cons, GetDetectorMaterial("ROHACELL_FOAM_51"),
                                                                         (boost::format("rohacellext_straight_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      // make curved sections for rohacell foam
      const double rh_phic_begin[2] = {-dphi_c, 0.0};
      const double rh_dphic[2] = {dphi_c, dphi_c};
      G4Tubs *rohacell_curve_cons[2];
      G4LogicalVolume *rohacell_curve_volume[2];
      G4Tubs *rohacellext_curve_cons[2];
      G4LogicalVolume *rohacellext_curve_volume[2];
      for (int i = 0; i < 2; i++)
      {
        rohacell_curve_cons[i] = new G4Tubs((boost::format("rohacell_curve_cons_%d_%d_%d") % inttlayer % itype % i).str(),
                                            0., Rcmin, stave_z / 2., rh_phic_begin[i], rh_dphic[i]);
        rohacell_curve_volume[i] = new G4LogicalVolume(rohacell_curve_cons[i], GetDetectorMaterial("ROHACELL_FOAM_51"),
                                                       (boost::format("rohacell_curve_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);
        rohacellext_curve_cons[i] = new G4Tubs((boost::format("rohacellext_curve_cons_%d_%d_%d") % inttlayer % itype % i).str(),
                                               0., Rcmin, hdiext_z / 2., rh_phic_begin[i], rh_dphic[i]);
        rohacellext_curve_volume[i] = new G4LogicalVolume(rohacellext_curve_cons[i], GetDetectorMaterial("ROHACELL_FOAM_51"),
                                                          (boost::format("rohacellext_curve_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);
      }

      // make trapezoidal sections for rohacell foam
      G4GenericTrap *rohacell_trap_cons[2];
      G4LogicalVolume *rohacell_trap_volume[2];
      G4GenericTrap *rohacellext_trap_cons[2];
      G4LogicalVolume *rohacellext_trap_volume[2];
      for (int i = 0; i < 2; i++)
      {
        double shift = 1.e-5;  // To mitigate fm order level overlaps reported by GEANT4...
        std::vector<G4TwoVector> rohatrap(8);
        if (i == 0)
        {
          rohatrap[0] = G4TwoVector(0. * cm, 0. * cm);
          rohatrap[1] = G4TwoVector(Rcmin * cos(dphi_c) - shift, -Rcmin * sin(dphi_c));
          rohatrap[2] = G4TwoVector(Rcmin * (1. - cos(dphi_c)) - shift, -stave_slant_cooler_y * cos(dphi_c) - Rcmin * sin(dphi_c));
          rohatrap[3] = G4TwoVector(0. * cm, -stave_slant_cooler_y * cos(dphi_c) - Rcmin * sin(dphi_c));
        }
        else
        {
          rohatrap[0] = G4TwoVector(0. * cm, +stave_slant_cooler_y * cos(dphi_c) + Rcmin * sin(dphi_c));
          rohatrap[1] = G4TwoVector(Rcmax * (1. - cos(dphi_c)) - shift, +stave_slant_cooler_y * cos(dphi_c) + Rcmin * sin(dphi_c));
          rohatrap[2] = G4TwoVector(Rcmin * cos(dphi_c) - shift, +Rcmin * sin(dphi_c));
          rohatrap[3] = G4TwoVector(0. * cm, 0. * cm);
        }
        rohatrap[4] = rohatrap[0];
        rohatrap[5] = rohatrap[1];
        rohatrap[6] = rohatrap[2];
        rohatrap[7] = rohatrap[3];

        rohacell_trap_cons[i] = new G4GenericTrap((boost::format("rohacell_trap_cons_%d_%d_%d") % inttlayer % itype % i).str(), stave_z / 2., rohatrap);
        rohacell_trap_volume[i] = new G4LogicalVolume(rohacell_trap_cons[i], GetDetectorMaterial("ROHACELL_FOAM_51"),
                                                      (boost::format("rohacell_trap_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);

        rohacellext_trap_cons[i] = new G4GenericTrap((boost::format("rohacellext_trap_cons_%d_%d_%d") % inttlayer % itype % i).str(), hdiext_z / 2., rohatrap);
        rohacellext_trap_volume[i] = new G4LogicalVolume(rohacellext_trap_cons[i], GetDetectorMaterial("ROHACELL_FOAM_51"),
                                                         (boost::format("rohacellext_trap_volume_%d_%d_%d") % inttlayer % itype % i).str(), 0, 0, 0);
      }

      m_DisplayAction->AddVolume(rohacell_straight_volume, "RohaCell");
      m_DisplayAction->AddVolume(rohacellext_straight_volume, "RohaCell");
      for (int i = 0; i < 2; i++)
      {
        m_DisplayAction->AddVolume(rohacell_curve_volume[i], "RohaCell");
        m_DisplayAction->AddVolume(rohacellext_curve_volume[i], "RohaCell");
        m_DisplayAction->AddVolume(rohacell_trap_volume[i], "RohaCell");
        m_DisplayAction->AddVolume(rohacellext_trap_volume[i], "RohaCell");
      }

      // Now we combine the elements of a stave defined above into a stave
      // Create a stave volume to install the stave sections into. The volume has to be big enouigh to contain the cooling tube
      double cooler_gap_x = 0.3 * cm;                      // id of cooling tube in cm
      double cooler_wall = stave_thickness;                // outer wall thickness of cooling tube
      double cooler_x = cooler_gap_x + 2.0 * cooler_wall;  // thickness of the formed sheet, the flat sheet, and the gap b/w the sheets
      stave_x = cooler_x;
      stave_y = hdi_y;

      // Make stave volume. Drop two corners in positive x to prevent ladder_volume overlapping
      // with neighbouring ladders because of small clearance in the latest configuration
      G4RotationMatrix *stv_rot_pos = new G4RotationMatrix();
      stv_rot_pos->rotateZ(-15. * M_PI / 180.);
      G4ThreeVector stvTranspos(stave_x / 2., stave_y / 2., 0.);

      G4RotationMatrix *stv_rot_neg = new G4RotationMatrix();
      stv_rot_neg->rotateZ(+15. * M_PI / 180.);
      G4ThreeVector stvTransneg(stave_x / 2., -stave_y / 2., 0.);

      G4VSolid *stave_basebox = new G4Box((boost::format("stave_basebox_%d_%d") % inttlayer % itype).str(), stave_x / 2., stave_y / 2., stave_z / 2.);
      G4VSolid *stave_subtbox = new G4Box((boost::format("stave_subtbox_%d_%d") % inttlayer % itype).str(), stave_x / 1.5, stave_y / 1.5, stave_z / 1.);  // has to be longer in z to avoid coincident surface

      G4VSolid *stave_box1 = new G4SubtractionSolid((boost::format("stave_box1_%d_%d") % inttlayer % itype).str(), stave_basebox, stave_subtbox, stv_rot_pos, stvTranspos);

      G4VSolid *stave_box = new G4SubtractionSolid((boost::format("stave_box_%d_%d") % inttlayer % itype).str(), stave_box1, stave_subtbox, stv_rot_neg, stvTransneg);

      stave_volume = new G4LogicalVolume(stave_box, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")),
                                         (boost::format("stave_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      G4VSolid *staveext_basebox = new G4Box((boost::format("staveext_basebox_%d_%d") % inttlayer % itype).str(), stave_x / 2., stave_y / 2., hdiext_z / 2.);
      G4VSolid *staveext_subtbox = new G4Box((boost::format("staveext_subtbox_%d_%d") % inttlayer % itype).str(), stave_x / 1.5, stave_y / 1.5, hdiext_z / 1.);  // has to be longer in z to avoid coincident surface

      G4VSolid *staveext_box1 = new G4SubtractionSolid((boost::format("staveext_box1_%d_%d") % inttlayer % itype).str(), staveext_basebox, staveext_subtbox, stv_rot_pos, stvTranspos);

      G4VSolid *staveext_box = new G4SubtractionSolid((boost::format("staveext_box_%d_%d") % inttlayer % itype).str(), staveext_box1, staveext_subtbox, stv_rot_neg, stvTransneg);

      staveext_volume = new G4LogicalVolume(staveext_box, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")),
                                            (boost::format("staveext_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      // the rotation matrices are just used by G4VSolid, ownership is not taken over
      delete stv_rot_pos;
      delete stv_rot_neg;

      m_DisplayAction->AddVolume(stave_volume, "StaveBox");
      m_DisplayAction->AddVolume(staveext_volume, "StaveBox");

      // Assemble the elements into the stave volume and the stave extension volume
      // They are place relative to the center of the stave box. Thus the offset of the center of the segment is relative to the center of the satev box.
      // But we want the segment to be located relative to the lowest x limit of the stave box.
      if (laddertype == PHG4InttDefs::SEGMENTATION_Z)  // Obsolete!!
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
                +stave_straight_cooler_y / 2. + 2. * curve_length_y + stave_straight_outer_y / 2.,
                -stave_straight_cooler_y / 2. - 2. * curve_length_y - stave_straight_outer_y / 2.};

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
                +stave_straight_cooler_y / 2.,
                +stave_straight_cooler_y / 2. + 2. * curve_length_y};

        for (int i = 0; i < 4; i++)
        {
          new G4PVPlacement(0, G4ThreeVector(x_off_cooler[i], y_off_cooler[i], 0.0), stave_curve_volume[i],
                            (boost::format("stave_curve_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector(x_off_cooler[i], y_off_cooler[i], 0.0), stave_curve_ext_volume[i],
                            (boost::format("stave_curve_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
        }
      }
      else if (laddertype == PHG4InttDefs::SEGMENTATION_PHI)  // The type PHG4InttDefs::SEGMENTATION_PHI ladder
      {
        // First place the straight sections, do the extension at the same time
        // we alternate  positive and negative y values here
        double x_off_str[6] =
            {
                (Rcmax + Rcmin) / 2. - stave_x / 2. + stave_thickness,  // inner straight section
                (Rcmax + Rcmax) / 4. - stave_x / 2. + stave_thickness,  // slant section
                (Rcmax + Rcmax) / 4. - stave_x / 2. + stave_thickness,  // slant section
                (Rcmax - Rcmin) / 2. - stave_x / 2. + stave_thickness,  // outer straight section
                (Rcmax - Rcmin) / 2. - stave_x / 2. + stave_thickness,  // outer straight section
                (Rcmax - Rcmin) / 2. - stave_x / 2.                     // bottom section
            };
        double y_off_str[6] =
            {
                0.0,                                                                                                                     // inner straight section
                -stave_straight_cooler_y / 2. - 1. * curve_length_y - cos(dphi_c) * stave_slant_cooler_y / 2.,                           // slant section
                +stave_straight_cooler_y / 2. + 1. * curve_length_y + cos(dphi_c) * stave_slant_cooler_y / 2.,                           // slant section
                -stave_straight_cooler_y / 2. - 2. * curve_length_y - cos(dphi_c) * stave_slant_cooler_y - stave_straight_outer_y / 2.,  // outer straight section
                +stave_straight_cooler_y / 2. + 2. * curve_length_y + cos(dphi_c) * stave_slant_cooler_y + stave_straight_outer_y / 2.,  // outer straight section
                0.0
                // bottom straight section
            };

        for (int i = 0; i < 6; i++)
        {
          if (i == 0)  // inner straight section
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_volume,
                              (boost::format("stave_straight_cooler_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_cooler_ext_volume,
                              (boost::format("stave_straight_cooler_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
          }
          else if (i == 1 || i == 2)  // slant section
          {
            G4RotationMatrix rotation;
            if (i == 1)
              rotation.rotateZ(-1. * dphi_c);
            else if (i == 2)
              rotation.rotateZ(dphi_c);
            new G4PVPlacement(G4Transform3D(rotation, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0)), stave_slant_cooler_volume,
                              (boost::format("stave_slant_cooler_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(G4Transform3D(rotation, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0)), stave_slant_cooler_ext_volume,
                              (boost::format("stave_slant_cooler_ext_%d_%d_%d") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
          }
          else if (i == 3 || i == 4)  // outer straight section
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_volume,
                              (boost::format("stave_straight_outer_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_straight_outer_ext_volume,
                              (boost::format("stave_straight_outer_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
          }
          else  // bottom
          {
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_bottom_cooler_volume,
                              (boost::format("stave_bottom_cooler_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
            new G4PVPlacement(0, G4ThreeVector(x_off_str[i], y_off_str[i], 0.0), stave_bottom_cooler_ext_volume,
                              (boost::format("stave_bottom_cooler_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
          }
        }

        //// Place the curved sections
        //// here we do in order of increasing y

        double x_off_curve[4] =
            {
                // increasing in y
                +Rcavge - cooler_gap_x / 2. + stave_thickness / 2.,
                -Rcavge + cooler_gap_x / 2. + stave_thickness / 2.,
                -Rcavge + cooler_gap_x / 2. + stave_thickness / 2.,
                +Rcavge - cooler_gap_x / 2. + stave_thickness / 2.};
        double y_off_curve[4] =
            {
                // increasing in y
                -stave_straight_cooler_y / 2. - 2. * curve_length_y - cos(dphi_c) * stave_slant_cooler_y,
                -stave_straight_cooler_y / 2.,
                +stave_straight_cooler_y / 2.,
                +stave_straight_cooler_y / 2. + 2. * curve_length_y + cos(dphi_c) * stave_slant_cooler_y};

        for (int i = 0; i < 4; i++)
        {
          new G4PVPlacement(0, G4ThreeVector(x_off_curve[i], y_off_curve[i], 0.0), stave_curve_volume[i], (boost::format("stave_curve_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector(x_off_curve[i], y_off_curve[i], 0.0), stave_curve_ext_volume[i], (boost::format("stave_curve_ext_%d_%d_%s") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
        }

        // Place the rohacell foam
        // straight box section
        double x_off_roha_str[2] =
            {
                // increasing in y
                -cooler_wall / 2. + stave_thickness / 2.,
                -cooler_wall / 2. + stave_thickness / 2.};
        double y_off_roha_str[2] =
            {
                // increasing in y
                -3. / 2. - rohacell_straight_y / 2.,
                +3. / 2. + rohacell_straight_y / 2.};

        for (int i = 0; i < 2; i++)
        {
          new G4PVPlacement(0, G4ThreeVector(x_off_roha_str[i], y_off_roha_str[i], 0.0), rohacell_straight_volume, (boost::format("rohacell_straight_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector(x_off_roha_str[i], y_off_roha_str[i], 0.0), rohacellext_straight_volume, (boost::format("rohacell_straight_ext_%d_%d_%d") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
        }

        //// curve section
        double x_off_roha_curve[2] =
            {
                // increasing in y
                -Rcavge + cooler_gap_x / 2. + stave_thickness / 2.,
                -Rcavge + cooler_gap_x / 2. + stave_thickness / 2.};
        double y_off_roha_curve[2] =
            {
                // increasing in y
                -3. / 2. - rohacell_straight_y,
                +3. / 2. + rohacell_straight_y};

        for (int i = 0; i < 2; i++)
        {
          new G4PVPlacement(0, G4ThreeVector(x_off_roha_curve[i], y_off_roha_curve[i], 0.0), rohacell_curve_volume[i], (boost::format("rohacell_curve_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector(x_off_roha_curve[i], y_off_roha_curve[i], 0.0), rohacellext_curve_volume[i], (boost::format("rohacell_curve_ext_%d_%d_%d") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
        }

        // trapezoidal section
        double x_off_roha_trap[2] =
            {
                // increasing in y
                -Rcmin - cooler_wall / 2. + cooler_gap_x / 2. + stave_thickness / 2.,
                -Rcmin - cooler_wall / 2. + cooler_gap_x / 2. + stave_thickness / 2.};
        double y_off_roha_trap[2] =
            {
                // increasing in y
                -3. / 2. - rohacell_straight_y,
                +3. / 2. + rohacell_straight_y};

        for (int i = 0; i < 2; i++)
        {
          new G4PVPlacement(0, G4ThreeVector(x_off_roha_trap[i], y_off_roha_trap[i], 0.0), rohacell_trap_volume[i], (boost::format("rohacell_trap_%d_%d_%d") % inttlayer % itype % i).str(), stave_volume, false, 0, OverlapCheck());
          new G4PVPlacement(0, G4ThreeVector(x_off_roha_trap[i], y_off_roha_trap[i], 0.0), rohacellext_trap_volume[i], (boost::format("rohacell_trap_ext_%d_%d_%d") % inttlayer % itype % i).str(), staveext_volume, false, 0, OverlapCheck());
        }

        // place glue box, cooling pipe and water inside
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), stave_pipe_volume, (boost::format("stave_pipe_%d_%d") % inttlayer % itype).str(), stave_glue_volume, false, 0, OverlapCheck());
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), staveext_pipe_volume, (boost::format("stave_pipe_ext_%d_%d") % inttlayer % itype).str(), staveext_glue_volume, false, 0, OverlapCheck());
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), stave_water_volume, (boost::format("stave_water_%d_%d") % inttlayer % itype).str(), stave_glue_volume, false, 0, OverlapCheck());
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), staveext_water_volume, (boost::format("stave_water_ext_%d_%d") % inttlayer % itype).str(), staveext_glue_volume, false, 0, OverlapCheck());

        // place of stave_glue_volume -cooler_wall / 2. + stave_thickness / 2. is actually 0. But I don't put 0 directry to make the origin of the value clear
        new G4PVPlacement(0, G4ThreeVector(-cooler_wall / 2. + stave_thickness / 2., 0.0, 0.0), stave_glue_volume, (boost::format("stave_glue_%d_%d") % inttlayer % itype).str(), stave_volume, false, 0, OverlapCheck());
        new G4PVPlacement(0, G4ThreeVector(-cooler_wall / 2. + stave_thickness / 2., 0.0, 0.0), staveext_glue_volume, (boost::format("stave_glue_ext_%d_%d") % inttlayer % itype).str(), staveext_volume, false, 0, OverlapCheck());
      }
      else
      {
        std::cout << PHWHERE << "invalid laddertype " << laddertype << std::endl;
        gSystem->Exit(1);
      }

      // ----- Step 2 ---------------------------------------------------
      // We place Si-sensor, FPHX, HDI, and stave in the ladder  volume.
      // ================================================================

      // Make the ladder volume first
      // We are still in the loop over inner or outer sensors. This is the ladder volume corresponding to this sensor.
      // But the thickness of the glue for FPHX is used since it's taller than the sensor in x.
      const double ladder_x = stave_x + hdi_kapton_x + hdi_copper_x + fphx_glue_x + fphx_x;
      double ladder_y = hdi_y;

      // Make ladder volume. Drop two corners in positive x as done for stave volume
      G4RotationMatrix *lad_box_rotpos = new G4RotationMatrix();
      lad_box_rotpos->rotateZ(-15. * M_PI / 180.);
      G4ThreeVector ladTranspos(ladder_x / 2., ladder_y / 2., 0.);

      G4RotationMatrix *lad_box_rotneg = new G4RotationMatrix();
      lad_box_rotneg->rotateZ(+15. * M_PI / 180.);
      G4ThreeVector ladTransneg(ladder_x / 2., -ladder_y / 2., 0.);

      G4VSolid *ladder_basebox = new G4Box((boost::format("ladder_basebox_%d_%d") % inttlayer % itype).str(), ladder_x / 2., ladder_y / 2., hdi_z / 2.);
      G4VSolid *ladder_subtbox = new G4Box((boost::format("ladder_subtbox_%d_%d") % inttlayer % itype).str(), stave_x / 1.5, ladder_y / 1.5, hdi_z / 1.);  // has to be longer in z to avoid coincident surface

      G4VSolid *ladder_box1 = new G4SubtractionSolid((boost::format("ladder_box1_%d_%d") % inttlayer % itype).str(), ladder_basebox, ladder_subtbox, lad_box_rotpos, ladTranspos);

      G4VSolid *ladder_box = new G4SubtractionSolid((boost::format("ladder_box_%d_%d") % inttlayer % itype).str(), ladder_box1, ladder_subtbox, lad_box_rotneg, ladTransneg);

      G4LogicalVolume *ladder_volume = new G4LogicalVolume(ladder_box, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")), (boost::format("ladder_%d_%d") % inttlayer % itype).str(), 0, 0, 0);

      G4VSolid *ladderext_basebox = new G4Box((boost::format("ladderext_basebox_%d_%d") % inttlayer % itype).str(), ladder_x / 2., ladder_y / 2., hdiext_z / 2.);
      G4VSolid *ladderext_subtbox = new G4Box((boost::format("ladderext_subtbox_%d_%d") % inttlayer % itype).str(), stave_x / 1.5, ladder_y / 1.5, hdiext_z / 1.);  // has to be longer in z to avoid coincident surface

      G4VSolid *ladderext_box1 = new G4SubtractionSolid((boost::format("ladderext_box1_%d_%d") % inttlayer % itype).str(), ladderext_basebox, ladderext_subtbox, lad_box_rotpos, ladTranspos);

      G4VSolid *ladderext_box = new G4SubtractionSolid((boost::format("ladderext_box_%d_%d") % inttlayer % itype).str(), ladderext_box1, ladderext_subtbox, lad_box_rotneg, ladTransneg);

      G4LogicalVolume *ladderext_volume = new G4LogicalVolume(ladderext_box, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")), (boost::format("ladderext_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      // the rotation matrices are just used by G4VSolid, ownership is not taken over
      delete lad_box_rotpos;
      delete lad_box_rotneg;
      m_DisplayAction->AddVolume(ladder_volume, "Ladder");
      m_DisplayAction->AddVolume(ladderext_volume, "Ladder");

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

      // HDI Kapton
      const double TVhdi_kapton_x = TVstave_x - stave_x / 2. - hdi_kapton_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVhdi_kapton_x, TVstave_y, 0.0), hdi_kapton_volume, (boost::format("hdikapton_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVhdi_kapton_x, TVstave_y, 0.0), hdiext_kapton_volume, (boost::format("hdiextkapton_%d_%s") % inttlayer % itype).str(), ladderext_volume, false, 0, OverlapCheck());

      // HDI copper
      const double TVhdi_copper_x = TVhdi_kapton_x - hdi_kapton_x / 2. - hdi_copper_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVhdi_copper_x, TVstave_y, 0.0), hdi_copper_volume, (boost::format("hdicopper_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVhdi_copper_x, TVstave_y, 0.0), hdiext_copper_volume, (boost::format("hdiextcopper_%d_%s") % inttlayer % itype).str(), ladderext_volume, false, 0, OverlapCheck());

      // Glue for Si-sensor
      const double TVsi_glue_x = TVhdi_copper_x - hdi_copper_x / 2. - si_glue_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVsi_glue_x, TVstave_y, 0.0), si_glue_volume, (boost::format("si_glue_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());

      // Si-sensor
      double TVSi_y = 0.0;
      // sensor is not centered in y in the ladder volume for the Z sensitive ladders
      if (laddertype == PHG4InttDefs::SEGMENTATION_Z)
        TVSi_y = +sensor_offset_y;
      //const double TVSi_x = TVhdi_copper_x - hdi_copper_x / 2. - siactive_x / 2.;
      const double TVSi_x = TVsi_glue_x - si_glue_x / 2. - siactive_x / 2.;
      new G4PVPlacement(0, G4ThreeVector(TVSi_x, TVSi_y, 0.0), siinactive_volume,
                        (boost::format("siinactive_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(TVSi_x, TVSi_y, 0.0), siactive_volume,
                        (boost::format("siactive_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());

      // FPHX glue
      const double TVfphx_glue_x = TVhdi_copper_x - hdi_copper_x / 2. - fphx_glue_x / 2.;
      double TVfphx_glue_y = sifull_y / 2. + gap_sensor_fphx + fphx_y / 2.;
      if (laddertype == PHG4InttDefs::SEGMENTATION_Z)
        TVfphx_glue_y -= sensor_offset_y;

      // laddertype PHG4InttDefs::SEGMENTATION_Z has only one FPHX, laddertype PHG4InttDefs::SEGMENTATION_PHI has two
      if (laddertype == PHG4InttDefs::SEGMENTATION_PHI)
      {
        new G4PVPlacement(0, G4ThreeVector(TVfphx_glue_x, +TVfphx_glue_y, 0.0), fphx_gluecontainer_volume, (boost::format("fphx_gluecontainerp_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      }
      new G4PVPlacement(0, G4ThreeVector(TVfphx_glue_x, -TVfphx_glue_y, 0.0), fphx_gluecontainer_volume, (boost::format("fphx_gluecontainerm_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());

      // FPHX container
      const double TVfphx_x = TVfphx_glue_x - fphx_glue_x / 2. - fphx_x / 2.;
      double TVfphx_y = sifull_y / 2. + gap_sensor_fphx + fphx_y / 2.;
      if (laddertype == PHG4InttDefs::SEGMENTATION_Z)
        TVfphx_y -= sensor_offset_y;

      // laddertype PHG4InttDefs::SEGMENTATION_Z has only one FPHX, laddertype PHG4InttDefs::SEGMENTATION_PHI has two
      if (laddertype == PHG4InttDefs::SEGMENTATION_PHI)
      {
        new G4PVPlacement(0, G4ThreeVector(TVfphx_x, +TVfphx_y, 0.0), fphxcontainer_volume, (boost::format("fphxcontainerp_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());
      }
      new G4PVPlacement(0, G4ThreeVector(TVfphx_x, -TVfphx_y, 0.0), fphxcontainer_volume, (boost::format("fphxcontainerm_%d_%d") % inttlayer % itype).str(), ladder_volume, false, 0, OverlapCheck());

      // ----- Step 3 -----------------------------------------------------------------------------------------------
      // We install the section of ladder for this sensor at all requested phi values and at positive and negative Z
      //=============================================================================================================

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
        if (laddertype == PHG4InttDefs::SEGMENTATION_Z)
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
          m_ActiveVolumeTuple.insert(std::make_pair(pointer_negz, std::make_tuple(inttlayer, itype, icopy, -1)));
          m_ActiveVolumeTuple.insert(std::make_pair(pointer_posz, std::make_tuple(inttlayer, itype, icopy, 1)));
        }

        // The net effect of the above manipulations for the Z sensitive ladders is that the center of the sensor is at dphi * icopy and at the requested radius
        // That us all that the geometry object needs to know, so no changes to that are necessary

        if (itype != 0)
        {
          // We have added the outer sensor above, now we add the HDI extension tab to the end of the outer sensor HDI
          const double posz_ext = (hdi_z_arr[inttlayer][0] + hdi_z) + hdiext_z / 2.;

          new G4PVPlacement(G4Transform3D(ladderrotation, G4ThreeVector(posx, posy, -posz_ext)), ladderext_volume,
                            (boost::format("ladderext_%d_%d_%d_negz") % inttlayer % itype % icopy).str(), trackerenvelope, false, 0, OverlapCheck());
          new G4PVPlacement(G4Transform3D(ladderrotation, G4ThreeVector(posx, posy, +posz_ext)), ladderext_volume,
                            (boost::format("ladderext_%d_%d_%d_posz") % inttlayer % itype % icopy).str(), trackerenvelope, false, 0, OverlapCheck());
        }

        if (Verbosity() > 100)
          std::cout << "   Ladder copy " << icopy << " radius " << radius << " phi " << phi << " itype " << itype << " posz " << m_PosZ[inttlayer][itype]
                    << " fRotate " << fRotate << " posx " << posx << " posy " << posy
                    << std::endl;

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
  const PHParameters *supportparams = m_ParamsContainer->GetParameters(PHG4InttDefs::SUPPORTPARAMS);

  // rails
  G4Tubs *rail_tube = new G4Tubs((boost::format("si_support_rail")).str(),
                                 supportparams->get_double_param("rail_inner_radius") * cm,
                                 supportparams->get_double_param("rail_outer_radius") * cm,
                                 supportparams->get_double_param("rail_length") * cm / 2.0,
                                 0, 2.0 * M_PI);
  G4LogicalVolume *rail_volume = new G4LogicalVolume(rail_tube, GetDetectorMaterial("CFRP_INTT"),
                                                     "rail_volume", 0, 0, 0);
  if (m_IsSupportActive > 0)
  {
    m_PassiveVolumeTuple.insert(std::make_pair(rail_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::SUPPORT_RAIL)));
  }
  m_DisplayAction->AddVolume(rail_volume, "Rail");

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
  // G4Tubs *outer_skin_cfcin_tube = new G4Tubs("si_outer_skin_cfcin",
  //                                            supportparams->get_double_param("outer_skin_cfcin_inner_radius") * cm,
  //                                            supportparams->get_double_param("outer_skin_cfcin_outer_radius") * cm,
  //                                            supportparams->get_double_param("outer_skin_cfcin_length") * cm / 2.,
  //                                            0, 2.0 * M_PI);
  // G4LogicalVolume *outer_skin_cfcin_volume = new G4LogicalVolume(outer_skin_cfcin_tube, GetDetectorMaterial("CFRP_INTT"),
  //                                                                "outer_skin_cfcin_volume", 0, 0, 0);

  // G4Tubs *outer_skin_foam_tube = new G4Tubs("si_outer_skin_foam",
  //                                           supportparams->get_double_param("outer_skin_foam_inner_radius") * cm,
  //                                           supportparams->get_double_param("outer_skin_foam_outer_radius") * cm,
  //                                           supportparams->get_double_param("outer_skin_foam_length") * cm / 2.,
  //                                           0, 2.0 * M_PI);
  // G4LogicalVolume *outer_skin_foam_volume = new G4LogicalVolume(outer_skin_foam_tube, GetDetectorMaterial("ROHACELL_FOAM_110"),
  //                                                               "outer_skin_foam_volume", 0, 0, 0);

  // G4Tubs *outer_skin_cfstd::cout_tube = new G4Tubs("si_outer_skin_cfstd::cout",
  //                                             supportparams->get_double_param("outer_skin_cfstd::cout_inner_radius") * cm,
  //                                             supportparams->get_double_param("outer_skin_cfstd::cout_outer_radius") * cm,
  //                                             supportparams->get_double_param("outer_skin_cfstd::cout_length") * cm / 2.,
  //                                             0, 2.0 * M_PI);
  // G4LogicalVolume *outer_skin_cfstd::cout_volume = new G4LogicalVolume(outer_skin_cfstd::cout_tube, GetDetectorMaterial("CFRP_INTT"),
  //                                                                 "outer_skin_cfstd::cout_volume", 0, 0, 0);

  G4Tubs *outer_skin_tube = new G4Tubs("si_outer_skin",
                                       supportparams->get_double_param("outer_skin_inner_radius") * cm,
                                       supportparams->get_double_param("outer_skin_outer_radius") * cm,
                                       supportparams->get_double_param("outer_skin_length") * cm / 2.,
                                       0, 2.0 * M_PI);
  G4LogicalVolume *outer_skin_volume = new G4LogicalVolume(outer_skin_tube, GetDetectorMaterial("CFRP_INTT"),
                                                           "outer_skin_volume", 0, 0, 0);

  if (m_IsSupportActive > 0)
  {
    // m_PassiveVolumeTuple.insert(std::make_pair(outer_skin_cfcin_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::INTT_OUTER_SKIN)));
    // m_PassiveVolumeTuple.insert(std::make_pair(outer_skin_foam_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::INTT_OUTER_SKIN)));
    // m_PassiveVolumeTuple.insert(std::make_pair(outer_skin_cfstd::cout_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::INTT_OUTER_SKIN)));
    m_PassiveVolumeTuple.insert(std::make_pair(outer_skin_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::INTT_OUTER_SKIN)));
  }
  // m_DisplayAction->AddVolume(outer_skin_cfcin_volume, "Skin");
  // m_DisplayAction->AddVolume(outer_skin_foam_volume, "Skin");
  // m_DisplayAction->AddVolume(outer_skin_cfstd::cout_volume, "Skin");
  // new G4PVPlacement(0, G4ThreeVector(0, 0.0, 0), outer_skin_cfcin_volume,
  //                   "si_support_outer_skin_cfcin", trackerenvelope, false, 0, OverlapCheck());
  // new G4PVPlacement(0, G4ThreeVector(0, 0.0), outer_skin_foam_volume,
  //                   "si_support_outer_skin_foam", trackerenvelope, false, 0, OverlapCheck());
  // new G4PVPlacement(0, G4ThreeVector(0, 0.0), outer_skin_cfstd::cout_volume,
  //                   "si_support_outer_skin_cfstd::cout", trackerenvelope, false, 0, OverlapCheck());
  m_DisplayAction->AddVolume(outer_skin_volume, "Skin");
  new G4PVPlacement(0, G4ThreeVector(0, 0.0, 0), outer_skin_volume,
                    "si_support_outer_skin_cfcin", trackerenvelope, false, 0, OverlapCheck());

  // Inner skin
  G4Tubs *inner_skin_tube = new G4Tubs("si_inner_skin",
                                       supportparams->get_double_param("inner_skin_inner_radius") * cm,
                                       supportparams->get_double_param("inner_skin_outer_radius") * cm,
                                       supportparams->get_double_param("inner_skin_length") * cm / 2.,
                                       0, 2.0 * M_PI);
  G4LogicalVolume *inner_skin_volume = new G4LogicalVolume(inner_skin_tube, GetDetectorMaterial("CFRP_INTT"),
                                                           "inner_skin_volume", 0, 0, 0);
  if (m_IsSupportActive > 0)
  {
    m_PassiveVolumeTuple.insert(std::make_pair(inner_skin_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::INTT_INNER_SKIN)));
  }
  m_DisplayAction->AddVolume(inner_skin_volume, "Skin");

  new G4PVPlacement(0, G4ThreeVector(0, 0.0), inner_skin_volume,
                    "si_support_inner_skin", trackerenvelope, false, 0, OverlapCheck());

  // Service barrel outer  ////////////////////////////////////////////////////////////////////////////////////
  G4Tubs *service_barrel_outer_tube = new G4Tubs("si_service_barrel_outer",
                                                 supportparams->get_double_param("service_barrel_outer_inner_radius") * cm,
                                                 supportparams->get_double_param("service_barrel_outer_outer_radius") * cm,
                                                 supportparams->get_double_param("service_barrel_outer_length") * cm / 2.,
                                                 0, 2.0 * M_PI);
  G4LogicalVolume *service_barrel_outer_volume = new G4LogicalVolume(service_barrel_outer_tube, GetDetectorMaterial("CFRP_INTT"),
                                                                     "service_barrel_outer_volume", 0, 0, 0);
  if (m_IsSupportActive > 0)
  {
    m_PassiveVolumeTuple.insert(std::make_pair(service_barrel_outer_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::SERVICE_BARREL_OUTER)));
  }
  m_DisplayAction->AddVolume(service_barrel_outer_volume, "Skin");

  new G4PVPlacement(0, G4ThreeVector(0, 0.0), service_barrel_outer_volume,
                    "si_support_service_barrel_outer", trackerenvelope, false, 0, OverlapCheck());

  // Support Tube  ////////////////////////////////////////////////////////////////////////////////////
  G4Tubs *support_tube_tube = new G4Tubs("si_support_tube",
                                         supportparams->get_double_param("support_tube_inner_radius") * cm,
                                         supportparams->get_double_param("support_tube_outer_radius") * cm,
                                         supportparams->get_double_param("support_tube_length") * cm / 2.,
                                         0, 2.0 * M_PI);
  G4LogicalVolume *support_tube_volume = new G4LogicalVolume(support_tube_tube, GetDetectorMaterial("CFRP_INTT"),
                                                             "support_tube_volume", 0, 0, 0);
  if (m_IsSupportActive > 0)
  {
    m_PassiveVolumeTuple.insert(std::make_pair(support_tube_volume, std::make_tuple(PHG4InttDefs::SUPPORT_DETID, PHG4InttDefs::SUPPORT_TUBE)));
  }
  m_DisplayAction->AddVolume(support_tube_volume, "Skin");

  new G4PVPlacement(0, G4ThreeVector(0, 0.0), support_tube_volume,
                    "si_support_support_tube", trackerenvelope, false, 0, OverlapCheck());

  // Endcap ring in simulations = Endcap rings + endcap staves
  int inttlayer = (m_LayerBeginEndIteratorPair.first)->second;
  const PHParameters *params1 = m_ParamsContainer->GetParameters(inttlayer);
  const int laddertype = params1->get_int_param("laddertype");
  const PHParameters *params = m_ParamsContainer->GetParameters(laddertype);

  // Carbon PEEK ring ////////////////////////////////////////////////////////////////////////////////////
  G4Tubs *endcap_CP_ring = new G4Tubs("endcap_CP_ring",
                                      supportparams->get_double_param("endcap_CPring_inner_radius") * cm,
                                      supportparams->get_double_param("endcap_CPring_outer_radius") * cm,
                                      supportparams->get_double_param("endcap_CPring_length") * cm / 2.,
                                      0, 2.0 * M_PI);

  G4LogicalVolume *endcap_CP_ring_volume = new G4LogicalVolume(endcap_CP_ring, GetDetectorMaterial("CF30_PEEK70"),
                                                               "endcap_CP_ring_volume", 0, 0, 0);
  m_DisplayAction->AddVolume(endcap_CP_ring_volume, "EndcapCPRing");

  // new Al-PEEK ring from Jan/2021    //////////////////////////////////////////////////////////////////////////
  // outermost part (Al)
  G4Tubs *endcap_AlPEEK_Alring_1 = new G4Tubs("endcap_AlPEEK_Alring_1",
                                              supportparams->get_double_param("endcap_AlPEEK_Cring_1_outer_radius") * cm,
                                              supportparams->get_double_param("endcap_AlPEEK_Alring_1_outer_radius") * cm,
                                              supportparams->get_double_param("endcap_AlPEEK_Alring_length") * cm / 2.,
                                              0, 2.0 * M_PI);

  G4LogicalVolume *endcap_AlPEEK_Alring_1_volume = new G4LogicalVolume(endcap_AlPEEK_Alring_1, GetDetectorMaterial("G4_Al"),
                                                                       "endcap_AlPEEK_Alring_1_volume", 0, 0, 0);
  m_DisplayAction->AddVolume(endcap_AlPEEK_Alring_1_volume, "EndcapAlPEEK_Al1");

  // 2nd outermost part (Carbon PEEK)
  G4Tubs *endcap_AlPEEK_Cring_1 = new G4Tubs("endcap_AlPEEK_Cring_1",
                                             supportparams->get_double_param("endcap_AlPEEK_Alring_2_outer_radius") * cm,
                                             supportparams->get_double_param("endcap_AlPEEK_Cring_1_outer_radius") * cm,
                                             supportparams->get_double_param("endcap_AlPEEK_Cring_length") * cm / 2.,
                                             0, 2.0 * M_PI);

  G4LogicalVolume *endcap_AlPEEK_Cring_1_volume = new G4LogicalVolume(endcap_AlPEEK_Cring_1, GetDetectorMaterial("CF30_PEEK70"),
                                                                      "endcap_AlPEEK_Cring_1_volume", 0, 0, 0);
  m_DisplayAction->AddVolume(endcap_AlPEEK_Cring_1_volume, "EndcapAlPEEK_C1");

  // 3rd outermost part (Al)
  G4Tubs *endcap_AlPEEK_Alring_2 = new G4Tubs("endcap_AlPEEK_Alring_2",
                                              supportparams->get_double_param("endcap_AlPEEK_Cring_2_outer_radius") * cm,
                                              supportparams->get_double_param("endcap_AlPEEK_Alring_2_outer_radius") * cm,
                                              supportparams->get_double_param("endcap_AlPEEK_Alring_length") * cm / 2.,
                                              0, 2.0 * M_PI);

  G4LogicalVolume *endcap_AlPEEK_Alring_2_volume = new G4LogicalVolume(endcap_AlPEEK_Alring_2, GetDetectorMaterial("G4_Al"),
                                                                       "endcap_AlPEEK_Alring_2_volume", 0, 0, 0);
  m_DisplayAction->AddVolume(endcap_AlPEEK_Alring_2_volume, "EndcapAlPEEK_Al2");

  // 4th outermost part (Carbon PEEK)
  G4Tubs *endcap_AlPEEK_Cring_2 = new G4Tubs("endcap_AlPEEK_Cring_2",
                                             supportparams->get_double_param("endcap_AlPEEK_Alring_3_outer_radius") * cm,
                                             supportparams->get_double_param("endcap_AlPEEK_Cring_2_outer_radius") * cm,
                                             supportparams->get_double_param("endcap_AlPEEK_Cring_length") * cm / 2.,
                                             0, 2.0 * M_PI);

  G4LogicalVolume *endcap_AlPEEK_Cring_2_volume = new G4LogicalVolume(endcap_AlPEEK_Cring_2, GetDetectorMaterial("CF30_PEEK70"),
                                                                      "endcap_AlPEEK_Cring_2_volume", 0, 0, 0);
  m_DisplayAction->AddVolume(endcap_AlPEEK_Cring_2_volume, "EndcapAlPEEK_C2");

  // 5th outermost part (innermost) (Al)
  G4Tubs *endcap_AlPEEK_Alring_3 = new G4Tubs("endcap_AlPEEK_Alring_3",
                                              supportparams->get_double_param("endcap_AlPEEK_Alring_3_inner_radius") * cm,
                                              supportparams->get_double_param("endcap_AlPEEK_Alring_3_outer_radius") * cm,
                                              supportparams->get_double_param("endcap_AlPEEK_Alring_length") * cm / 2.,
                                              0, 2.0 * M_PI);

  G4LogicalVolume *endcap_AlPEEK_Alring_3_volume = new G4LogicalVolume(endcap_AlPEEK_Alring_3, GetDetectorMaterial("G4_Al"),
                                                                       "endcap_AlPEEK_Alring_3_volume", 0, 0, 0);
  m_DisplayAction->AddVolume(endcap_AlPEEK_Alring_3_volume, "EndcapAlPEEK_Al3");

  if (m_IsEndcapActive)
  {
    double endcap_outer_edge_z = 0.0;                           // absolute z-coordinate of outer edge (furthest side from the origin) of the endcap, used for bus extender
    if (supportparams->get_int_param("endcap_ring_type") == 0)  // Place Al endcap rings
    {
      // Aluminum ring
      G4Tubs *endcap_Al_ring = new G4Tubs("endcap_Al_ring",
                                          supportparams->get_double_param("endcap_Alring_inner_radius") * cm,
                                          supportparams->get_double_param("endcap_Alring_outer_radius") * cm,
                                          supportparams->get_double_param("endcap_Alring_length") * cm / 2.,
                                          0, 2.0 * M_PI);

      G4LogicalVolume *endcap_Al_ring_volume = new G4LogicalVolume(endcap_Al_ring, GetDetectorMaterial("Al6061T6"),
                                                                   "endcap_Al_ring_volume", 0, 0, 0);

      // Stainlees steal ring
      G4Tubs *endcap_SS_ring = new G4Tubs("endcap_SS_ring",
                                          supportparams->get_double_param("endcap_SSring_inner_radius") * cm,
                                          supportparams->get_double_param("endcap_SSring_outer_radius") * cm,
                                          supportparams->get_double_param("endcap_SSring_length") * cm / 2.,
                                          0, 2.0 * M_PI);

      G4LogicalVolume *endcap_SS_ring_volume = new G4LogicalVolume(endcap_SS_ring, GetDetectorMaterial("SS316"),
                                                                   "endcap_SS_ring_volume", 0, 0, 0);

      // Water Glycol ring
      G4Tubs *endcap_WG_ring = new G4Tubs("endcap_WG_ring",
                                          supportparams->get_double_param("endcap_WGring_inner_radius") * cm,
                                          supportparams->get_double_param("endcap_WGring_outer_radius") * cm,
                                          supportparams->get_double_param("endcap_WGring_length") * cm / 2.,
                                          0, 2.0 * M_PI);

      G4LogicalVolume *endcap_WG_ring_volume = new G4LogicalVolume(endcap_WG_ring, GetDetectorMaterial("WaterGlycol_INTT"),
                                                                   "endcap_WG_ring_volume", 0, 0, 0);

      double endcap_ring_z = supportparams->get_double_param("endcap_ring_z") * cm;
      for (int i = 0; i < 2; i++)  // i=0 : positive z, i=1 negative z
      {
        endcap_ring_z = (i == 0) ? endcap_ring_z : -1.0 * endcap_ring_z;

        double width_WGring_z = supportparams->get_double_param("endcap_WGring_length") * cm;
        double width_SSring_z = supportparams->get_double_param("endcap_SSring_length") * cm;
        double width_Alring_z = supportparams->get_double_param("endcap_Alring_length") * cm;

        for (int j = 0; j < 2; j++)  // j=0 : positive side z, j=1 negative side z
        {
          width_WGring_z = (j == 0) ? width_WGring_z : -1.0 * width_WGring_z;
          width_SSring_z = (j == 0) ? width_SSring_z : -1.0 * width_SSring_z;
          width_Alring_z = (j == 0) ? width_Alring_z : -1.0 * width_Alring_z;

          double cent_WGring_z = endcap_ring_z + width_WGring_z / 2.;
          double cent_SSring_z = endcap_ring_z + width_WGring_z + width_SSring_z / 2.;
          double cent_Alring_z = endcap_ring_z + width_WGring_z + width_SSring_z + width_Alring_z / 2.;
          endcap_outer_edge_z = fabs(endcap_ring_z) + fabs(width_WGring_z + width_SSring_z + width_Alring_z / 2.);  // absolute distance from origin

          new G4PVPlacement(0, G4ThreeVector(0, 0, cent_WGring_z),
                            endcap_WG_ring_volume,
                            (boost::format("endcap_WG_ring_pv_%d_%d") % i % j).str(),
                            trackerenvelope, false, 0, OverlapCheck());

          new G4PVPlacement(0, G4ThreeVector(0, 0, cent_SSring_z),
                            endcap_SS_ring_volume,
                            (boost::format("endcap_SS_ring_pv_%d_%d") % i % j).str(),
                            trackerenvelope, false, 0, OverlapCheck());

          new G4PVPlacement(0, G4ThreeVector(0, 0, cent_Alring_z),
                            endcap_Al_ring_volume,
                            (boost::format("endcap_Al_ring_pv_%d_%d") % i % j).str(),
                            trackerenvelope, false, 0, OverlapCheck());
        }                                                            // end of the loop over positive/negative sides with j
      }                                                              // end of the loop over positive/negative sides with i
    }                                                                // end of endcap_ring_type == 0
    else if (supportparams->get_int_param("endcap_ring_type") == 1)  // Place CP endcap rings
    {
      double endcap_ring_z = supportparams->get_double_param("endcap_CPring_z") * cm;

      for (int i = 0; i < 2; i++)  // i=0 : positive z, i=1 negative z
      {
        endcap_ring_z = (i == 0) ? endcap_ring_z : -1.0 * endcap_ring_z;
        endcap_outer_edge_z = fabs(endcap_ring_z);

        new G4PVPlacement(0, G4ThreeVector(0, 0, endcap_ring_z),
                          endcap_CP_ring_volume,
                          (boost::format("endcap_CP_ring_pv_%d") % i).str(),
                          trackerenvelope, false, 0, OverlapCheck());

      }                                                              // end of the loop over positive/negative sides
    }                                                                // end of endcap_ring_type == 1
    else if (supportparams->get_int_param("endcap_ring_type") == 2)  // the new endcap introduced in Jan/2021
    {
      double si_0_width = params->get_double_param("strip_z_0") * params->get_int_param("nstrips_z_sensor_0") * cm + 2 * params->get_double_param("sensor_edge_z") * cm;  // length of the smaller cells
      double si_1_width = params->get_double_param("strip_z_1") * params->get_int_param("nstrips_z_sensor_1") * cm + 2 * params->get_double_param("sensor_edge_z") * cm;  // length of the larger cells
      double sifull_width = si_0_width + si_1_width;                                                                                                                      // length of the Si module
      double hdi_width = sifull_width + params->get_double_param("hdi_edge_z") * cm;
      double hdiext_width = params->get_double_param("halfladder_inside_z") * cm - sifull_width;
      double hdifull_width = hdi_width + hdiext_width;  // length of a half lader
      double endcap_ring_z = hdifull_width + supportparams->get_double_param("endcap_AlPEEK_Cring_length") / 2.0 * cm;

      for (int i = 0; i < 2; i++)  // i=0 : positive z, i=1 negative z
      {
        endcap_ring_z = (i == 0) ? endcap_ring_z : -1.0 * endcap_ring_z;
        endcap_outer_edge_z = fabs(endcap_ring_z) + supportparams->get_double_param("endcap_AlPEEK_Alring_length") * cm / 2.0;

        new G4PVPlacement(0, G4ThreeVector(0, 0, endcap_ring_z),
                          endcap_AlPEEK_Alring_1_volume,
                          (boost::format("endcap_AlPEEK_Alring_1_pv_%d") % i).str(),
                          trackerenvelope, false, 0, OverlapCheck());

        new G4PVPlacement(0, G4ThreeVector(0, 0, endcap_ring_z),
                          endcap_AlPEEK_Cring_1_volume,
                          (boost::format("endcap_AlPEEK_Cring_1_pv_%d") % i).str(),
                          trackerenvelope, false, 0, OverlapCheck());

        new G4PVPlacement(0, G4ThreeVector(0, 0, endcap_ring_z),
                          endcap_AlPEEK_Alring_2_volume,
                          (boost::format("endcap_AlPEEK_Alring_2_pv_%d") % i).str(),
                          trackerenvelope, false, 0, OverlapCheck());

        new G4PVPlacement(0, G4ThreeVector(0, 0, endcap_ring_z),
                          endcap_AlPEEK_Cring_2_volume,
                          (boost::format("endcap_AlPEEK_Cring_2_pv_%d") % i).str(),
                          trackerenvelope, false, 0, OverlapCheck());

        new G4PVPlacement(0, G4ThreeVector(0, 0, endcap_ring_z),
                          endcap_AlPEEK_Alring_3_volume,
                          (boost::format("endcap_AlPEEK_Alring_3_pv_%d") % i).str(),
                          trackerenvelope, false, 0, OverlapCheck());
      }
    }

    /////////////////////////////////////////////////////////////////////
    // Mimic cylinders for the bus extender (very simplified for the moment)
    double bus_extender_radius_innermost = supportparams->get_double_param("bus_extender_radius") * cm;
    double bus_extender_length = supportparams->get_double_param("bus_extender_length") * cm;
    double bus_extender_copper_x = supportparams->get_double_param("bus_extender_copper_x") * cm;
    double bus_extender_kapton_x = supportparams->get_double_param("bus_extender_kapton_x") * cm;

    // copper layer of the bus extender for the inner barrel
    double inner_radius = bus_extender_radius_innermost;
    G4Tubs *bus_extender_copper_inner = new G4Tubs("bus_extender_coppe_innerr",
                                                   inner_radius, inner_radius + bus_extender_copper_x,
                                                   bus_extender_length / 2.0, 0, 2.0 * M_PI);

    G4LogicalVolume *bus_extender_copper_inner_volume = new G4LogicalVolume(bus_extender_copper_inner, GetDetectorMaterial("G4_Cu"),
                                                                            "bus_extender_copper_inner_volume", 0, 0, 0);
    m_DisplayAction->AddVolume(bus_extender_copper_inner_volume, "BusExtenderCopperInner");

    // kapton layer of the bus extender for the inner barrel
    inner_radius += bus_extender_copper_x;
    G4Tubs *bus_extender_kapton_inner = new G4Tubs("bus_extender_kapton_inner",
                                                   inner_radius, inner_radius + bus_extender_kapton_x,
                                                   bus_extender_length / 2.0, 0, 2.0 * M_PI);

    G4LogicalVolume *bus_extender_kapton_inner_volume = new G4LogicalVolume(bus_extender_kapton_inner, GetDetectorMaterial("G4_KAPTON"),
                                                                            "bus_extender_kapton_inner_volume", 0, 0, 0);
    m_DisplayAction->AddVolume(bus_extender_kapton_inner_volume, "BusExtenderKaptonInner");

    // copper layer of the bus extender for the outer outerbarrel
    inner_radius += bus_extender_kapton_x;
    G4Tubs *bus_extender_copper_outer = new G4Tubs("bus_extender_copper_outer",
                                                   inner_radius, inner_radius + bus_extender_copper_x,
                                                   bus_extender_length / 2.0, 0, 2.0 * M_PI);

    G4LogicalVolume *bus_extender_copper_outer_volume = new G4LogicalVolume(bus_extender_copper_outer, GetDetectorMaterial("G4_Cu"),
                                                                            "bus_extender_copper_outer_volume", 0, 0, 0);
    m_DisplayAction->AddVolume(bus_extender_copper_outer_volume, "BusExtenderCopperOuter");

    // kapton layer of the bus extender for the outer barrel
    inner_radius += bus_extender_copper_x;
    G4Tubs *bus_extender_kapton_outer = new G4Tubs("bus_extender_kapton_outer",
                                                   inner_radius, inner_radius + bus_extender_kapton_x,
                                                   bus_extender_length / 2.0, 0, 2.0 * M_PI);

    G4LogicalVolume *bus_extender_kapton_outer_volume = new G4LogicalVolume(bus_extender_kapton_outer, GetDetectorMaterial("G4_KAPTON"),
                                                                            "bus_extender_kapton_outer_volume", 0, 0, 0);
    m_DisplayAction->AddVolume(bus_extender_kapton_outer_volume, "BusExtenderKaptonOuter");

    double bus_extender_z = endcap_outer_edge_z;
    for (int i = 0; i < 2; i++)  // loop over both positive and negative sides to put the cylinders,  i=0 : positive z, i=1 negative z
    {
      // copper layer of bus extender for the inner barrel
      double cent_bus_extender_z = bus_extender_z + bus_extender_length / 2.0;
      cent_bus_extender_z *= (i == 0) ? 1.0 : -1.0;
      new G4PVPlacement(0, G4ThreeVector(0, 0, cent_bus_extender_z),
                        bus_extender_copper_inner_volume,
                        (boost::format("bus_extender_copper_inner_layer_pv_%d") % i).str(),
                        trackerenvelope, false, 0, OverlapCheck());

      // kapton layer of bus extender for the inner barrel
      new G4PVPlacement(0, G4ThreeVector(0, 0, cent_bus_extender_z),
                        bus_extender_kapton_inner_volume,
                        (boost::format("bus_extender_kapton_inner_layer_pv_%d") % i).str(),
                        trackerenvelope, false, 0, OverlapCheck());

      // copper layer of bus extender for the outer barrel
      new G4PVPlacement(0, G4ThreeVector(0, 0, cent_bus_extender_z),
                        bus_extender_copper_outer_volume,
                        (boost::format("bus_extender_copper_outer_layer_pv_%d") % i).str(),
                        trackerenvelope, false, 0, OverlapCheck());

      // kapton layer of bus extender for the outer barrel
      new G4PVPlacement(0, G4ThreeVector(0, 0, cent_bus_extender_z),
                        bus_extender_kapton_outer_volume,
                        (boost::format("bus_extender_kapton_outer_layer_pv_%d") % i).str(),
                        trackerenvelope, false, 0, OverlapCheck());
    }
  }

  return 0;
}

void PHG4InttDetector::AddGeometryNode()
{
  int active = 0;
  //  std::map<int, int>::const_iterator iter;
  for (auto iter = m_IsActiveMap.begin(); iter != m_IsActiveMap.end(); ++iter)
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
      CylinderGeomIntt *mygeom = new CylinderGeomIntt(
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

std::map<G4VPhysicalVolume *, std::tuple<int, int, int, int>>::const_iterator
PHG4InttDetector::get_ActiveVolumeTuple(G4VPhysicalVolume *physvol) const
{
  auto iter = m_ActiveVolumeTuple.find(physvol);
  if (iter == m_ActiveVolumeTuple.end())
  {
    std::cout << PHWHERE << " Volume " << physvol->GetName() << " not in active volume tuple" << std::endl;
    gSystem->Exit(1);
  }
  return iter;
}

std::map<G4LogicalVolume *, std::tuple<int, int>>::const_iterator
PHG4InttDetector::get_PassiveVolumeTuple(G4LogicalVolume *logvol) const
{
  auto iter = m_PassiveVolumeTuple.find(logvol);
  if (iter == m_PassiveVolumeTuple.end())
  {
    std::cout << PHWHERE << " Volume " << logvol->GetName() << " not in passive volume tuple" << std::endl;
    gSystem->Exit(1);
  }
  return iter;
}
