// sPHENIX magnet chimney (valve box extension) and flux-return end doors for Geant4.
// Geometry from BNL 46010010 / 46010200 / 46010300 / 46010400 (P. Kovach, 2015, design review deck), SP00-000-006
// (Ruggiero 2014, overall dimensions) and the BaBar as-built (Ansaldo 620RM07241). Lab frame, cryostat centered at z = 0,
// chimney at the top (phi = 90 deg) on the south side (z < 0). All parameters in cm.
//
// Layout (SP00-000-006): the original chimney box ("camino") sits on the cryostat END FLANGE, outside the barrel, from the
// end face (z = -192.5) + 11.5 cm to 58.9 cm beyond it, straddling the outer radius (34 cm below, 38.5 cm above). The
// 8.00 in OD aluminum extension tube leaves its south face 16.45 cm above the box floor (r ~ 160 cm) and runs 122.5 cm
// south to the duplicate chimney box (46010200: 476 mm triangular footprint, 725 mm tall, 115 lb frame), which carries the
// riser (dia 275 mm port) up to the valve box. Inside the tube: 7.25 in OD heat shield (46010400), two helium lines and the
// copper lead stack (46010300/311/321: 48 x 55 mm with a dia 15 mm coolant hole) with G-10 clamps.
//
// Three parts, selected with the "part" parameter, so that volumes can go to the right mother:
//   0  walls of box 1 above the OHCal envelope inner radius (placed inside the OHCal envelope via SetMotherSubsystem)
//   1  everything in the world: box 1 below the envelope, tube and contents, box 2, riser, valve box
//   2  flux-return end caps (8 in plates from z = 330.81 cm as in the released drawings) with an opening in the south one: door_mode 1 = the real D-slot (205-0600-3036),
//      2 = vertical chase (BaBar) of half-width chase_halfwidth from r = chase_rlo outward, 3 = full south door (no opening)
//      moved outward by door_south_shift, i.e. beyond the second box and the riser
#include "PHG4ChimneyDetector.h"
#include "PHG4ChimneyDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VPhysicalVolume.hh>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

PHG4ChimneyDetector::PHG4ChimneyDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4ChimneyDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
{
}

int PHG4ChimneyDetector::IsInChimney(G4VPhysicalVolume *volume) const
{
  return (m_LogicalVolSet.find(volume->GetLogicalVolume()) != m_LogicalVolSet.end()) ? 1 : 0;
}

void PHG4ChimneyDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  G4Material *aluminum = GetDetectorMaterial("G4_Al");
  G4Material *copper = GetDetectorMaterial("G4_Cu");
  G4Material *g10 = GetDetectorMaterial("G10");
  G4Material *steel = GetDetectorMaterial("Steel_1006");
  G4Material *ssteel = GetDetectorMaterial("G4_STAINLESS-STEEL");

  const int part = m_Params->get_int_param("part");
  const double r_env = m_Params->get_double_param("r_env") * cm;  // OHCal envelope inner radius
  const double wall = m_Params->get_double_param("box_wall") * cm;
  // box 1 (camino on the end flange)
  const double b1_zmin = m_Params->get_double_param("box1_zmin") * cm, b1_zmax = m_Params->get_double_param("box1_zmax") * cm;
  const double b_rlo = m_Params->get_double_param("box_rlo") * cm, b_rhi = b_rlo + m_Params->get_double_param("box_dr") * cm;
  const double b1_hx = 0.5 * m_Params->get_double_param("box1_dx") * cm;
  // box 2 (duplicate chimney with the riser)
  const double b2_zmin = m_Params->get_double_param("box2_zmin") * cm, b2_zmax = m_Params->get_double_param("box2_zmax") * cm;
  const double b2_hx = 0.5 * m_Params->get_double_param("box2_dx") * cm;
  // tube and contents
  const double tube_od = m_Params->get_double_param("tube_od") * cm, tube_wall = m_Params->get_double_param("tube_wall") * cm;
  const double tube_r = m_Params->get_double_param("tube_r_axis") * cm;
  const double shield_od = m_Params->get_double_param("shield_od") * cm, shield_wall = m_Params->get_double_param("shield_wall") * cm;
  const double he_od = m_Params->get_double_param("heline_od") * cm, he_wall = m_Params->get_double_param("heline_wall") * cm;
  const double cu_dx = m_Params->get_double_param("cu_dx") * cm, cu_dy = m_Params->get_double_param("cu_dy") * cm;
  const double hole_d = m_Params->get_double_param("hole_d") * cm, leg_len = m_Params->get_double_param("leg_len") * cm;
  const double g10_w = m_Params->get_double_param("g10_w") * cm, g10_t = m_Params->get_double_param("g10_t") * cm;
  // riser and valve box
  const double riser_od = m_Params->get_double_param("riser_od") * cm, riser_wall = m_Params->get_double_param("riser_wall") * cm;
  const double vb_rlo = m_Params->get_double_param("vbox_rlo") * cm, vb_rhi = m_Params->get_double_param("vbox_rhi") * cm;
  const double vb_od = m_Params->get_double_param("vbox_od") * cm, vb_wall = m_Params->get_double_param("vbox_wall") * cm;
  // doors
  const int door_mode = m_Params->get_int_param("door_mode");
  const double door_z1 = m_Params->get_double_param("door_z1") * cm, door_z2 = m_Params->get_double_param("door_z2") * cm;
  const double door_rin = m_Params->get_double_param("door_rin") * cm, door_rout = m_Params->get_double_param("door_rout") * cm;
  const double chase_hw = m_Params->get_double_param("chase_halfwidth") * cm, chase_rlo = m_Params->get_double_param("chase_rlo") * cm;
  const double door_hole_d = m_Params->get_double_param("door_hole_d") * cm;
  const double slot_w = m_Params->get_double_param("door_slot_w") * cm, slot_h = m_Params->get_double_param("door_slot_h") * cm;
  const double slot_r = m_Params->get_double_param("door_slot_r") * cm;   // slot center radius
  const double xoff = m_Params->get_double_param("ext_xoff") * cm;        // transverse offset of the extension tube, box 2, riser and valve box
  const double door_south_shift = m_Params->get_double_param("door_south_shift") * cm;  // mode 3: full south door moved outward by this much

  // The OHCal envelope is placed with rot_y = 180 deg: envelope-local x = -x_lab, z = -z_lab.
  const bool flip = (std::string(logicWorld->GetName()) == "OHCal_envelope");
  auto P = [flip](double x, double y, double z) { return flip ? G4ThreeVector(-x, y, -z) : G4ThreeVector(x, y, z); };
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ChimneyDetector part " << part << ": mother volume " << logicWorld->GetName() << (flip ? " (rot_y 180: flipping x, z)" : "")
              << "; box 1 z " << b1_zmin / cm << ".." << b1_zmax / cm << " r " << b_rlo / cm << ".." << b_rhi / cm
              << "; tube axis r " << tube_r / cm << ", OD " << tube_od / cm << "; box 2 z " << b2_zmin / cm << ".." << b2_zmax / cm
              << "; door mode " << door_mode << std::endl;
  }
  int copy = 0;
  // Al plate from full lab extents, placed in the current mother
  auto plate = [&](const std::string &name, double x0, double x1, double y0, double y1, double z0, double z1) {
    if (x1 <= x0 || y1 <= y0 || z1 <= z0) return;
    G4Box *b = new G4Box(name + "_solid", 0.5 * (x1 - x0), 0.5 * (y1 - y0), 0.5 * (z1 - z0));
    G4LogicalVolume *lv = new G4LogicalVolume(b, aluminum, name);
    m_DisplayAction->AddVolume(lv, "CHIMNEYBOX");
    m_LogicalVolSet.insert(lv);
    new G4PVPlacement(nullptr, P(0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1)), lv, name, logicWorld, false, copy++, OverlapCheck());
  };
  // six walls of a box, restricted to y in [ylo, yhi] (the split at the envelope radius)
  auto boxwalls = [&](const std::string &tag, double hx, double zmin, double zmax, double ylo, double yhi, double xc = 0) {
    auto seg = [&](const std::string &n, double x0, double x1, double y0, double y1, double z0, double z1) {
      plate(tag + n, xc + x0, xc + x1, std::max(y0, ylo), std::min(y1, yhi), z0, z1);
    };
    seg("Floor", -hx, hx, b_rlo, b_rlo + wall, zmin, zmax);
    seg("Top", -hx, hx, b_rhi - wall, b_rhi, zmin, zmax);
    seg("SideL", -hx, -hx + wall, b_rlo + wall, b_rhi - wall, zmin, zmax);
    seg("SideR", hx - wall, hx, b_rlo + wall, b_rhi - wall, zmin, zmax);
    seg("EndS", -hx + wall, hx - wall, b_rlo + wall, b_rhi - wall, zmin, zmin + wall);
    seg("EndN", -hx + wall, hx - wall, b_rlo + wall, b_rhi - wall, zmax - wall, zmax);
  };
  // the envelope boundary is a cylinder: keep world pieces below y_lo_env and envelope pieces above r_env
  const double y_lo_env = std::sqrt(r_env * r_env - std::max(b1_hx, b2_hx) * std::max(b1_hx, b2_hx)) - 2 * mm;
  const double y_hi_env = r_env + 2 * mm;

  if (part == 0)
  {
    boxwalls("ChimneyBox1U", b1_hx, b1_zmin, b1_zmax, y_hi_env, 1e9);
    return;
  }

  if (part == 1)
  {
    boxwalls("ChimneyBox1L", b1_hx, b1_zmin, b1_zmax, -1e9, y_lo_env);
    boxwalls("ChimneyBox2", b2_hx, b2_zmin, b2_zmax, -1e9, 1e9, xoff);

    // tube between the boxes, with the heat shield inside
    const double tz0 = b2_zmax, tz1 = b1_zmin, tlen = tz1 - tz0, tzc = 0.5 * (tz0 + tz1);
    G4Tubs *tube = new G4Tubs("ChimneyTube", 0.5 * tube_od - tube_wall, 0.5 * tube_od, 0.5 * tlen, 0, 360 * deg);
    G4LogicalVolume *logTube = new G4LogicalVolume(tube, aluminum, "CHIMNEYTUBE");
    m_DisplayAction->AddVolume(logTube, "CHIMNEYTUBE");
    m_LogicalVolSet.insert(logTube);
    new G4PVPlacement(nullptr, P(xoff, tube_r, tzc), logTube, "CHIMNEYTUBE", logicWorld, false, 0, OverlapCheck());
    G4Tubs *shield = new G4Tubs("ChimneyShield", 0.5 * shield_od - shield_wall, 0.5 * shield_od, 0.5 * tlen - 1 * cm, 0, 360 * deg);
    G4LogicalVolume *logShield = new G4LogicalVolume(shield, aluminum, "CHIMNEYSHIELD");
    m_DisplayAction->AddVolume(logShield, "CHIMNEYSHIELD");
    m_LogicalVolSet.insert(logShield);
    new G4PVPlacement(nullptr, P(xoff, tube_r, tzc), logShield, "CHIMNEYSHIELD", logicWorld, false, 0, OverlapCheck());
    // two helium lines (stainless), upper half of the shield
    G4Tubs *he = new G4Tubs("ChimneyHeLine", 0.5 * he_od - he_wall, 0.5 * he_od, 0.5 * tlen - 2 * cm, 0, 360 * deg);
    G4LogicalVolume *logHe = new G4LogicalVolume(he, ssteel, "CHIMNEYHELINE");
    m_DisplayAction->AddVolume(logHe, "CHIMNEYHELINE");
    m_LogicalVolSet.insert(logHe);
    for (int s = -1; s <= 1; s += 2)
    {
      new G4PVPlacement(nullptr, P(xoff + s * 5.0 * cm, tube_r + 4.5 * cm, tzc), logHe, "CHIMNEYHELINE", logicWorld, false, (s > 0), OverlapCheck());
    }
    // copper lead stack along the tube (inside the shield), with the coolant hole, and the G-10 clamp bars
    const double cu_len = tlen - 2 * cm;
    G4Box *cuBox = new G4Box("ChimneyCuBox", 0.5 * cu_dx, 0.5 * cu_dy, 0.5 * cu_len);
    G4Tubs *hole = new G4Tubs("ChimneyCuHole", 0, 0.5 * hole_d, 0.5 * cu_len + 1 * mm, 0, 360 * deg);
    G4SubtractionSolid *cuBar = new G4SubtractionSolid("ChimneyCuBar", cuBox, hole);
    G4LogicalVolume *logCu = new G4LogicalVolume(cuBar, copper, "CHIMNEYBUS");
    m_DisplayAction->AddVolume(logCu, "CHIMNEYBUS");
    m_LogicalVolSet.insert(logCu);
    new G4PVPlacement(nullptr, P(xoff, tube_r, tzc), logCu, "CHIMNEYBUS", logicWorld, false, 0, OverlapCheck());
    G4Box *g10Box = new G4Box("ChimneyG10", 0.5 * g10_t, 0.5 * g10_w, 0.5 * cu_len);
    G4LogicalVolume *logG10 = new G4LogicalVolume(g10Box, g10, "CHIMNEYG10");
    m_DisplayAction->AddVolume(logG10, "CHIMNEYG10");
    m_LogicalVolSet.insert(logG10);
    for (int side = -1; side <= 1; side += 2)
    {
      new G4PVPlacement(nullptr, P(xoff + side * (0.5 * cu_dx + 0.5 * g10_t + 1 * mm), tube_r, tzc), logG10, "CHIMNEYG10", logicWorld, false, (side > 0), OverlapCheck());
    }
    // vertical copper legs: down to the floor of box 1 (toward the coil), up in box 2 (toward the riser)
    {
      const double y_top = tube_r + 0.5 * cu_dy, y_floor = b_rlo + wall + 1 * mm;
      const double len1 = std::min(leg_len, y_top - y_floor);
      G4Box *leg1 = new G4Box("ChimneyCuLeg1", 0.5 * cu_dx, 0.5 * len1, 0.5 * cu_dy);
      G4LogicalVolume *logLeg1 = new G4LogicalVolume(leg1, copper, "CHIMNEYBUSLEG1");
      m_DisplayAction->AddVolume(logLeg1, "CHIMNEYBUSLEG");
      m_LogicalVolSet.insert(logLeg1);
      new G4PVPlacement(nullptr, P(xoff, y_top - 0.5 * len1, b1_zmin + wall + 0.5 * cu_dy + 2 * mm), logLeg1, "CHIMNEYBUSLEG1", logicWorld, false, 0, OverlapCheck());
      const double y_bot = tube_r - 0.5 * cu_dy;
      const double len2 = std::min(leg_len, b_rhi - wall - 1 * mm - y_bot);
      G4Box *leg2 = new G4Box("ChimneyCuLeg2", 0.5 * cu_dx, 0.5 * len2, 0.5 * cu_dy);
      G4LogicalVolume *logLeg2 = new G4LogicalVolume(leg2, copper, "CHIMNEYBUSLEG2");
      m_DisplayAction->AddVolume(logLeg2, "CHIMNEYBUSLEG");
      m_LogicalVolSet.insert(logLeg2);
      new G4PVPlacement(nullptr, P(xoff, y_bot + 0.5 * len2, b2_zmax - wall - 0.5 * cu_dy - 2 * mm), logLeg2, "CHIMNEYBUSLEG2", logicWorld, false, 0, OverlapCheck());
    }
    // riser from the top of box 2 to the valve box, and the valve box (vertical cylinders, axis along y)
    G4RotationMatrix *rotX = new G4RotationMatrix();
    rotX->rotateX(90 * deg);
    const double zc2 = 0.5 * (b2_zmin + b2_zmax);
    const double rlen = vb_rlo - b_rhi;
    if (rlen > 0)
    {
      G4Tubs *riser = new G4Tubs("ChimneyRiser", 0.5 * riser_od - riser_wall, 0.5 * riser_od, 0.5 * rlen, 0, 360 * deg);
      G4LogicalVolume *logRiser = new G4LogicalVolume(riser, aluminum, "CHIMNEYRISER");
      m_DisplayAction->AddVolume(logRiser, "CHIMNEYRISER");
      m_LogicalVolSet.insert(logRiser);
      new G4PVPlacement(rotX, P(xoff, b_rhi + 0.5 * rlen, zc2), logRiser, "CHIMNEYRISER", logicWorld, false, 0, OverlapCheck());
    }
    if (vb_rhi > vb_rlo)
    {
      G4Tubs *vbs = new G4Tubs("ChimneyValveBoxShell", 0.5 * vb_od - vb_wall, 0.5 * vb_od, 0.5 * (vb_rhi - vb_rlo), 0, 360 * deg);
      G4LogicalVolume *logVb = new G4LogicalVolume(vbs, ssteel, "CHIMNEYVALVEBOX");
      m_DisplayAction->AddVolume(logVb, "CHIMNEYVALVEBOX");
      m_LogicalVolSet.insert(logVb);
      new G4PVPlacement(rotX, P(xoff, 0.5 * (vb_rlo + vb_rhi), zc2), logVb, "CHIMNEYVALVEBOX", logicWorld, false, 0, OverlapCheck());
      G4Tubs *vbe = new G4Tubs("ChimneyValveBoxEnd", 0, 0.5 * vb_od - vb_wall, 0.5 * vb_wall, 0, 360 * deg);
      G4LogicalVolume *logVbe = new G4LogicalVolume(vbe, ssteel, "CHIMNEYVALVEBOXEND");
      m_DisplayAction->AddVolume(logVbe, "CHIMNEYVALVEBOX");
      m_LogicalVolSet.insert(logVbe);
      new G4PVPlacement(rotX, P(xoff, vb_rlo + 0.5 * vb_wall, zc2), logVbe, "CHIMNEYVALVEBOXEND", logicWorld, false, 0, OverlapCheck());
      new G4PVPlacement(rotX, P(xoff, vb_rhi - 0.5 * vb_wall, zc2), logVbe, "CHIMNEYVALVEBOXEND", logicWorld, false, 1, OverlapCheck());
    }
    return;
  }

  if (part == 2 && door_mode > 0)
  {
    const double dlen = door_z2 - door_z1, dzc = 0.5 * (door_z1 + door_z2);
    G4Tubs *disk = new G4Tubs("ChimneyDoorDisk", door_rin, door_rout, 0.5 * dlen, 0, 360 * deg);
    // north door: full disk
    G4LogicalVolume *logN = new G4LogicalVolume(disk, steel, "CHIMNEYDOORN");
    m_DisplayAction->AddVolume(logN, "CHIMNEYDOOR");
    m_LogicalVolSet.insert(logN);
    new G4PVPlacement(nullptr, P(0, 0, dzc), logN, "CHIMNEYDOORN", logicWorld, false, 0, OverlapCheck());
    // south door with the opening for the chimney
    G4VSolid *south = disk;
    double zshift = 0;
    if (door_mode == 3)
    {
      zshift = door_south_shift;
    }
    else if (door_mode == 1)
    {
      // 205-0600-3036 (southwest half): D-slot at the center edge, slot_h tall, slot_w deep (straight part plus a semicircular end
      // of radius slot_h/2), centered at r = slot_r, on the +x side; the extension tube passes through it
      G4Box *sb = new G4Box("ChimneyDoorSlotBox", 0.5 * (slot_w - 0.5 * slot_h), 0.5 * slot_h, 0.5 * dlen + 1 * cm);
      G4Tubs *se = new G4Tubs("ChimneyDoorSlotEnd", 0, 0.5 * slot_h, 0.5 * dlen + 1 * cm, 0, 360 * deg);
      G4SubtractionSolid *s1 = new G4SubtractionSolid("ChimneyDoorSouth1", disk, sb, nullptr, G4ThreeVector(0.5 * (slot_w - 0.5 * slot_h), slot_r, 0));
      south = new G4SubtractionSolid("ChimneyDoorSouth", s1, se, nullptr, G4ThreeVector(slot_w - 0.5 * slot_h, slot_r, 0));
      (void) door_hole_d;
    }
    else
    {
      const double ytop = door_rout + 5 * cm;
      G4Box *slot = new G4Box("ChimneyDoorChase", chase_hw, 0.5 * (ytop - chase_rlo), 0.5 * dlen + 1 * cm);
      south = new G4SubtractionSolid("ChimneyDoorSouth", disk, slot, nullptr, G4ThreeVector(0, 0.5 * (ytop + chase_rlo), 0));
    }
    G4LogicalVolume *logS = new G4LogicalVolume(south, steel, "CHIMNEYDOORS");
    m_DisplayAction->AddVolume(logS, "CHIMNEYDOOR");
    m_LogicalVolSet.insert(logS);
    new G4PVPlacement(nullptr, P(0, 0, -dzc - zshift), logS, "CHIMNEYDOORS", logicWorld, false, 0, OverlapCheck());
    return;
  }
  return;
}
