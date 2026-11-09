#include "PHG4ChimneySubsystem.h"

#include "PHG4ChimneyDetector.h"
#include "PHG4ChimneyDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>

#include <phool/PHCompositeNode.h>

#include <iostream>

PHG4ChimneySubsystem::PHG4ChimneySubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

PHG4ChimneySubsystem::~PHG4ChimneySubsystem()
{
  delete m_DisplayAction;
}

int PHG4ChimneySubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  m_DisplayAction = new PHG4ChimneyDisplayAction(Name());
  m_Detector = new PHG4ChimneyDetector(this, topNode, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());
  return 0;
}

void PHG4ChimneySubsystem::Print(const std::string &what) const
{
  std::cout << "Chimney parameters:" << std::endl;
  GetParams()->Print();
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
}

PHG4Detector *PHG4ChimneySubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4ChimneySubsystem::SetDefaultParameters()
{
  // all lengths in cm. Sources: BNL 46010010/200/300/311/321/340/400 (Kovach 2015), SP00-000-006 (Ruggiero 2014), Ansaldo 620RM07241.
  set_default_int_param("part", 1);              // 0: box 1 above the OHCal envelope radius (mother = OHCal envelope), 1: world parts, 2: doors
  set_default_double_param("r_env", 177.42);     // OHCal envelope inner radius
  set_default_double_param("box_wall", 1.5);     // effective Al wall (80 kg camino, 52 kg frame + covers)
  set_default_double_param("box1_zmin", -251.4); // SP00-000-006: 23.18 in beyond the cryostat end (-192.5)
  set_default_double_param("box1_zmax", -203.0); // CAD: box 1 spans z -251..-198, but the PR 3555 BaBar "CONNECTOR" stub occupies z -202.5..-192.5 at this r: stop at -203
  set_default_double_param("box_rlo", 143.0);    // box bottom 34 cm below the cryostat outer radius (SP00-000-006 side view)
  set_default_double_param("box_dr", 72.5);      // 725 mm tall (46010200) -> top at 215.5 cm, below the OHCal steel (229)
  set_default_double_param("box1_dx", 51.1);     // camino width
  set_default_double_param("box2_zmin", -421.5); // box 2: 48.2 in of tube from box 1, 476 mm deep (46010200)
  set_default_double_param("box2_zmax", -373.9);
  set_default_double_param("box2_dx", 51.0);     // CAD: box 2 about 51 cm wide, centered at ext_xoff
  set_default_double_param("tube_od", 20.32);    // 8.00 in OD aluminum (46010010 note 3, 600ASA-800 flange)
  set_default_double_param("tube_wall", 0.32);   // 46010020: 8.00 in OD spool, 41.76 in, 19 lb with two flanges -> about 1/8 in wall
  set_default_double_param("tube_r_axis", 163.8);// CAD: busbar 46010321 at x 63.9-66.6 in, heat shield halves 61.5-67.4 in -> axis 64.5 in
  set_default_double_param("shield_od", 18.4);   // 7.25 in OD heat shield (46010400)
  set_default_double_param("shield_wall", 0.15); // ASSUMED
  set_default_double_param("heline_od", 1.6);    // ASSUMED 5/8 in helium lines
  set_default_double_param("heline_wall", 0.15);
  set_default_double_param("cu_dx", 4.8);
  set_default_double_param("cu_dy", 5.5);
  set_default_double_param("hole_d", 1.5);
  set_default_double_param("leg_len", 35.6);
  set_default_double_param("g10_w", 7.62);
  set_default_double_param("g10_t", 0.635);
  set_default_double_param("riser_od", 27.3);    // dia 275 mm port (46010200); 10.75 in tube ASSUMED
  set_default_double_param("riser_wall", 0.64);
  set_default_double_param("vbox_rlo", 334.0);   // valve box bottom (SP00-000-006: 112.12 in from the cryostat OD to the top with fittings)
  set_default_double_param("vbox_rhi", 434.0);   // ASSUMED 1 m tall body
  set_default_double_param("vbox_od", 68.7);     // 27.04 in body
  set_default_double_param("vbox_wall", 1.0);    // ASSUMED stainless
  set_default_int_param("door_mode", 0);         // 0 none, 1 real D-slot for the tube (released drawing), 2 BaBar-type chase, 3 full south door moved outward
  set_default_double_param("door_z1", 330.81);   // inner face: doors start at z = +-3308.1 mm (Dan), same as G4_PlugDoor.C
  set_default_double_param("door_z2", 351.13);   // 8.00 in thick plate (sPHENIX end cap procurement 2020: 212 in OD, 12 in ID, 8 in thick, 22.5 t per half)
  set_default_double_param("door_rin", 15.24);
  set_default_double_param("door_rout", 269.2);  // 212 in OD
  set_default_double_param("door_south_shift", 115.0);  // mode 3: south door inner face at -445.8, beyond box 2 (-421.5) and the valve box (-432)
  set_default_double_param("chase_halfwidth", 30.0);  // ASSUMED: clears the 51 cm box
  set_default_double_param("chase_rlo", 140.0);       // ASSUMED: down to the box floor
  set_default_double_param("door_hole_d", 24.0);      // unused (the slot is known)
  set_default_double_param("door_slot_w", 33.0);      // 205-0600-3036: D-slot 13 in deep from the center edge into the southwest half
  set_default_double_param("door_slot_h", 33.0);      // 13 in tall (r 149-182 cm), rounded end R 6.5 in
  set_default_double_param("door_slot_r", 165.1);     // slot center 65 in above the beam line
  set_default_double_param("ext_xoff", 15.0);         // CAD: tube, heat shield and box 2 centered about 6 in off the chimney plane, toward the slot
}
