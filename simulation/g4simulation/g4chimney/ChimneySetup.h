// Register the chimney model (and optionally the flux-return doors) with PHG4Reco. Call after G4Setup().
//   door_mode: 0 no doors, 1 doors with a round hole for the tube, 2 doors with a BaBar-type chase in the south door, 3 full south door moved beyond the chimney
//   with_chimney: 0 none, 1 box 1 + tube + contents + box 2 + riser + valve box
#ifndef CHIMNEYSETUP_H
#define CHIMNEYSETUP_H
#include "PHG4ChimneySubsystem.h"
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4Subsystem.h>
#include <iostream>
inline void RegisterChimney(PHG4Reco *g4, const int with_chimney, const int door_mode, const bool overlapcheck = false, const int verbosity = 0)
{
  if (with_chimney)
  {
    PHG4Subsystem *ohcal = g4->getSubsystem("HCALOUT_0");
    for (int part = 0; part <= 1; ++part)
    {
      PHG4ChimneySubsystem *c = new PHG4ChimneySubsystem("CHIMNEY", part);
      c->set_int_param("part", part);
      c->OverlapCheck(overlapcheck);
      c->Verbosity(verbosity);
      if (part == 0)
      {
        if (!ohcal) { std::cout << "RegisterChimney: no HCALOUT_0, skipping the envelope part" << std::endl; delete c; continue; }
        c->SetMotherSubsystem(ohcal);
      }
      g4->registerSubsystem(c);
    }
  }
  if (door_mode > 0)
  {
    PHG4ChimneySubsystem *d = new PHG4ChimneySubsystem("CHIMNEY", 2);
    d->set_int_param("part", 2);
    d->set_int_param("door_mode", door_mode);
    d->OverlapCheck(overlapcheck);
    d->Verbosity(verbosity);
    g4->registerSubsystem(d);
  }
  std::cout << "RegisterChimney: chimney " << with_chimney << ", door mode " << door_mode << std::endl;
}
#endif
