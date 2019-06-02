#include "PHG4Particlev3.h"
#include "PHG4Particle.h"

#include <Geant4/G4SystemOfUnits.hh>

#include <cmath>
#include <string>

using namespace std;

PHG4Particlev3::PHG4Particlev3()
  : PHG4Particlev2()
  , A(0)
  , Z(0)
  , ioncharge(NAN)
  , excitEnergy(NAN)
{
}

PHG4Particlev3::PHG4Particlev3(const PHG4Particle* in)
  : PHG4Particlev2(in)
  , A(in->get_A())
  , Z(in->get_Z())
  , ioncharge(in->get_IonCharge())
  , excitEnergy(in->get_ExcitEnergy())
{
}

void PHG4Particlev3::set_NumCharge(const int c)
{
  ioncharge = c * eplus;
}

void PHG4Particlev3::identify(std::ostream& os) const
{
  if (fname.size() > 0)
  {
    os << "PHG4Particlev3 name: " << fname << ", ";
  }
  else
  {
    os << "PHG4Particlev3 name: missing, ";
  }

  os << "track id: " << trkid
     << ", vtxid: " << vtxid
     << ", parent id: " << parentid
     << ", primary id: " << primaryid
     << ", pid: " << fpid
     << ", px: " << fpx
     << ", py: " << fpy
     << ", pz: " << fpz
     << ", e: " << fe
     << ", A: " << A
     << ", Z: " << Z
     << ", Eex: " << excitEnergy
     << ", ioncharge: " << ioncharge
     << endl;
  return;
}
