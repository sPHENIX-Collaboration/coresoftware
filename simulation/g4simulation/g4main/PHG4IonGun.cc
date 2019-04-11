#include "PHG4IonGun.h"

#include <Geant4/G4Event.hh>
#include <Geant4/G4IonTable.hh>
#include <Geant4/G4PrimaryParticle.hh>
#include <Geant4/G4PrimaryVertex.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>

using namespace std;

PHG4IonGun::PHG4IonGun()
  : A(0)
  , Z(0)
  , ioncharge(0)
  , excitEnergy(0)
{
  fill(begin(mom),end(mom),NAN);
}

void PHG4IonGun::SetCharge(const int c)
{
  ioncharge = c * eplus;
}
void PHG4IonGun::SetMom(const double px, const double py, const double pz)
{
  mom[0] = px * GeV;
  mom[1] = py * GeV;
  mom[2] = pz * GeV;
}

void PHG4IonGun::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
  G4ThreeVector position(0 * cm, 0 * cm, 0 * cm);
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, 0 * s);
  G4PrimaryParticle* g4part = new G4PrimaryParticle(ion);
  g4part->SetCharge(ioncharge);
  g4part->SetMomentum(mom[0], mom[1], mom[2]);
  vertex->SetPrimary(g4part);
  anEvent->AddPrimaryVertex(vertex);
  return;
}

void PHG4IonGun::Print(const string &what) const
{
  cout << "PHG4IonGun, using ions of" << endl;
  cout << "A: " << A << ", Z: " << Z << ", charge: " << ioncharge 
       << ", excitation Energy: " << excitEnergy << endl;
  cout << "px: " << mom[0]/GeV << " GeV, py: " << mom[1]/GeV << " GeV, pz: " << mom[2]/GeV << " GeV" << endl;
}
