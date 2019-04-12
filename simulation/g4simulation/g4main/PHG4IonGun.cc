#include "PHG4IonGun.h"

#include "PHG4Particlev3.h"
#include "PHG4InEvent.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

#include <Geant4/G4Event.hh>
#include <Geant4/G4IonTable.hh>
#include <Geant4/G4PrimaryParticle.hh>
#include <Geant4/G4PrimaryVertex.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>

using namespace std;

PHG4IonGun::PHG4IonGun(const string &name)
  : PHG4ParticleGeneratorBase(name)
  , ion(new PHG4Particlev3())
  , A(0)
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
  mom[0] = px;
  mom[1] = py;
  mom[2] = pz;
}

int PHG4IonGun::process_event(PHCompositeNode *topNode)
{
// get pdg code for ions
  UpdateParticle();
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");
  ReuseExistingVertex(topNode); // checks if we should reuse existing vertex
  int vtxindex = ineve->AddVtx(vtx_x,vtx_y,vtx_z,t0);
//   G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
   //G4PrimaryParticle* g4part = new G4PrimaryParticle(ion);
//   cout << "name: " << ion->GetParticleName() << ", pdgcode: " << ion->GetPDGEncoding() << endl;
      PHG4Particle *particle = new PHG4Particlev3(ion);
      SetParticleId(particle,ineve);
      ineve->AddParticle(vtxindex, particle);
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4IonGun::UpdateParticle()
{
  ion->set_A(A);
  ion->set_Z(Z);
  ion->set_NumCharge(ioncharge);
  ion->set_ExcitEnergy(excitEnergy);
  ion->set_px(mom[0]);
  ion->set_py(mom[1]);
  ion->set_pz(mom[2]);
  G4ParticleDefinition* iondef = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
  ion->set_name(iondef->GetParticleName());
  ion->set_pid(iondef->GetPDGEncoding());
  return;
}

void PHG4IonGun::Print(const string &what) const
{
  cout << "PHG4IonGun, using ions of" << endl;
  cout << "A: " << A << ", Z: " << Z << ", charge: " << ioncharge 
       << ", excitation Energy: " << excitEnergy << endl;
  cout << "px: " << mom[0]/GeV << " GeV, py: " << mom[1]/GeV << " GeV, pz: " << mom[2]/GeV << " GeV" << endl;
}
