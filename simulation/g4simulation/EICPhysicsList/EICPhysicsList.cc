// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        EICPhysicsList.hh                                            //
//  Description: Geant4 physics list for Electron Ion Collider detectors      //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        21 June 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "EICPhysicsList.hh"

#include "ProtonPhysics.hh"
#include "NeutronPhysics.hh"
#include "PionPhysics.hh"
#include "KaonPhysics.hh"
#include "HyperonPhysics.hh"
#include "AntiBaryonPhysics.hh"
#include "IonPhysics.hh"
#include "GammaLeptoNuclearPhysics.hh"

#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4UnitsTable.hh>

#include <Geant4/G4EmStandardPhysics.hh>
#include <Geant4/G4EmExtraPhysics.hh>
#include <Geant4/G4OpticalPhysics.hh>
#include <Geant4/G4EmParameters.hh>
#include <Geant4/G4DecayPhysics.hh>
#include <Geant4/G4NuclideTable.hh>
// #include <Geant4/G4RadioactiveDecayPhysics.hh>

// particles

#include <Geant4/G4BosonConstructor.hh>
#include <Geant4/G4LeptonConstructor.hh>
#include <Geant4/G4MesonConstructor.hh>
#include <Geant4/G4BaryonConstructor.hh>
#include <Geant4/G4IonConstructor.hh>
#include <Geant4/G4ShortLivedConstructor.hh>


EICPhysicsList::EICPhysicsList()
 :G4VModularPhysicsList()
{
  G4int verb = 1;
  SetVerboseLevel(verb);
  
  //add new units for radioActive decays
  //
  new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
  // 
  const G4double minute = 60*second;
  const G4double hour   = 60*minute;
  const G4double day    = 24*hour;
  const G4double year   = 365*day;
  new G4UnitDefinition("minute", "min", "Time", minute);
  new G4UnitDefinition("hour",   "h",   "Time", hour);
  new G4UnitDefinition("day",    "d",   "Time", day);
  new G4UnitDefinition("year",   "y",   "Time", year);

  // Mandatory for G4NuclideTable
  // Half-life threshold must be set small or many short-lived isomers 
  // will not be assigned life times (default to 0) 
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);
          
  // EM physics
  RegisterPhysics(new G4EmStandardPhysics());
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetAugerCascade(true);
#if G4VERSION_NUMBER >= 1004
   param->SetStepFunction(1., 1*CLHEP::mm);
   param->SetStepFunctionMuHad(1., 1*CLHEP::mm);
#endif

  // Decay
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
//  RegisterPhysics(new G4RadioactiveDecayPhysics());
            
  // Hadronic physics
  RegisterPhysics(new ProtonPhysics() );
  RegisterPhysics(new NeutronPhysics() );
  RegisterPhysics(new PionPhysics() );
  RegisterPhysics(new KaonPhysics() );
  RegisterPhysics(new HyperonPhysics() );
  RegisterPhysics(new AntiBaryonPhysics() );
  RegisterPhysics(new IonPhysics() );
  RegisterPhysics(new GammaLeptoNuclearPhysics() );

  // Gamma-Nuclear Physics
//  G4EmExtraPhysics* gnuc = new G4EmExtraPhysics(verb);
//  gnuc->ElectroNuclear(false);
//  gnuc->MuonNuclear(false);
//  RegisterPhysics(gnuc);

  // Optical physics
  RegisterPhysics(new G4OpticalPhysics() );
}


EICPhysicsList::~EICPhysicsList()
{}


void EICPhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}


void EICPhysicsList::SetCuts()
{
  SetCutValue(0.7*mm, "proton");
  SetCutValue(0.7*mm, "e-");
  SetCutValue(0.7*mm, "e+");
  SetCutValue(0.7*mm, "gamma");      
}

