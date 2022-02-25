////////////////////////////////////////////////////////////////////////////////
//
//  eASTPhysicsList.cc
//  Geant4 physics list for Electron Ion Collider detector simulation
//
//  History
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.05.2021 : migration to eAST - Makoto Asai (SLAC)
//    Dec.22.2021 : migration to Geant4 version 11.0 - Makoto Asai (JLab)
//
////////////////////////////////////////////////////////////////////////////////

#include "eASTPhysicsList.hh"
#include "eASTPhysicsListMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessTable.hh"
#include "G4RegionStore.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4HadronicParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4NuclideTable.hh"

#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterPhysics.hh"

#include "eASTProtonPhysics.hh"
#include "eASTNeutronPhysics.hh"
#include "eASTPionPhysics.hh"
#include "eASTKaonPhysics.hh"
#include "eASTHyperonPhysics.hh"
#include "eASTAntiBaryonPhysics.hh"
#include "eASTIonPhysics.hh"
#include "eASTGammaLeptoNuclearPhysics.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4Version.hh"

eASTPhysicsList::eASTPhysicsList(G4int verb)
 :G4VModularPhysicsList()
{
  SetVerboseLevel(verb);
  
  //add new units for radioActive decays
  //
  new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
  // 
#if G4VERSION_NUMBER < 1100
  const G4double minute = 60*second;
  const G4double hour   = 60*minute;
  const G4double day    = 24*hour;
  const G4double year   = 365*day;
  new G4UnitDefinition("minute", "min", "Time", minute);
  new G4UnitDefinition("hour",   "h",   "Time", hour);
  new G4UnitDefinition("day",    "d",   "Time", day);
  new G4UnitDefinition("year",   "y",   "Time", year);
#endif

  // Mandatory for G4NuclideTable
  // Half-life threshold must be set small or many short-lived isomers 
  // will not be assigned life times (default to 0) 
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

  globalCuts[2] = 0.7*mm; //gamma
  globalCuts[0] = 0.7*mm; //e-
  globalCuts[1] = 0.7*mm; //e+
  globalCuts[3] = 0.7*mm; //proton

  pMessenger = new eASTPhysicsListMessenger(this);
}


eASTPhysicsList::~eASTPhysicsList()
{ delete pMessenger; }

void eASTPhysicsList::ConstructParticle()
{
  if(!processesAreRegistered) SetupProcesses();

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

void eASTPhysicsList::SetupProcesses()
{
  // EM physics
  RegisterPhysics(new G4EmStandardPhysics());
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetAugerCascade(true);
  param->SetStepFunction(1., 1.*CLHEP::mm);
  param->SetStepFunctionMuHad(1., 1.*CLHEP::mm);

  // Decay
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
  if(addRDM)
  { RegisterPhysics(new G4RadioactiveDecayPhysics()); }

  // Hadronic physics
  RegisterPhysics(new eASTProtonPhysics() );
  RegisterPhysics(new eASTNeutronPhysics() );
  RegisterPhysics(new eASTPionPhysics() );
  RegisterPhysics(new eASTKaonPhysics() );
  RegisterPhysics(new eASTHyperonPhysics() );
  RegisterPhysics(new eASTAntiBaryonPhysics() );
  RegisterPhysics(new eASTIonPhysics() );
  RegisterPhysics(new eASTGammaLeptoNuclearPhysics() );

  // Optical physics
  if(addOptical)
  { RegisterPhysics(new G4OpticalPhysics() ); }

  // Step limiter
  if(stepLimit_opt>=0)
  { RegisterPhysics(new G4StepLimiterPhysics()); }

  processesAreRegistered = true;
}

void eASTPhysicsList::ConstructProcess()
{
  auto verbose = G4ProcessTable::GetProcessTable()->GetVerboseLevel();
  SetVerboseLevel(verbose);
  G4EmParameters::Instance()->SetVerbose(verbose);

#if G4VERSION_NUMBER >= 1100
  G4HadronicParameters::Instance()->SetVerboseLevel(verbose);
#endif
  G4VModularPhysicsList::ConstructProcess();
}

void eASTPhysicsList::SetCuts()
{
  SetCutValue(globalCuts[2],"gamma"); // gamma should be defined first!
  SetCutValue(globalCuts[0],"e-");
  SetCutValue(globalCuts[1],"e+");
  SetCutValue(globalCuts[3],"proton");
}

G4Region* eASTPhysicsList::FindRegion(const G4String& reg) const
{
  auto store = G4RegionStore::GetInstance();
  return store->GetRegion(reg);
}

void eASTPhysicsList::SetGlobalCuts(G4double val)
{
  for(G4int i=0; i<4; i++)
  { SetGlobalCut(i,val); }
  SetCuts();
}

void eASTPhysicsList::SetGlobalCut(G4int i, G4double val)
{
  globalCuts[i] = val;
  SetCuts();
}

G4Region* eASTPhysicsList::SetLocalCut(const G4String& reg,G4int i,G4double val)
{
  auto regPtr = FindRegion(reg);
  if(!regPtr) return regPtr;

  auto cuts = regPtr->GetProductionCuts();
  if(!cuts)
  {
    cuts = new G4ProductionCuts();
    regPtr->SetProductionCuts(cuts);
  }

  cuts->SetProductionCut(val,i);
  return regPtr;
}

G4double eASTPhysicsList::GetLocalCut(const G4String& reg,G4int i) const
{
  auto regPtr = FindRegion(reg);
  G4double val = -1.0;
  if(regPtr)
  {
    auto cuts = regPtr->GetProductionCuts();
    if(cuts) val = cuts->GetProductionCut(i);
  }
  return val;
}

#include "G4UserLimits.hh"

G4Region* eASTPhysicsList::SetLocalStepLimit(const G4String& reg,G4double val)
{
  auto regPtr = FindRegion(reg);
  if(!regPtr) return regPtr;

  auto uLim = regPtr->GetUserLimits();
  if(!uLim)
  {
    uLim = new G4UserLimits(val);
    regPtr->SetUserLimits(uLim);
  }
  else
  { uLim->SetMaxAllowedStep(val); }
  return regPtr;
}

#include "G4Track.hh"
G4double eASTPhysicsList::GetLocalStepLimit(const G4String& reg) const
{
  static G4Track dummyTrack;
  auto regPtr = FindRegion(reg);
  G4double val = -1.0;
  if(regPtr)
  {
    auto uLim = regPtr->GetUserLimits();
    if(uLim) val = uLim->GetMaxAllowedStep(dummyTrack);
  }
  return val;
}

void eASTPhysicsList::SetGlobalStepLimit(G4double val)
{ SetLocalStepLimit("DefaultRegionForTheWorld",val); }

G4double eASTPhysicsList::GetGlobalStepLimit() const
{ return GetLocalStepLimit("DefaultRegionForTheWorld"); }

