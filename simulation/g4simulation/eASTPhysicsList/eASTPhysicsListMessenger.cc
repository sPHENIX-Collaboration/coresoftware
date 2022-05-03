// ********************************************************************
//
// eASTPhysicsListMessenger.cc
//   A messenger class that handles physics list options.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "eASTPhysicsListMessenger.hh"

#include "eASTPhysicsList.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

eASTPhysicsListMessenger::eASTPhysicsListMessenger(eASTPhysicsList* pl)
: pPL(pl)
{
  G4UIparameter* param = nullptr;

  physDir = new G4UIdirectory("/eAST/physics/");
  physDir->SetGuidance("eAST physics selection");

  //addHPCmd = new G4UIcmdWithoutParameter("/eAST/physics/addHP",this);
  //addHPCmd->AvailableForStates(G4State_PreInit);
  //addHPCmd->SetToBeBroadcasted(false);
  //addHPCmd->SetGuidance("Add High-Precision neutron model.");
  //addHPCmd->SetGuidance(" Note: Shielding option has already had HP. This command does not make effect to Shielding option.");

  addRDMCmd = new G4UIcmdWithoutParameter("/eAST/physics/addRDM",this);
  addRDMCmd->AvailableForStates(G4State_PreInit);
  addRDMCmd->SetToBeBroadcasted(false);
  addRDMCmd->SetGuidance("Add Radioactive Decay model.");
  addRDMCmd->SetGuidance(" Note: Shielding option has already had RDM. This command does not make effect to Shielding option.");

  addOpticalCmd = new G4UIcmdWithoutParameter("/eAST/physics/addOptical",this);
  addOpticalCmd->AvailableForStates(G4State_PreInit);
  addOpticalCmd->SetToBeBroadcasted(false);
  addOpticalCmd->SetGuidance("Add Optical physics");

  addStepLimitCmd = new G4UIcmdWithAString("/eAST/physics/addStepLimit",this);
  addStepLimitCmd->AvailableForStates(G4State_PreInit);
  addStepLimitCmd->SetToBeBroadcasted(false);
  addStepLimitCmd->SetGuidance("Add step-limiter process to artificially limit step length.");
  addStepLimitCmd->SetGuidance("Specify particle types to be applied.");
  addStepLimitCmd->SetGuidance("  charged (default) : applied only to the charged particles");
  addStepLimitCmd->SetGuidance("  neutral : applied only to the neutral particles");
  addStepLimitCmd->SetGuidance("  all : applied to all particle types");
  addStepLimitCmd->SetGuidance("  e+/- : applied only to e+/e-");
  addStepLimitCmd->SetGuidance(" Note: In addition to this command, you need to specify the limitation value by");
  addStepLimitCmd->SetGuidance("       /eAST/physics/limit/stepLimit or /eAST/physics/limit/localStepLimt command.");
  addStepLimitCmd->SetParameterName("particle",true);
  addStepLimitCmd->SetDefaultValue("charged");
  addStepLimitCmd->SetCandidates("charged neutral all e+/-");

  physLimitDir = new G4UIdirectory("/eAST/physics/limit/");
  physLimitDir->SetGuidance("Specify step limitation");

  setStepLimitCmd = new G4UIcmdWithADoubleAndUnit("/eAST/physics/limit/stepLimit",this);
  setStepLimitCmd->AvailableForStates(G4State_Idle);
  setStepLimitCmd->SetToBeBroadcasted(false);
  setStepLimitCmd->SetParameterName("length",false);
  setStepLimitCmd->SetDefaultUnit("mm");
  setStepLimitCmd->SetGuidance("Define the limitation of the step length");
  setStepLimitCmd->SetGuidance("This limitation is applied to the entire geometry except regions that has its dedicated limit.");

  setRegionStepLimitCmd = new G4UIcommand("/eAST/physics/limit/regionStepLimit",this);
  setRegionStepLimitCmd->AvailableForStates(G4State_Idle);
  setRegionStepLimitCmd->SetToBeBroadcasted(false);
  setRegionStepLimitCmd->SetGuidance("Define the limitation of the step length for the specified region");
  setRegionStepLimitCmd->SetGuidance("   [usage] /eAST/physics/limit/regionStepLimit region length [unit]");
  setRegionStepLimitCmd->SetGuidance("      region (string) : region name");
  setRegionStepLimitCmd->SetGuidance(" Note: Region has to be defined in advance to this command.");
  setRegionStepLimitCmd->SetGuidance("       If new region is necessary, use /eAST/geometry/createRegion to create it.");
  param = new G4UIparameter("region",'s',false);
  setRegionStepLimitCmd->SetParameter(param);
  param = new G4UIparameter("length",'d',false);
  setRegionStepLimitCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setRegionStepLimitCmd->SetParameter(param);

  physCutDir = new G4UIdirectory("/eAST/physics/cuts/");
  physCutDir->SetGuidance("Specify production thresholds (a.k.a. cuts)");

  setCutCmd = new G4UIcmdWithADoubleAndUnit("/eAST/physics/cuts/setCuts",this);
  setCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  setCutCmd->SetToBeBroadcasted(false);
  setCutCmd->SetParameterName("length",false);
  setCutCmd->SetDefaultUnit("mm");
  setCutCmd->SetGuidance("Specify production thresholds (a.k.a. cuts) that is applied to the entire geometry");
  setCutCmd->SetGuidance("This threshold is applied to all of e-, e+, gamma and proton.");
  setCutCmd->SetGuidance("Threshold of each particle can be overwitted by /eAST/physics/cuts/setParticleCut command");

  setCutParticleCmd = new G4UIcommand("/eAST/physics/cuts/setParticleCut",this);
  setCutParticleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  setCutParticleCmd->SetToBeBroadcasted(false);
  setCutParticleCmd->SetGuidance("Specify production threshold (a.k.a. cut) for the specified particle that is applied to the entire geometry");
  setCutParticleCmd->SetGuidance("  [usage] /eAST/physics/setParticleCut particle cut unit");
  param = new G4UIparameter("particle",'s',false);
  param->SetParameterCandidates("e- e+ gamma proton");
  setCutParticleCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  setCutParticleCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setCutParticleCmd->SetParameter(param);

  setCutRegionCmd = new G4UIcommand("/eAST/physics/cuts/setRegionCut",this);
  setCutRegionCmd->AvailableForStates(G4State_Idle);
  setCutRegionCmd->SetToBeBroadcasted(false);
  setCutRegionCmd->SetGuidance("Specify production threshold (a.k.a. cut) that is applied to the specified region");
  setCutRegionCmd->SetGuidance("  [usage] /eAST/physics/setRegionCut region cut unit");
  setCutRegionCmd->SetGuidance("This threshold is applied to all of e-, e+, gamma and proton.");
  setCutRegionCmd->SetGuidance("Threshold of each particle can be overwitted by /eAST/physics/cuts/setRegionParticleCut command");
  setCutRegionCmd->SetGuidance(" Note: Region has to be defined in advance to this command.");
  setCutRegionCmd->SetGuidance("       If new region is necessary, use /eAST/geometry/createRegion to create it.");
  param = new G4UIparameter("region",'s',false);
  setCutRegionCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  setCutRegionCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setCutRegionCmd->SetParameter(param);

  setCutRegionParticleCmd = new G4UIcommand("/eAST/physics/cuts/setRegionParticleCut",this);
  setCutRegionParticleCmd->AvailableForStates(G4State_Idle);
  setCutRegionParticleCmd->SetToBeBroadcasted(false);
  setCutRegionParticleCmd->SetGuidance("Specify production threshold (a.k.a. cut) that is applied to the specified region");
  setCutRegionParticleCmd->SetGuidance("  [usage] /eAST/physics/setRegionParticleCut region particle cut unit");
  setCutRegionParticleCmd->SetGuidance(" Note: Region has to be defined in advance to this command.");
  setCutRegionParticleCmd->SetGuidance("       If new region is necessary, use /eAST/geometry/createRegion to create it.");
  param = new G4UIparameter("region",'s',false);
  setCutRegionParticleCmd->SetParameter(param);
  param = new G4UIparameter("particle",'s',false);
  param->SetParameterCandidates("e- e+ gamma proton");
  setCutRegionParticleCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  setCutRegionParticleCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setCutRegionParticleCmd->SetParameter(param);

}

eASTPhysicsListMessenger::~eASTPhysicsListMessenger()
{
  //delete addHPCmd;
  delete addRDMCmd;
  delete addOpticalCmd;
  delete addStepLimitCmd;
  delete setStepLimitCmd;
  delete setRegionStepLimitCmd;
  delete setCutCmd;
  delete setCutParticleCmd;
  delete setCutRegionCmd;
  delete setCutRegionParticleCmd;

  delete physLimitDir;
  delete physCutDir;
  delete physDir;
}

#include "G4Tokenizer.hh"

void eASTPhysicsListMessenger::SetNewValue(G4UIcommand* cmd, G4String val)
{
  //if(cmd==addHPCmd)
  //{ pPL->AddHP(); }
  //else 
  if(cmd==addRDMCmd)
  { pPL->AddRDM(); }
  else if(cmd==addOpticalCmd)
  { pPL->AddOptical(); }
  else if(cmd==addStepLimitCmd)
  {
    G4int opt = 0;
    if(val=="neutral") opt = 1; 
    else if(val=="all") opt = 2; 
    else if(val=="e+/-") opt = 3; 
    pPL->AddStepLimit(opt);
  }
  else if(cmd==setStepLimitCmd)
  { pPL->SetGlobalStepLimit(setStepLimitCmd->GetNewDoubleValue(val)); }
  else if(cmd==setRegionStepLimitCmd)
  {
    G4Tokenizer next(val);
    G4String reg = next();
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    auto regPtr = pPL->SetLocalStepLimit(reg,setRegionStepLimitCmd->ConvertToDimensionedDouble(newVal));
    if(!regPtr)
    {
      G4ExceptionDescription ed;
      ed << "Region <" << reg << "> is not defined.";
      setRegionStepLimitCmd->CommandFailed(ed);
    }
  }
  else if(cmd==setCutCmd)
  { pPL->SetGlobalCuts(setCutCmd->GetNewDoubleValue(val)); }
  else if(cmd==setCutParticleCmd)
  {
    G4Tokenizer next(val);
    G4String pat = next();
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    G4int i = 0;
    if(pat=="e-") i = 0; 
    else if(pat=="e+") i = 1; 
    else if(pat=="gamma") i = 2; 
    else if(pat=="proton") i = 3; 
    pPL->SetGlobalCut(i,setCutParticleCmd->ConvertToDimensionedDouble(newVal));
  }
  else if(cmd==setCutRegionCmd)
  {
    G4Tokenizer next(val);
    G4String reg = next();
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    auto regPtr = pPL->SetLocalCuts(reg,setCutRegionCmd->ConvertToDimensionedDouble(newVal));
    if(!regPtr)
    {
      G4ExceptionDescription ed;
      ed << "Region <" << reg << "> is not defined.";
      setRegionStepLimitCmd->CommandFailed(ed);
    }
  }
  else if(cmd==setCutRegionParticleCmd)
  {
    G4Tokenizer next(val);
    G4String reg = next();
    G4String pat = next();
    G4int i = 0;
    if(pat=="e-") i = 0; 
    else if(pat=="e+") i = 1; 
    else if(pat=="gamma") i = 2; 
    else if(pat=="proton") i = 3; 
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    auto regPtr = pPL->SetLocalCut(reg,i,setCutRegionParticleCmd->ConvertToDimensionedDouble(newVal));
    if(!regPtr)
    {
      G4ExceptionDescription ed;
      ed << "Region <" << reg << "> is not defined.";
      setRegionStepLimitCmd->CommandFailed(ed);
    }
  }
  
}

G4String eASTPhysicsListMessenger::GetCurrentValue(G4UIcommand* cmd)
{
  G4String val("");

  //if(cmd==addHPCmd)
  //{ val = cmd->ConvertToString(pPL->IfHP()); }
  //else 
  if(cmd==addRDMCmd)
  { val = cmd->ConvertToString(pPL->IfRDM()); }
  else if(cmd==addOpticalCmd)
  { val = cmd->ConvertToString(pPL->IfOptical()); }
  else if(cmd==addStepLimitCmd)
  {
    auto opt = pPL->IfStepLimit();
    switch(opt)
    {
      case 0: val =  "charged"; break;
      case 1: val =  "neutral"; break;
      case 2: val =  "all"; break;
      case 3: val =  "e+/-"; break;
      default : val = "undefined"; break;
    }
  }
  return val;
}


