//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4TBFieldMessenger.cc,v 1.5 2012/07/10 16:48:20 pinkenbu Exp $
// GEANT4 tag $Name:  $
//
//

#include "G4TBFieldMessenger.hh"

#include "G4TBMagneticFieldSetup.hh"

#include <Geant4/G4ApplicationState.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4UIcmdWithADoubleAndUnit.hh>
#include <Geant4/G4UIcmdWithAString.hh>
#include <Geant4/G4UIcmdWithAnInteger.hh>
#include <Geant4/G4UIcmdWithoutParameter.hh>
#include <Geant4/G4UIdirectory.hh>

class G4UIcommand;

//////////////////////////////////////////////////////////////////////////////

G4TBFieldMessenger::G4TBFieldMessenger(G4TBMagneticFieldSetup* pEMfield)
  : fEFieldSetup(pEMfield)
{
  G4TBdetDir = new G4UIdirectory("/field/");
  G4TBdetDir->SetGuidance("G4TB field tracking control.");

  StepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType", this);
  StepperCmd->SetGuidance("Select stepper type for electric field");
  StepperCmd->SetParameterName("choice", true);
  StepperCmd->SetDefaultValue(4);
  StepperCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/field/update", this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

  ElFieldCmd = new G4UIcmdWithADoubleAndUnit("/field/setFieldZ", this);
  ElFieldCmd->SetGuidance("Define uniform Electric field.");
  ElFieldCmd->SetGuidance("Electric field will be in Z direction.");
  ElFieldCmd->SetGuidance("Value of Electric field has to be given in volt/m");
  ElFieldCmd->SetParameterName("Ez", false, false);
  ElFieldCmd->SetDefaultUnit("volt/m");
  ElFieldCmd->AvailableForStates(G4State_Idle);

  MinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep", this);
  MinStepCmd->SetGuidance("Define minimal step");
  MinStepCmd->SetParameterName("min step", false, false);
  MinStepCmd->SetDefaultUnit("mm");
  MinStepCmd->AvailableForStates(G4State_Idle);

  AbsMaterCmd = new G4UIcmdWithAString("/field/setAbsMat", this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice", true);
  AbsMaterCmd->SetDefaultValue("Xe");
  AbsMaterCmd->AvailableForStates(G4State_Idle);
}

///////////////////////////////////////////////////////////////////////////////

G4TBFieldMessenger::~G4TBFieldMessenger()
{
  delete StepperCmd;
  delete ElFieldCmd;
  delete MinStepCmd;
  delete G4TBdetDir;
  delete UpdateCmd;

  delete AbsMaterCmd;
}

////////////////////////////////////////////////////////////////////////////
//
//

void G4TBFieldMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == StepperCmd)
  {
    fEFieldSetup->SetStepperType(StepperCmd->GetNewIntValue(newValue));
  }
  if (command == UpdateCmd)
  {
    fEFieldSetup->UpdateField();
  }
  if (command == ElFieldCmd)
  {
    fEFieldSetup->SetFieldValue(ElFieldCmd->GetNewDoubleValue(newValue));
  }
  if (command == MinStepCmd)
  {
    fEFieldSetup->SetMinStep(MinStepCmd->GetNewDoubleValue(newValue));
  }
}

//
//
/////////////////////////////////////////////////////////////////////////
