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
// $Id: G4TBMagneticFieldSetup.cc,v 1.24 2014/11/14 21:47:38 mccumber Exp $
// GEANT4 tag $Name:  $
//
//
//   User Field class implementation.
//

#include "G4TBMagneticFieldSetup.hh"
#include "G4TBFieldMessenger.hh"
#include "PHG4MagneticField.h"

#include <Geant4/G4CashKarpRKF45.hh>
#include <Geant4/G4ChordFinder.hh>
#include <Geant4/G4ClassicalRK4.hh>
#include <Geant4/G4ExplicitEuler.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4ImplicitEuler.hh>
#include <Geant4/G4MagIntegratorDriver.hh>
#include <Geant4/G4MagIntegratorStepper.hh>
#include <Geant4/G4Mag_UsualEqRhs.hh>
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4SimpleHeum.hh>
#include <Geant4/G4SimpleRunge.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4Types.hh>  // for G4double, G4int
#include <Geant4/G4UniformMagField.hh>

#include <cassert>
#include <cstdlib>  // for exit, size_t
#include <iostream>
#include <string>

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

G4TBMagneticFieldSetup::G4TBMagneticFieldSetup(PHField* phfield)
{
  assert(phfield);

  fEMfield = new PHG4MagneticField(phfield);
  fFieldMessenger = new G4TBFieldMessenger(this);
  fEquation = new G4Mag_UsualEqRhs(fEMfield);
  fMinStep = 0.005 * mm;  // minimal step of 10 microns
  fStepperType = 4;       // ClassicalRK4 -- the default stepper

  fFieldManager = GetGlobalFieldManager();
  UpdateField();
  double point[4] = {0, 0, 0, 0};
  fEMfield->GetFieldValue(&point[0], &magfield_at_000[0]);
  for (size_t i = 0; i < sizeof(magfield_at_000) / sizeof(double); i++)
  {
    magfield_at_000[i] = magfield_at_000[i] / tesla;
  }
  if (verbosity > 0)
  {
    std::cout << "field: x" << magfield_at_000[0]
         << ", y: " << magfield_at_000[1]
         << ", z: " << magfield_at_000[2]
         << std::endl;
  }
}

//G4TBMagneticFieldSetup::G4TBMagneticFieldSetup(const float magfield)
//  : verbosity(0), fChordFinder(0), fStepper(0), fIntgrDriver(0)
//{
//  //solenoidal field along the axis of the cyclinders?
//  fEMfield = new G4UniformMagField(G4ThreeVector(0.0, 0.0, magfield*tesla));
//  fFieldMessenger = new G4TBFieldMessenger(this) ;
//  fEquation = new  G4Mag_UsualEqRhs(fEMfield);
//  fMinStep     = 0.005 * mm ; // minimal step of 10 microns
//  fStepperType = 4 ;        // ClassicalRK4 -- the default stepper
//
//  fFieldManager = GetGlobalFieldManager();
//  UpdateField();
//
//  double point[4] = {0,0,0,0};
//  fEMfield->GetFieldValue(&point[0],&magfield_at_000[0]);
//  for (size_t i=0; i<sizeof(magfield_at_000)/sizeof(double);i++)
//    {
//      magfield_at_000[i] = magfield_at_000[i]/tesla;
//    }
//}
//
///////////////////////////////////////////////////////////////////////////////////
//
//G4TBMagneticFieldSetup::G4TBMagneticFieldSetup(const string &fieldmapname, const int dim, const float magfield_rescale)
//  : verbosity(0),
//    fChordFinder(0),
//    fStepper(0),
//    fIntgrDriver(0)
//{
//
//  switch(dim)
//    {
//    case 1:
//      fEMfield = new PHG4FieldsPHENIX(fieldmapname,magfield_rescale);
//      break;
//    case 2:
//      fEMfield = new PHG4Field2D(fieldmapname,0,magfield_rescale);
//      break;
//    case 3:
//      fEMfield = new PHG4Field3D(fieldmapname,0,magfield_rescale);
//      break;
//    default:
//      std::cout << "Invalid dimension, valid is 2 for 2D, 3 for 3D" << std::endl;
//      exit(1);
//    }
//  fFieldMessenger = new G4TBFieldMessenger(this) ;
//  fEquation = new  G4Mag_UsualEqRhs(fEMfield);
//  fMinStep     = 0.005*mm ; // minimal step of 10 microns
//  fStepperType = 4 ;        // ClassicalRK4 -- the default stepper
//
//  fFieldManager = GetGlobalFieldManager();
//  UpdateField();
//  double point[4] = {0,0,0,0};
//  fEMfield->GetFieldValue(&point[0],&magfield_at_000[0]);
//  for (size_t i=0; i<sizeof(magfield_at_000)/sizeof(double);i++)
//    {
//      magfield_at_000[i] = magfield_at_000[i]/tesla;
//    }
//  if (verbosity > 0)
//    {
//      std::cout << "field: x" << magfield_at_000[0]
//	   << ", y: " << magfield_at_000[1]
//	   << ", z: " << magfield_at_000[2]
//	   << std::endl;
//    }
//}

////////////////////////////////////////////////////////////////////////////////

G4TBMagneticFieldSetup::~G4TBMagneticFieldSetup()
{
  delete fChordFinder;
  delete fStepper;
  delete fFieldMessenger;
  delete fEquation;
  delete fEMfield;
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void G4TBMagneticFieldSetup::UpdateField()
{
  SetStepper();

  fFieldManager->SetDetectorField(fEMfield);

  delete fChordFinder;

  fIntgrDriver = new G4MagInt_Driver(fMinStep,
                                     fStepper,
                                     fStepper->GetNumberOfVariables());

  fChordFinder = new G4ChordFinder(fIntgrDriver);

  fFieldManager->SetChordFinder(fChordFinder);
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void G4TBMagneticFieldSetup::SetStepper()
{
  G4int nvar = 8;

  delete fStepper;

  std::stringstream message;

  switch (fStepperType)
  {
  case 0:
    fStepper = new G4ExplicitEuler(fEquation, nvar);
    message << "Stepper in use: G4ExplicitEuler";
    break;
  case 1:
    fStepper = new G4ImplicitEuler(fEquation, nvar);
    message << "Stepper in use: G4ImplicitEuler";
    break;
  case 2:
    fStepper = new G4SimpleRunge(fEquation, nvar);
    message << "Stepper in use: G4SimpleRunge";
    break;
  case 3:
    fStepper = new G4SimpleHeum(fEquation, nvar);
    message << "Stepper in use: G4SimpleHeum";
    break;
  case 4:
    fStepper = new G4ClassicalRK4(fEquation, nvar);
    message << "Stepper in use: G4ClassicalRK4 (default)";
    break;
  case 5:
    fStepper = new G4CashKarpRKF45(fEquation, nvar);
    message << "Stepper in use: G4CashKarpRKF45";
    break;
  case 6:
    fStepper = nullptr;  // new G4RKG3_Stepper( fEquation, nvar );
    message << "G4RKG3_Stepper is not currently working for Magnetic Field";
    break;
  case 7:
    fStepper = nullptr;  // new G4HelixExplicitEuler( fEquation );
    message << "G4HelixExplicitEuler is not valid for Magnetic Field";
    break;
  case 8:
    fStepper = nullptr;  // new G4HelixImplicitEuler( fEquation );
    message << "G4HelixImplicitEuler is not valid for Magnetic Field";
    break;
  case 9:
    fStepper = nullptr;  // new G4HelixSimpleRunge( fEquation );
    message << "G4HelixSimpleRunge is not valid for Magnetic Field";
    break;
  default:
    fStepper = nullptr;
  }

  if (verbosity > 0)
  {
    std::cout << " ---------- G4TBMagneticFieldSetup::SetStepper() -----------" << std::endl;
    std::cout << "  " << message.str() << std::endl;
    std::cout << "  Minimum step size: " << fMinStep / mm << " mm" << std::endl;
    std::cout << " -----------------------------------------------------------" << std::endl;
  }

  if (!fStepper)
  {
    std::cout << "no stepper set, edxiting now" << std::endl;
    exit(1);
  }

  return;
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void G4TBMagneticFieldSetup::SetFieldValue(const G4double fieldValue)
{
  G4ThreeVector fieldVector(0.0, 0.0, fieldValue);

  SetFieldValue(fieldVector);
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field value to fieldVector
//

void G4TBMagneticFieldSetup::SetFieldValue(const G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr = GetGlobalFieldManager();

  if (fieldVector != G4ThreeVector(0., 0., 0.))
  {
    if (fEMfield) delete fEMfield;
    fEMfield = new G4UniformMagField(fieldVector);

    fEquation->SetFieldObj(fEMfield);  // must now point to the new field

    // UpdateField();

    fieldMgr->SetDetectorField(fEMfield);
  }
  else
  {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    delete fEMfield;
    fEMfield = nullptr;
    fEquation->SetFieldObj(fEMfield);  // As a double check ...
    fieldMgr->SetDetectorField(fEMfield);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager* G4TBMagneticFieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
