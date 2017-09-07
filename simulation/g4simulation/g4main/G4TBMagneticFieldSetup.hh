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
// $Id: G4TBMagneticFieldSetup.hh,v 1.10 2014/11/14 21:47:38 mccumber Exp $
// GEANT4 tag $Name:  $
//
//    A class for control of the Electric Field of the detector.
//  The field for this case is uniform.
//
//  It is simply a 'setup' class that creates the field and necessary other parts
//


#ifndef G4TBElectricFieldSetup_H
#define G4TBElectricFieldSetup_H

#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4UniformMagField.hh>

#include <string>

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4MagInt_Driver; 
class G4TBFieldMessenger;
class PHField;

class G4TBMagneticFieldSetup 
{
public:


  G4TBMagneticFieldSetup(PHField * phfield) ;
//  G4TBMagneticFieldSetup(const float magfield) ;
//  G4TBMagneticFieldSetup(const std::string &fieldmapfile, const int mapdim, const float magfield_rescale = 1.0) ;

  virtual ~G4TBMagneticFieldSetup() ;  

  void Verbosity(const int verb) {verbosity = verb;}
  
  void SetStepperType( const G4int i) { fStepperType = i ; }

  void SetStepper();

  void SetMinStep(const G4double s) { fMinStep = s ; }

  void UpdateField();

  void SetFieldValue(const G4ThreeVector fieldVector);
  void SetFieldValue(const G4double      fieldValue);

  double get_magfield_at_000(const int i) const {return magfield_at_000[i];}

protected:

      // Find the global Field Manager
  G4FieldManager*         GetGlobalFieldManager() ;

private:
  int verbosity;
  
  G4FieldManager*         fFieldManager ;

  G4ChordFinder*          fChordFinder ;

 G4Mag_UsualEqRhs*   fEquation ;

  G4MagneticField*        fEMfield;
 
  G4ThreeVector           fElFieldValue ; 

  G4MagIntegratorStepper* fStepper ;
  G4MagInt_Driver*        fIntgrDriver;

  G4int                   fStepperType ;

  G4double                fMinStep ;
 
  G4TBFieldMessenger*      fFieldMessenger;

  double magfield_at_000[3];
};

#endif
