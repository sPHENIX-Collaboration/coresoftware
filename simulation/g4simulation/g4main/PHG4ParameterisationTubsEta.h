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
// $Id: PHG4ParameterisationTubsEta.h,v 1.3 2012/07/10 16:48:20 pinkenbu Exp $
// GEANT4 tag $Name:  $
//
// classes G4ParameterisationTubsRho
//         G4ParameterisationTubsPhi
//         G4ParameterisationTubsZ
//
// Class description:
//
// This class represents the parameterised positioning equivalent to 
// dividing a G4Tubs along one of each axis Rho, Phi, Z.

// History:
// -------
// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------
#ifndef PHG4ParameterisationTubsETA_H
#define PHG4ParameterisationTubsETA_H 1

#include <vector>
#include <Geant4/G4ParameterisationTubs.hh>

class G4VPhysicalVolume;

class PHG4ParameterisationTubsEta : public G4VParameterisationTubs
{ 
public:  // with description

  PHG4ParameterisationTubsEta( EAxis axis, G4int nCopies,
			       G4double offset, G4double step,
			       G4VSolid* motherSolid, DivisionType divType=DivNDIVandWIDTH );
  virtual ~PHG4ParameterisationTubsEta();
  
  virtual G4double GetMaxParameter() const;
  
  virtual void ComputeTransformation(const G4int copyNo,
				     G4VPhysicalVolume* physVol) const;
  void ComputeDimensions(G4Tubs& tubs, const G4int copyNo,
			 const G4VPhysicalVolume* physVol) const;

private:
  std::vector<double> _zpos; // Z positions of the rings
  std::vector<double> _zhalf; // Z half-widths of the rings

};

#endif
