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
// $Id: EvtGenExtDecayerPhysics.hh,v 1.3 2015/01/05 23:49:39 mccumber Exp $
//
/// \file eventgenerator/pythia/decayer6/include/EvtGenExtDecayerPhysics.hh
/// \brief Definition of the EvtGenExtDecayerPhysics class
///
/// \author I. Hrivnacova; IPN Orsay

#ifndef EVT_EXT_DECAYER_PHYSICS_H
#define EVT_EXT_DECAYER_PHYSICS_H

#include <Geant4/G4String.hh>
#include <Geant4/G4VPhysicsConstructor.hh>

/// The builder for external decayer.
///
/// The external decayer is added to all instantiated decay
/// processes
///
/// \author I. Hrivnacova; IPN Orsay

class EvtGenExtDecayerPhysics : public G4VPhysicsConstructor
{
 public:
  EvtGenExtDecayerPhysics(const G4String& name = "ExtDecayer");
  virtual ~EvtGenExtDecayerPhysics();

 protected:
  // methods
  // construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();

 private:
  /// Not implemented
  EvtGenExtDecayerPhysics(const EvtGenExtDecayerPhysics& right);
  /// Not implemented
  EvtGenExtDecayerPhysics& operator=(const EvtGenExtDecayerPhysics& right);
};

#endif  // P6D_EXT_DECAYER_PHYSICS_H
