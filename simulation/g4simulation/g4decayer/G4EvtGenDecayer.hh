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
// $Id: G4EvtGenDecayer.hh,v 1.2 2014/10/07 03:06:54 mccumber Exp $
//
/// \file eventgenerator/pythia/decayer6/include/G4EvtGenDecayer.hh
/// \brief Definition of the G4EvtGenDecayer class
//

//#include "G4EvtGenDecayerMessenger.hh"


#include <Geant4/G4ThreeVector.hh>              // for G4ThreeVector
#include <Geant4/G4Types.hh>                    // for G4int, G4bool
#include <Geant4/G4VExtDecayer.hh>
#include <TParticle.h>   // for operator<<
#include <TClonesArray.h>   // for operator<<
#include "PHEvtGenDecayer.h"


class G4DecayProducts;
class G4DynamicParticle;
class G4ParticleDefinition;
class G4Track;

namespace CLHEP { class HepLorentzVector; }

/// Pythia6 decayer
///
/// Implements the G4VExtDecayer abstract class using the Pythia6 interface.
/// According to TPythia6Decayer class in Root:
/// http://root.cern.ch/
/// see http://root.cern.ch/root/License.html

class G4EvtGenDecayer : public G4VExtDecayer
{
  public:

    G4EvtGenDecayer();
    virtual ~G4EvtGenDecayer();

    virtual G4DecayProducts* ImportDecayProducts(const G4Track& track);
    
    void SetVerboseLevel(G4int verboseLevel) { fVerboseLevel =  verboseLevel; }
    
  private:

    /// Not implemented
    G4EvtGenDecayer(const G4EvtGenDecayer& right);
    /// Not implemented
    G4EvtGenDecayer& operator=(const G4EvtGenDecayer& right);
	PHEvtGenDecayer* myEvtGenDecayer = NULL;
	EvtGen *myGenerator;
  
    G4ParticleDefinition*
    GetParticleDefinition(const TParticle* p,G4bool warn = true) const;
    
 
    bool IsG4Detectable(const TParticle* p,G4bool warn = true) const;

    



    G4int            fVerboseLevel;        ///< verbose level
  //  ParticleVector*  fDecayProductsArray ; ///< array of decay products
	TLorentzVector FourMom;
	
};

// ----------------------------------------------------------------------------


