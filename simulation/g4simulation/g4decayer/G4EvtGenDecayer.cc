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
// $Id: G4EvtGenDecayer.cc,v 1.2 2014/10/07 03:06:54 mccumber Exp $
//
/// \file eventgenerator/pythia/decayer6/src/G4EvtGenDecayer.cc
/// \brief Implementation of the G4EvtGenDecayer class

// ----------------------------------------------------------------------------
// According to TPythia6Decayer class in Root:
// http://root.cern.ch/
// see http://root.cern.ch/root/License.html
// ----------------------------------------------------------------------------

#include "G4EvtGenDecayer.hh"
#include "PHEvtGenRandomEngine.hh"  // for PHEvtGenRandomEngine

#include <EvtGen/EvtGen.hh>                      // for EvtGen
#include <EvtGenBase/EvtAbsRadCorr.hh>           // for EvtAbsRadCorr
#include <EvtGenBase/EvtHepMCEvent.hh>           // for EvtHepMCEvent
#include <EvtGenBase/EvtId.hh>                   // for EvtId
#include <EvtGenBase/EvtPDL.hh>                  // for EvtPDL
#include <EvtGenBase/EvtParticle.hh>             // for EvtParticle
#include <EvtGenBase/EvtParticleFactory.hh>      // for EvtParticleFactory
#include <EvtGenBase/EvtRandom.hh>               // for EvtRandom
#include <EvtGenBase/EvtVector4R.hh>             // for EvtVector4R
#include <EvtGenExternal/EvtExternalGenList.hh>  // for EvtExternalGenList

#include <HepMC3/Attribute.h>        // for string
#include <HepMC3/FourVector.h>       // for FourVector
#include <HepMC3/GenEvent.h>         // for GenEvent
#include <HepMC3/GenParticle.h>      // for GenParticle
#include <HepMC3/GenParticle_fwd.h>  // for GenParticlePtr
#include <HepMC3/GenVertex.h>        // for GenVertex
#include <HepMC3/GenVertex_fwd.h>    // for GenVertexPtr
#include <HepMC3/Print.h>            // for Print

#include <Geant4/G4DecayProducts.hh>
#include <Geant4/G4DynamicParticle.hh>
#include <Geant4/G4ParticleDefinition.hh>  // for G4ParticleDefinition
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Track.hh>
#include <Geant4/G4Types.hh>        // for G4int, G4double
#include <Geant4/G4VExtDecayer.hh>  // for G4VExtDecayer

#include <algorithm>  // for copy, max
#include <cstdlib>    // for abs
#include <iostream>   // for operator<<, basic_ostream:...
#include <memory>     // for __shared_ptr_access
#include <stack>      // for operator<<
#include <string>     // for operator<<
#include <vector>     // for vector

class EvtDecayBase;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EvtGenDecayer::G4EvtGenDecayer()
  : G4VExtDecayer("G4EvtGenDecayer")
  , fVerboseLevel(0)
{
  mEvtGenRandomEngine = new PHEvtGenRandomEngine();

  EvtRandom::setRandomEngine(static_cast<EvtRandomEngine*>(mEvtGenRandomEngine));

  radCorrEngine = genList.getPhotosModel();

  extraModels = genList.getListOfModels();

  // the hardcoded paths are temporary
  const char* offline_main = getenv("OFFLINE_MAIN");

  if (!offline_main)
  {
    exit(1);
  }

  string decay = string(offline_main) + "/share/EvtGen/DECAY.DEC";  // Using PDG 2019 reference as the input for now
  string evt = string(offline_main) + "/share/EvtGen/evt.pdl";

  mEvtGen = new EvtGen(decay, evt, static_cast<EvtRandomEngine*>(mEvtGenRandomEngine), radCorrEngine, &extraModels);
  extraModels.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EvtGenDecayer::~G4EvtGenDecayer()
{
  delete mEvtGen;

  delete radCorrEngine;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4EvtGenDecayer::
    GetParticleDefinition(int ParPDGID) const
{
  /// Return G4 particle definition for given TParticle

  // get particle definition from G4ParticleTable
  G4int pdgEncoding = ParPDGID;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particleDefinition = nullptr;
  if (pdgEncoding != 0)
  {
    particleDefinition = particleTable->FindParticle(pdgEncoding);
  }

  return particleDefinition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DecayProducts* G4EvtGenDecayer::ImportDecayProducts(const G4Track& track)
{
  /// Import decay products

  G4ThreeVector momentum = track.GetMomentum();
  G4double etot = track.GetDynamicParticle()->GetTotalEnergy();

  G4ParticleDefinition* particleDef = track.GetDefinition();
  G4int pdgEncoding = particleDef->GetPDGEncoding();

  EvtVector4R p_init(etot / GeV, momentum.x() / GeV, momentum.y() / GeV, momentum.z() / GeV);

  EvtId parentID = EvtPDL::evtIdFromLundKC(pdgEncoding);

  auto mParticle = EvtParticleFactory::particleFactory(parentID, p_init);
  // Decay the particle
  mEvtGen->generateDecay(mParticle);

  EvtHepMCEvent theEvent;

  theEvent.constructEvent(mParticle);

  HepMC3::GenEvent* evt = theEvent.getEvent();  // Directly use HepMC3 records

  G4DecayProducts* decayProducts = new G4DecayProducts(*(track.GetDynamicParticle()));

  std::stack<HepMC3::GenParticlePtr> stack;

  auto part = evt->particles()[0];

  if (part->pdg_id() != pdgEncoding)
  {
    std::cout << "Issue Found: The first particle in the decay chain is NOT the incoming particle decayed by EvtGen" << std::endl;
  }

  stack.push(part);

  while (!stack.empty())
  {
    auto particle = stack.top();
    stack.pop();

    for (const auto& children : particle->children())
    {
      float EvtWidth = 9999;
      int pdg = children->pdg_id();

      G4ParticleDefinition* particleDefinition = GetParticleDefinition(pdg);
      if (EvtPDL::evtIdFromLundKC(pdg).getId() != -1)
      {
        EvtWidth = EvtPDL::getWidth(EvtPDL::evtIdFromLundKC(pdg)) * 1000000;  // Decay Width Unit: GeV -> keV to define stable particles in EvtGen -> G4
      }

      if (particleDefinition && EvtWidth < WidthThreshold)
      {
        bool SameVtxPos = true;

        if (children->production_vertex()->position().x() != particle->end_vertex()->position().x() || children->production_vertex()->position().y() != particle->end_vertex()->position().y() || children->production_vertex()->position().z() != particle->end_vertex()->position().z() || children->production_vertex()->position().t() != particle->end_vertex()->position().t())
        {
          SameVtxPos = false;
        }
        if (!SameVtxPos)
        {
          std::cout << "Issue Found: in this vertex, particles pushed to GEANT have different production positions, need to check and understand why!!!" << std::endl;
          std::cout << "This is due to the particle PDGID: " << children->pdg_id() << "  from the decay vertex of the parent particle:" << particle->pdg_id() << std::endl;
          std::cout << "Here is the Event Content" << std::endl;
          const HepMC3::GenEvent& CheckEvent = *evt;
          HepMC3::Print::content(CheckEvent);
        }

        G4ThreeVector G4Mom = G4ThreeVector(children->momentum().px() * GeV, children->momentum().py() * GeV, children->momentum().pz() * GeV);
        G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDefinition, G4Mom);

        if (dynamicParticle)
        {
          if (fVerboseLevel > 0)
          {
            std::cout << "  G4 particle name: "
                      << dynamicParticle->GetDefinition()->GetParticleName()
                      << std::endl;
          }

          decayProducts->PushProducts(dynamicParticle);
        }
        else
        {
          std::cout << "Dynamical particle to be pushed from EvtGen to GEANT 4 is not properly defined: need to check why this happens!" << std::endl;
          exit(1);
        }
      }
      else
      {
        stack.push(children);
      }
    }
  }

  evt->clear();
  mParticle->deleteTree();
  return decayProducts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EvtGenDecayer::SetDecayTable(const string& decayTable, bool useXml)
{
  mEvtGen->readUDecay(decayTable, useXml);
}
