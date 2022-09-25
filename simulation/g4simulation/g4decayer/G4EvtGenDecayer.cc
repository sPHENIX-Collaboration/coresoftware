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

#include <Geant4/G4DecayProducts.hh>
#include <Geant4/G4DynamicParticle.hh>
#include <Geant4/G4ParticleDefinition.hh>  // for G4ParticleDefinition
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4VExtDecayer.hh>  // for G4VExtDecayer

#include <CLHEP/Vector/LorentzVector.h>

#include <cstdlib>   // for abs
#include <iostream>  // for operator<<, basic_ostream:...
#include <string>    // for operator<<

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EvtGenDecayer::G4EvtGenDecayer()
  : G4VExtDecayer("G4EvtGenDecayer")
  , fVerboseLevel(0)
{
  //bool WilluseXml =false;
  EvtRandomEngine* eng = 0;
  eng = new EvtSimpleRandomEngine();
  EvtRandom::setRandomEngine((EvtRandomEngine*) eng);
  PHEvtGenDecayer();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EvtGenDecayer::~G4EvtGenDecayer()
{
  /// Destructor
  //std::cout << __PRETTY_FUNCTION__ << std::endl;

  //	delete fDecayProductsArray;
}

//
// private methods
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4EvtGenDecayer::
    GetParticleDefinition(int ParPDGID, G4bool warn) const
{
  /// Return G4 particle definition for given TParticle

  // get particle definition from G4ParticleTable
  G4int pdgEncoding = ParPDGID;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particleDefinition = 0;
  if (pdgEncoding != 0)
    particleDefinition = particleTable->FindParticle(pdgEncoding);

  if (particleDefinition == 0 && warn)
  {
    std::cerr
        << "G4EvtGenDecayer: GetParticleDefinition: " << std::endl
        << "G4ParticleTable::FindParticle() for particle with PDG = "
        << pdgEncoding
        << " failed." << std::endl;
  }

  return particleDefinition;
}

bool G4EvtGenDecayer::
    IsG4Detectable(int ParPDGID, G4bool warn) const
{
  bool DetectorParticle = false;
  G4int pdgEncoding = ParPDGID;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particleDefinition = 0;
  if (pdgEncoding != 0) particleDefinition = particleTable->FindParticle(pdgEncoding);

  if (particleDefinition == 0 && warn)
  {
    DetectorParticle = false;
  }
  else
  {
    DetectorParticle = true;
  }

  return DetectorParticle;
}

//
// public methods
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DecayProducts* G4EvtGenDecayer::ImportDecayProducts(const G4Track& track)
{
  /// Import decay products
  // get particle momentum

  G4ThreeVector momentum = track.GetMomentum();
  G4double etot = track.GetDynamicParticle()->GetTotalEnergy();
 
  CLHEP::HepLorentzVector p;
  p[0] = momentum.x() / GeV;
  p[1] = momentum.y() / GeV;
  p[2] = momentum.z() / GeV;
  p[3] = etot / GeV;

  // get particle PDG
  // ask G4EvtGenDecayer to get PDG encoding
  // (in order to get PDG from extended TDatabasePDG
  // in case the standard PDG code is not defined)
  G4ParticleDefinition* particleDef = track.GetDefinition();
  G4int pdgEncoding = particleDef->GetPDGEncoding();
  //EvtGen Declaration//

  G4LorentzVector* ParMom = new G4LorentzVector;

  ParMom->setPx(p[0]);
  ParMom->setPy(p[1]);
  ParMom->setPz(p[2]);
  ParMom->setE(p[3]);

  std::vector<int> DecayPDGID;
  std::vector<int> DecayStatus;
  std::vector<G4LorentzVector> DecayMom;
  std::vector<G4LorentzVector> DecayVtx;

  Decay(pdgEncoding, ParMom);
  ImportParticles(DecayPDGID, DecayStatus, DecayMom, DecayVtx);

  int nTrk = DecayPDGID.size();

  // convert decay products EvtGen Products type
  // to G4DecayProducts
  G4DecayProducts* decayProducts = new G4DecayProducts(*(track.GetDynamicParticle()));

  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {

    const G4ParticleDefinition* particleDefinition = GetParticleDefinition(DecayPDGID[iTrk]);

    G4ThreeVector G4Mom = G4ThreeVector(DecayMom[iTrk].px() * GeV, DecayMom[iTrk].py() * GeV, DecayMom[iTrk].pz() * GeV);

//    int pdg = DecayPDGID[iTrk];
  //  int Status = DecayStatus[iTrk];

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
  }

  return decayProducts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EvtGenDecayer::PHEvtGenDecayer()
{
  //if(mEvtGen) return; // trust that mEvtGen is properly initialized by the user
  mOwner = true;

  //#ifdef EVTGEN_CPP11
  // Use the Mersenne-Twister generator (C++11 only)
  //	mEvtGenRandomEngine = new EvtMTRandomEngine(stime);
  //#else
  mEvtGenRandomEngine = new EvtSimpleRandomEngine();
  //#endif
  EvtRandom::setRandomEngine((EvtRandomEngine*) mEvtGenRandomEngine);
  EvtExternalGenList genList;
  EvtAbsRadCorr* radCorrEngine = genList.getPhotosModel();
  std::list<EvtDecayBase*> extraModels = genList.getListOfModels();

  // the hardcoded paths are temporary
  const string decay = 	  string(getenv("OFFLINE_MAIN"))  + "/share/EvtGen/DECAY.DEC"; // Using PDG 2019 reference as the input for now
  const string evt   =  string(getenv("OFFLINE_MAIN"))  + "/share/EvtGen/evt.pdl";

  std::ifstream in(decay);

  /*
	//	in.Open(decay);
	if (!in.good()) { 
	decay = "$(PHENIX)/"+decay;  decay = gSystem->ExpandPathName( decay.Data() );
	evt   = "$(PHENIX)/"+evt;    evt   = gSystem->ExpandPathName( evt.Data()   );
	}
	*/

  mEvtGen = new EvtGen(decay,
                       evt,
                       (EvtRandomEngine*) mEvtGenRandomEngine,
                       radCorrEngine,
                       &extraModels);

  //theEvent=new EvtHepMCEvent();
}

void G4EvtGenDecayer::Decay(int pdgId, G4LorentzVector* _p)
{
  // Clear the event from the last run
  ClearEvent();

  // Add the particle to the pythia stack
  AppendParticle(pdgId, _p);

  // Decay the particle
  mEvtGen->generateDecay(mParticle);
}

void G4EvtGenDecayer::ClearEvent()
{
  mVertex = new EvtVector4R(0, 0, 0, 0);
}

void G4EvtGenDecayer::AppendParticle(int pdg, G4LorentzVector* _p)
{
  // Append a particle to the stack to be decayed
  EvtVector4R p_init(_p->e(), _p->px(), _p->py(), _p->pz());
  EvtId parentID = EvtPDL::evtIdFromLundKC(pdg);
  //std::cout << "PDGID = " << pdg << "     p_init.Mag() = " << p_init.d3mag() << std::endl;
  mParticle = EvtParticleFactory::particleFactory(parentID, p_init);
  mParticle->setDiagonalSpinDensity();
}

int G4EvtGenDecayer::ImportParticles(std::vector<int>& DecayPDGID, std::vector<int>& DecayStatus, std::vector<G4LorentzVector>& DecayMom, std::vector<G4LorentzVector>& DecayVtx)
{
  // Save the decay products



	DecayPDGID.clear();
	DecayStatus.clear();
	DecayMom.clear();
	DecayVtx.clear();


	EvtHepMCEvent theEvent;


	theEvent.constructEvent(mParticle);



	HepMC3::GenEvent * evt = theEvent.getEvent(); //Directly use HepMC3 records
												  //HepMC::GenEvent * evt = ConvertHepMCGenEvent_3to2( *evt3); // Convert HepMC3 to HepMC2 Format


	G4LorentzVector FourMom;
	G4LorentzVector FourPos;


	
	int nparts = 0;
	int NUndefined = 0;

	std::vector<int> UndefinedParticles;
		

	for (auto vtx: evt->vertices()) {

		for (auto part: vtx->particles_in()) {


			int pdgidwithin = part->pdg_id();
			bool Detectable = IsG4Detectable(pdgidwithin);

			if(!Detectable){

	//			std::cout << "Within Vertex PDGID = " << pdgidwithin << "  is NOT Detectable by G4 " << std::endl;
				UndefinedParticles.push_back(pdgidwithin);

			}else{

				//std::cout << "Particle in the Chain: " << ptl0->GetName() << "   PDGID = " << pdg << "  IS Detectable by G4 " << std::endl;
		//		std::cout << "Within Vertex PDGID = " << pdgidwithin << "  is Detectable by G4 " << std::endl;	
			}

			

		}

	}

	NUndefined = UndefinedParticles.size();


	int TotalVertexProcessed = 0;
	int TotalVertexToProcess = 1;
	int TotalParticleToProcess = NUndefined;
	int TotalParticleProcessed = 0;

	for (auto vtx: evt->vertices()) {



			

		if(TotalVertexProcessed >= TotalVertexToProcess && TotalParticleProcessed >= TotalParticleToProcess){
		

			break;

		}
		bool IsUndefined = 0;

		for (auto part: vtx->particles_in()) {



			int pdgidwithin = part->pdg_id();

			
			for(int q = 0; q < NUndefined; q++ ){

				if(pdgidwithin == UndefinedParticles[q]){

			
					IsUndefined = 1;


				}
				
			}

			FourMom.setPx(part->momentum().px());
			FourMom.setPy(part->momentum().py());
			FourMom.setPz(part->momentum().pz());
			FourMom.setE(part->momentum().e());




			FourPos.setX(mVertex->get(1)+part->production_vertex()->position().x());
			FourPos.setY(mVertex->get(2)+part->production_vertex()->position().y());
			FourPos.setZ(mVertex->get(3)+part->production_vertex()->position().z());
			FourPos.setT(mVertex->get(0)+part->production_vertex()->position().t());



			nparts++;

		}


		for (auto part: vtx->particles_out()) {
	
			if(TotalVertexProcessed >= TotalVertexToProcess && IsUndefined == 0){
				std::cout << "Process Done, Skip" << std::endl;
				continue;
			}



			//Do Within the Particle Importation//
			int pdgidwithin = part->pdg_id();
			bool Detectable = IsG4Detectable(pdgidwithin);


			if(!Detectable){

				continue;
			}else{

			}



			DecayPDGID.push_back(part->pdg_id());			
			if(nparts == 0){
				DecayStatus.push_back(part->status()==1?11:-11);

			}else{
				DecayStatus.push_back(part->status()==1?91:-91);
			}


			FourMom.setPx(part->momentum().px());
			FourMom.setPy(part->momentum().py());
			FourMom.setPz(part->momentum().pz());
			FourMom.setE(part->momentum().e());


			DecayMom.push_back(FourMom);

			FourPos.setX(mVertex->get(1)+part->production_vertex()->position().x());
			FourPos.setY(mVertex->get(2)+part->production_vertex()->position().y());
			FourPos.setZ(mVertex->get(3)+part->production_vertex()->position().z());
			FourPos.setT(mVertex->get(0)+part->production_vertex()->position().t());




			DecayVtx.push_back(FourPos);

			nparts++;

		}
	
		TotalVertexProcessed = TotalVertexProcessed + 1;

	}




	return nparts;

}

void G4EvtGenDecayer::SetDecayTable(const string decayTable, bool useXml)
{
  mEvtGen->readUDecay(decayTable, useXml);
}

void G4EvtGenDecayer::SetVertex(G4LorentzVector* r)
{
  mVertex->set(r->t(), r->x(), r->y(), r->z());
}
