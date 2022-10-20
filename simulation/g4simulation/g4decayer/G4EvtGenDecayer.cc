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
#include <stack>    // for operator<<
#include <queue>    // for operator<<

#include <phool/PHRandomSeed.h>  //Introduce Random Number with Dedicated PHENIX Random Seed


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EvtGenDecayer::G4EvtGenDecayer()
	: G4VExtDecayer("G4EvtGenDecayer")
	  , fVerboseLevel(0)
{

	//  mOwner = true;

	//#ifdef EVTGEN_CPP11
	// Use the Mersenne-Twister generator (C++11 only)
	//	mEvtGenRandomEngine = new EvtMTRandomEngine(stime);
	//#else
	mEvtGenRandomEngine = new PHEvtGenRandomEngine();
	//#endif
	EvtRandom::setRandomEngine((EvtRandomEngine*) mEvtGenRandomEngine);

	radCorrEngine = genList.getPhotosModel();

	//	std::list<EvtDecayBase*> extraModels = genList.getListOfModels();

	extraModels = genList.getListOfModels();
	//auto extraModels = genList.getListOfModels();


	// the hardcoded paths are temporary
	const string decay = 	  string(getenv("OFFLINE_MAIN"))  + "/share/EvtGen/DECAY.DEC"; // Using PDG 2019 reference as the input for now
	const string evt   =  string(getenv("OFFLINE_MAIN"))  + "/share/EvtGen/evt.pdl";

	std::ifstream in(decay);

	/*
	   mEvtGen = new EvtGen(decay,
	   evt,
	   (EvtRandomEngine*) mEvtGenRandomEngine,
	   radCorrEngine,
	   &extraModels);
	   */
	mEvtGen = new EvtGen(decay,  evt, (EvtRandomEngine*) mEvtGenRandomEngine,radCorrEngine,&extraModels);
	extraModels.clear();

	// delete mEvtGen;	

	//bool WilluseXml =false;
	//SetDecayTable("EvtGenDecayFiles/Bc.D0STARDs.DEC",WilluseXml); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EvtGenDecayer::~G4EvtGenDecayer()
{

	delete mEvtGen;		
	//delete evt;		
	delete radCorrEngine;
	delete mEvtGenRandomEngine;
	delete mParticle;
	//mParticle->deleteTree();

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

//	std::cout << "Indeed Using EvtGen Decayer - Final" << std::endl;

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

	std::unique_ptr <G4LorentzVector> ParMom (new G4LorentzVector);

	ParMom->setPx(p[0]);
	ParMom->setPy(p[1]);
	ParMom->setPz(p[2]);
	ParMom->setE(p[3]);


	//	int DecayPDGID;
	//	int DecayStatus;	
	//int MotherPDG;


	//  Decay(pdgEncoding, ParMom);
	//  ImportParticles(DecayPDGID, DecayStatus, DecayMom, DecayVtx);

	//   EvtVector4R p_init(ParMom->e(), ParMom->px(), ParMom->py(), ParMom->pz());
	EvtVector4R p_init(p[3],p[0],p[1],p[2]);

	EvtId parentID = EvtPDL::evtIdFromLundKC(pdgEncoding);
	//std::unique_ptr <EvtParticle*> mParticle = std::unique_ptr<EvtParticle*>(new EvtParticle());
	mParticle = EvtParticleFactory::particleFactory(parentID, p_init);
	//   mParticle->setDiagonalSpinDensity();


	// Decay the particle
	mEvtGen->generateDecay(mParticle);


	EvtHepMCEvent theEvent;


	theEvent.constructEvent(mParticle);



	HepMC3::GenEvent * evt = theEvent.getEvent(); //Directly use HepMC3 records
												  //HepMC::GenEvent * evt = ConvertHepMCGenEvent_3to2( *evt3); // Convert HepMC3 to HepMC2 Format


	G4LorentzVector FourMom;
	G4LorentzVector FourPos;
	G4LorentzVector MotherPos;


	//std::vector<int> CheckWhat;
	std::queue<int> stack;

	G4DecayProducts* decayProducts = new G4DecayProducts(*(track.GetDynamicParticle()));

	int MotherToProcess = 0; 
	bool MotherMatched;

	for (auto vtx: evt->vertices()) {

		MotherMatched = false;
		if(!stack.empty()) MotherToProcess = stack.front();
		if(stack.empty()) MotherMatched = true;


		for (auto part: vtx->particles_in()) {

			int pdgidwithin = part->pdg_id();
			bool Detectable = IsG4Detectable(pdgidwithin);


			//std::cout << "ParIn: pdgidwithin = " << pdgidwithin << std::endl;
		 
			if(!stack.empty()){
				if(pdgidwithin == MotherToProcess) MotherMatched = true;
				if(pdgidwithin != MotherToProcess){
		
				//	std::cout << "Mother PDG to process  does NOT match with the incoming " << std::endl;
				}
			}
		
			if(Detectable && stack.empty()) stack.push(pdgidwithin);

		}


		if(!MotherMatched){
			//	std::cout << "Mother DO NOT Match... Skip the Vertex!" << std::endl;
				//std::cout << "These daughters are from Mother: " << MotherPDG << "   They are NOT the remaining mother " << MotherToProcess << std::endl;
				continue;
		}

		if(stack.empty()) continue;

		//int MotherToProcess = stack.top();
		//int MotherToProcess = stack.front();

		stack.pop();

		for (auto part: vtx->particles_out()) {


			//		stack.top();
	 
			int pdgidwithin = part->pdg_id();
			bool Detectable = IsG4Detectable(pdgidwithin);
			//std::cout << "ParOut: pdgidwithin = " << pdgidwithin << std::endl;



			bool SameVtxPos = 1;


			if(Detectable){

				for (auto ParentPart: part->parents()) {



					//		MotherPDG = ParentPart->pdg_id();



					//std::cout << "PDGID = " << part->pdg_id()  << "   Parengt Size = " <<  part->parents().size() <<  "   part->production_vertex()->position().x()  = " << part->production_vertex()->position().x()   << "   MotherPDGID = " << MotherPDG  << "    Mother X Pos = " << ParentPart->end_vertex()->position().x() << std::endl;


					if(part->production_vertex()->position().x() != ParentPart->end_vertex()->position().x() || part->production_vertex()->position().y() != ParentPart->end_vertex()->position().y() || part->production_vertex()->position().z() != ParentPart->end_vertex()->position().z() || part->production_vertex()->position().t() != ParentPart->end_vertex()->position().t() ){
						SameVtxPos = 0;
					}

				}


				const G4ParticleDefinition* particleDefinition = GetParticleDefinition(pdgidwithin);
				G4ThreeVector G4Mom = G4ThreeVector(part->momentum().px() * GeV, part->momentum().py() * GeV, part->momentum().pz() * GeV);
				G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDefinition, G4Mom);
				if (dynamicParticle)
				{
					if (fVerboseLevel > 0)
					{
						std::cout << "  G4 particle name: "
							<< dynamicParticle->GetDefinition()->GetParticleName()
							<< std::endl;
					}


					if(SameVtxPos){

						//		std::cout << "Particle: " << pdgidwithin << " is registered!" << std::endl;
						decayProducts->PushProducts(dynamicParticle); 

					}
					if(!SameVtxPos){

						std::cout << "Issue Found: in this vertex, particles pushed to GEANT have different production positions, need to check and understand why!!!" << std::endl;
						std::cout << "Here is the Event Content" << std::endl;
						const HepMC3::GenEvent& CheckEvent =  *evt;
						HepMC3::Print::content(CheckEvent);
						//delete CheckEvent;
					}
				}


			}


//			if(Detectable) CheckWhat.push_back(pdgidwithin);

		//	std::cout << "pdgidwithin = " << pdgidwithin << "   SIZE of the STACK NOW = " << stack.size() << std::endl;
			if(!Detectable){
				//		std::cout << "pdg = " << pdgidwithin << "  Is NOT detectible! Moving on to the next tree" << std::endl;
				stack.push(pdgidwithin);
			}
		//	std::cout << "pdgidwithin = " << pdgidwithin << "   SIZE of the STACK POST = " << stack.size() << std::endl;

		}

		if(stack.empty()){

		//	std::cout << "Vertex Assembly DONE" << std::endl;
			break;

		}

	}

/*
	int SizeVec = CheckWhat.size();

	std::cout << "------------------------------------------------------------------------------------------"<< std::endl;

	for(int q = 0; q < SizeVec; q++){

		std::cout << "q = " << q << "   What we have pushed:  " << CheckWhat[q] << std::endl;

	}



	std::cout << "------------------------------------------------------------------------------------------"<< std::endl;

	std::cout << "G4 Decay Products: " << decayProducts->entries() << std::endl;
*/

	evt->clear();

//	delete mEvtGen;
//	delete mParticle;
//	delete (*theEvent);	
//	delete evt;
	mParticle->deleteTree();
	return decayProducts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EvtGenDecayer::SetDecayTable(const string decayTable, bool useXml)
{
	mEvtGen->readUDecay(decayTable, useXml);
}


