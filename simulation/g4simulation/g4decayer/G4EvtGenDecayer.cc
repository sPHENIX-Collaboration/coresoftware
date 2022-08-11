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
#include "TFile.h"
#include "TTree.h"


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
#include "../Share.hh"

TFile * CheckDecay;
TTree * DecayTree;
int EventID;
int PDGID;
float px;
float py;
float pz;
float pt;
float E;
float vx;
float vy;
float vz;

bool useXml;

//const EDecayType G4EvtGenDecayer::fgkDefaultDecayType = kAll;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EvtGenDecayer::G4EvtGenDecayer()
  : G4VExtDecayer("G4EvtGenDecayer")
  , fVerboseLevel(0)
{

	useXml = false;




	EventID = 0;


	EvtRandomEngine* eng = 0;
	eng = new EvtSimpleRandomEngine();
	EvtRandom::setRandomEngine((EvtRandomEngine*)eng);
	EvtAbsRadCorr* radCorrEngine = 0;
	std::list<EvtDecayBase*> extraModels;
	std::cout << "Pass 8" << std::endl;

	EvtExternalGenList genList;
	std::cout << "Pass 9" << std::endl;
	radCorrEngine = genList.getPhotosModel();


	extraModels = genList.getListOfModels();

	std::string Decay_DEC = "InputDECAYFiles/DECAY.DEC";
	std::string Evt_pdl = "InputDECAYFiles/evt.pdl";

	myGenerator = new EvtGen(Decay_DEC, Evt_pdl, (EvtRandomEngine*)eng, radCorrEngine, &extraModels);
	myEvtGenDecayer = new PHEvtGenDecayer(myGenerator);
	std::cout << "Setting Decay Table" << std::endl;

	myEvtGenDecayer->SetDecayTable("InputDECAYFiles/Bs.KKPiPi.DEC",useXml);

	
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
GetParticleDefinition(const TParticle* particle, G4bool warn) const
{
	/// Return G4 particle definition for given TParticle

	// get particle definition from G4ParticleTable
	G4int pdgEncoding = particle->GetPdgCode();
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particleDefinition = 0;
	if (pdgEncoding != 0)
		particleDefinition = particleTable->FindParticle(pdgEncoding);
/*
	std::cout << "Just Check" << std::endl;

	for(int i = 0; i < 500; i ++){
		particleDefinition = particleTable->FindParticle(i);

		if (particleDefinition == 0 && warn){

			std::cout << "PDGID " << i << "   is NOT in the " << std::endl;

		}
		else{

			std::cout << "PDGID " << i << "   IS inside the "  << std::endl;

		}
		


	}
*/	


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
IsG4Detectable(const TParticle* particle, G4bool warn) const
{
	bool DetectorParticle = false;
	G4int pdgEncoding = particle->GetPdgCode();
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particleDefinition = 0;
	if (pdgEncoding != 0) particleDefinition = particleTable->FindParticle(pdgEncoding);

	if (particleDefinition == 0 && warn){
		DetectorParticle = false;
	}
	else{
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
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
  // get particle momentum

    std::cout << "OK New" << std::endl;	

	G4ThreeVector momentum = track.GetMomentum();
	G4double etot = track.GetDynamicParticle()->GetTotalEnergy();
	;
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
	ParentID = pdgEncoding;
    ParentMass = p.m();
    ParentE = p[3];
    ParentPx = p[0];
    ParentPy = p[1];
    ParentPz = p[2];


	


	TClonesArray daughters("TParticle", 10);
	TLorentzVector* ParMom= new TLorentzVector;
	ParMom->SetPxPyPzE(p[0],p[1],p[2],p[3]);

	myEvtGenDecayer->Decay(pdgEncoding, ParMom);

	myEvtGenDecayer->ImportParticles(&daughters);



	int nTrk = daughters.GetEntriesFast();






	// convert decay products EvtGen Products type
	// to G4DecayProducts
    TLorentzVector * FourMomentum = new TLorentzVector;
    TLorentzVector * TotalFourMom = new TLorentzVector;
    TotalFourMom->SetXYZM(0,0,0,0);

    float ParMass = 0;
    TotalParticle = 0;

	G4DecayProducts* decayProducts = new G4DecayProducts(*(track.GetDynamicParticle()));

	
	for (int iTrk = 0; iTrk < nTrk; ++iTrk)
	{

		TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
		int pdg = ptl0->GetPdgCode();
		int Status = ptl0->GetStatusCode();
		//int MomID = ptl0->GetFirstMother();
		ptl0->Momentum(FourMom);
	

		
		TVector3 rcPosSig(ptl0->Vx() * cm, ptl0->Vy() * cm, ptl0->Vz() * cm);
		TVector3 Mom = FourMom.Vect(); 

		bool Detectable = IsG4Detectable(ptl0);

		if(!Detectable){
		
	//		std::cout << "Particle in the Chain: " << ptl0->GetName() << "   PDGID = " << pdg << "  is NOT Detectable by G4 " << std::endl;
			continue;
		}else{

	//		std::cout << "Particle in the Chain: " << ptl0->GetName() << "   PDGID = " << pdg << "  IS Detectable by G4 " << std::endl;

		}
		

		const G4ParticleDefinition* particleDefinition = GetParticleDefinition(ptl0);

		G4ThreeVector G4Mom = G4ThreeVector(Mom[0] * GeV,Mom[1] * GeV,Mom[2] * GeV);
			
		if(abs(pdg) == 130 || abs(pdg) == 310 || abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 22 || abs(pdg) == 11 || abs(pdg) == 2212 || abs(pdg) == 13 ||  abs(pdg) == 12 ||  abs(pdg) == 14 || abs(pdg) == 16){

			if(Status != 91) std::cout<<"Now PDGID: "  <<  pdg  << "   Status = " << Status << endl;

		}

		if(Status != -11 && Status != -91  && (abs(pdg) == 130 || abs(pdg) == 310 || abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 22 || abs(pdg) == 11 || abs(pdg) == 2212 || abs(pdg) == 13 ||  abs(pdg) == 12 ||  abs(pdg) == 14 || abs(pdg) == 16 || abs(pdg) == 2112)){

		TotalPx = TotalPx + Mom[0];
		TotalPy = TotalPy + Mom[1];
		TotalPz = TotalPz + Mom[2];



		TotalParticle = TotalParticle + 1;
	

		if(pdg == 211) TotalPiP = TotalPiP + 1;
		if(pdg == 321) TotalKP = TotalKP + 1;
		if(pdg == 2212) TotalPP = TotalPP + 1;
		if(pdg == 11) TotalElecM = TotalElecM + 1;
		if(pdg == 13) TotalMuM = TotalMuM + 1;
		if(pdg == -211) TotalPiM = TotalPiM + 1;
		if(pdg == -321) TotalKM = TotalKM + 1;
		if(pdg == -2212) TotalPM = TotalPM + 1;
	
		if(pdg == -11) TotalElecP = TotalElecP + 1;
		if(pdg == -13) TotalMuP = TotalMuP + 1;
		if(pdg == 130) TotalKs = TotalKs + 1;
		if(pdg == 310) TotalKL = TotalKL + 1;
		if(pdg == 22) TotalGamma = TotalGamma + 1;
		if(abs(pdg) == 12 || abs(pdg) == 14 ||  abs(pdg) == 16) TotalNu = TotalNu + 1;
		if(abs(pdg) == 2112) TotalN = TotalN + 1;


		if(abs(pdg) == 211) ParMass = 0.13957039;
		if(abs(pdg) == 321) ParMass = 0.493677;
		if(abs(pdg) == 2212) ParMass = 0.9382720881;
		if(abs(pdg) == 11) ParMass =  0.000511;
		if(abs(pdg) == 13) ParMass = 0.1056583755;
		if(abs(pdg) == 22) ParMass = 0.0;
		if(abs(pdg) == 12 || abs(pdg) == 14 ||  abs(pdg) == 16) ParMass =0;
		if(abs(pdg) == 130) ParMass = 0.497611;
		if(abs(pdg) == 310) ParMass = 0.497611;
		if(abs(pdg) == 2112) ParMass = 0.939565420;

		if(ParMass == -1) std::cout << "Bad PDGID " << pdg << "     Status:  " <<  Status  <<  std::endl;
		FourMomentum->SetXYZM(Mom[0],Mom[1],Mom[2],ParMass);


		ParMass = -1;

//		if(abs(pdg)!= 22) *TotalFourMom = *TotalFourMom + *FourMom;
		*TotalFourMom = *TotalFourMom + *FourMomentum;
		TotalE = FourMomentum->E() + TotalE;

		}

		G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDefinition, G4Mom);
		if (dynamicParticle)
		{
			if (fVerboseLevel > 0)
			{
				std::cout << "  G4 particle name: "
					<< dynamicParticle->GetDefinition()->GetParticleName()
					<< std::endl;
			}

			// add dynamicParticle to decayProducts
	//		if(id != pdgEncoding) decayProducts->PushProducts(dynamicParticle); 
			//if(abs(id) != 531 && abs(id) != 443 && abs(id) != 333) decayProducts->PushProducts(dynamicParticle); 
//			if(abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 11 || abs(pdg) == 2212 || abs(pdg) == 13 ) decayProducts->PushProducts(dynamicParticle); 
	//	 	if(Status != -11) decayProducts->PushProducts(dynamicParticle); 

		//	if(abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 11 || abs(pdg) == 2212 || abs(pdg) == 13  || abs(pdg) == 22) decayProducts->PushProducts(dynamicParticle); 
			if( Status != -11 && Status != -91 && (abs(pdg) == 130 || abs(pdg) == 310 || abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 22 || abs(pdg) == 11 || abs(pdg) == 2212 || abs(pdg) == 13 ||  abs(pdg) == 12 ||  abs(pdg) == 14 || abs(pdg) == 16 || abs(pdg) == 2112))  decayProducts->PushProducts(dynamicParticle); 

		}
	}
    TotalMass = TotalFourMom->M();

	return decayProducts;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

