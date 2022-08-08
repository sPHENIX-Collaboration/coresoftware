#include "PHEvtGenDecayer.h"

#include <HepMC3/HepMC3.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenParticle.h>


#include "HepMC3/interfaces/HepMCCompatibility/include/HepMCCompatibility.h"

#include <EvtGenBase/EvtParticle.hh>
#include <EvtGenBase/EvtParticleFactory.hh>
#include <EvtGenBase/EvtPDL.hh>
#include <EvtGenBase/EvtRandom.hh>
//#include <EvtGenBase/EvtHepMCEvent.hh>
#include <EvtGenBase/EvtSimpleRandomEngine.hh>
#include <EvtGenBase/EvtRandomEngine.hh>

#include <EvtGenBase/EvtMTRandomEngine.hh>
#include <EvtGenBase/EvtAbsRadCorr.hh>
#include <EvtGenBase/EvtDecayBase.hh>
#include <EvtGenExternal/EvtExternalGenList.hh>
#include "../evtgen/EvtGenBase/EvtHepMCEvent.hh"
//#include <EvtGenBase/EvtHepMCEvent.hh>
//#include "EvtHepMCEvent.cpp" //Overwriting EvtGenPackage to Access HEPMC3 Events Records
//#include "EvtHepMCEvent.hh"

#include <list>
#include <TParticle.h>
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TSystem.h"
//#include "StMessMgr.h"
#include <iostream>
#include <fstream>




using namespace HepMC;
//using namespace HepMC3;
using namespace std;



PHEvtGenDecayer::PHEvtGenDecayer(EvtGen* evtGen): mEvtGen(evtGen), mEvtGenRandomEngine(NULL), mParticle(NULL), mVertex(NULL), mOwner(false), mDebug(0)
{
	if(mEvtGen) return; // trust that mEvtGen is properly initialized by the user
	mOwner = true;

//#ifdef EVTGEN_CPP11
	// Use the Mersenne-Twister generator (C++11 only)
//	mEvtGenRandomEngine = new EvtMTRandomEngine(stime);
//#else
	mEvtGenRandomEngine = new EvtSimpleRandomEngine();
//#endif
	EvtRandom::setRandomEngine((EvtRandomEngine*)mEvtGenRandomEngine);
	EvtExternalGenList genList;
	EvtAbsRadCorr* radCorrEngine = genList.getPhotosModel();
	std::list<EvtDecayBase*> extraModels = genList.getListOfModels();

	// the hardcoded paths are temporary
//	TString decay = "InputDECAYFiles/DECAY.DEC";
//	TString evt   = "InputDECAYFiles/evt.pdl";
 
	const string decay = "InputDECAYFiles/DECAY.DEC";
	const string evt   = "InputDECAYFiles/evt.pdl";
 


	std::ifstream in(decay);

/*
//	in.Open(decay);
	if (!in.good()) { 
		decay = "$(PHENIX)/"+decay;  decay = gSystem->ExpandPathName( decay.Data() );
		evt   = "$(PHENIX)/"+evt;    evt   = gSystem->ExpandPathName( evt.Data()   );
	}
*/
	
	mEvtGen = new EvtGen( decay, 
			evt,
			(EvtRandomEngine*)mEvtGenRandomEngine,
			radCorrEngine, 
			&extraModels);

	//theEvent=new EvtHepMCEvent();
}

PHEvtGenDecayer::~PHEvtGenDecayer()
{
	if(mOwner) {
		delete mEvtGen;
		delete mEvtGenRandomEngine;
	}
	mParticle->deleteTree();
	delete mVertex;
	//delete theEvent;
}

void PHEvtGenDecayer::Init()
{
	//cout << " Init Done" << endl;
}

void PHEvtGenDecayer::Decay(int pdgId, TLorentzVector* _p)
{
	//cout << "Decay pdgid=" << pdgId << endl;
	
	//cout << "SUCK BRI ZS" << endl;
	// Clear the event from the last run
	ClearEvent();

	//cout << "Pass Clear Event" << endl;

	// Add the particle to the pythia stack
	AppendParticle(pdgId, _p );
	//cout << "Pass AppendParticle" << endl;

	// Decay the particle
	mEvtGen->generateDecay(mParticle);

	//cout << "Pass Generated Decay" << endl;
	
}
void PHEvtGenDecayer::ClearEvent()
{

	//cout << "Now Clear" << endl;

	//mVertex= 
	//(0,0,0,0);// Reset the postion for mother particle,Default is(0,0,0,0)
	//mVertex->SetXYZT(0,0,0,0);PHEvtGenDecayer.cxx
	mVertex=new EvtVector4R(0,0,0,0);
	//cout << "After mVertex->mass() = " << mVertex->mass() << endl;


	////cout << "mVertex " << mVertex->M() << endl;

	if(mParticle) mParticle->deleteTree(); // this deletes the daughter and mParticle itself (the object commits suicide)
}

void PHEvtGenDecayer::AppendParticle(Int_t pdg, TLorentzVector* _p)
{
	// Append a particle to the stack to be decayed
	EvtVector4R p_init(_p->E(), _p->Px(), _p->Py(), _p->Pz());
	EvtId parentID = EvtPDL::evtIdFromLundKC(pdg);
	//cout << "parentID = " << parentID << "     p_init.Mag() = " << p_init.d3mag() << endl;
	mParticle = EvtParticleFactory::particleFactory(parentID, p_init);
	mParticle->setDiagonalSpinDensity();
}

Int_t PHEvtGenDecayer::ImportParticles(TClonesArray* _array)
{
	// Save the decay products
	
	cout<< "Pass ImportParticles = 1" <<  endl;

	assert(_array);
	//cout<< "Pass ImportParticles = 2" <<  endl;
	
	TClonesArray &array = *_array;
	//cout<< "Pass ImportParticles = 3" <<  endl;

	array.Clear();
	//cout<< "Pass ImportParticles = 4" <<  endl;

	EvtHepMCEvent theEvent;

	//cout<< "Pass ImportParticles = 5" <<  endl;
	theEvent.constructEvent(mParticle);
	// Print list of EvtGen lines on debug
/*
	if (mDebug) {
		theEvent.getEvent()->print(std::cout);
	}
*/
	//cout<< "Pass ImportParticles = 6" <<  endl;

	//cout << "DEBUG --- NEW EVTGEN VERSION" << endl;


//	theEvent.getEvent()->print(std::cout);

	//Convert HEPMC3 to HEPMC2


	HepMC3::GenEvent * evt3 = theEvent.getEvent();
	cout << "Now Convert" << endl;

	HepMC::GenEvent * evt = ConvertHepMCGenEvent_3to2( *evt3);

	cout << "Now Try" << endl;

//	evt->print(std::cout);


	Int_t nparts = 0;
	Int_t particle_barcode[100];
	//cout<< "Pass ImportParticles = 7" <<  endl;
	////cout << " theEvent.getEvent()->size() = " << theEvent.getEvent()->size() << endl;

	for ( GenEvent::vertex_const_iterator vtx = evt->vertices_begin();
			vtx != evt->vertices_end(); ++vtx ) {
		    //cout<< "Pass ImportParticles = 7.1" <<  endl;
		
		if(vtx==evt->vertices_begin()) {
			//cout << "Pass ImportParticles = 7.14" << endl;
			//cout << "Pass ImportParticles = 7.16" << endl;
		//	GenVertex::particles_in_const_iterator part2 = (*vtx)->particles_in_const_begin();
			//cout << "(*vtx)->particles_in_size()" << (*vtx)->particles_in_size() << endl;
		//	assert(part2);
		//	if(!part2) //cout << "(*vtx)->particles_in_const_begin() DOES NOT Exist, need to check why" << endl;
			//cout << "Pass ImportParticles = 7.18" << endl;

			//OLD HEPMC Codes
			for (GenVertex::particles_in_const_iterator part = (*vtx)->particles_in_const_begin();part != (*vtx)->particles_in_const_end(); part++ ) {
				//cout << "nparts = " << nparts << endl;
				particle_barcode[nparts]=(*part)->barcode();
				nparts++;
			}
	    //cout<< "Pass ImportParticles = 7.2" <<  endl;
		
		}
		for ( GenVertex::particles_out_const_iterator part = (*vtx)->particles_out_const_begin();part != (*vtx)->particles_out_const_end(); part++ ) {
			particle_barcode[nparts]=(*part)->barcode();
			nparts++;
	    //cout<< "Pass ImportParticles = 7.3" <<  endl;
		}
	
/*		//New HEPMC3 Codes

		for (auto part: (*vtx)->particles_in()) {
				//cout << "nparts = " << nparts << endl;
				particle_barcode[nparts]=(*part)->barcode();
				nparts++;
			}
	    //cout<< "Pass ImportParticles = 7.2" <<  endl;
		
		}
		for ( auto part: (*vtx)->particles_out() ) {
			particle_barcode[nparts]=(*part)->barcode();
			nparts++;
	    //cout<< "Pass ImportParticles = 7.3" <<  endl;
		}
*/

	}
	//cout<< "Pass ImportParticles = 8" <<  endl;

	GenParticle* p;
	for(Int_t np=0;np<nparts;np++) {
		p=evt->barcode_to_particle(particle_barcode[np]);
		int firstmother=0, firstdaughter=0, lastdaughter=0;
		int firstmother_barcode=-1, firstdaughter_barcode=-1;
		if(np>0) firstmother_barcode=(*(p->production_vertex()->particles_in_const_begin()))->barcode();
		if(p->status()==2&&p->end_vertex()) {
			firstdaughter_barcode=(*(p->end_vertex()->particles_out_const_begin()))->barcode();
		}
		for(int i=0;i<nparts;i++) {
			if(particle_barcode[i]==firstmother_barcode) firstmother=i;
			if(particle_barcode[i]==firstdaughter_barcode) {
				firstdaughter=i;
				lastdaughter=i;
				if(p->status()==2&&p->end_vertex()) lastdaughter=i+p->end_vertex()->particles_out_size()-1;
			}
		}
		if(np==0) {
			new(array[np]) TParticle(p->pdg_id(), p->status()==1?11:-11,
					firstmother, 0,
					firstdaughter, lastdaughter,
					p->momentum().px(),
					p->momentum().py(),
					p->momentum().pz(),
					p->momentum().e(),
					mVertex->get(1),
					mVertex->get(2),
					mVertex->get(3),
					mVertex->get(0));
		} else {
			new(array[np]) TParticle(p->pdg_id(), p->status()==1?91:-91,
					firstmother, 0,
					firstdaughter, lastdaughter,
					p->momentum().px(),
					p->momentum().py(),
					p->momentum().pz(),
					p->momentum().e(),
					mVertex->get(1)+p->production_vertex()->position().x(),
					mVertex->get(2)+p->production_vertex()->position().y(),
					mVertex->get(3)+p->production_vertex()->position().z(),
					mVertex->get(0)+p->production_vertex()->position().t());
		}
		
	//	cout << "--------------------------  Decay Particle Info  -----------------------------------" << endl;
	//	cout<<"HEMPMC - PDG ID =  "<<p->pdg_id() <<  "  px =  "  << p->momentum().x() << "  py =  " << p->momentum().y() << "  pz =  " << p->momentum().z() <<  "   E  = "<< p->momentum().t() << "   Mass = "  <<  sqrt(p->momentum().t() *  p->momentum().t() -   p->momentum().x()  *  p->momentum().x()  -  p->momentum().y() *  p->momentum().y() -  p->momentum().z()  *  p->momentum().z())   <<  endl;

		if (mDebug) {
			cout<<np<<" "<<((TParticle*)array[np])->GetStatusCode()
				<<" "<<((TParticle*)array[np])->GetFirstDaughter()
				<<" "<<((TParticle*)array[np])->GetLastDaughter()
				<<" "<<((TParticle*)array[np])->Vx()
				<<" "<<((TParticle*)array[np])->T()<<" ";
			((TParticle*)array[np])->Print();
		}
		//		((TParticle*)array[np])->Print() ;
		
	}
	//cout<< "Pass ImportParticles = 10" <<  endl;

	return nparts;
}


void PHEvtGenDecayer::SetForceDecay(Int_t type)
{
	cout << "PHEvtGenDecayer::SetForceDecay method is not implemented in this class" << type <<endl;
}

void PHEvtGenDecayer::ForceDecay()
{
	//cout << "PHEvtGenDecayer::ForceDecay method is not implemented in this class" <<endl;
}

Float_t PHEvtGenDecayer::GetPartialBranchingRatio(Int_t ipart)
{
	cout << "PHEvtGenDecayer::GetPartialBranchingRatio method is not implemented in this class" << ipart <<endl;
	return 1.0;
}

void PHEvtGenDecayer::ReadDecayTable()
{
	//cout << "PHEvtGenDecayer::ReadDecayTable method is not implemented in this class" <<endl;
}

Float_t PHEvtGenDecayer::GetLifetime(Int_t pdg)
{
	return (EvtPDL::getctau(EvtPDL::evtIdFromLundKC(pdg)) * 3.3333e-12) ;
}

void PHEvtGenDecayer::SetDecayTable(const string decayTable, bool useXml)
{
	mEvtGen->readUDecay(decayTable,useXml);
}

void PHEvtGenDecayer::SetVertex(TLorentzVector* r)
{
	mVertex->set(r->T(),r->X(),r->Y(),r->Z());
}
