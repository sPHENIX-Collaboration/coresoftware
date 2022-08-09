/**
  \class PHEvtGenDecayer

  \brief PH wrapper for EvtGen Decayer

Authors: Xiaozhi Bai (xiaozhi@uic.edu),
Mustafa Mustafa (mmustafa@lbl.gov)
Zhaozhong Shi for Modification
*/

#ifndef PHEVTGENDECAYER_H
#define PHEVTGENDECAYER_H



#include <cstddef>

#include <TVirtualMCDecayer.h>
#include <TString.h>
#include <EvtGen/EvtGen.hh>
#include <TLorentzVector.h>
#include <EvtGenBase/EvtSimpleRandomEngine.hh>
#include <EvtGenBase/EvtRandomEngine.hh>
#include <EvtGenExternal/EvtExternalGenList.hh>
#include <EvtGenBase/EvtRandom.hh>

#include <HepMC3/HepMC3.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenParticle.h>


class TLorentzVector;
class TClonesArray;
class EvtStdlibRandomEngine;
class EvtParticle;
class EvtRandomEngine;

using namespace std;

class PHEvtGenDecayer : public TVirtualMCDecayer
{
	public:
		PHEvtGenDecayer(EvtGen* evtGen = NULL);
		~PHEvtGenDecayer();
		void Init();
		void Decay(Int_t pdgId, TLorentzVector* p);
		Int_t ImportParticles(TClonesArray* particles);
		void SetForceDecay(Int_t type);
		void ForceDecay();
		void ReadDecayTable();
		Float_t GetPartialBranchingRatio(Int_t ipart);
		Float_t GetLifetime(Int_t pdgid);

		void setVertex(TLorentzVector* r);
		void setDecayTable(const string decayTable, bool useXml);
		void ClearEvent();
		void AppendParticle(Int_t pdg, TLorentzVector* _p);	
		void SetDecayTable(const string decayTable, bool useXml);
		void SetVertex(TLorentzVector* r);	


	private:
		// EvtStdlibRandomEngine* mEvtGenRandomEngine;
		EvtGen* mEvtGen = NULL;	
		EvtRandomEngine  * mEvtGenRandomEngine = NULL;
		EvtParticle* mParticle;
	//	TLorentzVector * mVertex = new TLorentzVector;  
		EvtVector4R * mVertex;
		bool mOwner;
		bool mDebug = false; 
		TLorentzVector * r_init = new TLorentzVector;   
};



inline void PHEvtGenDecayer::setDecayTable(const string decayTable, bool useXml) { mEvtGen->readUDecay(decayTable, useXml); }
inline void PHEvtGenDecayer::setVertex(TLorentzVector* r) {r_init->SetXYZT(r->X(),r->Y(),r->Z(),r->T()); }

#endif  /* PHEVTGENDECAYER_H */

