#include "QAG4Decayer.h"

#include "QAHistManagerDef.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <decayfinder/DecayFinder.h>
#include <decayfinder/DecayFinderContainerBase.h>    // for DecayFinderContainerBase::Iter
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <CLHEP/Vector/LorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <TDatabasePDG.h>

#include <g4eval/SvtxEvalStack.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLatex.h>

#include <TVector3.h>
#include <TLorentzVector.h>

const int NHFQA = 13;
int QAVtxPDGID[NHFQA] = {411,421,431,4122,511,521,531,553,443,100443,200443,100553,200553};



/*
 *  QA module to check decay branching ratio, decay lifetime, and momentum conservation for inclusive heavy flavor hadron decay, which is handle by EvtGen as default
 *  Authors: Zhaozhong Shi
 *  Date: November 2022
 */

QAG4Decayer::QAG4Decayer(const std::string &name) 
	  : SubsysReco(name)
	  , m_write_nTuple(false) 
	  , m_write_QAHists(true)
	  , m_SaveFiles(false)		  
{

	std::cout << "New QA Clean Up for Valgrind Test" << std::endl;

}

QAG4Decayer::~QAG4Decayer(){

}

int QAG4Decayer::Init(PHCompositeNode *topNode){

	Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
	assert(hm);

	TH1 *h(nullptr);

	

	h = new TH1F("SVtoPVDistance","",100,0,0.2);
	hm->registerHisto(h);
	

	h = new TH1F("TotalVtxStat","",20,-0.5,19.5);
	hm->registerHisto(h);

	h = new TH1F("BadVtxStat","",20,-0.5,19.5);
	hm->registerHisto(h);
	
	h = new TH1F("BadVtxPercent","",20,-0.5,19.5);
	hm->registerHisto(h);
	


	h = new TH1F("QAPx","",100,-1,1);
	hm->registerHisto(h);

	h = new TH1F("QAPy","",100,-1,1);
	hm->registerHisto(h);
	
	h = new TH1F("QAPz","",100,-1,1);
	hm->registerHisto(h);

	h = new TH1F("QAE","",100,-1,1);
	hm->registerHisto(h);

	h = new TH1F("QACosTheta","",120,-1.2,1.2);
	hm->registerHisto(h);

	h = new TH1F("MassHis","",100,1.6,2.1);
	hm->registerHisto(h);
	
	h = new TH1F("NPartHis","",10,-0.5,9.5);
	hm->registerHisto(h);


	for(int i = 0; i < NHFQA; i++){

		h = new TH1F(Form("BR1DHis_%d",i),"",10,-0.5,9.5);
		hm->registerHisto(h);

		h = new TH1F(Form("ProperLifeTime_%d",i),"",100,0,0.05);
		hm->registerHisto(h);

	}

	h = new TH1F("HFHadronStat","",14,-0.5,13.5);
	hm->registerHisto(h);


	h = new TH1F("OtherHis","",1,-0.5,0.5);
	hm->registerHisto(h);



	NParticles = 0;

	gStyle->SetOptStat(0);


	assert(topNode);

		
	EvtID = 0;
	LifeTime = 0;
	NParticles = 0;
	PVDaughtersPDGID.clear();


	if(m_SaveFiles){
		fout = new TFile("MyQAFile.root","RECREATE");
		QATree = new TTree("QATree","QATree");
		QATree->Branch("EvtID",&EvtID,"EvtID/I");
		QATree->Branch("NParticles",&NParticles,"NParticles/I");	
		QATree->Branch("PVDaughtersPDGID",&PVDaughtersPDGID);		

		QATree->Branch("LifeTime",&LifeTime,"LifeTime/F");

	}

	return 0;

}

int QAG4Decayer::process_event(PHCompositeNode *topNode)
{




	NParticles = 0;

	Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
	assert(hm);

	TH1F * SVtoPVDistance = dynamic_cast<TH1F *>(hm->getHisto("SVtoPVDistance"));
	assert(SVtoPVDistance);



	TH1F * TotalVtxStat = dynamic_cast<TH1F *>(hm->getHisto("TotalVtxStat"));
	assert(TotalVtxStat);

	TH1F * BadVtxStat = dynamic_cast<TH1F *>(hm->getHisto("BadVtxStat"));
	assert(BadVtxStat);


	TH1F * BadVtxPercent = dynamic_cast<TH1F *>(hm->getHisto("BadVtxPercent"));
	assert(BadVtxPercent);


	TH1F * QAPx = dynamic_cast<TH1F *>(hm->getHisto("QAPx"));
	assert(QAPx);


	TH1F * QAPy = dynamic_cast<TH1F *>(hm->getHisto("QAPy"));
	assert(QAPy);

	TH1F * QAPz = dynamic_cast<TH1F *>(hm->getHisto("QAPz"));
	assert(QAPz);

	TH1F * QAE = dynamic_cast<TH1F *>(hm->getHisto("QAE"));
	assert(QAE);


	TH1F * QACosTheta = dynamic_cast<TH1F *>(hm->getHisto("QACosTheta"));
	assert(QACosTheta);


	TH1F * MassHis = dynamic_cast<TH1F *>(hm->getHisto("MassHis"));
	assert(MassHis);



	TH1F * BR1DHis[NHFQA];
	TH1F * ProperLifeTime[NHFQA];

	for(int i = 0; i < NHFQA; i++){

		BR1DHis[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("BR1DHis_%d",i)));
		assert(BR1DHis[i]);


		ProperLifeTime[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("ProperLifeTime_%d",i)));
		assert(ProperLifeTime[i]);

	}

	TH1F * NPartHis = dynamic_cast<TH1F *>(hm->getHisto("NPartHis"));
	assert(NPartHis);
	
	TH1F * HFHadronStat = dynamic_cast<TH1F *>(hm->getHisto("HFHadronStat"));
	assert(HFHadronStat);

 

	if (!_svtxevalstack)
	{
		_svtxevalstack = new SvtxEvalStack(topNode);
		_svtxevalstack->set_strict(false);
		_svtxevalstack->set_verbosity(0);
		_svtxevalstack->set_use_initial_vertex(1);
		_svtxevalstack->set_use_genfit_vertex(1);
		_svtxevalstack->next_event(topNode);
	}
	else
	{
		_svtxevalstack->next_event(topNode);
	}


	m_truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

	bool VtxToQA = false;



	int BadVtx = 0;
	float SVtoPVDis= 0;
	float SVtoPVTau = 0;


	float DevPx;
	float DevPy;
	float DevPz;
	float DevE;



	std::vector<int> ParentTrkInfo;
	std::vector<int> ParentEInfo;
	std::vector<double> DiffEPerVertex;

	std::vector<int> ParentPxInfo;
	std::vector<double> DiffPxPerVertex;

	std::vector<int> ParentPyInfo;
	std::vector<double> DiffPyPerVertex;

	std::vector<int> ParentPzInfo;
	std::vector<double> DiffPzPerVertex;


	std::vector<std::vector<int>> DaughterInfo;
	std::vector<int> VertexInfo;
	std::vector<int> HFIndexInfo;


	float CosTheta = -2;







	PHG4TruthInfoContainer::ConstRange range = m_truth_info->GetParticleRange();
	for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
			iter != range.second; ++iter) 
	{


		PHG4Particle* g4particle = iter->second;

	
		int gflavor = g4particle->get_pid();

		NParticles = NParticles + 1;


		int ParentPDGID = -1;
		int GrandParentPDGID = -1;
		int ParentTrkId = 0;

		double ParentE = 0;
		double ParentPx = 0;
		double ParentPy = 0;
		double ParentPz = 0;

	

		PHG4Particle* mother = NULL;



		TVector3 HFParMom(0,0,0);
		TVector3 HFProdVtx(0,0,0);
		TVector3 HFDecayVtx(0,0,0);

		TLorentzVector HFParFourMom(0,0,0,0);


		if (g4particle->get_parent_id() == 0)
		{
			ParentPDGID = 0; 


		}else{

			mother = m_truth_info->GetParticle(g4particle->get_parent_id());
			ParentPDGID = mother->get_pid();
			ParentTrkId = mother->get_track_id();
			ParentE = mother->get_e();
			ParentPx = mother->get_px();
			ParentPy = mother->get_py();
			ParentPz = mother->get_pz();

			HFParMom.SetXYZ(mother->get_px(),mother->get_py(),mother->get_pz());
			HFParFourMom.SetXYZT(mother->get_px(),mother->get_py(),mother->get_pz(),mother->get_e());

			if (mother->get_parent_id() == 0)	GrandParentPDGID = 0;



		}

		int NDig  = (int) log10(abs(gflavor));
		int firstDigit = (int)(abs(gflavor) / pow(10, NDig)); 	
		if((firstDigit == 4 || firstDigit == 5) && ParentPDGID == 0){

			int HFFillIndex = -99;

			for(int q = 0; q < NHFQA; q++){

				if(abs(gflavor) == QAVtxPDGID[q])  HFFillIndex = q;
			}

			HFHadronStat->Fill(HFFillIndex);

		}



		int VtxSize = ParentTrkInfo.size(); 
		bool NewVtx = true;
		int Index = -1;
		int HFIndex = -1;


		for(int i = 0; i < VtxSize; i++){

			if(ParentTrkId != 0 && ParentTrkId == ParentTrkInfo[i]){
				NewVtx = false;
				Index = i;

			}

		}



		for(int p = 0; p < NHFQA; p++){

			if(abs(ParentPDGID) == QAVtxPDGID[p]){

				VtxToQA = true;
				HFIndex = p;
			}

		}




		if((ParentTrkId > 0 || abs(gflavor) == abs(ParentPDGID)) && VtxToQA == true){   
			if(NewVtx){
				ParentTrkInfo.push_back(ParentTrkId);
				ParentEInfo.push_back(ParentE);
				DiffEPerVertex.push_back(ParentE - g4particle->get_e());
				ParentPxInfo.push_back(ParentPx);
				DiffPxPerVertex.push_back(ParentPx - g4particle->get_px());
				ParentPyInfo.push_back(ParentPy);
				DiffPyPerVertex.push_back(ParentPy - g4particle->get_py());
				ParentPzInfo.push_back(ParentPz);
				DiffPzPerVertex.push_back(ParentPz - g4particle->get_pz());


				VertexInfo.push_back(ParentPDGID);
				HFIndexInfo.push_back(HFIndex);

				std::vector<int> Daughters;

				Daughters.push_back(gflavor);
				DaughterInfo.push_back(Daughters);

			}
			if(!NewVtx){
				DiffEPerVertex[Index] = DiffEPerVertex[Index]  - g4particle->get_e();
				DiffPxPerVertex[Index] = DiffPxPerVertex[Index]  - g4particle->get_px();
				DiffPyPerVertex[Index] = DiffPyPerVertex[Index]  - g4particle->get_py();
				DiffPzPerVertex[Index] = DiffPzPerVertex[Index]  - g4particle->get_pz();

				DaughterInfo[Index].push_back(gflavor);



			}

		}


		SvtxTruthEval* trutheval = _svtxevalstack->get_truth_eval();
		PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);			

		PHG4VtxPoint* ParentVtx = trutheval->get_vertex(mother);			


		if(GrandParentPDGID == 0  && VtxToQA){

			float ParMass = sqrt(mother->get_e() * mother->get_e() - mother->get_px()  * mother->get_px()  - mother->get_py()  * mother->get_py() - mother->get_pz()  * mother->get_pz());
			float ParP = sqrt(mother->get_px()  * mother->get_px()  +  mother->get_py()  * mother->get_py() + mother->get_pz()  * mother->get_pz());

			HFProdVtx.SetXYZ(ParentVtx->get_x(),ParentVtx->get_y(),ParentVtx->get_z());			
			HFDecayVtx.SetXYZ(vtx->get_x(),vtx->get_y(),vtx->get_z());

			SVtoPVDis = (HFDecayVtx - HFProdVtx).Mag();
			SVtoPVTau = SVtoPVDis/ParP * ParMass;


			CosTheta = ((HFDecayVtx - HFProdVtx).Dot(HFParMom))/((HFDecayVtx - HFProdVtx).Mag() * HFParMom.Mag());



			QACosTheta->Fill(CosTheta);
			MassHis->Fill(ParMass);

			if(HFIndex > -1) ProperLifeTime[HFIndex]->Fill(SVtoPVTau); 


		}






	}


	int VtxSizeFinal = DiffEPerVertex.size();

	for(int q = 0; q < VtxSizeFinal; q++){


		if(abs(DiffEPerVertex[q]) > 0.004){

			BadVtx = BadVtx + 1;

		}

		DevPx = DiffPxPerVertex[q]/ParentPxInfo[q];
		DevPy = DiffPyPerVertex[q]/ParentPyInfo[q];
		DevPz = DiffPzPerVertex[q]/ParentPzInfo[q];
		DevE = DiffEPerVertex[q]/ParentEInfo[q];

		QAPx->Fill(DevPx);
		QAPy->Fill(DevPy);
		QAPz->Fill(DevPz);
		QAE->Fill(DevE);




	}
	//BR Working here

	for(int q = 0; q < VtxSizeFinal; q++){
		int HFIndexToFill = HFIndexInfo[q];
	
		std::vector<int> ChannelID = Channel(VertexInfo[q],DaughterInfo[q]);
		int ChannelSize = ChannelID.size();



		if(HFIndexToFill < 0){
	
			continue;
		}

		for(int r = 0 ; r <  ChannelSize; r++){


	
			BR1DHis[HFIndexToFill]->Fill(ChannelID[r]);
			
		}



		//		std::cout << "------------------------------------------------------------------------------"  << std::endl;
		ProperLifeTime[HFIndexToFill]->Fill(SVtoPVTau);

	}



	SVtoPVDistance->Fill(SVtoPVDis);



	TotalVtxStat->Fill(VtxSizeFinal);

	BadVtxStat->Fill(BadVtx);

	BadVtxPercent->Fill(BadVtx);




	LifeTime = SVtoPVTau;

	if(m_write_nTuple) QATree->Fill();

	EvtID = EvtID + 1;



	NPartHis->Fill(NParticles);



	
	return Fun4AllReturnCodes::EVENT_OK;


}



std::vector<int> QAG4Decayer::Channel(int pdgid, std::vector<int> Daughter){

	Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
	assert(hm);

	TH1F * OtherHis = dynamic_cast<TH1F *>(hm->getHisto("OtherHis"));
	assert(OtherHis);


	std::vector<int> Channel;

	int DaughterSize = Daughter.size();

	int TotalPiP = 0;
	int TotalKP = 0;
	int TotalPP = 0;
	int TotalElecP = 0;
	int TotalMuP = 0;
	int TotalPiM = 0;
	int TotalKM = 0;		
	int TotalPM = 0;
	int TotalElecM = 0;
	int TotalMuM = 0;
	int TotaleNu = 0;
	int TotalmuNu = 0;
	int TotaltauNu = 0;
	int TotaleNuBar = 0;
	int TotalmuNuBar = 0;
	int TotaltauNuBar = 0;
	int TotalN = 0;
	int TotalNBar = 0;
	int TotalKs = 0;
	int TotalKL = 0;
	int TotalGamma = 0;
	int TotalPiZ = 0;
	int TotalX = 0;

	//Intermediate State Quicky Decaying: Light Flavors

	int TotalRhoZ = 0;
	int TotalRhoP = 0;
	int TotalRhoM = 0;

	int TotalPhi = 0;
//	int TotalKStarZ = 0;
	int TotalKStarP = 0;
	int TotalKStarM = 0;
//	int TotalKStarZBar = 0;

	int TotalEta = 0;
	int TotalEtaPrime = 0;
	int Totalomega = 0;

	//Intermediate State Quicky Decaying: Heavy Flavors

	int TotalDZ = 0;
	int TotalDZBar = 0;	

	int TotalDP = 0;
	int TotalDM = 0;

	int TotalDsP = 0;
	int TotalDsM = 0;	
	int TotalJpsi = 0;

	int TotalOther = 0;
	int TotalKZ = 0;
	int TotalKStarZ = 0;

	//Extra Intermediate States
	int TotalKStar1430 = 0;
	int TotalKStar1430Bar = 0;


	for(int s = 0; s < DaughterSize; s++){

		//std::cout << "Daughter[s] = " << Daughter[s] << std::endl;

		if(Daughter[s] == 211) TotalPiP++;
		else if(Daughter[s] == 321) TotalKP++;
		else if(Daughter[s] == 2212) TotalPP++;
		else if(Daughter[s] == -11) TotalElecP++;
		else if(Daughter[s] == -13) TotalMuP++;
		else if(Daughter[s] == -211) TotalPiM++;
		else if(Daughter[s] == -321) TotalKM++;
		else if(Daughter[s] == -2212) TotalPM++;
		else if(Daughter[s] == 11) TotalElecM++;
		else if(Daughter[s] == 13) TotalMuM++;
		else if(Daughter[s] == 12) TotaleNu++;
		else if(Daughter[s] == 14) TotalmuNu++;
		else if(Daughter[s] == 16) TotaltauNu++;
		else if(Daughter[s] == -12) TotaleNuBar++;
		else if(Daughter[s] == -14) TotalmuNuBar++;
		else if(Daughter[s] == -16) TotaltauNuBar++;
		else if(Daughter[s] == 2112) TotalN++;
		else if(Daughter[s] == -2112) TotalNBar++;
		else if(Daughter[s] == 310) TotalKs++;
		else if(Daughter[s] == 130) TotalKL++;
		else if(Daughter[s] == 22) TotalGamma++;
		else if(Daughter[s] == 111) TotalPiZ++;
		else if(Daughter[s] == 421) TotalDZ++;
		else if(Daughter[s] == -421) TotalDZBar++;
		else if(Daughter[s] == 411) TotalDP++;
		else if(Daughter[s] == -411) TotalDM++;
		else if(Daughter[s] == 431) TotalDsP++;
		else if(Daughter[s] == -431) TotalDsM++;
		else if(Daughter[s] == 443) TotalJpsi++;
		else if(abs(Daughter[s]) == 20443) TotalX++;
		else if(Daughter[s] == 333) TotalPhi++;
		else if(abs(Daughter[s]) == 311) TotalKZ++;
		else if(abs(Daughter[s]) == 313) TotalKStarZ++;
		else if(Daughter[s] == 221) TotalEta++;
		else if(Daughter[s] == 10311) TotalKStar1430++;
		else if(Daughter[s] == -10311) TotalKStar1430Bar++;

		//	else if(Daughter[s] == -313) TotalKStarZ++;		
		else{
			//std::cout << "Daughter[s]  = " << Daughter[s]  << std::endl;
			TotalOther++;

		}

	}



	if(TotalOther > 0){
		OtherHis->Fill(0);
		Channel.push_back(-1);
		return Channel;
	}

//	std::cout << "Pass 2" << std::endl;

	if(pdgid == 411){  //D+

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 2 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 		

	//	if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 1 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 1 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 1 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 1 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); 		

		if(TotalDZ == 1) Channel.push_back(7); 		
		if(TotalPhi == 1) Channel.push_back(8);  //D0 -> K-e+ve
		if(TotalElecM == 1) Channel.push_back(9);  //D0 -> K-mu+vmu
	}

//	std::cout << "Pass 3" << std::endl;


	if(pdgid == 421){



		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 1 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 3 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 2 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); 		

		if(TotalKs > 0) Channel.push_back(7); 		
		if(TotalPiZ > 0) Channel.push_back(8); 				

		if(TotalElecP == 1) Channel.push_back(9);  //D0 -> K-mu+vmu
	}

	if(pdgid == 431){

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 1 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 1 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		 //Ds -> mu+ vmu phi

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 1 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 		 //Ds -> mu+ vmu


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 1 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 1 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 1 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); 		


		if(TotalDsM == 1) Channel.push_back(7); 		
		if(TotalPhi == 1) Channel.push_back(8); 				
		if(TotalElecP == 1) Channel.push_back(9);  //D0 -> K-mu+vmu

	}

	if(pdgid == 4122){   //Lambda_c

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		 //Ds -> mu+ vmu phi

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 		 //Ds -> mu+ vmu


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 1 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 2 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 1 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); 		


		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); 		

		if(TotalKStarZ==0&&TotalKZ==0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); 		

		if(TotalPP == 1) Channel.push_back(7); 		
		if(TotalN == 1) Channel.push_back(8); 				
		if(TotalElecP == 1) Channel.push_back(9);  //D0 -> K-mu+vmu

	}
//	std::cout << "Pass 4" << std::endl;

	if(pdgid == 511){

		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); //Here	
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); //	
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); //	
		//if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 1 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); //	
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 1 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0  && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 1 &&  TotalDsM == 0) Channel.push_back(6); //Leptonic Mode
		if(TotalX > 0) Channel.push_back(7); //Leptonic Mode
		if(TotalDZ == 1 || TotalDZBar == 1 ||  TotalDP == 1 ||  TotalDM == 1 ||  TotalDsP == 1||  TotalDsM == 1) Channel.push_back(8);  
		if(TotalJpsi == 1) Channel.push_back(9); 


	}


	if(pdgid == 521){

		//Decay B+ -> K+ K- pi+ without any intermediate resonances 
		if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 0 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		
		
		//Decay B+ -> K+ K- K+ without any resonances 		
		if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 0 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 2 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 
		
		//if(TotalX == 1 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); //Leptonic Mode	
		//if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 0 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 2 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); 
		
		//Decay B+ -> K+ phi(1020) (phi likely to decay into K+ K-)	
		if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 0 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 1 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); 

		//Decay B+ -> K+ Jpsi	
		if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 0 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); 
	//	if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); //Leptonic Mode

		//Decay B+ -> K+ DBar	
		if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 0 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 1 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); //Leptonic Mode
	
		//Decay B+ -> DBar Ds+	
		
		if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 0 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 1 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 1 &&  TotalDsM == 0) Channel.push_back(5); //Leptonic Mode
	//	if(TotalX == 1 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); //Leptonic Mode		

		//Decay B+ -> K+ K*0Bar(1430) (K*0Bar(1430) -> K- pi+)

		if(TotalKStar1430 == 0 &&  TotalKStar1430Bar == 1 && TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); //Leptonic Mode

	//	if(TotalPP == 1  && TotalPM == 1) Channel.push_back(6); //3 prong p p bar + anyhadron (p and p bar identifiable) 
		//Decay B+ -> e ve X
		if(TotalElecP == 1  && TotaleNu == 1) Channel.push_back(7); //Leptonic Mode

		//Decay B+ -> D + X
		if(TotalDZ == 1 || TotalDZBar == 1 ||  TotalDP == 1 ||  TotalDM == 1 ||  TotalDsP == 1||  TotalDsM == 1) Channel.push_back(8);  

		//Decay B+ -> Jpsi + X
		if(TotalJpsi == 1) Channel.push_back(9); 

	}


	if(pdgid == 531){



		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 1 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 3 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 1 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); //Leptonic Mode	
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 2 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 1 && TotalDM == 1 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 1 &&  TotalDsM == 1) Channel.push_back(6); //Leptonic Mode
		if(TotalElecP == 1  && TotaleNu == 1) Channel.push_back(7); //Leptonic Mode
		if(TotalDsM == 1) Channel.push_back(8);  
		if(TotalJpsi == 1) Channel.push_back(9); 




	}




	if(pdgid == 443){


		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); //Leptonic Mode	
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); //Leptonic Mode
		if(TotalGamma == 1) Channel.push_back(7); //Leptonic Mode
		if(TotalPhi == 1) Channel.push_back(8);  
		if(TotalDZ == 1 || TotalDZBar == 1 ||  TotalDP == 1 ||  TotalDM == 1 ||  TotalDsP == 1||  TotalDsM == 1) Channel.push_back(9); 


	}

	//std::cout << "Pass 5" << std::endl;


	if(pdgid == 553){


		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(0); 		
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(1); 
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(2); //Leptonic Mode	
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(3); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(4); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(5); //Leptonic Mode
		if(TotalKStarZ==0&&TotalKZ==0&&TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0  && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0&& TotalDsP == 0 &&  TotalDsM == 0) Channel.push_back(6); //Leptonic Mode
		if(TotalGamma == 1) Channel.push_back(7); //Leptonic Mode
		if(TotalJpsi == 1) Channel.push_back(8);  
		if(TotalDZ == 1 || TotalDZBar == 1 ||  TotalDP == 1 ||  TotalDM == 1 ||  TotalDsP == 1||  TotalDsM == 1) Channel.push_back(9); 


	}


	return Channel;


}


int QAG4Decayer::End(PHCompositeNode *topNode) 
{
	assert(topNode);
	
	if(m_SaveFiles){
		fout->cd();
		if(m_write_nTuple) QATree->Write();
		fout->Close();
	}
	
	return Fun4AllReturnCodes::EVENT_OK;
}


