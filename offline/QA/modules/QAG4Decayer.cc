#include "QAG4Decayer.h"

#include "QAHistManagerDef.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <decayfinder/DecayFinder.h>
#include <decayfinder/DecayFinderContainerBase.h>  // for DecayFinderContainerBase::Iter
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <CLHEP/Vector/LorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <TDatabasePDG.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <g4eval/SvtxEvalStack.h>

#include <TLorentzVector.h>
#include <TVector3.h>

int BadParticle;

TCanvas* c;

/*
 *  QA module to check that decayers obey laws of physics
 *  Can also measure geometrical acceptance in eta and pT
 *  (Useful for cross-section measurements and similar)
 *  Authors: Cameron Dean and Thomas Marshall
 *  Date: July 2022
 */

QAG4Decayer::QAG4Decayer()
{
}

QAG4Decayer::~QAG4Decayer()
{
}

int QAG4Decayer::Init(PHCompositeNode* topNode)
{
  BadParticle = 0;
  NParticles = 0;

  gStyle->SetOptStat(0);
  //	std::cout << "Using Updated Version of QA Module - Fixed - GrandMom - FIT B Meson - More BR Here - With Stat" << std::endl;

  c = new TCanvas("c", "c", 600, 600);

  assert(topNode);

  //Declear All Histograms//

  SVtoPVDistance = new TH1D("SVtoPVDistance", "", 100, 0, 0.2);
  SVtoPVDistance->GetXaxis()->SetTitle("SV to PV Distance (cm)");
  SVtoPVDistance->GetYaxis()->SetTitle("Number of Events");
  SVtoPVDistance->GetXaxis()->CenterTitle();
  SVtoPVDistance->GetYaxis()->CenterTitle();
  SVtoPVDistance->GetYaxis()->SetTitleOffset(1.4);
  SVtoPVDistance->SetTitle("D^{0} Meson SV to PV Distance Distribution");

  ProperLifeTime = new TH1D("ProperLifeTime", "", 100, 0, 0.05);
  ProperLifeTime->GetXaxis()->SetTitle("Decay Life Time c #tau (cm)");
  ProperLifeTime->GetYaxis()->SetTitle("Number of Events");
  ProperLifeTime->GetXaxis()->CenterTitle();
  ProperLifeTime->GetYaxis()->CenterTitle();
  ProperLifeTime->GetYaxis()->SetTitleOffset(1.4);
  ProperLifeTime->SetTitle("#Omega_{c}^{0} Baryon Proper Lifetime #tau Distribution");

  DecayBR = new TH2D("DecayBR", "", 10, -0.5, 9.5, 10, -0.5, 9.5);
  DecayBR->GetYaxis()->SetTitle("Decay Channel");
  DecayBR->GetXaxis()->SetTitle("Particle Type`");
  DecayBR->GetXaxis()->CenterTitle();
  DecayBR->GetYaxis()->CenterTitle();
  DecayBR->GetYaxis()->SetTitleOffset(1.4);
  DecayBR->SetTitle("Decay Branchong Ratio Statistics");

  TotalVtxStat = new TH1D("TotalVtxStat", "", 20, -0.5, 19.5);
  TotalVtxStat->GetXaxis()->SetTitle("Number of Vertex per Event");
  TotalVtxStat->GetYaxis()->SetTitle("Number of Events");
  TotalVtxStat->GetXaxis()->CenterTitle();
  TotalVtxStat->GetYaxis()->CenterTitle();
  TotalVtxStat->GetYaxis()->SetTitleOffset(1.4);
  TotalVtxStat->GetXaxis()->SetTitle("D^{0} Event Decay Vertices");

  BadVtxStat = new TH1D("BadVtxStat", "", 20, -0.5, 19.5);
  BadVtxStat->GetXaxis()->SetTitle("Number of Bad Vertex per Event");
  BadVtxStat->GetYaxis()->SetTitle("Number of Events");
  BadVtxStat->GetXaxis()->CenterTitle();
  BadVtxStat->GetYaxis()->CenterTitle();
  BadVtxStat->GetYaxis()->SetTitleOffset(1.4);
  BadVtxStat->SetTitle("D^{0} Event Decay With Missing Vertices");

  BadVtxPercent = new TH1D("BadVtxPercent", "", 20, -0.5, 19.5);
  BadVtxPercent->GetXaxis()->SetTitle("Number of Bad Vertex per Event");
  BadVtxPercent->GetYaxis()->SetTitle("Percentage (%)");
  BadVtxPercent->GetXaxis()->CenterTitle();
  BadVtxPercent->GetYaxis()->CenterTitle();
  BadVtxPercent->GetYaxis()->SetTitleOffset(1.4);
  BadVtxPercent->GetXaxis()->SetTitle("D^{0} Event Decay with Missing Vertices");

  QAPx = new TH1D("QAPx", "", 100, -1, 1);
  QAPx->GetXaxis()->SetTitle("( Parent p_{x} - Children Summed p_{x})/Parent p_{x}");
  QAPx->GetYaxis()->SetTitle("Number of Vertices");
  QAPx->GetXaxis()->CenterTitle();
  QAPx->GetYaxis()->CenterTitle();
  QAPx->GetYaxis()->SetTitleOffset(1.4);
  QAPx->SetTitle("Children Deviation from Parent: Px");

  QAPy = new TH1D("QAPy", "", 100, -1, 1);
  QAPy->GetXaxis()->SetTitle("(Parent p_{y} - Children Summed p_{y})/Parent p_{y}");
  QAPy->GetYaxis()->SetTitle("Number of Vertices");
  QAPy->GetXaxis()->CenterTitle();
  QAPy->GetYaxis()->CenterTitle();
  QAPy->GetYaxis()->SetTitleOffset(1.4);
  QAPy->SetTitle("Children Deviation from Parent: Py");

  QAPz = new TH1D("QAPz", "", 100, -1, 1);
  QAPz->GetXaxis()->SetTitle("(Parent p_{z} - Children Summed p_{z} )/Parent p_{z}");
  QAPz->GetYaxis()->SetTitle("Number of Vertices");
  QAPz->GetXaxis()->CenterTitle();
  QAPz->GetYaxis()->CenterTitle();
  QAPz->GetYaxis()->SetTitleOffset(1.4);
  QAPz->SetTitle("Children Deviation from Parent: Pz");

  QAE = new TH1D("QAE", "", 100, -1, 1);
  QAE->GetXaxis()->SetTitle("(Parent E - Children Summed E)/Parent E");
  QAE->GetYaxis()->SetTitle("Number of Vertices");
  QAE->GetXaxis()->CenterTitle();
  QAE->GetYaxis()->CenterTitle();
  QAE->GetYaxis()->SetTitleOffset(1.4);
  QAE->SetTitle("Children Deviation from Parent: E");

  QACosTheta = new TH1D("QACosTheta", "", 120, -1.2, 1.2);
  QACosTheta->GetXaxis()->SetTitle("Cos(#theta)");
  QACosTheta->GetYaxis()->SetTitle("Number of Event");
  QACosTheta->GetXaxis()->CenterTitle();
  QACosTheta->GetYaxis()->CenterTitle();
  QACosTheta->GetYaxis()->SetTitleOffset(1.4);
  QACosTheta->SetTitle("Cos(#theta) Distribution for Sanity Check");

  MassHis = new TH1D("MassHis", "", 100, 1.6, 2.1);
  MassHis->GetXaxis()->SetTitle("D^{0} Invariant Mass (GeV/c^{2})");
  MassHis->GetYaxis()->SetTitle("Number of Event");
  MassHis->GetXaxis()->CenterTitle();
  MassHis->GetYaxis()->CenterTitle();
  MassHis->GetYaxis()->SetTitleOffset(1.4);
  MassHis->SetTitle("Reconstructed D^{0} Invariant Mass from GEANT in PHG4TruthInfoContainer");

  NPartHis = new TH1D("NPartHis", "", 10, -0.5, 9.5);
  NPartHis->GetXaxis()->SetTitle("Number of Particles Per Event");
  NPartHis->GetYaxis()->SetTitle("Number of Event");
  NPartHis->GetXaxis()->CenterTitle();
  NPartHis->GetYaxis()->CenterTitle();
  NPartHis->GetYaxis()->SetTitleOffset(1.4);
  NPartHis->SetTitle("Number of Particles registered in D^{0} -> Kpi in PHG4TruthInfoContainer");

  BR1DHis = new TH1D("BR1DHis", "", 10, -0.5, 9.5);
  BR1DHis->GetXaxis()->SetTitle("Decay Channel");
  BR1DHis->GetYaxis()->SetTitle("Branching Ratio");
  BR1DHis->GetXaxis()->CenterTitle();
  BR1DHis->GetYaxis()->CenterTitle();
  BR1DHis->GetYaxis()->SetTitleOffset(1.4);
  BR1DHis->SetTitle("Decay Branching Ratio Statisticsr");

  HFHadronStat = new TH1D("HFHadronStat", "", 14, -0.5, 13.5);
  HFHadronStat->GetXaxis()->SetTitle("Heavy Flavor Particle ID");
  HFHadronStat->GetYaxis()->SetTitle("Counts");
  HFHadronStat->GetXaxis()->CenterTitle();
  HFHadronStat->GetYaxis()->CenterTitle();
  HFHadronStat->GetYaxis()->SetTitleOffset(1.4);
  HFHadronStat->SetTitle("Promptly Produced Heavy Flavor Particle Statistics");

  //    _svtxevalstack = new SvtxEvalStack(topNode);

  //	fout->cd();
  fout = new TFile("MyQAFile.root", "RECREATE");
  QATree = new TTree("QATree", "QATree");
  EvtID = 0;
  LifeTime = 0;

  NParticles = 0;
  D0Px = 0;
  D0Py = 0;
  D0Pz = 0;
  D0E = 0;

  PVDaughtersPDGID.clear();

  QATree->Branch("EvtID", &EvtID, "EvtID/I");
  QATree->Branch("NParticles", &NParticles, "NParticles/I");
  QATree->Branch("D0Px", &D0Px, "D0Px/F");
  QATree->Branch("D0Py", &D0Py, "D0Py/F");
  QATree->Branch("D0Pz", &D0Pz, "D0Pz/F");
  QATree->Branch("D0E", &D0E, "D0E/F");
  QATree->Branch("PVDaughtersPDGID", &PVDaughtersPDGID);

  QATree->Branch("LifeTime", &LifeTime, "LifeTime/F");

  return 0;
}

int QAG4Decayer::process_event(PHCompositeNode* topNode)
{
  NParticles = 0;
  PVDaughtersPDGID.clear();

  /*

	   std::cout << "Process QADecayer GEANT 4 External Decayer Shits" << std::endl;

	   std::cout << "Now Print The Nodes We have for QADecayer.cc" <<	std::endl;
	   std::cout << "------------------------------------------------------------------------------------" <<	std::endl;
	   topNode->print();
	   std::cout << "------------------------------------------------------------------------------------" <<	std::endl;


	   std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
	   */
  //  auto m_truth_info2 = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // assert(m_truth_info2);

  //  std::shared_ptr<SvtxEvalStack>  _svtxevalstack;

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

  //	std::cout << "------------------------------------------ Asserted ----------------------------------------------" << std::endl;

  _svtxevalstack = new SvtxEvalStack(topNode);

  m_truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  /*

	   auto _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	   assert(_vertexmap);
	   
	   */
  //m_truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  //	const int NHFQA = 9;
  //	int QAVtxPDGID[NHFQA] = {411,421,431,443,511,521,531,553,4122};

  const int NHFQA = 13;
  int QAVtxPDGID[NHFQA] = {411, 421, 431, 511, 521, 531, 553, 4122, 443, 100443, 200443, 100553, 200553};
  bool VtxToQA = false;

  int BadVtx = 0;
  float SVtoPVDis = 0;
  float SVtoPVTau = 0;
  //	float SVtoPVTau2 = 0;

  float DevPx;
  float DevPy;
  float DevPz;
  float DevE;

  /*	
		float ParentPx;
		float TotalPx;
		float ParentPy;
		float TotalPy;
		float ParentPz;
		float TotalPz;

		float ParentE;
		float TotalE;
		*/

  //	float ParentE;

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

  //	float ReduceE;

  float CosTheta = -2;

  //std::cout << "PHG4TruthInfoContainer?  Anything "<< std::endl;

  PHG4TruthInfoContainer::ConstRange range = m_truth_info->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second; ++iter)
  {
    //m_g4particle = iter->second;
    PHG4Particle* g4particle = iter->second;

    //	int gtrackID = g4particle->get_track_id();
    int gflavor = g4particle->get_pid();

    //std::cout << "QA: gflavor = " << gflavor << std::endl;
    //if((abs(gflavor) > 500 && abs(gflavor) < 550) || (abs(gflavor) == 5)) std::cout << "beauty hadron detected " << std::endl;
    //if(abs(gflavor) == 10441 || abs(gflavor) ==  20443 || abs(gflavor) == 4332) std::cout << "Bad Particle Registered " << std::endl;

    NParticles = NParticles + 1;

    int ParentPDGID = -1;
    int GrandParentPDGID = -1;
    int ParentTrkId = 0;

    double ParentE = 0;
    double ParentPx = 0;
    double ParentPy = 0;
    double ParentPz = 0;

    //std::cout << "Before: gtrackID = " << gtrackID << "   gflavor = " << gflavor << std::endl;

    PHG4Particle* mother = NULL;

    TVector3 HFParMom(0, 0, 0);
    TVector3 HFProdVtx(0, 0, 0);
    TVector3 HFDecayVtx(0, 0, 0);

    TLorentzVector HFParFourMom(0, 0, 0, 0);

    if (g4particle->get_parent_id() == 0)
    {
      ParentPDGID = 0;
      //		mother = m_truth_info->GetParticle(ParentPDGID);
    }
    else
    {
      mother = m_truth_info->GetParticle(g4particle->get_parent_id());
      ParentPDGID = mother->get_pid();
      ParentTrkId = mother->get_track_id();
      ParentE = mother->get_e();
      ParentPx = mother->get_px();
      ParentPy = mother->get_py();
      ParentPz = mother->get_pz();

      HFParMom.SetXYZ(mother->get_px(), mother->get_py(), mother->get_pz());
      HFParFourMom.SetXYZT(mother->get_px(), mother->get_py(), mother->get_pz(), mother->get_e());

      if (mother->get_parent_id() == 0) GrandParentPDGID = 0;
    }

    int NDig = (int) log10(abs(gflavor));
    int firstDigit = (int) (abs(gflavor) / pow(10, NDig));
    //std::cout << "firstDigit = " << firstDigit << std::endl;
    if ((firstDigit == 4 || firstDigit == 5) && ParentPDGID == 0)
    {
      int HFFillIndex = -99;

      for (int q = 0; q < NHFQA; q++)
      {
        if (abs(gflavor) == QAVtxPDGID[q]) HFFillIndex = q;
      }

      HFHadronStat->Fill(HFFillIndex);
    }

    int VtxSize = ParentTrkInfo.size();
    bool NewVtx = true;
    int Index = -1;
    int HFIndex = -1;

    for (int i = 0; i < VtxSize; i++)
    {
      if (ParentTrkId != 0 && ParentTrkId == ParentTrkInfo[i])
      {
        NewVtx = false;
        Index = i;
      }
    }

    //std::cout << "gtrackID = " << gtrackID << "   gflavor = " << gflavor << "   Energy = " <<   g4particle->get_e() << "   ParentTrkId = " << ParentTrkId << "  Parent Energy = " << ParentE << std::endl;

    //		std::cout << "gtrackID = " << gtrackID << "   gflavor = " << gflavor << "   ParentTrkId = " << ParentTrkId << std::endl;

    //		std::cout << "ParentPDGID  = " << ParentPDGID << "  NParticles = " << NParticles << std::endl;
    /*
		if(abs(ParentPDGID) == 4132){

			NParticles = NParticles + 1;
			PVDaughtersPDGID.push_back(gflavor);

		}
	*/

    if (gflavor == 521)
    {
      D0Px = g4particle->get_px();
      D0Py = g4particle->get_py();
      D0Pz = g4particle->get_pz();
      D0E = g4particle->get_e();
    }

    for (int p = 0; p < NHFQA; p++)
    {
      if (abs(ParentPDGID) == QAVtxPDGID[p])
      {
        VtxToQA = true;
        HFIndex = p;
      }
    }

    if ((ParentTrkId > 0 || abs(gflavor) == abs(ParentPDGID)) && VtxToQA == true)
    {  //abs(gflavor) == abs(ParentPDGID) HF oscillation for instance Bs -> anti-Bs
      if (NewVtx)
      {
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
        //	Daughters.clear();
        //				std::cout << "New Vertex: " << "  Parent Track ID = " << ParentTrkId << "   Parent Energy = " << ParentE << "    g4particle->get_e() = " << g4particle->get_e() << std::endl;
      }
      if (!NewVtx)
      {
        DiffEPerVertex[Index] = DiffEPerVertex[Index] - g4particle->get_e();
        DiffPxPerVertex[Index] = DiffPxPerVertex[Index] - g4particle->get_px();
        DiffPyPerVertex[Index] = DiffPyPerVertex[Index] - g4particle->get_py();
        DiffPzPerVertex[Index] = DiffPzPerVertex[Index] - g4particle->get_pz();

        DaughterInfo[Index].push_back(gflavor);

        //				std::cout << "No New Vertex: " << "  Parent Track ID = " << ParentTrkId << "   Parent Energy = " << ParentE << "    g4particle->get_e() = " << g4particle->get_e() << std::endl;
      }
    }

    SvtxTruthEval* trutheval = _svtxevalstack->get_truth_eval();
    PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);

    PHG4VtxPoint* ParentVtx = trutheval->get_vertex(mother);

    if (GrandParentPDGID == 0)
    {
      float ParMass = sqrt(mother->get_e() * mother->get_e() - mother->get_px() * mother->get_px() - mother->get_py() * mother->get_py() - mother->get_pz() * mother->get_pz());
      float ParP = sqrt(mother->get_px() * mother->get_px() + mother->get_py() * mother->get_py() + mother->get_pz() * mother->get_pz());

      HFProdVtx.SetXYZ(ParentVtx->get_x(), ParentVtx->get_y(), ParentVtx->get_z());
      HFDecayVtx.SetXYZ(vtx->get_x(), vtx->get_y(), vtx->get_z());

      //		std::cout << "HFDecayVtx.Mag()" << HFDecayVtx.Mag() << "    HFProdVtx.Mag() = " << HFProdVtx.Mag() << "    (HFDecayVtx - HFProdVtx).Mag() = " << (HFDecayVtx - HFProdVtx).Mag() << std::endl;

      SVtoPVDis = (HFDecayVtx - HFProdVtx).Mag();
      SVtoPVTau = SVtoPVDis / ParP * ParMass;

      CosTheta = ((HFDecayVtx - HFProdVtx).Dot(HFParMom)) / ((HFDecayVtx - HFProdVtx).Mag() * HFParMom.Mag());

      QACosTheta->Fill(CosTheta);
      MassHis->Fill(ParMass);
    }
  }

  int VtxSizeFinal = DiffEPerVertex.size();

  for (int q = 0; q < VtxSizeFinal; q++)
  {
    //			std::cout << "q = " << q << "    ParentTrkInfo = " << ParentTrkInfo[q] << "   DiffEPerVertex = " << DiffEPerVertex[q] << std::endl;
    if (abs(DiffEPerVertex[q]) > 0.004)
    {
      //	std::cout << "Unconserving Vtx Detected: q = " << q << "    ParentTrkInfo = " << ParentTrkInfo[q] << "   DiffEPerVertex = " << DiffEPerVertex[q] << std::endl;
      BadVtx = BadVtx + 1;
    }

    DevPx = DiffPxPerVertex[q] / ParentPxInfo[q];
    DevPy = DiffPyPerVertex[q] / ParentPyInfo[q];
    DevPz = DiffPzPerVertex[q] / ParentPzInfo[q];
    DevE = DiffEPerVertex[q] / ParentEInfo[q];

    QAPx->Fill(DevPx);
    QAPy->Fill(DevPy);
    QAPz->Fill(DevPz);
    QAE->Fill(DevE);
  }
  //BR Working here

  //	std::cout << "Now Work on Branching Ratio Studies" << std::endl;

  for (int q = 0; q < VtxSizeFinal; q++)
  {
    /*
		int DaughterSize = DaughterInfo[q].size();
		//int MotherID = VertexInfo[q];

		for(int s = 0; s < DaughterSize; s++){

			//			std::cout << "s = " << s << "   Vertex Child: " <<  DaughterInfo[q][s] << std::endl;


		}
*/

    int ChannelID = Channel(VertexInfo[q], DaughterInfo[q]);

    if (ChannelID > -1)
    {
      DecayBR->Fill(HFIndexInfo[q], ChannelID, 1);
    }
    //std::cout << "q = " << q << "   Vertex Parent: " << VertexInfo[q] << "   ChannelID = " << ChannelID << std::endl;

    if (VertexInfo[q] == 521) BR1DHis->Fill(ChannelID);

    //		std::cout << "------------------------------------------------------------------------------"  << std::endl;
  }

  //	std::cout << "DONE on Branching Ratio Studies" << std::endl;

  SVtoPVDistance->Fill(SVtoPVDis);

  ProperLifeTime->Fill(SVtoPVTau);

  TotalVtxStat->Fill(VtxSizeFinal);

  BadVtxStat->Fill(BadVtx);

  BadVtxPercent->Fill(BadVtx);

  //	DecayBR->Fill(Particle,BR);

  LifeTime = SVtoPVTau;

  QATree->Fill();

  EvtID = EvtID + 1;

  //std::cout << "NParticles" << std::endl;

  NPartHis->Fill(NParticles);

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4Decayer::Channel(int pdgid, std::vector<int> Daughter)
{
  int Channel = -1;

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
  int TotalKStarZ = 0;
  int TotalKStarP = 0;
  int TotalKStarM = 0;

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

  for (int s = 0; s < DaughterSize; s++)
  {
    if (Daughter[s] == 211) TotalPiP++;
    if (Daughter[s] == 321) TotalKP++;
    if (Daughter[s] == 2212) TotalPP++;
    if (Daughter[s] == -11) TotalElecP++;
    if (Daughter[s] == -13) TotalMuP++;

    if (Daughter[s] == -211) TotalPiM++;
    if (Daughter[s] == -321) TotalKM++;
    if (Daughter[s] == -2212) TotalPM++;
    if (Daughter[s] == 11) TotalElecM++;
    if (Daughter[s] == 13) TotalMuM++;

    if (Daughter[s] == 12) TotaleNu++;
    if (Daughter[s] == 14) TotalmuNu++;
    if (Daughter[s] == 16) TotaltauNu++;

    if (Daughter[s] == -12) TotaleNuBar++;
    if (Daughter[s] == -14) TotalmuNuBar++;
    if (Daughter[s] == -16) TotaltauNuBar++;

    if (Daughter[s] == 2112) TotalN++;
    if (Daughter[s] == -2112) TotalNBar++;

    if (Daughter[s] == 310) TotalKs++;
    if (Daughter[s] == 130) TotalKL++;

    if (Daughter[s] == 22) TotalGamma++;
    if (Daughter[s] == 111) TotalPiZ++;

    //Intermediate Particles
    if (Daughter[s] == 113) TotalRhoZ++;
    if (Daughter[s] == 213) TotalRhoP++;
    if (Daughter[s] == -213) TotalRhoM++;

    if (Daughter[s] == 313) TotalKStarZ++;
    if (Daughter[s] == 323) TotalKStarP++;
    if (Daughter[s] == -323) TotalKStarM++;
    if (Daughter[s] == 333) TotalPhi++;
    if (Daughter[s] == 221) TotalEta++;
    if (Daughter[s] == 331) TotalEtaPrime++;
    if (Daughter[s] == 223) Totalomega++;

    if (Daughter[s] == 421) TotalDZ++;
    if (Daughter[s] == -421) TotalDZBar++;

    if (Daughter[s] == 411) TotalDP++;
    if (Daughter[s] == -411) TotalDM++;

    if (Daughter[s] == 431) TotalDsP++;
    if (Daughter[s] == -431) TotalDsM++;

    if (Daughter[s] == 443) TotalJpsi++;
    if (abs(Daughter[s]) == 20443) TotalX++;
  }

  if (pdgid == 411)
  {  //D+

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;

    if (TotalPiP == 2 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 1 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 1 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 1 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 1 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;

    if (TotalDZ == 1) Channel = 7;
    if (TotalPhi == 1) Channel = 8;    //D0 -> K-e+ve
    if (TotalElecM == 1) Channel = 9;  //D0 -> K-mu+vmu
  }

  if (pdgid == 421)
  {
    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 1 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;

    if (TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;

    if (TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 3 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 2 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;

    if (TotalKs > 0) Channel = 7;
    if (TotalPiZ > 0) Channel = 8;

    if (TotalElecP == 1) Channel = 9;  //D0 -> K-mu+vmu
  }

  if (pdgid == 431)
  {
    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 1 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 1 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;  //Ds -> mu+ vmu phi

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 1 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;  //Ds -> mu+ vmu

    if (TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;

    if (TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 1 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 1 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 1 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 1 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;

    if (TotalDsM == 1) Channel = 7;
    if (TotalPhi == 1) Channel = 8;
    if (TotalElecP == 1) Channel = 9;  //D0 -> K-mu+vmu
  }

  if (pdgid == 4122)
  {  //Lambda_c

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;  //Ds -> mu+ vmu phi

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;  //Ds -> mu+ vmu

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 1 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;

    if (TotalPiP == 1 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 1 && TotalKL == 0 && TotalPiZ == 2 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;

    if (TotalPiP == 0 && TotalKP == 0 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 1 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;

    if (TotalPiP == 0 && TotalKP == 1 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;

    if (TotalPiP == 0 && TotalKP == 1 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;

    if (TotalPP == 1) Channel = 7;
    if (TotalN == 1) Channel = 8;
    if (TotalElecP == 1) Channel = 9;  //D0 -> K-mu+vmu
  }

  if (pdgid == 511)
  {
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Here
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 1 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 1 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 1 && TotalDsM == 0) Channel = 6;  //Leptonic Mode
    if (TotalX > 0) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           //Leptonic Mode
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 8;
    if (TotalJpsi == 1) Channel = 9;
  }

  if (pdgid == 521)
  {
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 1 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;    //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 1 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;    //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && (TotalKs == 1 || TotalKL == 1) && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 1 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 1 && TotalDsM == 0) Channel = 5;    //Leptonic Mode
    if (TotalX == 1 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;    //Leptonic Mode
    if (TotalElecP == 1 && TotaleNu == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       //Leptonic Mode
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 8;
    if (TotalJpsi == 1) Channel = 9;
  }

  if (pdgid == 531)
  {
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 1 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 1 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 3 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 1 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 2 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 1 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 1 && TotalDM == 1 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 1 && TotalDsM == 1) Channel = 6;  //Leptonic Mode
    if (TotalElecP == 1 && TotaleNu == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     //Leptonic Mode
    if (TotalDsM == 1) Channel = 8;
    if (TotalJpsi == 1) Channel = 9;
  }

  if (pdgid == 443)
  {
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;  //Leptonic Mode
    if (TotalGamma == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      //Leptonic Mode
    if (TotalPhi == 1) Channel = 8;
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 9;
  }

  if (pdgid == 100443)
  {
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;  //Leptonic Mode
    if (TotalGamma == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      //Leptonic Mode
    if (TotalPhi == 1) Channel = 8;
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 9;
  }

  if (pdgid == 200443)
  {
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;  //Leptonic Mode
    if (TotalGamma == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      //Leptonic Mode
    if (TotalPhi == 1) Channel = 8;
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 9;
  }

  if (pdgid == 553)
  {
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;  //Leptonic Mode
    if (TotalGamma == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      //Leptonic Mode
    if (TotalJpsi == 1) Channel = 8;
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 9;
  }

  if (pdgid == 100553)
  {
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;  //Leptonic Mode
    if (TotalGamma == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      //Leptonic Mode
    if (TotalJpsi == 1) Channel = 8;
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 9;
  }

  if (pdgid == 200553)
  {
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 0;
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 3 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 1 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 1;
    if (TotalX == 0 && TotalPiP == 1 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 1 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 2;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 1 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 1 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 3;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 4;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 1 && TotalMuP == 0 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 1 && TotalMuM == 0 && TotalGamma == 1 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 5;  //Leptonic Mode
    if (TotalX == 0 && TotalPiP == 0 && TotalKP == 0 && TotalPP == 0 && TotalElecP == 0 && TotalMuP == 1 && TotalPiM == 0 && TotalKM == 0 && TotalPM == 0 && TotalElecM == 0 && TotalMuM == 1 && TotalGamma == 0 && TotaleNu == 0 && TotalmuNu == 0 && TotaltauNu == 0 && TotaleNuBar == 0 && TotalmuNuBar == 0 && TotaltauNuBar == 0 && TotalN == 0 && TotalNBar == 0 && TotalKs == 0 && TotalKL == 0 && TotalPiZ == 0 && TotalRhoZ == 0 && TotalRhoP == 0 && TotalRhoM == 0 && TotalKStarZ == 0 && TotalKStarP == 0 && TotalKStarM == 0 && TotalPhi == 0 && TotalEta == 0 && TotalEtaPrime == 0 && Totalomega == 0 && TotalDZ == 0 && TotalDZBar == 0 && TotalDP == 0 && TotalDM == 0 && TotalJpsi == 0 && TotalDsP == 0 && TotalDsM == 0) Channel = 6;  //Leptonic Mode
    if (TotalGamma == 1) Channel = 7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      //Leptonic Mode
    if (TotalJpsi == 1) Channel = 8;
    if (TotalDZ == 1 || TotalDZBar == 1 || TotalDP == 1 || TotalDM == 1 || TotalDsP == 1 || TotalDsM == 1) Channel = 9;
  }

  return Channel;
}

int QAG4Decayer::End(PHCompositeNode* topNode)
{
  assert(topNode);
  /*
	c->SetRightMargin(0.15);

	c->cd();

	SVtoPVDistance->Draw();
	c->SaveAs("Plots/SVtoPVDistance.png");

	ProperLifeTime->Draw();
	c->SaveAs("Plots/ProperLifeTime.png");

	TF1 * func = new TF1("Func","[0] * TMath::Exp(-[1] * x)",0.005,0.05);
	func->SetParLimits(1,10,1000);
	ProperLifeTime->Fit(func,"R");

	float tau = 1/func->GetParameter(1) * 10000;
	float tauError = func->GetParError(1)/(func->GetParameter(1) * func->GetParameter(1)) * 10000;

	//std::cout << "Lifetime tau = " << tau << "   Error = " << tauError << std::endl;



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.040);
	lat->DrawLatex(0.35,0.80,Form("c #tau = %.1f #pm %.1f (#mu m)",tau,tauError));	

	c->SaveAs("Plots/ProperLifeTimeFitted.png");

	TotalVtxStat->Draw("hist");
	c->SaveAs("Plots/TotalVtxStat.png");

	BadVtxStat->Draw("hist");
	c->SaveAs("Plots/BadVtxStat.png");

	BadVtxPercent->Scale(1.0/BadVtxPercent->Integral());
	c->SaveAs("Plots/BadVtxPercent.png");



	QAPx->Draw();
	c->SaveAs("Plots/QAPx.png");
	c->SaveAs("Plots/QAE.png");

	QAPy->Draw();
	c->SaveAs("Plots/QAPy.png");

	QAPz->Draw();
	c->SaveAs("Plots/QAPz.png");

	QAE->Draw();
	c->SaveAs("Plots/QAE.png");

	QACosTheta->Draw();
	c->SaveAs("Plots/QACosTheta.png");

	DecayBR->Draw("COLZ");
	c->SaveAs("Plots/DecayBR.png");

	MassHis->Draw("COLZ");
	c->SaveAs("Plots/MassHis.png");

	NPartHis->Draw();
	c->SaveAs("Plots/NPartHis.png");


	//BR1DHis->Scale(1.0/(EvtID + 1));
	BR1DHis->Draw();
	c->SaveAs("Plots/BR1DHis.png");

	HFHadronStat->Draw();
	c->SaveAs("Plots/HFHadronStat.png");

*/

  fout->cd();
  //	QATree->Write();
  HFHadronStat->Write();
  BR1DHis->Write();
  fout->Close();

  std::cout << "BadParticle = " << BadParticle << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
