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

//#include <g4eval/SvtxEvalStack.h>
#include <TF1.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <map>

const int NHFQA = 9;
int QAVtxPDGID[NHFQA] = {411, 421, 431, 4122, 511, 521, 531, 443, 553};

std::multimap<std::vector<int>, int> decaymap[NHFQA];

/*
 *  QA module to check decay branching ratio, decay lifetime, and momentum conservation for inclusive heavy flavor hadron decay, which is handle by EvtGen as default
 *  Authors: Zhaozhong Shi
 *  Date: November 2022
 */

QAG4Decayer::QAG4Decayer(const std::string &name)
  : SubsysReco(name)
  , m_write_nTuple(false)
  //  , m_write_QAHists(true)
  , m_SaveFiles(false)
{
}

QAG4Decayer::~QAG4Decayer()
{
}

int QAG4Decayer::Init(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  for (int i = 0; i < NHFQA; i++)
  {
    h = new TH1F(Form("QAPx_%d", i), "", 200, -1, 1);
    hm->registerHisto(h);

    h = new TH1F(Form("QAPy_%d", i), "", 200, -1, 1);
    hm->registerHisto(h);

    h = new TH1F(Form("QAPz_%d", i), "", 200, -1, 1);
    hm->registerHisto(h);

    h = new TH1F(Form("QAE_%d", i), "", 200, -1, 1);
    hm->registerHisto(h);

    h = new TH1F(Form("InvMass_%d", i), "", 200, 1, 6);
    hm->registerHisto(h);

    h = new TH1F(Form("QACosTheta_%d", i), "", 120, -1.2, 1.2);
    hm->registerHisto(h);

    h = new TH1F(Form("BR1DHis_%d", i), "", 10, -0.5, 9.5);
    hm->registerHisto(h);

    h = new TH1F(Form("ProperLifeTime_%d", i), "", 100, 0.0001, 0.05);
    hm->registerHisto(h);
  }

  h = new TH1F("HFHadronStat", "", 9, -0.5, 8.5);
  hm->registerHisto(h);

  NParticles = 0;

  gStyle->SetOptStat(0);

  assert(topNode);

  EvtID = 0;
  LifeTime = 0;
  NParticles = 0;

  // D+
  decaymap[0].insert({{111, 211}, 0});
  decaymap[0].insert({{-211, 211, 211}, 1});
  decaymap[0].insert({{111, 211, 310}, 2});
  decaymap[0].insert({{-321, 211, 211}, 3});
  decaymap[0].insert({{-211, -211, 211, 211, 211}, 4});
  decaymap[0].insert({{211, 333}, 5});
  decaymap[0].insert({{-13, 14}, 6});
  decaymap[0].insert({{310, 321}, 7});

  // D0
  decaymap[1].insert({{-321, 211}, 0});
  decaymap[1].insert({{-321, 111, 211}, 1});
  decaymap[1].insert({{-211, 321}, 2});
  decaymap[1].insert({{-321, 310, 321}, 3});
  decaymap[1].insert({{310, 310, 310}, 4});
  decaymap[1].insert({{111, 111}, 5});
  decaymap[1].insert({{-211, 211}, 6});

  // Ds
  decaymap[2].insert({{-13, 14, 333}, 0});
  decaymap[2].insert({{-13, 14}, 1});
  decaymap[2].insert({{-321, 211, 321}, 2});
  decaymap[2].insert({{310, 321}, 3});
  decaymap[2].insert({{111, 111, 211}, 4});
  decaymap[2].insert({{-2112, 2212}, 5});
  decaymap[2].insert({{211, 333}, 6});

  // LambdaC
  decaymap[3].insert({{-321, 211, 2212}, 0});
  decaymap[3].insert({{-311, 2212}, 1});
  decaymap[3].insert({{333, 2212}, 2});
  decaymap[3].insert({{-321, 111, 211, 2212}, 3});
  decaymap[3].insert({{-311, -211, 211, 2212}, 4});
  decaymap[3].insert({{-211, -211, 211, 211, 2112}, 5});
  decaymap[3].insert({{-211, 211, 2212}, 6});

  // B0
  decaymap[4].insert({{-211, 321}, 0});
  decaymap[4].insert({{-211, 111, 321}, 1});
  decaymap[4].insert({{-211, 211}, 2});
  decaymap[4].insert({{111, 111}, 3});
  decaymap[4].insert({{-411, 211}, 4});
  decaymap[4].insert({{-411, 431}, 5});
  decaymap[4].insert({{-211, 321, 443}, 6});
  decaymap[4].insert({{-211, 211, 443}, 7});

  // B+
  decaymap[5].insert({{-321, 211, 321}, 0});
  decaymap[5].insert({{-321, 321, 321}, 1});
  decaymap[5].insert({{321, 333}, 2});
  decaymap[5].insert({{321, 443}, 3});
  decaymap[5].insert({{-421, 321}, 4});
  decaymap[5].insert({{-421, 431}, 5});
  decaymap[5].insert({{-10311, 321}, 6});

  // Bs
  decaymap[6].insert({{-321, 211}, 0});
  decaymap[6].insert({{-321, 321}, 1});
  decaymap[6].insert({{333, 333}, 2});
  decaymap[6].insert({{22, 333}, 3});
  decaymap[6].insert({{-431, 211}, 4});
  decaymap[6].insert({{-431, 431}, 5});
  decaymap[6].insert({{-321, 321, 443}, 6});
  decaymap[6].insert({{333, 443}, 7});

  // Jpsi
  decaymap[7].insert({{-11, 11}, 0});
  decaymap[7].insert({{-13, 13}, 1});
  decaymap[7].insert({{-211, 211}, 2});
  decaymap[7].insert({{-211, 111, 211}, 3});
  decaymap[7].insert({{-321, 111, 321}, 4});
  decaymap[7].insert({{-321, -211, 211, 321}, 5});
  decaymap[7].insert({{-2212, 2212}, 6});
  decaymap[7].insert({{-2112, 2112}, 7});
  decaymap[7].insert({{22, 22, 22}, 8});
  decaymap[7].insert({{22, 111}, 9});

  // Upsilon (1S)
  decaymap[8].insert({{-11, 11}, 0});
  decaymap[8].insert({{-13, 13}, 1});
  decaymap[8].insert({{-211, 22, 211}, 2});
  decaymap[8].insert({{-211, 111, 211}, 3});
  decaymap[8].insert({{-321, 22, 321}, 4});
  decaymap[8].insert({{22, 111, 111}, 5});
  decaymap[8].insert({{-321, 321, 333}, 6});
  decaymap[8].insert({{-211, 111, 111, 211}, 7});
  decaymap[8].insert({{-2212, -211, 22, 211, 2212}, 8});

  // Finish Construction

  if (m_SaveFiles)
  {
    fout = new TFile("MyQAFile.root", "RECREATE");
    QATree = new TTree("QATree", "QATree");
    QATree->Branch("EvtID", &EvtID, "EvtID/I");
    QATree->Branch("NParticles", &NParticles, "NParticles/I");
    QATree->Branch("PVDaughtersPDGID", &PVDaughtersPDGID);

    QATree->Branch("LifeTime", &LifeTime, "LifeTime/F");
  }

  return 0;
}

int QAG4Decayer::process_event(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  // Individual QA for every single heavy flavor particle//

  TH1F *QAPx[NHFQA];
  TH1F *QAPy[NHFQA];
  TH1F *QAPz[NHFQA];
  TH1F *QAE[NHFQA];
  TH1F *QACosTheta[NHFQA];
  TH1F *BR1DHis[NHFQA];
  TH1F *ProperLifeTime[NHFQA];
  TH1F *InvMass[NHFQA];

  for (int i = 0; i < NHFQA; i++)
  {
    QAPx[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("QAPx_%d", i)));
    assert(QAPx[i]);

    QAPy[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("QAPy_%d", i)));
    assert(QAPy[i]);

    QAPz[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("QAPz_%d", i)));
    assert(QAPz[i]);

    QAE[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("QAE_%d", i)));
    assert(QAE[i]);

    QACosTheta[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("QAPx_%d", i)));
    assert(QAPx[i]);

    QAPx[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("QAPx_%d", i)));
    assert(QAPx[i]);

    BR1DHis[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("BR1DHis_%d", i)));
    assert(BR1DHis[i]);

    InvMass[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("InvMass_%d", i)));
    assert(InvMass[i]);

    ProperLifeTime[i] = dynamic_cast<TH1F *>(hm->getHisto(Form("ProperLifeTime_%d", i)));
    assert(ProperLifeTime[i]);
  }

  TH1F *HFHadronStat = dynamic_cast<TH1F *>(hm->getHisto("HFHadronStat"));
  assert(HFHadronStat);

  m_truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  bool VtxToQA = false;

  float SVtoPVDis = 0;
  float SVtoPVTau = 0;

  float DevPx;
  float DevPy;
  float DevPz;
  float DevE;

  std::vector<int> ParentTrkInfo;
  std::vector<double> ParentEInfo;
  std::vector<double> DiffEPerVertex;

  std::vector<double> ParentPxInfo;
  std::vector<double> DiffPxPerVertex;

  std::vector<double> ParentPyInfo;
  std::vector<double> DiffPyPerVertex;

  std::vector<double> ParentPzInfo;
  std::vector<double> DiffPzPerVertex;

  std::vector<std::vector<int>> DaughterInfo;
  std::vector<std::vector<double>> DaughterEInfo;
  std::vector<int> VertexInfo;
  std::vector<int> HFIndexInfo;

  float CosTheta = -2;

  PHG4TruthInfoContainer::ConstRange range = m_truth_info->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second; ++iter)
  {
    PHG4Particle *g4particle = iter->second;

    int gflavor = g4particle->get_pid();

    int ParentPDGID = -1;
    int GrandParentPDGID = -1;
    int ParentTrkId = 0;

    double ParentE = 0;
    double ParentPx = 0;
    double ParentPy = 0;
    double ParentPz = 0;

    PHG4Particle *mother = nullptr;

    TVector3 HFParMom(0, 0, 0);
    TVector3 HFProdVtx(0, 0, 0);
    TVector3 HFDecayVtx(0, 0, 0);

    TLorentzVector HFParFourMom(0, 0, 0, 0);

    if (g4particle->get_parent_id() == 0)
    {
      ParentPDGID = 0;
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

    VtxToQA = false;

    for (int p = 0; p < NHFQA; p++)
    {
      if (abs(ParentPDGID) == QAVtxPDGID[p])
      {
        VtxToQA = true;
        HFIndex = p;
      }
    }

    if ((ParentTrkId > 0 || abs(gflavor) == abs(ParentPDGID)) && VtxToQA == true)
    {
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

        // Check Energy
        std::vector<double> DaughtersE;
        DaughtersE.push_back(g4particle->get_e());
        DaughterEInfo.push_back(DaughtersE);
      }
      if (!NewVtx)
      {
        DiffEPerVertex[Index] = DiffEPerVertex[Index] - g4particle->get_e();
        DiffPxPerVertex[Index] = DiffPxPerVertex[Index] - g4particle->get_px();
        DiffPyPerVertex[Index] = DiffPyPerVertex[Index] - g4particle->get_py();
        DiffPzPerVertex[Index] = DiffPzPerVertex[Index] - g4particle->get_pz();

        DaughterInfo[Index].push_back(gflavor);
        DaughterEInfo[Index].push_back(g4particle->get_e());
      }
    }

    PHG4VtxPoint *vtx = m_truth_info->GetVtx(g4particle->get_vtx_id());
    PHG4VtxPoint *ParentVtx = nullptr;
    if (mother) ParentVtx = m_truth_info->GetVtx(mother->get_vtx_id());

    if (GrandParentPDGID == 0 && ParentVtx && VtxToQA)
    {
      float ParMass = sqrt(mother->get_e() * mother->get_e() - mother->get_px() * mother->get_px() - mother->get_py() * mother->get_py() - mother->get_pz() * mother->get_pz());
      float ParP = sqrt(mother->get_px() * mother->get_px() + mother->get_py() * mother->get_py() + mother->get_pz() * mother->get_pz());

      HFProdVtx.SetXYZ(ParentVtx->get_x(), ParentVtx->get_y(), ParentVtx->get_z());
      HFDecayVtx.SetXYZ(vtx->get_x(), vtx->get_y(), vtx->get_z());

      SVtoPVDis = (HFDecayVtx - HFProdVtx).Mag();
      SVtoPVTau = SVtoPVDis / ParP * ParMass;

      if (SVtoPVDis > 0) CosTheta = ((HFDecayVtx - HFProdVtx).Dot(HFParMom)) / ((HFDecayVtx - HFProdVtx).Mag() * HFParMom.Mag());

      // QACosTheta->Fill(CosTheta);
      // MassHis->Fill(ParMass);

      if (HFIndex > -1)
      {
        ProperLifeTime[HFIndex]->Fill(SVtoPVTau);
        QACosTheta[HFIndex]->Fill(CosTheta);
        InvMass[HFIndex]->Fill(ParMass);
      }
    }
  }

  // BR Working here
  int VtxSizeFinal = DiffEPerVertex.size();

  for (int q = 0; q < VtxSizeFinal; q++)
  {
    int HFIndexToFill = HFIndexInfo[q];

    DevPx = DiffPxPerVertex[q] / ParentPxInfo[q];
    DevPy = DiffPyPerVertex[q] / ParentPyInfo[q];
    DevPz = DiffPzPerVertex[q] / ParentPzInfo[q];
    DevE = DiffEPerVertex[q] / ParentEInfo[q];

    QAPx[HFIndexToFill]->Fill(DevPx);
    QAPy[HFIndexToFill]->Fill(DevPy);
    QAPz[HFIndexToFill]->Fill(DevPz);
    QAE[HFIndexToFill]->Fill(DevE);

    sort(DaughterInfo[q].begin(), DaughterInfo[q].end());

    if (HFIndexToFill < 0)
    {
      continue;
    }

    int key = -1;
    std::vector<int> ChannelID;

    if (decaymap[HFIndexToFill].find({DaughterInfo[q]}) != decaymap[HFIndexToFill].end()) key = decaymap[HFIndexToFill].find({DaughterInfo[q]})->second;
    ChannelID.push_back(key);

    if (HFIndexToFill == 0)
    {
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -11) != DaughterInfo[q].end()) ChannelID.push_back(8);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 333) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    if (HFIndexToFill == 1)
    {
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 310) != DaughterInfo[q].end()) ChannelID.push_back(7);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 111) != DaughterInfo[q].end()) ChannelID.push_back(8);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -11) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    if (HFIndexToFill == 2)
    {
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 310) != DaughterInfo[q].end()) ChannelID.push_back(7);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 333) != DaughterInfo[q].end()) ChannelID.push_back(8);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -11) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    if (HFIndexToFill == 3)
    {
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -311) != DaughterInfo[q].end()) ChannelID.push_back(7);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 2212) != DaughterInfo[q].end()) ChannelID.push_back(8);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -11) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    if (HFIndexToFill == 4)
    {
      if ((std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 411) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -411) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 421) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -421) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 431) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -431) != DaughterInfo[q].end())) ChannelID.push_back(8);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 443) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    if (HFIndexToFill == 5)
    {
      if ((std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 411) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -411) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 421) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -421) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 431) != DaughterInfo[q].end()) || (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -431) != DaughterInfo[q].end())) ChannelID.push_back(8);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 443) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    if (HFIndexToFill == 6)
    {
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 431) != DaughterInfo[q].end() || std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), -431) != DaughterInfo[q].end()) ChannelID.push_back(8);
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 443) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    if (HFIndexToFill == 8)
    {
      if (std::find(DaughterInfo[q].begin(), DaughterInfo[q].end(), 443) != DaughterInfo[q].end()) ChannelID.push_back(9);
    }

    int ChannelSize = ChannelID.size();

    for (int r = 0; r < ChannelSize; r++)
    {
      BR1DHis[HFIndexToFill]->Fill(ChannelID[r]);
    }
  }

  LifeTime = SVtoPVTau;

  if (m_write_nTuple) QATree->Fill();

  EvtID = EvtID + 1;

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4Decayer::End(PHCompositeNode *topNode)
{
  assert(topNode);

  if (m_SaveFiles)
  {
    fout->cd();
    if (m_write_nTuple) QATree->Write();
    fout->Close();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
