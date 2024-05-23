//____________________________________________________________________________..


#include "JetKinematicCheck.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetContainer.h>
#include <jetbase/Jetv2.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TPad.h>
#include <TLegend.h>
#include <cmath>
#include <string>
#include <vector>


//____________________________________________________________________________..

JetKinematicCheck::JetKinematicCheck(const std::string& recojetnameR02, const std::string& recojetnameR03, const std::string& recojetnameR04):
SubsysReco("JetKinematicCheck")
  , m_recoJetNameR02(recojetnameR02)
  , m_recoJetNameR03(recojetnameR03)
  , m_recoJetNameR04(recojetnameR04)
  , m_etaRange(-1.1, 1.1)
  , m_ptRange(10, 100)


{
  std::cout << "JetKinematicCheck::JetKinematicCheck(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
JetKinematicCheck::~JetKinematicCheck()
{
  std::cout << "JetKinematicCheck::~JetKinematicCheck() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int JetKinematicCheck::Init(PHCompositeNode *topNode)
{


  hm = QAHistManagerDef::getHistoManager();
  assert(hm);


  // initialize histograms

  jet_spectra_r02 = new TH1D("h_spectra_r02", "", 19, 5, 100);
  jet_spectra_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  jet_spectra_r03 = new TH1D("h_spectra_r03", "", 19, 5, 100);
  jet_spectra_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  jet_spectra_r04 = new TH1D("h_spectra_r04", "", 19, 5, 100);
  jet_spectra_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  jet_eta_phi_r02 = new TH2D("h_eta_phi_r02","", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  jet_eta_phi_r02->GetXaxis()->SetTitle("#eta");
  jet_eta_phi_r02->GetYaxis()->SetTitle("#Phi");

  jet_eta_phi_r03 = new TH2D("h_eta_phi_r03","", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  jet_eta_phi_r03->GetXaxis()->SetTitle("#eta");
  jet_eta_phi_r03->GetYaxis()->SetTitle("#Phi");

  jet_eta_phi_r04 = new TH2D("h_eta_phi_r04","", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  jet_eta_phi_r04->GetXaxis()->SetTitle("#eta");
  jet_eta_phi_r04->GetYaxis()->SetTitle("#Phi");

  jet_mass_pt_r02 = new TH2D("h_jet_mass_pt_r02","", 19, 5, 100, 15, 0, 15);
  jet_mass_pt_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_r02->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_pt_r03 = new TH2D("h_jet_mass_pt_r03","", 19, 5, 100, 15, 0, 15);
  jet_mass_pt_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_r03->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_pt_r04 = new TH2D("h_jet_mass_pt_r04","", 19, 5, 100, 15, 0, 15);
  jet_mass_pt_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_r04->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");


  jet_mass_eta_r02 = new TH2D("h_jet_mass_eta_r02","", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_r02->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_r02->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_eta_r03 = new TH2D("h_jet_mass_eta_r03","", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_r03->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_r03->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_eta_r04 = new TH2D("h_jet_mass_eta_r04","", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_r04->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_r04->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");


  jet_mass_pt_1D_r02 = new TH1D("h_jet_mass_pt_1D_r02", "", 19, 5, 100);
  jet_mass_pt_1D_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_pt_1D_r03 = new TH1D("h_jet_mass_pt_1D_r03", "", 19, 5, 100);
  jet_mass_pt_1D_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_pt_1D_r04 = new TH1D("h_jet_mass_pt_1D_r04", "", 19, 5, 100);
  jet_mass_pt_1D_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");


  jet_mass_eta_1D_r02 = new TH2D("h_jet_mass_eta_1D_r02","", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_1D_r02->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_eta_1D_r03 = new TH2D("h_jet_mass_eta_1D_r03","", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_1D_r03->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_eta_1D_r04 = new TH2D("h_jet_mass_eta_1D_r04","", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_1D_r04->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");




  std::cout << "JetKinematicCheck::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int JetKinematicCheck::InitRun(PHCompositeNode *topNode)
{
  std::cout << "JetKinematicCheck::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}




//____________________________________________________________________________..
int JetKinematicCheck::process_event(PHCompositeNode *topNode)
{


  std::vector<std::string> m_recoJetName_array = {m_recoJetNameR02, m_recoJetNameR03, m_recoJetNameR04};
  m_radii = {0.2, 0.3, 0.4};
  int n_radii = m_radii.size();


  // Loop over each reco jet radii from array
  for(int i = 0; i < n_radii; i++){

    std::string recoJetName = m_recoJetName_array[i];

    JetContainer* jets = findNode::getClass<JetContainer>(topNode, recoJetName);
    if (!jets)
      {
  	std::cout
  	  << "MyJetAnalysis::process_event - Error can not find DST Reco JetContainer node "
  	  << recoJetName << std::endl;
  	exit(-1);
      }


    //loop over jets
    for (auto jet : *jets)
      {
  	bool eta_cut = (jet->get_eta() >= m_etaRange.first) and (jet->get_eta() <= m_etaRange.second);
  	bool pt_cut = (jet->get_pt() >= m_ptRange.first) and (jet->get_pt() <= m_ptRange.second);
  	if ((not eta_cut) or (not pt_cut)) continue;
  	if(jet->get_pt() < 1) continue; // to remove noise jets

        if(i == 0){
          jet_spectra_r02->Fill(jet->get_pt());
          jet_eta_phi_r02->Fill(jet->get_eta(), jet->get_phi());
	  jet_mass_pt_r02->Fill(jet->get_pt(), jet->get_mass());
	  jet_mass_eta_r02->Fill(jet->get_eta(), jet->get_mass());
        }

        else if(i == 1){
	  jet_spectra_r03->Fill(jet->get_pt());
  	  jet_eta_phi_r03->Fill(jet->get_eta(), jet->get_phi());
	  jet_mass_pt_r03->Fill(jet->get_pt(), jet->get_mass());
	  jet_mass_eta_r03->Fill(jet->get_eta(), jet->get_mass());

  	}

  	else if(i == 2){
  	  jet_spectra_r04->Fill(jet->get_pt());
  	  jet_eta_phi_r04->Fill(jet->get_eta(), jet->get_phi());
	  jet_mass_pt_r04->Fill(jet->get_pt(), jet->get_mass());
	  jet_mass_eta_r04->Fill(jet->get_eta(), jet->get_mass());
  	}
      }
  }



  std::cout << "JetKinematicCheck::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}



//____________________________________________________________________________..
int JetKinematicCheck::EndRun(const int runnumber)
{
  std::cout << "JetKinematicCheck::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetKinematicCheck::End(PHCompositeNode *topNode)
{

  std::cout << "JetKinematicCheck::End(PHCompositeNode *topNode) Entering the end" << std::endl;



  //for jet spectra [R02]
  TLegend *leg1 = new TLegend(.7,.9,.9,1);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.06);
  leg1->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg1->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_spectra_r02->SetMarkerStyle(8);
  jet_spectra_r02->SetMarkerColor(1);
  jet_spectra_r02->SetLineColor(1);
  jet_spectra_r02->SetTitle("Jet Spectra [R = 0.2]");
  jet_spectra_r02->GetYaxis()->SetTitle("Counts");
  jet_spectra_r02->GetListOfFunctions()->Add(leg1);
  jet_spectra_r02->SetStats(0);


  //for jet spectra [R03]
  TLegend *leg2 = new TLegend(.7,.9,.9,1);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.06);
  leg2->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg2->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_spectra_r03->SetMarkerStyle(8);
  jet_spectra_r03->SetMarkerColor(1);
  jet_spectra_r03->SetLineColor(1);
  jet_spectra_r03->SetTitle("Jet Spectra [R = 0.3]");
  jet_spectra_r03->GetYaxis()->SetTitle("Counts");
  jet_spectra_r03->GetListOfFunctions()->Add(leg2);
  jet_spectra_r03->SetStats(0);

  //for jet spectra [R04]
  TLegend *leg3 = new TLegend(.7,.9,.9,1);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.06);
  leg3->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg3->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_spectra_r04->SetMarkerStyle(8);
  jet_spectra_r04->SetMarkerColor(1);
  jet_spectra_r04->SetLineColor(1);
  jet_spectra_r04->SetTitle("Jet Spectra [R = 0.4]");
  jet_spectra_r04->GetYaxis()->SetTitle("Counts");
  jet_spectra_r04->GetListOfFunctions()->Add(leg3);
  jet_spectra_r04->SetStats(0);


  //for jet eta-phi [R02]
  TLegend *leg4 = new TLegend(.7,.9,.9,1);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.06);
  leg4->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg4->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_eta_phi_r02->SetStats(0);
  jet_eta_phi_r02->SetTitle("Jet Eta-Phi [R = 0.2]");
  jet_eta_phi_r02->GetListOfFunctions()->Add(leg4);


  //for jet eta-phi [R03]
  TLegend *leg5 = new TLegend(.7,.9,.9,1);
  leg5->SetFillStyle(0);
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.06);
  leg5->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg5->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_eta_phi_r03->SetStats(0);
  jet_eta_phi_r03->SetTitle("Jet Eta-Phi [R = 0.3]");
  jet_eta_phi_r03->GetListOfFunctions()->Add(leg5);

  //for jet eta-phi [R04]
  TLegend *leg6 = new TLegend(.7,.9,.9,1);
  leg6->SetFillStyle(0);
  leg6->SetBorderSize(0);
  leg6->SetTextSize(0.06);
  leg6->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg6->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_eta_phi_r04->SetStats(0);
  jet_eta_phi_r04->SetTitle("Jet Eta-Phi [R = 0.4]");
  jet_eta_phi_r04->GetListOfFunctions()->Add(leg6);

  //for jet mass vs pt [R02]
  TLegend *leg7 = new TLegend(.7,.9,.9,1);
  leg7->SetFillStyle(0);
  leg7->SetBorderSize(0);
  leg7->SetTextSize(0.06);
  leg7->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg7->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_pt_r02->SetStats(0);
  jet_mass_pt_r02->SetTitle("Jet Mass vs p_{T} [R = 0.2]");
  jet_mass_pt_r02->GetListOfFunctions()->Add(leg7);

  //for average jet mass vs pt [R02]
  TLegend *leg8 = new TLegend(.7,.9,.9,1);
  leg8->SetFillStyle(0);
  leg8->SetBorderSize(0);
  leg8->SetTextSize(0.06);
  leg8->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg8->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_pt_1D_r02 = (TH1D*)jet_mass_pt_r02->ProfileX();
  jet_mass_pt_1D_r02->SetStats(0);
  jet_mass_pt_1D_r02->SetTitle("Average Jet Mass vs p_{T} [R = 0.2]");
  jet_mass_pt_1D_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_pt_1D_r02->SetMarkerStyle(8);
  jet_mass_pt_1D_r02->SetMarkerColor(1);
  jet_mass_pt_1D_r02->SetLineColor(1);
  jet_mass_pt_1D_r02->GetListOfFunctions()->Add(leg8);


  //for jet mass vs pt [R03]
  TLegend *leg9 = new TLegend(.7,.9,.9,1);
  leg9->SetFillStyle(0);
  leg9->SetBorderSize(0);
  leg9->SetTextSize(0.06);
  leg9->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg9->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_pt_r03->SetStats(0);
  jet_mass_pt_r03->SetTitle("Jet Mass vs p_{T} [R = 0.3]");
  jet_mass_pt_r03->GetListOfFunctions()->Add(leg9);


  //for average jet mass vs pt [R03]
  TLegend *leg10 = new TLegend(.7,.9,.9,1);
  leg10->SetFillStyle(0);
  leg10->SetBorderSize(0);
  leg10->SetTextSize(0.06);
  leg10->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg10->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_pt_1D_r03 = (TH1D*)jet_mass_pt_r03->ProfileX();
  jet_mass_pt_1D_r03->SetStats(0);
  jet_mass_pt_1D_r03->SetTitle("Average Jet Mass vs p_{T} [R = 0.3]");
  jet_mass_pt_1D_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_pt_1D_r03->SetMarkerStyle(8);
  jet_mass_pt_1D_r03->SetMarkerColor(1);
  jet_mass_pt_1D_r03->SetLineColor(1);
  jet_mass_pt_1D_r03->GetListOfFunctions()->Add(leg10);


  //for jet mass vs pt [R04]
  TLegend *leg11 = new TLegend(.7,.9,.9,1);
  leg11->SetFillStyle(0);
  leg11->SetBorderSize(0);
  leg11->SetTextSize(0.06);
  leg11->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg11->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_pt_r04->SetStats(0);
  jet_mass_pt_r04->SetTitle("Jet Mass vs p_{T} [R = 0.4]");
  jet_mass_pt_r04->GetListOfFunctions()->Add(leg11);

  //for average jet mass vs pt [R04]
  TLegend *leg12 = new TLegend(.7,.9,.9,1);
  leg12->SetFillStyle(0);
  leg12->SetBorderSize(0);
  leg12->SetTextSize(0.06);
  leg12->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg12->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_pt_1D_r04 = (TH1D*)jet_mass_pt_r04->ProfileX();
  jet_mass_pt_1D_r04->SetStats(0);
  jet_mass_pt_1D_r04->SetTitle("Average Jet Mass vs p_{T} [R = 0.4]");
  jet_mass_pt_1D_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_pt_1D_r04->SetMarkerStyle(8);
  jet_mass_pt_1D_r04->SetMarkerColor(1);
  jet_mass_pt_1D_r04->SetLineColor(1);
  jet_mass_pt_1D_r04->GetListOfFunctions()->Add(leg12);


  //for jet mass vs eta [R02]
  TLegend *leg13 = new TLegend(.7,.9,.9,1);
  leg13->SetFillStyle(0);
  leg13->SetBorderSize(0);
  leg13->SetTextSize(0.06);
  leg13->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg13->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_eta_r02->SetStats(0);
  jet_mass_eta_r02->SetTitle("Jet Mass vs #eta [R = 0.2]");
  jet_mass_eta_r02->GetListOfFunctions()->Add(leg13);

  //for average jet mass vs eta [R02]
  TLegend *leg14 = new TLegend(.7,.9,.9,1);
  leg14->SetFillStyle(0);
  leg14->SetBorderSize(0);
  leg14->SetTextSize(0.06);
  leg14->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg14->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_eta_1D_r02 = (TH1D*)jet_mass_eta_r02->ProfileX();
  jet_mass_eta_1D_r02->SetStats(0);
  jet_mass_eta_1D_r02->SetTitle("Average Jet Mass vs #eta [R = 0.2]");
  jet_mass_eta_1D_r02->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_eta_1D_r02->SetMarkerStyle(8);
  jet_mass_eta_1D_r02->SetMarkerColor(1);
  jet_mass_eta_1D_r02->SetLineColor(1);
  jet_mass_eta_1D_r02->GetListOfFunctions()->Add(leg14);

  //for jet mass vs eta [R03]
  TLegend *leg15 = new TLegend(.7,.9,.9,1);
  leg15->SetFillStyle(0);
  leg15->SetBorderSize(0);
  leg15->SetTextSize(0.06);
  leg15->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg15->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_eta_r03->SetStats(0);
  jet_mass_eta_r03->SetTitle("Jet Mass vs #eta [R = 0.3]");
  jet_mass_eta_r03->GetListOfFunctions()->Add(leg15);

  //for average jet mass vs eta [R03]
  TLegend *leg16 = new TLegend(.7,.9,.9,1);
  leg16->SetFillStyle(0);
  leg16->SetBorderSize(0);
  leg16->SetTextSize(0.06);
  leg16->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg16->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_eta_1D_r03 = (TH1D*)jet_mass_eta_r03->ProfileX();
  jet_mass_eta_1D_r03->SetStats(0);
  jet_mass_eta_1D_r03->SetTitle("Average Jet Mass vs #eta [R = 0.3]");
  jet_mass_eta_1D_r03->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_eta_1D_r03->SetMarkerStyle(8);
  jet_mass_eta_1D_r03->SetMarkerColor(1);
  jet_mass_eta_1D_r03->SetLineColor(1);
  jet_mass_eta_1D_r03->GetListOfFunctions()->Add(leg16);


  //for jet mass vs eta [R04]
  TLegend *leg17 = new TLegend(.7,.9,.9,1);
  leg17->SetFillStyle(0);
  leg17->SetBorderSize(0);
  leg17->SetTextSize(0.06);
  leg17->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg17->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_eta_r04->SetStats(0);
  jet_mass_eta_r04->SetTitle("Jet Mass vs #eta [R = 0.4]");
  jet_mass_eta_r04->GetListOfFunctions()->Add(leg17);

  //for average jet mass vs eta [R04]
  TLegend *leg18 = new TLegend(.7,.9,.9,1);
  leg18->SetFillStyle(0);
  leg18->SetBorderSize(0);
  leg18->SetTextSize(0.06);
  leg18->AddEntry((TObject*)0, Form("%2.0f < p_{T} < %2.0f [GeV/c]", m_ptRange.first, m_ptRange.second),"");
  leg18->AddEntry((TObject*)0, Form("%1.1f < #eta < %1.1f", m_etaRange.first, m_etaRange.second),"");
  jet_mass_eta_1D_r04 = (TH1D*)jet_mass_eta_r04->ProfileX();
  jet_mass_eta_1D_r04->SetStats(0);
  jet_mass_eta_1D_r04->SetTitle("Average Jet Mass vs #eta [R = 0.4]");
  jet_mass_eta_1D_r04->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_eta_1D_r04->SetMarkerStyle(8);
  jet_mass_eta_1D_r04->SetMarkerColor(1);
  jet_mass_eta_1D_r04->SetLineColor(1);
  jet_mass_eta_1D_r04->GetListOfFunctions()->Add(leg18);



  hm->registerHisto(jet_spectra_r02);
  hm->registerHisto(jet_spectra_r03);
  hm->registerHisto(jet_spectra_r04);
  hm->registerHisto(jet_eta_phi_r02);
  hm->registerHisto(jet_eta_phi_r03);
  hm->registerHisto(jet_eta_phi_r04);
  hm->registerHisto(jet_mass_pt_r02);
  hm->registerHisto(jet_mass_pt_r03);
  hm->registerHisto(jet_mass_pt_r04);
  hm->registerHisto(jet_mass_pt_1D_r02);
  hm->registerHisto(jet_mass_pt_1D_r03);
  hm->registerHisto(jet_mass_pt_1D_r04);
  hm->registerHisto(jet_mass_eta_r02);
  hm->registerHisto(jet_mass_eta_r03);
  hm->registerHisto(jet_mass_eta_r04);
  hm->registerHisto(jet_mass_eta_1D_r02);
  hm->registerHisto(jet_mass_eta_1D_r03);
  hm->registerHisto(jet_mass_eta_1D_r04);




  std::cout << "JetKinematicCheck::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetKinematicCheck::Reset(PHCompositeNode *topNode)
{
 std::cout << "JetKinematicCheck::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void JetKinematicCheck::Print(const std::string &what) const
{
  std::cout << "JetKinematicCheck::Print(const std::string &what) const Printing info for " << what << std::endl;
}
