//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in LiteCaloEval.h.
//
// LiteCaloEval(const std::string &name = "LiteCaloEval")
// everything is keyed to LiteCaloEval, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// LiteCaloEval::~LiteCaloEval()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int LiteCaloEval::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int LiteCaloEval::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int LiteCaloEval::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int LiteCaloEval::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int LiteCaloEval::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int LiteCaloEval::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int LiteCaloEval::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void LiteCaloEval::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "LiteCaloEval.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>

#include <TF1.h>
#include <TGraph.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <cstdlib>
#include <iostream>

using namespace std;
//____________________________________________________________________________..
LiteCaloEval::LiteCaloEval(const string &name, const string& caloname, const string& filename):
 SubsysReco(name)
 , user_choice(0)
 , _caloname(caloname)
 , _filename(filename)
{

//_caloname = caloname;
  cout << "LiteCaloEval::LiteCaloEval(const string &name) Calling ctor" << endl;

}

//____________________________________________________________________________..
LiteCaloEval::~LiteCaloEval()
{
  cout << "LiteCaloEval::~LiteCaloEval() Calling dtor" << endl;
}

//____________________________________________________________________________..
int LiteCaloEval::Init(PHCompositeNode *topNode)
{
  cout << "LiteCaloEval::Init(PHCompositeNode *topNode) Initializing" << endl;
  
  _ievent = 0;

  cal_output = new TFile(_filename.c_str(),"create");
  
  if(user_choice == 'i' || user_choice == 'I')
   {
     
     hcalin_energy_eta = new TH2F("hcalin_energy_eta", "hcalin energy eta", 1000,0,10,240,-1.1,1.1);
     hcalin_e_eta_phi = new TH3F("hcalin_e_eta_phi", "hcalin e eta phi",50,0,10,24,-1.1,1.1,64,-3.14159,3.14159);
     for (int i = 0; i<24; i++)
       {
	 for (int j = 0; j<64; j++)
	   {
	     TString i1;
	     TString j1;
	     i1.Form("%d",i);
	     j1.Form("%d",j);
	      TString hist_name = "hcal_in_eta_" + i1 + "_phi_" + j1;
	      
	      hcal_in_eta_phi[i][j] = new TH1F(hist_name.Data(),"Hcal_in_energy",1000,0,10);
	   }
       }

     for(int i = 0; i < 24; i++)
       {
	 TString i1;
	 i1.Form("%d",i);
	 TString hist_name = "hcalin_eta_" + i1;
	 
	 hcalin_eta[i] = new TH1F (hist_name.Data(), "hcalin eta's", 1000,0,10);
	}
     
   }
  else   if(user_choice == 'o' || user_choice == 'O')  
    {
      hcalout_energy_eta = new TH2F("hcalout_energy_eta", "hcalout energy eta", 10,0,10,24000,-1.1,1.1); 
      hcalout_e_eta_phi = new TH3F("hcalout_e_eta_phi", "hcalout e eta phi",50,0,10,24,-1.1,1.1,64,-3.14159,3.14159);
      for(int i = 0; i < 24; i++)
	{
	  for(int j = 0; j < 64; j++)
	    {
	      TString i1;
	      TString j1;
	      i1.Form("%d",i);
	      j1.Form("%d",j);
	      TString hist_name = "hcal_out_eta_" + i1 + "_phi_" + j1;
	     
	      hcal_out_eta_phi[i][j] = new TH1F (hist_name.Data(),"Hcal_out energy",1000,0,10);
	    }
	}
      
      for(int i = 0; i < 24; i++)
	{
	  TString i1;
	  i1.Form("%d",i);
	  TString hist_name = "hcalout_eta_" + i1;
	  
	  hcalout_eta[i] = new TH1F (hist_name.Data(), "hcalout eta's",1000,0,10);
	}
     
    }
  else   if (user_choice == 'e' || user_choice == 'E') 
   {
     
     for (int i = 0; i<96; i++)
       {
	 for (int j = 0; j<258; j++)
	   {
	     TString i1;
	     TString j1;
	     i1.Form("%d",i);
	     j1.Form("%d",j);
	     TString hist_name = "emc_ieta" + i1 + "_phi"+ j1;
	     
	     cemc_hist_eta_phi[i][j] = new TH1F(hist_name.Data(),"Hist_ieta_phi_leaf(e)",1000,0,10);
	   }
       }
     
     
     //create the eta 1d histos
     for (int i = 0; i < 96; i++)
       {
	 gStyle->SetOptFit(1);
	 TString a;
	 a.Form("%d",i);
	 TString b = "eta_" + a;
	 
	 eta_hist[i] = new TH1F (b.Data(),"eta and all phi's",1000,0,10);
	 
       }
     
     
     //make 2d histo
     energy_eta_hist = new TH2F ("energy_eta_hist", "energy eta and all phi",10,0,10,9600,-1,1);
     
     //make 3d histo
     e_eta_phi = new TH3F ("e_eta_phi","e v eta v phi",50,0,10,100,-1,1,256,-3.14159,3.14159);
     
    }
     
  

 //  _tfile = new TFile(_filename.c_str(), "RECREATE");
  //  _ntp_tower = new TNtuple("ntp_tower", "tower => max truth primary",
  //  "event:towerID:ieta:iphi:eta:phi:e:x:y:z:");
			   
 /*                                               "gparticleID:gflavor:gnhits:"
                                               "geta:gphi:ge:gpt:gvx:gvy:gvz:"
                                               "gembed:gedep:"
                                               "efromtruth");
 */

  return Fun4AllReturnCodes::EVENT_OK;
  
}

//____________________________________________________________________________..
int LiteCaloEval::InitRun(PHCompositeNode *topNode)
{
  cout << "LiteCaloEval::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::process_event(PHCompositeNode *topNode)
{
  if (_ievent % 100 == 0)  cout << "LiteCaloEval::process_event(PHCompositeNode *topNode) Processing Event " << _ievent << endl;

  string towernode = "TOWER_CALIB_" + _caloname;
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, towernode.c_str());
    if (!towers)
    {
      cerr << PHWHERE << " ERROR: Can't find " << towernode << endl;
      exit(-1);
    }


    string towergeomnode = "TOWERGEOM_" + _caloname;
    RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode.c_str());
    if (!towergeom)
    {
      cerr << PHWHERE << " ERROR: Can't find " << towergeomnode << endl;
      exit(-1);
    }

    RawTowerContainer::ConstRange begin_end = towers->getTowers();
    RawTowerContainer::ConstIterator rtiter;
    for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
    {
      RawTower* tower = rtiter->second;

      if (tower->get_energy() < 0.0005) continue;

      RawTowerGeom* tower_geom = towergeom->get_tower_geometry(tower->get_id());
      if (!tower_geom)
      {
        cerr << PHWHERE << " ERROR: Can't find tower geometry for this tower hit: ";
        tower->identify();
        exit(-1);
      }

      const float towerid = tower->get_id();
      const float ieta = tower->get_bineta();
      const float iphi = tower->get_binphi();
      const float eta = tower_geom->get_eta();
      const float phi = tower_geom->get_phi();
      const float e = tower->get_energy();
      //      const float x = tower_geom->get_center_x();
      //const float y = tower_geom->get_center_y();
      //const float z = tower_geom->get_center_z();

      /*
      float tower_data[] = {(float) _ievent,
			    towerid,
			    ieta,
			    iphi,
			    eta,
			    phi,
			    e,
			    x,
			    y,
			    z

      };
      _ntp_tower->Fill(tower_data);
      */
     

      if (user_choice == 'e' || user_choice == 'E')
	{
	  //create ii variables for indexing the histos below
	  int iiphi = (int) iphi;
	  int iieta = (int) ieta;

	  if (towerid < 0) 
	    {
	      cout << "a towerid was less than 0 " << endl;
	      break;
	    }
	  
	  //fill the hist with energy data
	  cemc_hist_eta_phi[iieta][iiphi]->Fill(e);
	  	  
	  
	  //fill and fit the 1d hist eta and all phi
	  eta_hist[iieta]->Fill(e);
     
	  
	  //fill the 2d histo eta, energy and all phi
	  energy_eta_hist->Fill(e,eta);
	  
	  //fill 3d histo e_eta_phid
	  e_eta_phi->Fill(e,eta,phi);
	  
	}
      else if (user_choice == 'o' || user_choice == 'O')
	{
	  
	  //create ii variables for indexing the histos below
	  int iiphi = (int) iphi;
	  int iieta = (int) ieta;
	  
	  //fill the hist with energy data
	  //cout << iieta << " " <<  iiphi  << endl;
	  
	  hcal_out_eta_phi[iieta][iiphi]->Fill(e);
	  
	  hcalout_eta[iieta]->Fill(e);
	  
	  hcalout_energy_eta->Fill(e,eta);
	  
	  hcalout_e_eta_phi->Fill(e,eta,phi);
	  
	}
      else   if(user_choice == 'i' || user_choice == 'I')
	{
	  //create ii variables for indexing the histos below
	  int iiphi = (int) iphi;
	  int iieta = (int) ieta;
	  
	  //fill the hist with energy data
	  
	  
	  hcal_in_eta_phi[iieta][iiphi]->Fill(e);
	  
	  hcalin_eta[iieta]->Fill(e);
	  
	  hcalin_energy_eta->Fill(e,eta);
	  
	  hcalin_e_eta_phi->Fill(e,eta,phi);
	  
	}

    }
      
      
    _ievent++;

    return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::ResetEvent(PHCompositeNode *topNode)
{
  //  cout << "LiteCaloEval::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::EndRun(const int runnumber)
{
  cout << "LiteCaloEval::EndRun(const int runnumber) Ending Run for Run " << runnumber << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::End(PHCompositeNode *topNode)
{
  //  _tfile->cd();

  if(user_choice == 'i' || user_choice == 'I')
    {
  
      double par_value [24];
      double eta_value [24];
      double par_err [24];
      double rel_err [24];
      
      //same as above but for the second fit performed below in code
      double par_value2[24];
      double par_err2[24];
      double rel_err2[24];
      
      
      for (int i =0; i < 24; i++)
	{
	  //a,b,c,d hold the par value, error, par2 value, par2 error respectivley. The 2 refers to the second fit below
	  double a;
	  double b;
	  double c;
	  double d;
	  
	  //make functions to fit the eta histos
	  TF1 *f1 = new TF1("f1","expo",0.02,0.1);
	  TF1 *f2 = new TF1("f2","expo",0.1,1);
	  f2->SetLineColor(7);

	  hcalin_eta[i]->Fit("f1","R");
	  hcalin_eta[i]->Fit("f2","R+");
	  
	  a = f1->GetParameter(1);
	  b = f1->GetParError(1);

	  c = f2->GetParameter(1);
	  d = f2->GetParError(1);
       
	  //assign retreived parameter values
	  par_value[i] = abs(a);
	  par_err[i] = b;
	  eta_value[i] = i;
	  rel_err[i] = par_err[i] / par_value[i];

	  par_value2[i] = abs(c);
	  par_err2[i] = d;
	  rel_err2[i] = par_err2[i] / par_value2[i];

	}
   
  
      //create graphs
      TGraph g1(24,eta_value,par_value);
      g1.SetTitle("HCal In (0.02-0.1 GeV); eta; p1");
      g1.SetMarkerStyle(20);
      g1.Draw("ap");
      g1.SetName("Fit1_hcalin");
      g1.Write();
  

      TGraph g2(24,eta_value,rel_err);
      g2.SetTitle("HCal In Error (0.02-0.1 GeV); eta; p1 rel error");
      g2.SetMarkerStyle(20);
      g2.Draw("ap");
      g2.SetName("Fit1_err_hcalin");
      g2.Write();
   
      TGraph g3(24,eta_value,par_value2);
      g3.SetTitle("HCal In (0.1-1 GeV); eta; p1");
      g3.SetMarkerStyle(20);
      g3.Draw("ap");
      g3.SetName("Fit2_hcalin");
      g3.Write();
  

      TGraph g4(24, eta_value, rel_err2);
      g4.SetTitle("HCal In Error (0.1-1 GeV); eta; p1 rel error");
      g4.SetMarkerStyle(20);
      g4.Draw("ap");
      g4.SetName("Fit2_err_hcalin");
      g4.Write();

      
      cal_output->Write();
      
    }
  else if(user_choice == 'o' || user_choice == 'O')
    {

      double par_value [24];
      double eta_value [24];
      double par_err [24];
      double rel_err [24];

      //same as above but for the second/third fit performed below in code
      double par_value2[24];
      double par_err2[24];
      double rel_err2[24];

      double par_value3[24];
      double par_err3[24];
      double rel_err3[24];

      
      for (int i =0; i < 24; i++)
	{
	  //a,b hold the p1 value, p1 error, and then repeats for the second fit (c,d) and third fit (e,f)
	  double a;
	  double b;
	  double c;
	  double d;
	  double e;
	  double f;

	  //make functions to fit the eta histos
	  TF1 *f1 = new TF1("f1","expo",0.05,0.2);
	  TF1 *f2 = new TF1("f2","expo",0.2,1);
	  TF1 *f3 = new TF1("f3","expo",1,2);
	  f2->SetLineColor(7);
	  f3->SetLineColor(1);
	  
	  hcalout_eta[i]->Fit("f1","R");
	  hcalout_eta[i]->Fit("f2","R+");
	  hcalout_eta[i]->Fit("f3","R+");

	  a = f1->GetParameter(1);
	  b = f1->GetParError(1);
	  
	  c = f2->GetParameter(1);
	  d = f2->GetParError(1);

	  e = f3->GetParameter(1);
	  f = f3->GetParError(1);
	  
	  //assign retreived parameter values
	  par_value[i] = abs(a);
	  par_err[i] = b;
	  eta_value[i] = i;
	  rel_err[i] = par_err[i] / par_value[i];
	  
	  par_value2[i] = abs(c);
	  par_err2[i] = d;
	  rel_err2[i] = par_err2[i] / par_value2[i];
	  
	  par_value3[i] = abs(e);
	  par_err3[i] = f;
	  rel_err3[i] = par_err3[i] / par_value3[i];
	  
	}
      
      
      //create graphs
      TGraph g1(24,eta_value,par_value);
      g1.SetTitle("HCal Out (0.05-0.2 GeV); eta; p1");
      g1.SetMarkerStyle(20);
      g1.Draw("ap");
      g1.SetName("Fit1_hcalout");
      g1.Write();
  

      TGraph g2(24,eta_value,rel_err);
      g2.SetTitle("HCal Out Error (0.05-0.2 GeV); eta; p1 rel err");
      g2.SetMarkerStyle(20);
      g2.Draw("ap");
      g2.SetName("Fit1_err_hcalout");
      g2.Write();
   
      TGraph g3(24,eta_value,par_value2);
      g3.SetTitle("HCal Out (0.2-1 GeV); eta; p1");
      g3.SetMarkerStyle(20);
      g3.Draw("ap");
      g3.SetName("Fit2_hcalout");
      g3.Write();
  

      TGraph g4(24, eta_value, rel_err2);
      g4.SetTitle("HCal Out Error (0.2-1 GeV); eta; p1 rel err");
      g4.SetMarkerStyle(20);
      g4.Draw("ap");
      g4.SetName("Fit2_err_hcalout");
      g4.Write();

      TGraph g5(24, eta_value, par_value3);
      g5.SetTitle("HCal Out (1-2 GeV); eta; p1");
      g5.SetMarkerStyle(20);
      g5.Draw("ap");
      g5.SetName("Fit3_hcalout");
      g5.Write();

      TGraph g6(24,eta_value, rel_err3);
      g6.SetTitle("HCal Out Error (1-2 GeV); eta; p1 rel err");
      g6.SetMarkerStyle(20);
      g6.Draw("ap");
      g6.SetName("Fit3_err_hcalout");
      g6.Write();


      cal_output->Write();

      
    }
  else if(user_choice == 'e' || user_choice == 'E')
    {

      //create arrays that holds parameter (p1) value and error
      double par_value [96];
      double eta_value [96];
      double par_err [96];
      double rel_err [96];
      
      //same as above but for the second fit performed below in code
      double par_value2[96];
      double par_err2[96];
      double rel_err2[96];
      
      double par_value3[96];
      double par_err3[96];
      double rel_err3[96];
      
      
            //create graphs for parameter (p1) vs eta
      for (int i =0; i < 96; i++)
	{
	  //a,b,c,d hold the par value, error, par2 value, par2 error respectivley. The 2 refers to the second fit below
	  double a;
	  double b;
	  double c;
	  double d;
	  double e;
	  double f;

	  //make functions to fit the eta histos
	  TF1 *f1 = new TF1("f1","expo",0.04,0.1);
	  TF1 *f2 = new TF1("f2","expo",0.1,0.4);
	  TF1 *f3 = new TF1("f3","expo",0.4,2);
	  f2->SetLineColor(7);
	  f3->SetLineColor(1);

	  eta_hist[i]->Fit("f1","R");
	  eta_hist[i]->Fit("f2","R+");
	  eta_hist[i]->Fit("f3","R+");

	  a = f1->GetParameter(1);
	  b = f1->GetParError(1);

	  c = f2->GetParameter(1);
	  d = f2->GetParError(1);

	  e = f3->GetParameter(1);
	  f = f3->GetParError(1);
       
	  //assign retreived parameter values
	  par_value[i] = abs(a);
	  par_err[i] = b;
	  eta_value[i] = i;
	  rel_err[i] = par_err[i] / par_value[i];

	  par_value2[i] = abs(c);
	  par_err2[i] = d;
	  rel_err2[i] = par_err2[i] / par_value2[i];

	  par_value3[i] = abs(e);
	  par_err3[i] = f;
	  rel_err3[i] = par_err3[i] / par_value3[i];

	}
   
  
      //create graphs
      TGraph g1(96,eta_value,par_value);
      g1.SetTitle("EMCal (0.04-0.1 GeV); eta; p1");
      g1.SetMarkerStyle(20);
      g1.Draw("ap");
      g1.SetName("Fit1_emc");
      g1.Write();
  

      TGraph g2(96,eta_value,rel_err);
      g2.SetTitle("EMCal Error (0.04-0.1 GeV); eta; p1 rel error");
      g2.SetMarkerStyle(20);
      g2.Draw("ap");
      g2.SetName("Fit1_err_emc");
      g2.Write();
   
      TGraph g3(96,eta_value,par_value2);
      g3.SetTitle("EMCal (0.1-0.4 GeV); eta; p1");
      g3.SetMarkerStyle(20);
      g3.Draw("ap");
      g3.SetName("Fit2_emc");
      g3.Write();
  

      TGraph g4(96, eta_value, rel_err2);
      g4.SetTitle("EMCal Error (0.1-0.4 GeV); eta; p1 rel error");
      g4.SetMarkerStyle(20);
      g4.Draw("ap");
      g4.SetName("Fit2_err_emc");
      g4.Write();

      TGraph g5(96, eta_value, par_value3);
      g5.SetTitle("EMCal (0.4-2 GeV); eta; p1");
      g5.SetMarkerStyle(20);
      g5.Draw("ap");
      g5.SetName("Fit3_emc");
      g5.Write();

      TGraph g6(96,eta_value,rel_err3);
      g6.SetTitle("EMCal Error (0.4-2 GeV); eta; p1 rel err");
      g6.SetMarkerStyle(20);
      g6.Draw("ap");
      g6.SetName("Fit3_err_emc");
      g6.Write();

      cal_output->Write();


      
    }



  //  _ntp_tower->Write();
  cout << "LiteCaloEval::End(PHCompositeNode *topNode) This is the End..." << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::Reset(PHCompositeNode *topNode)
{
  // cout << "LiteCaloEval::Reset(PHCompositeNode *topNode) being Reset" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void LiteCaloEval::Print(const string &what) const
{
  cout << "LiteCaloEval::Print(const string &what) const Printing info for " << what << endl;
}
