#include <cmath>
#include <TFile.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TCut.h>
#include <Eigen/Dense>

void plot_residuals()
{
  TFile fin("residuals.root");

  TNtuple *ntuple;
  fin.GetObject("ntp_residuals",ntuple);

  float seed_id;
  float layer;
  float dphi;
  float dz;
  float x;
  float y;
  float z;
  float pt;
  float px;
  float py;
  float pz;
  float crossing;
  float isSilicon;
  float isTpc;

  //take seedid count entries() whihc is nclusters) with seedid and switches 
  //count how many entries in ntuple have same seedid possibly thousands

  ntuple->SetBranchAddress("seed_id",&seed_id);
  ntuple->SetBranchAddress("layer",&layer);
  ntuple->SetBranchAddress("dphi",&dphi);
  ntuple->SetBranchAddress("dz",&dz);
  ntuple->SetBranchAddress("x",&x);
  ntuple->SetBranchAddress("y",&y);
  ntuple->SetBranchAddress("z",&z);
  ntuple->SetBranchAddress("pt",&pt);
  ntuple->SetBranchAddress("px",&px);
  ntuple->SetBranchAddress("py",&py);
  ntuple->SetBranchAddress("pz",&pz);
  ntuple->SetBranchAddress("crossing",&crossing);
  ntuple->SetBranchAddress("isSilicon",&isSilicon);
  ntuple->SetBranchAddress("isTpc",&isTpc);


  int entries = ntuple->GetEntries();
  
  TH2D *xy         = new TH2D("xy","cluster x vs. y",5000,-80,80,5000,-80,80);
  TH2D *tpc_xy     = new TH2D("tpc_xy","tpc cluster x vs. y",5000,-80,80,5000,-80,80);
  TH2D *si_xy      = new TH2D("si_xy","si cluster x vs. y",5000,-10,10,5000,-10,10);
  TH2D *dphiLayer  = new TH2D("dphiLayer","dphi vs. layer",5000,-0.02,0.02,5000,0,56);
  TH2D *dzLayer    = new TH2D("dzLayer","dz vs. layer",5000,-1,1,5000,0,56);
  TH2D *resLayer   = new TH2D("resLayer","residual vs. layer",5000,0,1,5000,0,56);
  TH2D *dphiPt     = new TH2D("dphiPt","dphi vs. pt",5000,-0.02,0.02,5000,0,7);
  TH2D *dphiPz     = new TH2D("dphiPz","dphi vs. pz",5000,-0.02,0.02,5000,0,5);
  

  //make seedid array
  //loop over entries
  // if seedid is allready in array add 1
  // if seed id not in array add to array then add 1
  //deltaphi v nhits

  for(int i=0; i<entries; ++i)
    {
      ntuple->GetEntry(i);
     

      float r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
      float residual = sqrt(pow(r*dphi,2)+pow(dz,2));

      resLayer->Fill(residual,layer);
      dphiLayer->Fill(dphi,layer);
      dzLayer->Fill(dz,layer);
      xy->Fill(x,y);
      dphiPt->Fill(dphi,pt);
      dphiPz->Fill(dphi,pz);

      if(isSilicon==0)
	{
	  tpc_xy->Fill(x,y);
	}
      if(isTpc==0)
	{ 
	  si_xy->Fill(x,y);
	}
    } 

  TCanvas *c1 = new TCanvas("c1","",10,10,800,800);
  xy->DrawCopy();


  TCanvas *c2 = new TCanvas("c2","",10,10,800,800);
  tpc_xy->DrawCopy();

  TCanvas *c3 = new TCanvas("c3","",10,10,600,600); 
  si_xy->DrawCopy();

  TCanvas *c4 = new TCanvas("c4","",10,10,800,800); 
  dphiLayer->DrawCopy();

  TCanvas *c5 = new TCanvas("c5","",10,10,800,800); 
  dzLayer->DrawCopy();

  TCanvas *c6 = new TCanvas("c6","",10,10,800,800); 
  resLayer->DrawCopy();


  TCanvas *c7 = new TCanvas("c7","",10,10,800,800); 
  dphiPt->DrawCopy();


  TCanvas *c8 = new TCanvas("c8","",10,10,800,800); 
  dphiPz->DrawCopy();

 


}
