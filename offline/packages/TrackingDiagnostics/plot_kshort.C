#include <cmath>
#include <TFile.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TCut.h>
#include <Eigen/Dense>

void plot_kshort()
{
  TFile fin("unidentified_eval_output/combined_ntp_mass_out.root");
  //TFile fin("hijing_eval_output/combined_ntp_mass_out.root");
  //TFile fin("old_combined_ntp_mass_out/feb15.root");
  

  TNtuple *ntuple;
  fin.GetObject("ntp_reco_info",ntuple);


  /* Below are values in massRecoAnalysis ntuple 3/6/23
"ntp_reco_info","decay_pairs","x1:y1:z1:px1:py1:pz1:dca3dxy1:dca3dz1:phi1:pca_rel1_x:pca_rel1_y:pca_rel1_z:eta1:charge1:tpcClusters_1:x2:y2:z2:px2:py2:pz2:dca3dxy2:dca3dz2:phi2:pca_rel2_x:pca_rel2_y:pca_rel2_z:eta2:charge2:tpcClusters_2:vertex_x:vertex_y:vertex_z:pair_dca:invariant_mass:invariant_pt:pathlength_x:pathlength_y:pathlength_z:pathlength:rapidity:pseudorapidity:projected_pos1_x:projected_pos1_y:projected_pos1_z:projected_pos2_x:projected_pos2_y:projected_pos2_z:projected_mom1_x:projected_mom1_y:projected_mom1_z:projected_mom2_x:projected_mom2_y:projected_mom2_z:projected_pca_rel1_x:projected_pca_rel1_y:projected_pca_rel1_z:projected_pca_rel2_x:projected_pca_rel2_y:projected_pca_rel2_z:projected_pair_dca:projected_pathlength_x:projected_pathlength_y:projected_pathlength_z:projected_pathlength:quality1:quality2:cosThetaReco"
  */



  float invariant_mass;
  float projected_pathlength_x;
  float projected_pathlength_y;
  float projected_pathlength_z;
  float projected_pathlength;
  float vertex_x;
  float vertex_y;
  float vertex_z;
  float projected_mom1_x;
  float projected_mom1_y;
  float projected_mom1_z;
  float projected_mom2_x;
  float projected_mom2_y;
  float projected_mom2_z;
  float projected_pair_dca;
  float invariant_pt;
  float quality1;
  float quality2;

  float dca3dxy1;
  float dca3dxy2;
  float dca3dz1;
  float dca3dz2;


  ntuple->SetBranchAddress("invariant_mass",&invariant_mass);
  ntuple->SetBranchAddress("projected_pathlength_x",&projected_pathlength_x);
  ntuple->SetBranchAddress("projected_pathlength_y",&projected_pathlength_y);
  ntuple->SetBranchAddress("projected_pathlength_z",&projected_pathlength_z);
  ntuple->SetBranchAddress("projected_pathlength",&projected_pathlength);
  ntuple->SetBranchAddress("vertex_x",&vertex_x);
  ntuple->SetBranchAddress("vertex_y",&vertex_y);
  ntuple->SetBranchAddress("vertex_z",&vertex_z);
  ntuple->SetBranchAddress("projected_mom1_x",&projected_mom1_x);
  ntuple->SetBranchAddress("projected_mom1_y",&projected_mom1_y);
  ntuple->SetBranchAddress("projected_mom1_z",&projected_mom1_z);
  ntuple->SetBranchAddress("projected_mom2_x",&projected_mom2_x);
  ntuple->SetBranchAddress("projected_mom2_y",&projected_mom2_y);
  ntuple->SetBranchAddress("projected_mom2_z",&projected_mom2_z);
  ntuple->SetBranchAddress("projected_pair_dca",&projected_pair_dca);
  ntuple->SetBranchAddress("invariant_pt",&invariant_pt);
  ntuple->SetBranchAddress("quality1",&quality1);
  ntuple->SetBranchAddress("quality2",&quality2);
  ntuple->SetBranchAddress("dca3dxy1",&dca3dxy1);
  ntuple->SetBranchAddress("dca3dxy2",&dca3dxy2);
  ntuple->SetBranchAddress("dca3dz1",&dca3dz1);
  ntuple->SetBranchAddress("dca3dz2",&dca3dz2);



  //TH2D *decay_length = new TH2D("decay_length","",5000,-1,1,5000,0,10);
  //kfparticle comaprison cuts
  //TCut cut = "K0_mass>0.1&&track_1_IP_xy>0.02&&track_2_IP_xy>0.02&&abs(track_1_track_2_DCA)<0.05&&K0_pT>0.1&&abs(K0_x)>0.1&&abs(K0_y)>0.1"; 
  //TCut cut = "abs(track_1_IP_xy)>0.02&&abs(track_2_IP_xy)>0.02&&abs(K0_decayLength)>0.2";
  //TCut cut = "invariant_mass>0.1&&abs(dca3dxy1)>0.02&&abs(dca3dxy2)>0.02&&abs(dca3dz1)>0.02&&abs(dca3dz2)>0.02&&abs(pair_dca)<0.05&&invariant_pt>0.1&&abs(pathlength_x)>0.1&&abs(pathlength_y)>0.1"; 
  //TCut cut = "invariant_mass>0.1&&abs(dca3dxy1)>0.02&&abs(dca3dxy2)>0.02&&abs(pair_dca)<0.05&&invariant_pt>0.1&&abs(projected_pathlength_x)>0.1&&abs(projected_pathlength_y)>0.1&&abs(projected_pathlength_z)>0.1"; 
  //TCut cut = "invariant_mass>0.1&&abs(dca3dxy1)>0.025&&abs(dca3dxy2)>0.025&&abs(pair_dca)<0.05&&invariant_pt>0.1&&abs(projected_pathlength_x)>0.1&&abs(projected_pathlength_y)>0.1"; 
  //TCut cut = "invariant_mass>0.1&&abs(dca3dxy1)>0.02&&abs(dca3dxy2)>0.02&&abs(dca3dz1)>0.02&&abs(dca3dz2)>0.02&&abs(pair_dca)<0.05&&invariant_pt>0.1&&abs(pathlength)>0.1"; 
  //TCut cut = "invariant_mass>0.1&&dca3dxy1>0.02&&dca3dxy2>0.02&&dca3dz1>0.02&&dca3dz2>0.02&&abs(pair_dca)<0.05&&invariant_pt>0.1&& sqrt(pow(pathlength_x,2) + pow(pathlength_y,2))>0.4";
  //TCut cut = "invariant_mass>0.1&&dca3dxy1>0.02&&dca3dxy2>0.02&&dca3dz1>0.02&&dca3dz2>0.02&&abs(pair_dca)<0.05&&invariant_pt>0.1"; 
  //TCut cut2 = "tpcClusters_1>40&&tpcClusters_2>40";
  //TCut kshort_cut = "invariant_mass > 0.485 && invariant_mass<0.515";


  int entries = ntuple->GetEntries();

  TH1D *invariant_mass_hist = new TH1D("invariant_mass_hist","Invariant Mass",1000,0.35,0.65);
  invariant_mass_hist->SetMinimum(0);
  TH2D *pathlengthpathlengthz = new TH2D("pathlengthpathlengthz","path v. pathlength_z",5000,-20.0,20.0,5000,0.0,10.0);

  TH2D *pathMass = new TH2D("pathMass","Invariant Mass vs. Path Length",5000,0.0,1.0,5000,0.0,10.0);


  int qual_cut = 10;
  float pathlength_cut = 0.1;
  double track_dca_cut = 0.0101;

  for(int i=0; i<entries; ++i)
    {
      ntuple->GetEntry(i);

      Eigen::Vector3f mom1(projected_mom1_x,projected_mom1_y,projected_mom1_z);
      Eigen::Vector3f mom2(projected_mom2_x,projected_mom2_y,projected_mom2_z);
      Eigen::Vector3f projected_momentum = mom1 + mom2; 
      Eigen::Vector3f pathLength (projected_pathlength_x,projected_pathlength_y,projected_pathlength_z);
      float pathLengthMag = pathLength.norm();

      Eigen::Vector3f vertex(vertex_x, vertex_y, vertex_z);
      float costheta = pathLength.dot(projected_momentum)/(projected_momentum.norm()*pathLength.norm());
      //float radius = sqrt(pow(projected_pathlength_x,2) + pow(projected_pathlength_y,2))
      
      if(costheta < 0.9995) continue;
      //if(abs(projected_pathlength_x) < 0.1 || abs(projected_pathlength_y) < 0.1) continue;
      if(abs(projected_pathlength_x) < pathlength_cut || abs(projected_pathlength_y) < pathlength_cut) continue;
      if(abs(projected_pair_dca) > 0.035) continue;
      if(invariant_pt < 0.1) continue;
      if(quality1 > qual_cut || quality2 > qual_cut) continue;
      if(dca3dxy1 < track_dca_cut || dca3dxy2 < track_dca_cut || dca3dz1 < track_dca_cut || dca3dz2 < track_dca_cut) continue;

      invariant_mass_hist->Fill(invariant_mass);
      pathlengthpathlengthz->Fill(projected_pathlength_z,pathLengthMag);
      pathMass->Fill(invariant_mass,pathLengthMag);
    } 

  TCanvas *c1 = new TCanvas("c1","",10,10,800,800);
 
  TF1* f = new TF1("f","gaus(0)+[3]+[4]*x",0.45,0.55);
  f->SetParameter(0,100);
  f->SetParameter(1,0.497);
  f->SetParameter(2,0.010);
  f->SetParameter(3,0.0);
  f->SetParameter(4,0.0);
  
  invariant_mass_hist->Fit("f","R");
  //invariant_mass_hist->DrawCopy();

  int binLow      = invariant_mass_hist->FindBin(0.48); 
  int binHigh     = invariant_mass_hist->FindBin(0.52);
  double integral = invariant_mass_hist->Integral(binLow,binHigh);
  std::cout << " Integral: "<< integral << std::endl;

  TF1* fbackground  = new TF1("fbackground","[0]+[1]*x");
  fbackground->SetParameter(0,f->GetParameter(3));
  fbackground->SetParameter(1,f->GetParameter(4));

  double backgroundIntegral    = fbackground->Integral(0.485,0.51);
  double foregroundIntegral    = f->Integral(0.485,0.51);
  double signal                = (foregroundIntegral - backgroundIntegral);
  double signalBackgroundRatio = signal/backgroundIntegral;
  std::cout<< "signal: " << signal<<"  Signal to background ratio: " << signalBackgroundRatio << std::endl;


  TF1* fsignal = new TF1("fsignal","gaus(0)");
  fsignal->SetParameter(0,f->GetParameter(0));
  fsignal->SetParameter(1,f->GetParameter(1));
  fsignal->SetParameter(2,f->GetParameter(2));
  double signal_integral = fsignal->Integral(0.35,0.65);
  double binwidth = invariant_mass_hist->GetBinWidth(1);
  std::cout << " Integral: "<< signal_integral<< " yield: "<< signal_integral/binwidth<<std::endl;

  
  TCanvas *c2 = new TCanvas("c2","",10,10,600,600); // path v path_z
  //pathlengthpathlengthz->DrawCopy();

  TCanvas *c3 = new TCanvas("c3","",10,10,600,600);
  //pathMass->DrawCopy();

  
  

  //increase dca cut may reduce background more than signal

  // TH1D *invariant_mass = new TH1D("invariant_mass","Invariant Mass",200,0.47,0.53);
  // ntuple->Draw("invariant_mass>>invariant_mass",cut);
  // invariant_mass->DrawCopy();

  // TCanvas *c2 = new TCanvas("c2","",10,10,600,600);
  // TH2D *invariants = new TH2D("invariants","Invariant mass v. Invariant pt",5000,0,7,5000,0,5);
  // ntuple->Draw("invariant_mass:invariant_pt>>invariants",cut);
  // invariants->DrawCopy();

  // TH2D *pseudoInv= new TH2D("pseudoInv","Invariant Mass and pseudorapidity",5000,0,3,5000,0,7);
  // ntuple->Draw("invariant_mass:pseudorapidity>>pseudoInv",cut);
  // pseudoInv->DrawCopy()    ;

  // TH2D *pairRapidity = new TH2D("pairRapidity","Pair DCA and pseudorapidity",5000,-2,2,5000,-0.06,0.06);
  // ntuple->Draw("pair_dca:pseudorapidity>>pairRapidity",cut);
  // pairRapidity->DrawCopy();


  /*
  TCanvas *c3 = new TCanvas("c3","",10,10,600,600);
  TH2D *pathMass = new TH2D("pathMass","Invariant Mass vs. Path Length",5000,0.0,5.0,5000,0.4,0.6);
  ntuple->Draw("invariant_mass:pathlength>>pathMass",cut+cut2);
  pathMass->DrawCopy();

  TCanvas *c4 = new TCanvas("c4","",10,10,600,600);
  TH2D *radiusMass = new TH2D("radiusMass","Invariant Mass vs. Radius",5000,0.0,5.0,5000,0.4,0.6);
  ntuple->Draw("invariant_mass:sqrt(pow(pathlength_x,2) + pow(pathlength_y,2))>>radiusMass",cut+cut2);
  radiusMass->DrawCopy();

  TCanvas *c5 = new TCanvas("c5","",10,10,600,600);
  TH2D *pathlengthpathlengthz = new TH2D("pathlengthpathlengthz","path v. pathlength_z",5000,-20.0,20.0,5000,0.0,10.0);
  ntuple->Draw("pathlength:pathlength_z>>pathlengthpathlengthz",cut+kshort_cut+cut2);
  pathlengthpathlengthz->DrawCopy();
  */


}
