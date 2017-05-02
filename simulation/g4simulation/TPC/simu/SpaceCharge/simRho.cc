#include <iostream>
#include <fstream>
#include "sHelix.h"
#include "sChargeMap.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

const int npid = 12+1;
int apid[npid] =     {11, 13, // 2
                      211, // 1
                      321, // 1
                      411, 431, // 2
                      511, // 1
                      2212, // 1
                      3222, 3112, 3312, 3334, // 6
                      0};
int acha[npid] = {-1,-1,+1,+1,+1,+1,+1,+1,+1,-1,-1,-1};
TString spid[npid] = {"e^{-}",
                      "#mu^{-}",
                      "#pi^{+}",
                      "K^{+}",
                      "D^{+}",
                      "D^{+}_{s}",
                      "B^{+}",
                      "p",
                      "#Sigma^{+}",
                      "#Sigma^{-}",
                      "#Xi^{-}",
                      "#Omega^{-}",

                      "Neutral"};
TH1F *hpid      = new TH1F("hpid","hpid",npid,-0.5,npid-0.5);
TH1F *hrpass    = new TH1F("hrpass","hrpass; r [cm]",100,-15,+15);
TH1F *hzpass    = new TH1F("hzpass","hzpass; z [cm]",100,-15,+15);
TH1F *hpidpass  = new TH1F("hpidpass","hpidpass",npid,-0.5,npid-0.5);
TH1F *hptpass   = new TH1F("hptpass","hptpass; pt [GeV]",100,0,10);
TH1F *hetapass  = new TH1F("hetapass","hetapass; eta",100,-3,+3);
TH1F *hlengthpass=new TH1F("hlengthpass","hlengthpass; length [cm]",100,0,150);
TH1F *hbxing    = new TH1F("hbxing","hbxing",150,-0.5,300-0.5);
TH3F *tracking  = new TH3F("tracking","Generated tracks: tracking;xi;yi;zi",120,-100,+100,120,-100,+100,120,-100,+100);

int main() {
  for(int n=0; n!=npid; ++n) hpid->GetXaxis()->SetBinLabel(1+n,spid[n].Data());
  for(int n=0; n!=npid; ++n) hpidpass->GetXaxis()->SetBinLabel(1+n,spid[n].Data());

  gStyle->SetOptStat(111111111);
  //================
  // Gas parameters
  //float gs_wf_me = 35; // [eV] mean energy for ionisation
  //float gs_density = 1e-3; // [gr/cm^3]
  //float gs_np = 14.35; // [cm^-1]
  float gs_nt = 47.80; // [cm^-1]
  float gs_ion_mobility = 4; // [cm^2/(V.s)]

  //================
  // TPC running conditions
  float tpc_magnetic_field = 1.5; // [Tesla]
  float tpc_electric_field = 400; // [V/cm]
  float tpc_drift_velocity_e = 8; // [cm/us]
  float tpc_drift_velocity_i = tpc_electric_field*gs_ion_mobility*1e-3; // [cm/ms]
  int lumi = 189; //  1/(189*106e-9) => 49915.1 Hz

  // input file
  std::ifstream inputf("/Users/cperez//ampt/ana/ampt.dat");

  sChargeMap *map = new sChargeMap(100,30,80,
				   1,0,TMath::TwoPi(),
				   4*int(80/(tpc_drift_velocity_i*(lumi*106e-6)))+1,80,
				   tpc_drift_velocity_e,tpc_drift_velocity_i);
  // main routine
  int ncoll;
  for(ncoll=0; ncoll!=3; ++ncoll) {
    // =========
    // COLLISION
    if(ncoll%100==0) std::cout << ncoll << std::endl;
    // event wise
    int eventno, testno, nhadrons, npartA, npartB, npartElaA, npartIneA, npartElaB, npartIneB;
    float impactPar;
    // hadron wise
    int pid;
    float px, py, pz, mass, x, y, z, t;
    inputf >> eventno;
    if(!inputf.good()) {
      inputf.close();
      inputf.open("ampt.dat");
      inputf >> eventno;
    }
    inputf >> testno >> nhadrons >> impactPar >> npartA >> npartB >> npartElaA >> npartIneA >> npartElaB >> npartIneB;
    for(int np=0; np!=nhadrons; ++np) {
      inputf >> pid >> px >> py >> pz >> mass >> x >> y >> z >> t;
      if( TMath::AreEqualAbs(t,0,1e-3) ) continue;
      x *= 1e-13;
      y *= 1e-13;
      z *= 1e-13;
      int fill = npid-1;
      for(int n=0; n!=npid-1; ++n) if(apid[n]==TMath::Abs(pid)) fill = n;
      hpid->Fill( fill );
      if(fill==npid-1) continue;
      int q;
      for(int n=0; n!=npid-1; ++n) if(apid[n]==TMath::Abs(pid)) q = acha[n];
      sHelix a(x,y,z,px,py,pz,q,tpc_magnetic_field);
      float t1 = a.findFirstInterceptTo(30,80);
      if(t1>999) continue;
      float t2 = a.findFirstInterceptTo(80,80);
      if(t2>999) continue;
      float length = a.s(t1,t2);
      if(length<1) continue;
      hrpass->Fill( sqrt(x*x+y*y) );
      hzpass->Fill( z );
      hpidpass->Fill( fill );
      hptpass->Fill( sqrt(px*px+py*py) );
      float pmom = sqrt(px*px+py*py+pz*pz);
      float eta = 1.e30;
      if (pmom != TMath::Abs(pz)) eta = 0.5*TMath::Log((pmom+pz)/(pmom-pz));
      hetapass->Fill( eta );
      hlengthpass->Fill( length );
      float nprim = gRandom->Poisson(gs_nt*length);
      float track[100][3];
      a.breakIntoPieces(t1,t2,track);
      for(int i=0; i!=100; ++i) {
	float ri = TMath::Sqrt( track[i][0]*track[i][0] + track[i][1]*track[i][1] );
	if(ri<30-0.1) break;
	if(ri>80+0.1) break;
	float pi = 1;
        map->Fill(ri, pi, track[i][2], -nprim/100);
        map->Fill(ri, pi, track[i][2], nprim/100);
	//if(ncoll==1&&
	//   (
	//    np==622||
	//    np==935
	//    )) {
	tracking->Fill(track[i][0],track[i][1],track[i][2]);
	//  if(i==0) {
	//  std::cout << "pmom " << pmom << "| eta " << eta;
	//  std::cout << "| x y z " << x << " " << y << " " << z << " " << px << " " << py << " " << pz << std::endl;
	//  std::cout << "t " << t1 << "|" << t2 << std::endl;
	//}
	//std::cout << ri << " ";
	//if(i==99)
	//  std::cout << std::endl;
	//}
      }
    }

    // ==========
    // DRIFT TIME
    int steps = gRandom->Poisson( lumi );
    hbxing->Fill(steps);
    float lapse = 106*steps*1e-6; // [ns] -> [ms]
    map->Propagate(lapse);
  }
  map->ScreenShot("map","root",ncoll);
  inputf.close();
  map->SaveIonMap("ionmap.root");
  //map->SaveRho("mark1_0.root",100,1,1600);
  
  TFile *qa = new TFile("qa.root","RECREATE");
  qa->WriteTObject(hpid,hpid->GetName());
  qa->WriteTObject(hrpass,hrpass->GetName());
  qa->WriteTObject(hzpass,hzpass->GetName());
  qa->WriteTObject(hpidpass,hpidpass->GetName());
  qa->WriteTObject(hptpass,hptpass->GetName());
  qa->WriteTObject(hetapass,hetapass->GetName());
  qa->WriteTObject(hlengthpass,hlengthpass->GetName());
  qa->Close();
  tracking->SaveAs("tracking.root","root");

  return 0;
}
