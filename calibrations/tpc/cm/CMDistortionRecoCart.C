// step 2
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include <TTime.h>

using namespace std;

int CMDistortionRecoCart(int nMaxEvents = -1) {
  int nbins = 35; 
  double low = -80.0;
  double high = 80.0;
  double deltaX, deltaY, deltaZ, deltaR, deltaPhi;
  int nEvents;
    
  //take in events
  const char * inputpattern="/sphenix/u/skurdi/CMCalibration/cmDistHitsTree_Event*.root"; 
  
  //find all files that match the input string (includes wildcards)
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  TString sourcefilename;
  
  //how many events
  if (nMaxEvents<0){
    nEvents=filelist->GetNFiles();
  } else if(nMaxEvents<filelist->GetNFiles()){
    nEvents=nMaxEvents;
  } else {
    nEvents= filelist->GetNFiles();
  }
  
  //canvas for checking data
  //TCanvas *canvas1=new TCanvas("canvas1","CMDistortionRecoCart1",1200,800);
  //canvas1->Divide(3,2);

  //canvas for time plot
  TCanvas *canvas=new TCanvas("canvas","CMDistortionRecoCart2",400,400);
  
  TVector3 *position, *newposition;
  position = new TVector3(1.,1.,1.);
  newposition = new TVector3(1.,1.,1.);

  //histogram to compare times
  TH1F *hTimePerEvent = new TH1F("hTimePerEvent","Time Per Event; time (ms)",20,0,10000);
    
  for (int ifile=0;ifile < nEvents;ifile++){
    //call to TTime before opening ttree
    TTime now;
    now=gSystem->Now();
    unsigned long before = now;
    
    //get data from ttree
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();
    
    char const *treename="cmDistHitsTree";
    TFile *input=TFile::Open(sourcefilename, "READ");
    TTree *inTree=(TTree*)input->Get("tree");
    
    inTree->SetBranchAddress("position",&position);
    inTree->SetBranchAddress("newposition",&newposition);

    //for forward only
   
    TH2F *hStripesPerBin = new TH2F("hStripesPerBin","CM Stripes Per Bin (z in stripes); x (cm); y (cm)",nbins,low,high,nbins,low,high);

    TH2F *hCartesianForward[3];
    hCartesianForward[0] = new TH2F("hForwardX","X Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[1] = new TH2F("hForwardY","Y Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[2] = new TH2F("hForwardZ","Z Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
 
    for (int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);

      hStripesPerBin->Fill(position->X(),position->Y(),1);
      
      deltaX = (newposition->X() - position->X())*(1e4); //convert from cm to micron 
      deltaY = (newposition->Y() - position->Y())*(1e4);
      deltaZ = (newposition->Z() - position->Z())*(1e4);

      deltaR = (newposition->Perp() - position->Perp())*(1e4);
      deltaPhi = newposition->DeltaPhi(*position);

      hCartesianForward[0]->Fill(position->X(),position->Y(),deltaX);
      hCartesianForward[1]->Fill(position->X(),position->Y(),deltaY);
      hCartesianForward[2]->Fill(position->X(),position->Y(),deltaZ);
    }

    TH2F *hCartesianAveShift[3];
    hCartesianAveShift[0] = new TH2F("AveShiftX","Average of CM Model X over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[1] = new TH2F("AveShiftY","Average of CM Model Y over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[2] = new TH2F("AveShiftZ","Average of CM Model Z over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    TH2F *hCylindricalAveShift[2];
     hCylindricalAveShift[0] = new TH2F("AveShiftRCart","Average of CM Model R over Stripes per Bin from Cartesian (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[1] = new TH2F("AveShiftPhiCart","Average of CM Model Phi over Stripes per Bin from Cartesian (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    
    for (int i = 0; i < 3; i ++){
      hCartesianAveShift[i]->Divide(hCartesianForward[i],hStripesPerBin);
    }

    //cyl models from cart coords
    for(int i = 0; i < nbins; i++){
      double x = low + ((high - low)/(1.0*nbins))*(i+0.5); //center of bin
      
      for(int j = 0; j < nbins; j++){
	double y = low + ((high - low)/(1.0*nbins))*(j+0.5); //center of bin
		
	int xbin = hCartesianAveShift[0]->FindBin(x,y);
	int ybin = hCartesianAveShift[1]->FindBin(x,y);
	double xaveshift = (hCartesianAveShift[0]->GetBinContent(xbin))*(1e-4); // converts  microns to cm 
	double yaveshift = (hCartesianAveShift[1]->GetBinContent(ybin))*(1e-4);
	
	TVector3 shifted, original;
	original.SetX(x);
	original.SetY(y);
	shifted.SetX(x+xaveshift);
	shifted.SetY(y+yaveshift);
	//x n y above for orig
	//shifted is orig + ave shift
	
	double raveshift = (shifted.Perp() - original.Perp())*(1e4);
	double phiaveshift = shifted.DeltaPhi(original);

	//fill with r from x n y
	hCylindricalAveShift[0]->Fill(x,y,raveshift);
	hCylindricalAveShift[1]->Fill(x,y,phiaveshift);
      } 
    } 
  
    //same range and bins for each coordinate, binned in cm
    //hardcoded numbers from average distortion file's hIntDistortionPosX
    int nphi = 82;
    int nr = 54;
    int nz = 82;
    
    double minphi = -0.078539819;
    double minr = 18.884615;
    //double minz = 5.0;
    double minz = -1.3187500;
    
    double maxphi = 6.3617253;
    double maxr = 79.115387;
    double maxz = 106.81875;

    TH3F *hCartesianCMModel[3];
    hCartesianCMModel[0]=new TH3F("hCMModelX", "CM Model: X Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz); //rad, cm, cm
    hCartesianCMModel[1]=new TH3F("hCMModelY", "CM Model: Y Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCartesianCMModel[2]=new TH3F("hCMModelZ", "CM Model: Z Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

    TH3F *hCylindricalCMModel[2];
    hCylindricalCMModel[0]=new TH3F("hCMModelRCart", "CM Model: Radial Shift Forward of Stripe Centers from Cartesian", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCylindricalCMModel[1]=new TH3F("hCMModelPhiCart", "CM Model: Phi Shift Forward of Stripe Centers from Cartesian", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
      
    double xshift, yshift, zshift, rshiftcart, phishiftcart;
      
    for(int i = 0; i < nphi; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin

      for(int j = 0; j < nr; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin

	double x = r*cos(phi); //cm
	double y = r*sin(phi);

	for(int k = 0; k < nz; k++){
	  double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	  xshift=(hCartesianAveShift[0]->Interpolate(x,y))*(1e-4);//coordinate of your stripe
	  yshift=(hCartesianAveShift[1]->Interpolate(x,y))*(1e-4);//convert micron to cm
	  zshift=(hCartesianAveShift[2]->Interpolate(x,y))*(1e-4);

	  rshiftcart=(hCylindricalAveShift[0]->Interpolate(x,y))*(1e-4);
	  phishiftcart=hCylindricalAveShift[1]->Interpolate(x,y);
	  
	  hCartesianCMModel[0]->Fill(phi,r,z,xshift*(1-z/105.5));
	  hCartesianCMModel[1]->Fill(phi,r,z,yshift*(1-z/105.5));
	  hCartesianCMModel[2]->Fill(phi,r,z,zshift*(1-z/105.5));

	  hCylindricalCMModel[0]->Fill(phi,r,z,rshiftcart*(1-z/105.5));
	  hCylindricalCMModel[1]->Fill(phi,r,z,phishiftcart*(1-z/105.5));
	}
      }
    }
    
    TFile *plots;

    plots=TFile::Open(Form("CMModelsCart_Event%d.root",ifile),"RECREATE");
    hStripesPerBin->Write(); 

    for(int i = 0; i < 3; i++){
      hCartesianForward[i]->Write();
      hCartesianAveShift[i]->Write();
      hCartesianCMModel[i]->Write();
    }

    for(int i = 0; i < 2; i++){
      hCylindricalAveShift[i]->Write();
      hCylindricalCMModel[i]->Write();
    }
    
    plots->Close();

    
    //call to TTime after outputting TH3Fs
    now=gSystem->Now();
    unsigned long after = now;
    
    hTimePerEvent->Fill(after-before);
    
    /*
    //to check histograms
    for (int i = 0; i < 3; i++){
      hCartesianForward[i]->SetStats(0);
    }

    hStripesPerBin->SetStats(0);
    canvas1->cd(1);
    hCartesianForward[0]->Draw("colz");
    canvas1->cd(2);
    hCartesianForward[1]->Draw("colz");
    canvas1->cd(3);
    hCartesianForward[2]->Draw("colz");
    canvas1->cd(4)->Clear();
    canvas1->cd(5)->Clear();
    canvas1->cd(6)->Clear();
      
    if(ifile == 0){ 
      canvas1->Print("CMDistortionRecoCart1.pdf(","pdf");
    } else if (ifile == nEvents - 1){
      canvas1->Print("CMDistortionRecoCart1.pdf)","pdf");
    } else {
      canvas1->Print("CMDistortionRecoCart1.pdf","pdf");
    }*/
    
    
  }

  canvas->cd();
  hTimePerEvent->Draw();
  canvas->Print("CMDistortionRecoCart2.pdf","pdf");

  return 0;
}
