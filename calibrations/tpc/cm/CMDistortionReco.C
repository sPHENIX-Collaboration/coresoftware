// step 2 with phi,r coords
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include <TTime.h>

using namespace std;

int CMDistortionReco(int nMaxEvents = -1) {
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
  //  nEvents = 2;
  
  TCanvas *canvas1=new TCanvas("canvas1","CMDistortionReco1",1200,800);
  canvas1->Divide(3,2);

  //canvas for time plot
  TCanvas *canvas=new TCanvas("canvas","CMDistortionReco2",400,400);
  
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

    //hardcoded numbers from average distortion file's hIntDistortionPosX
    int nbinsphi = 30; //when using 35, blank spots at around r = 22 cm, phi just above n below pi
    double lowphi = -0.078539819;
    double highphi = 6.3617253;
    int nbinsr = 30; // when using 35, blank stripe around r = 58 cm
    double lowr = 0.0;
    double highr = 90.0;
    
    //for forward only
   
    TH2F *hStripesPerBin = new TH2F("hStripesPerBin","CM Stripes Per Bin (z in stripes); x (cm); y (cm)",nbins,low,high,nbins,low,high);

    //phi,r binning
    TH2F *hStripesPerBinPhiR = new TH2F("hStripesPerBinPhiR","CM Stripes Per Bin (z in stripes); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);

    TH2F *hCartesianForward[3];
    hCartesianForward[0] = new TH2F("hForwardX","X Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[1] = new TH2F("hForwardY","Y Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[2] = new TH2F("hForwardZ","Z Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);

    //phi,r binning
    TH2F *hCartesianForwardPhiR[3];
    hCartesianForwardPhiR[0] = new TH2F("hForwardX_PhiR","X Shift Forward of Stripe Centers, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCartesianForwardPhiR[1] = new TH2F("hForwardY_PhiR","Y Shift Forward of Stripe Centers, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCartesianForwardPhiR[2] = new TH2F("hForwardZ_PhiR","Z Shift Forward of Stripe Centers, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    
     TH2F *hCylindricalForwardPhiR[2];
    hCylindricalForwardPhiR[0] = new TH2F("hForwardR_PhiR","Radial Shift Forward of Stripe Centers, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCylindricalForwardPhiR[1] = new TH2F("hForwardPhi_PhiR","Phi Shift Forward of Stripe Centers, Phi,R binning (rad); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    
    for (int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);

      hStripesPerBin->Fill(position->X(),position->Y(),1);

      //for phi,r binning
      double r = position->Perp();
      double phi = position->Phi();

      if(position->Phi() < 0.0){
	phi = position->Phi() + 2.0*TMath::Pi(); 
      }
      
      hStripesPerBinPhiR->Fill(phi,r,1);
      
      deltaX = (newposition->X() - position->X())*(1e4); //convert from cm to micron 
      deltaY = (newposition->Y() - position->Y())*(1e4);
      deltaZ = (newposition->Z() - position->Z())*(1e4);

      deltaR = (newposition->Perp() - position->Perp())*(1e4);
      deltaPhi = newposition->DeltaPhi(*position);

      hCartesianForward[0]->Fill(position->X(),position->Y(),deltaX);
      hCartesianForward[1]->Fill(position->X(),position->Y(),deltaY);
      hCartesianForward[2]->Fill(position->X(),position->Y(),deltaZ);

      // phi,r binning
      hCartesianForwardPhiR[0]->Fill(phi,r,deltaX);
      hCartesianForwardPhiR[1]->Fill(phi,r,deltaY);
      hCartesianForwardPhiR[2]->Fill(phi,r,deltaZ);

      hCylindricalForwardPhiR[0]->Fill(phi,r,deltaR);
      hCylindricalForwardPhiR[1]->Fill(phi,r,deltaPhi);
    }

    TH2F *hCartesianAveShift[3];
    hCartesianAveShift[0] = new TH2F("AveShiftX","Average of CM Model X over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[1] = new TH2F("AveShiftY","Average of CM Model Y over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[2] = new TH2F("AveShiftZ","Average of CM Model Z over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    TH2F *hCylindricalAveShift[2];
     hCylindricalAveShift[0] = new TH2F("AveShiftRCart","Average of CM Model R over Stripes per Bin from Cartesian (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[1] = new TH2F("AveShiftPhiCart","Average of CM Model Phi over Stripes per Bin from Cartesian (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    // phi,r binning
     TH2F *hCartesianAveShiftPhiR[3];
    hCartesianAveShiftPhiR[0] = new TH2F("AveShiftX_PhiR","Average of CM Model X over Stripes per Bin, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCartesianAveShiftPhiR[1] = new TH2F("AveShiftY_PhiR","Average of CM Model Y over Stripes per Bin, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCartesianAveShiftPhiR[2] = new TH2F("AveShiftZ_PhiR","Average of CM Model Z over Stripes per Bin, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 

    TH2F *hCylindricalAveShiftPhiR[2];
     hCylindricalAveShiftPhiR[0] = new TH2F("AveShiftR_PhiR","Average of CM Model R over Stripes per Bin, Phi,R binning (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCylindricalAveShiftPhiR[1] = new TH2F("AveShiftPhi_PhiR","Average of CM Model Phi over Stripes per Bin, Phi,R binning (rad); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    
    for (int i = 0; i < 3; i ++){
      hCartesianAveShift[i]->Divide(hCartesianForward[i],hStripesPerBin);
      hCartesianAveShiftPhiR[i]->Divide(hCartesianForwardPhiR[i],hStripesPerBinPhiR);
    }

    hCylindricalAveShiftPhiR[0]->Divide(hCylindricalForwardPhiR[0],hStripesPerBinPhiR);
    hCylindricalAveShiftPhiR[1]->Divide(hCylindricalForwardPhiR[1],hStripesPerBinPhiR);
    
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

    //phi,r binning
    TH3F *hCartesianCMModelPhiR[3];
    hCartesianCMModelPhiR[0]=new TH3F("hCMModelX_PhiR", "CM Model: X Shift Forward of Stripe Centers, Phi,R binning", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz); //rad, cm, cm
    hCartesianCMModelPhiR[1]=new TH3F("hCMModelY_PhiR", "CM Model: Y Shift Forward of Stripe Centers, Phi,R binning", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCartesianCMModelPhiR[2]=new TH3F("hCMModelZ_PhiR", "CM Model: Z Shift Forward of Stripe Centers, Phi,R binning", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

    TH3F *hCylindricalCMModelPhiR[2];
    hCylindricalCMModelPhiR[0]=new TH3F("hCMModelR_PhiR", "CM Model: Radial Shift Forward of Stripe Centers, Phi,R binning", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCylindricalCMModelPhiR[1]=new TH3F("hCMModelPhi_PhiR", "CM Model: Phi Shift Forward of Stripe Centers, Phi,R binning", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);  
      
    double xshift, yshift, zshift, rshiftcart, phishiftcart;
    double xshiftPhiR, yshiftPhiR, zshiftPhiR, rshiftcartPhiR, phishiftcartPhiR;
  
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

	  //phi,r binning
	  xshiftPhiR=(hCartesianAveShiftPhiR[0]->Interpolate(phi,r))*(1e-4);//coordinate of your stripe
	  yshiftPhiR=(hCartesianAveShiftPhiR[1]->Interpolate(phi,r))*(1e-4);//convert micron to cm
	  zshiftPhiR=(hCartesianAveShiftPhiR[2]->Interpolate(phi,r))*(1e-4);

	  rshiftcartPhiR=(hCylindricalAveShiftPhiR[0]->Interpolate(phi,r))*(1e-4);
	  phishiftcartPhiR=hCylindricalAveShiftPhiR[1]->Interpolate(phi,r);
	  
	  hCartesianCMModelPhiR[0]->Fill(phi,r,z,xshiftPhiR*(1-z/105.5));
	  hCartesianCMModelPhiR[1]->Fill(phi,r,z,yshiftPhiR*(1-z/105.5));
	  hCartesianCMModelPhiR[2]->Fill(phi,r,z,zshiftPhiR*(1-z/105.5));

	  hCylindricalCMModelPhiR[0]->Fill(phi,r,z,rshiftcartPhiR*(1-z/105.5));
	  hCylindricalCMModelPhiR[1]->Fill(phi,r,z,phishiftcartPhiR*(1-z/105.5));
	}
      }
    }
    
    TFile *plots;

    plots=TFile::Open(Form("CMModels_Event%d.root",ifile),"RECREATE");
    hStripesPerBin->Write(); 

    for(int i = 0; i < 3; i++){
      hCartesianForward[i]->Write();
      hCartesianAveShift[i]->Write();
      hCartesianCMModel[i]->Write();

      hCartesianForwardPhiR[i]->Write();
      hCartesianAveShiftPhiR[i]->Write();
      hCartesianCMModelPhiR[i]->Write();
    }

    for(int i = 0; i < 2; i++){
      hCylindricalAveShift[i]->Write();
      hCylindricalCMModel[i]->Write();

      hCylindricalForwardPhiR[i]->Write();
      hCylindricalAveShiftPhiR[i]->Write();
      hCylindricalCMModelPhiR[i]->Write();
    }
    
    plots->Close();

    
    //call to TTime after outputting TH3Fs
    now=gSystem->Now();
    unsigned long after = now;
    
    hTimePerEvent->Fill(after-before);
    
    //to check histograms
    for (int i = 0; i < 3; i++){
      hCartesianForward[i]->SetStats(0);
      hCartesianForwardPhiR[i]->SetStats(0);
    }

    hStripesPerBin->SetStats(0);
    canvas1->cd(1);
    hStripesPerBinPhiR->Draw("colz");
    canvas1->cd(2);
    hCartesianForward[1]->Draw("colz");
    canvas1->cd(3);
    hCartesianForward[2]->Draw("colz");
    canvas1->cd(4);
    hCartesianForwardPhiR[0]->Draw("colz");
    canvas1->cd(5);
    hCartesianForwardPhiR[1]->Draw("colz");
    canvas1->cd(6);
    hCartesianForwardPhiR[2]->Draw("colz");
  
    if(ifile == 0){ 
      canvas1->Print("CMDistortionReco1.pdf(","pdf");
    } else if (ifile == nEvents - 1){
      canvas1->Print("CMDistortionReco1.pdf)","pdf");
    } else {
      canvas1->Print("CMDistortionReco1.pdf","pdf");
    }
    
  }

  canvas->cd();
  hTimePerEvent->Draw();
  canvas->Print("CMDistortionReco2.pdf","pdf");

  return 0;
}
