//step 3 with phi,r coords
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"

using namespace std;

class Shifter {
public:
Shifter(TString sourcefilename);
  TFile *forward, *average;
  TH3F *hX, *hY, *hZ, *hR, *hPhi, *hXave, *hYave, *hZave, *hRave, *hPhiave;  
};

Shifter::Shifter(TString sourcefilename){
  //single event distortion file
  forward=TFile::Open(sourcefilename,"READ"); 

  hX=(TH3F*)forward->Get("hIntDistortionPosX");
  hY=(TH3F*)forward->Get("hIntDistortionPosY");
  hZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hR=(TH3F*)forward->Get("hIntDistortionPosR");
  hPhi=(TH3F*)forward->Get("hIntDistortionPosP");

  //average distortion file
  average=TFile::Open("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 
  
  hXave=(TH3F*)average->Get("hIntDistortionPosX");
  hYave=(TH3F*)average->Get("hIntDistortionPosY");
  hZave=(TH3F*)average->Get("hIntDistortionPosZ");
  
  hRave=(TH3F*)average->Get("hIntDistortionPosR");
  hPhiave=(TH3F*)average->Get("hIntDistortionPosP");
 
  //subtract average from total distortions to study fluctuations
  hX->Add(hXave,-1);
  hY->Add(hYave,-1);
  hZ->Add(hZave,-1);
  
  hR->Add(hRave,-1);
  hPhi->Add(hPhiave,-1);
}

int CMDistortionAnalysisPhiR(int nMaxEvents = -1) {
  Shifter *shifter;
  int nbins = 35; 
  double x, y, z;
  double low = -80.0;
  double high = 80.0;
  double deltaX, deltaY, deltaZ, deltaR, deltaPhi;
  int nEvents; 
  
  TCanvas *canvas=new TCanvas("canvas","CMDistortionAnalysisPhiR",2000,3000);

  int nsumbins = 20;
  int minsum = -10;
  int maxsum = 10;
  
  //set up summary plots
  TH1F *hDifferenceMeanR = new TH1F("hDifferenceMeanR", "Average Difference between R Model and True of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    TH1F *hDifferenceStdDevR = new TH1F("hDifferenceStdDevR", "Std Dev of Difference between R Model and True of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hTrueMeanR = new TH1F("hTrueMeanR", "Mean True R Distortion Model of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    TH1F *hTrueStdDevR = new TH1F("hTrueStdDevR", "Std Dev of True R Distortion Model of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hDifferenceMeanPhi = new TH1F("hDifferenceMeanPhi", "Average Difference between Phi Model and True of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    TH1F *hDifferenceStdDevPhi = new TH1F("hDifferenceStdDevPhi", "Std Dev of Difference between Phi Model and True of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hTrueMeanPhi = new TH1F("hTrueMeanPhi", "Mean True Phi Distortion Model of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    TH1F *hTrueStdDevPhi = new TH1F("hTrueStdDevPhi", "Std Dev of True Phi Distortion Model of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);

    const char * inputpattern="/sphenix/user/rcorliss/distortion_maps/2021.04/*h_Charge_*.root"; //updated
    
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

  for (int ifile=0;ifile < nEvents;ifile++){
    //for each file, find all histograms in that file.
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();

    //create shifter
    shifter = new Shifter(sourcefilename);
    
    TFile *plots;

    plots=TFile::Open(Form("CMModelsPhiR_Event%d.root",ifile),"READ");

    TH3F *hCartCMModelPhiR[3];
    hCartCMModelPhiR[0]=(TH3F*)plots->Get("hCMModelX_PhiR");
    hCartCMModelPhiR[1]=(TH3F*)plots->Get("hCMModelY_PhiR");
    hCartCMModelPhiR[2]=(TH3F*)plots->Get("hCMModelZ_PhiR");

    TH3F *hCylCMModelPhiR[2];
    hCylCMModelPhiR[0]=(TH3F*)plots->Get("hCMModelR_PhiR");
    hCylCMModelPhiR[1]=(TH3F*)plots->Get("hCMModelPhi_PhiR");
    
    //for forward only

    //same range and bins for each coordinate, binned in cm
    //hardcoded numbers from average distortion file's hIntDistortionPosX
    int nphi = 82;
    int nr = 54;
    int nz = 82;
    
    double minphi = -0.078539819;
    double minr = 18.884615;
    double minz = -1.3187500;
    
    double maxphi = 6.3617253;
    double maxr = 79.115387;
    double maxz = 106.81875;

    double rshiftcart, phishiftcart;

    int ndiff = 300;
    int mindiff = -20;
    int maxdiff = 20;
    
     TH1F *hCartesianShiftDifferencePhiR[3];
    hCartesianShiftDifferencePhiR[0] = new TH1F("hShiftDifferenceX_PhiR", "Difference between CM Model X and True, Phi,R binning (R > 30); #Delta X (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifferencePhiR[1] = new TH1F("hShiftDifferenceY_PhiR", "Difference between CM Model Y and True, Phi,R binning (R > 30); #Delta Y (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifferencePhiR[2] = new TH1F("hShiftDifferenceZ_PhiR", "Difference between CM Model Z and True, Phi,R binning (R > 30); #Delta Z (#mum)", ndiff, mindiff, maxdiff);
    
    TH1F *hCylindricalShiftDifferencePhiR[2];
    hCylindricalShiftDifferencePhiR[0] = new TH1F("hShiftDifferenceR_PhiR", "Difference between CM Model R and True, Phi,R binning (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    hCylindricalShiftDifferencePhiR[1] = new TH1F("hShiftDifferencePhi_PhiR", "Difference between CM Model Phi and True, Phi,R binning (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);

    TH1F *hRShiftTrue = new TH1F("hRShiftTrue", "True R Distortion Model (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    TH1F *hPhiShiftTrue = new TH1F("hPhiShiftTrue", "True Phi Distortion Model (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);
  
     TH2F *hCartesianDiffPhiR[6];
    hCartesianDiffPhiR[0] = new TH2F("hDiffXYX_PhiR", "Difference in PhiR for CM Model X, Phi,R binning; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiR[1] = new TH2F("hDiffRZX_PhiR", "Difference in RZ for CM Model X, Phi,R binning; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianDiffPhiR[2] = new TH2F("hDiffXYY_PhiR", "Difference in PhiR for CM Model Y, Phi,R binning; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiR[3] = new TH2F("hDiffRZY_PhiR", "Difference in RZ for CM Model Y, Phi,R binning; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianDiffPhiR[4] = new TH2F("hDiffXYZ_PhiR", "Difference in PhiR for CM Model Z, Phi,R binning; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiR[5] = new TH2F("hDiffRZZ_PhiR", "Difference in RZ for CM Model Z, Phi,R binning; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    
    TH2F *hCylindricalDiffPhiR[4];
    hCylindricalDiffPhiR[0] = new TH2F("hDiffXYR_PhiR", "Difference in PhiR for CM Model R, Phi,R binning; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalDiffPhiR[1] = new TH2F("hDiffRZR_PhiR", "Difference in RZ for CM Model R, Phi,R binning; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
    hCylindricalDiffPhiR[2] = new TH2F("hDiffXYPhi_PhiR", "Difference in PhiR for CM Model Phi, Phi,R binning; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalDiffPhiR[3] = new TH2F("hDiffRZPhi_PhiR", "Difference in RZ for CM Model Phi, Phi,R binning; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hCartesianAveDiffPhiR[6];
    hCartesianAveDiffPhiR[0] = new TH2F("hAveDiffXYX_PhiR", "X Model - Truth Averaged Over z, Phi,R binning (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiR[1] = new TH2F("hAveDiffRZX_PhiR", "X Model - Truth Averaged Over phi, Phi,R binning (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianAveDiffPhiR[2] = new TH2F("hAveDiffXYY_PhiR", "Y Model - Truth Averaged Over z, Phi,R binning (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiR[3] = new TH2F("hAveDiffRZY_PhiR", "Y Model - Truth Averaged Over phi, Phi,R binning (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianAveDiffPhiR[4] = new TH2F("hAveDiffXYZ_PhiR", "Z Model - Truth Averaged Over z, Phi,R binning (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiR[5] = new TH2F("hAveDiffRZZ_PhiR", "Z Model - Truth Averaged Over phi, Phi,R binning (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    
     TH2F *hCylindricalAveDiffPhiR[4];
    hCylindricalAveDiffPhiR[0] = new TH2F("hAveDiffXYR_PhiR", "R Model - Truth Averaged Over z, Phi,R binning (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalAveDiffPhiR[1] = new TH2F("hAveDiffRZR_PhiR", "R Model - Truth Averaged Over phi, Phi,R binning (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCylindricalAveDiffPhiR[2] = new TH2F("hAveDiffXYPhi_PhiR", "Phi Model - Truth Averaged Over z, Phi,R binning (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalAveDiffPhiR[3] = new TH2F("hAveDiffRZPhi_PhiR", "Phi Model - Truth Averaged Over phi, Phi,R binning (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

    TH2F *hSamplePerBinRZ = new TH2F("hSamplePerBinRZ", "Filling each rz bin; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

    TH2F *hSamplePerBinPhiR = new TH2F("hSamplePerBinPhiR", "Filling each PhiR bin; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);

    TH2F *hCompareRTrue_PhiR = new TH2F("hCompareRTrue_PhiR", "Compare Difference from R Model and True, Phi,R binning (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);
    TH2F *hComparePhiTrue_PhiR = new TH2F("hComparePhiTrue_PhiR", "Compare Difference from Phi Model and True, Phi,R binning (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);

    TH2F *hRDiffvR_PhiR = new TH2F("hRDiffvR_PhiR", "Difference between R Model and True vs. r, Phi,R binning (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvZ_PhiR = new TH2F("hRDiffvZ_PhiR", "Difference between R Model and True vs. z, Phi,R binning (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvPhi_PhiR = new TH2F("hRDiffvPhi_PhiR", "Difference between R Model and True vs. phi, Phi,R binning (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    TH2F *hPhiDiffvR_PhiR = new TH2F("hPhiDiffvR_PhiR", "Difference between Phi Model and True vs. r, Phi,R binning (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvZ_PhiR = new TH2F("hPhiDiffvZ_PhiR", "Difference between Phi Model and True vs. z, Phi,R binning (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvPhi_PhiR = new TH2F("hPhiDiffvPhi_PhiR", "Difference between Phi Model and True vs. phi, Phi,R binning (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);
    
    for(int i = 1; i < nphi - 1; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 1; j < nr - 1; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
	for(int k = 1; k < nz - 1; k++){
	  double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	  double shifttrueCart[3];
	  double shifttrueCyl[2];

	  double shiftrecoCartPhiR[3];
	  double differenceCartPhiR[3];

	  double shiftrecoCylPhiR[2];
	  double differenceCylPhiR[2];

	  double differenceR_PhiR, differencePhi_PhiR;	  

	  int binPhiR = hCartCMModelPhiR[0]->FindBin(phi,r,z);

	  if((r > 30.0) && (r < 76.0)){
	    //x y and z
	    shifttrueCart[0] = (shifter->hX->Interpolate(phi,r,z))*(1e4); //convert from cm to micron
	    shifttrueCart[1] = (shifter->hY->Interpolate(phi,r,z))*(1e4); //convert from cm to micron 
	    shifttrueCart[2] = (shifter->hZ->Interpolate(phi,r,z))*(1e4); //convert from cm to micron
	    //r and phi
	    shifttrueCyl[0] = (shifter->hR->Interpolate(phi,r,z))*(1e4); //convert from cm to micron
	    shifttrueCyl[1] = (shifter->hPhi->Interpolate(phi,r,z))*(1e4);
	    hRShiftTrue->Fill(shifttrueCyl[0]);
	    hPhiShiftTrue->Fill(shifttrueCyl[1]);
	    
	    for(int l = 0; l < 3; l ++){
	      shiftrecoCartPhiR[l] =  (hCartCMModelPhiR[l]->GetBinContent(binPhiR))*(1e4);
	  
	      differenceCartPhiR[l] = shiftrecoCartPhiR[l] - shifttrueCart[l]; 

	      hCartesianShiftDifferencePhiR[l]->Fill(differenceCartPhiR[l]);
	    }

	    //r
	    shiftrecoCylPhiR[0] =  (hCylCMModelPhiR[0]->GetBinContent(binPhiR))*(1e4);
	    differenceCylPhiR[0] = shiftrecoCylPhiR[0] - shifttrueCyl[0]; 
	    hCylindricalShiftDifferencePhiR[0]->Fill(differenceCylPhiR[0]);
	      
	    //phi
	    shiftrecoCylPhiR[1] = r*(1e4)*(hCylCMModelPhiR[1]->GetBinContent(binPhiR));
	    differenceCylPhiR[1] = (shiftrecoCylPhiR[1] - shifttrueCyl[1]); 
	    hCylindricalShiftDifferencePhiR[1]->Fill(differenceCylPhiR[1]);

	    //x
	    hCartesianDiffPhiR[0]->Fill(phi,r, differenceCartPhiR[0]);
	    hCartesianDiffPhiR[1]->Fill(z,r, differenceCartPhiR[0]);
	    //y
	    hCartesianDiffPhiR[2]->Fill(phi,r, differenceCartPhiR[1]);	  
	    hCartesianDiffPhiR[3]->Fill(z,r, differenceCartPhiR[1]);
	    //z
	    hCartesianDiffPhiR[4]->Fill(phi,r, differenceCartPhiR[2]);
	    hCartesianDiffPhiR[5]->Fill(z,r, differenceCartPhiR[2]);

	    //r
	    hCylindricalDiffPhiR[0]->Fill(phi,r, differenceCylPhiR[0]);
	    hCylindricalDiffPhiR[1]->Fill(z,r, differenceCylPhiR[0]);

	    hCompareRTrue_PhiR->Fill(shiftrecoCylPhiR[0],shifttrueCyl[0]);

	    hRDiffvR_PhiR->Fill(r,differenceCylPhiR[0],1);
	    hRDiffvPhi_PhiR->Fill(phi,differenceCylPhiR[0],1);
	    hRDiffvZ_PhiR->Fill(z,differenceCylPhiR[0],1);

	    //phi 
	    hCylindricalDiffPhiR[2]->Fill(phi,r, differenceCylPhiR[1]);
	    hCylindricalDiffPhiR[3]->Fill(z,r, differenceCylPhiR[1]);
	    	    
	    hComparePhiTrue_PhiR->Fill(shiftrecoCylPhiR[1],shifttrueCyl[1]);
	    
	    hPhiDiffvR_PhiR->Fill(r,differenceCylPhiR[1],1);
	    hPhiDiffvPhi_PhiR->Fill(phi,differenceCylPhiR[1],1);
	    hPhiDiffvZ_PhiR->Fill(z,differenceCylPhiR[1],1);

	    hSamplePerBinRZ->Fill(z,r,1);
	    hSamplePerBinPhiR->Fill(phi,r,1);
	  }
	}
      }
    }

    //average over z
    for (int m = 0; m < 6; m = m+2){
      hCartesianAveDiffPhiR[m]->Divide(hCartesianDiffPhiR[m],hSamplePerBinPhiR);
    }
    for (int m = 0; m < 4; m = m+2){
      hCylindricalAveDiffPhiR[m]->Divide(hCylindricalDiffPhiR[m],hSamplePerBinPhiR);
    }
    
    //average over phi
    for (int m = 1; m < 6; m = m+2){
      hCartesianAveDiffPhiR[m]->Divide(hCartesianDiffPhiR[m],hSamplePerBinRZ);
    }
    for (int m = 1; m < 4; m = m+2){
      hCylindricalAveDiffPhiR[m]->Divide(hCylindricalDiffPhiR[m],hSamplePerBinRZ);
    }

    //summary plots
    hDifferenceMeanR->Fill(hCylindricalShiftDifferencePhiR[0]->GetMean(1));
    hDifferenceStdDevR->Fill(hCylindricalShiftDifferencePhiR[0]->GetStdDev(1));

    hTrueMeanR->Fill(hRShiftTrue->GetMean(1));
    hTrueStdDevR->Fill(hRShiftTrue->GetStdDev(1));
    
    hDifferenceMeanPhi->Fill(hCylindricalShiftDifferencePhiR[1]->GetMean(1));
    hDifferenceStdDevPhi->Fill(hCylindricalShiftDifferencePhiR[1]->GetStdDev(1));

    hTrueMeanPhi->Fill(hPhiShiftTrue->GetMean(1));
    hTrueStdDevPhi->Fill(hPhiShiftTrue->GetStdDev(1));

    for (int m = 0; m < 6; m++){
      hCartesianAveDiffPhiR[m]->SetStats(0);
    }
    for (int m = 0; m < 4; m++){
      hCylindricalAveDiffPhiR[m]->SetStats(0);
    }
  

    hCompareRTrue_PhiR->SetStats(0);
    hComparePhiTrue_PhiR->SetStats(0);

    hRDiffvR_PhiR->SetStats(0);
    hRDiffvZ_PhiR->SetStats(0);
    hRDiffvPhi_PhiR->SetStats(0);
  
    hPhiDiffvR_PhiR->SetStats(0);
    hPhiDiffvZ_PhiR->SetStats(0);
    hPhiDiffvPhi_PhiR->SetStats(0);
    
    TPad *c1=new TPad("c1","",0.0,0.8,1.0,0.93); //can i do an array of pads?
    TPad *c2=new TPad("c2","",0.0,0.64,1.0,0.77);
    TPad *c3=new TPad("c3","",0.0,0.48,1.0,0.61);
    TPad *c4=new TPad("c4","",0.0,0.32,1.0,0.45);
    TPad *c5=new TPad("c5","",0.0,0.16,1.0,0.29);
    TPad *c6=new TPad("c6","",0.0,0.0,1.0,0.13);
    
    TPad *titlepad=new TPad("titlepad","",0.0,0.96,1.0,1.0);

    TPad *stitlepad1=new TPad("stitlepad1","",0.0,0.93,1.0,0.96);
    TPad *stitlepad2=new TPad("stitlepad2","",0.0,0.77,1.0,0.8);
    TPad *stitlepad3=new TPad("stitlepad3","",0.0,0.61,1.0,0.64);
    TPad *stitlepad4=new TPad("stitlepad4","",0.0,0.45,1.0,0.48);
    TPad *stitlepad5=new TPad("stitlepad5","",0.0,0.29,1.0,0.32);
    TPad *stitlepad6=new TPad("stitlepad6","",0.0,0.13,1.0,0.16);
    
    TLatex * title = new TLatex(0.0,0.0,"");

    TLatex * stitle1 = new TLatex(0.0,0.0,""); //array?
    TLatex * stitle2 = new TLatex(0.0,0.0,"");
    TLatex * stitle3 = new TLatex(0.0,0.0,"");
    TLatex * stitle4 = new TLatex(0.0,0.0,"");
    TLatex * stitle5 = new TLatex(0.0,0.0,"");
    TLatex * stitle6 = new TLatex(0.0,0.0,"");
    
    title->SetNDC();
    stitle1->SetNDC();
    stitle2->SetNDC();
    stitle3->SetNDC();
    stitle4->SetNDC();
    stitle5->SetNDC();
    stitle6->SetNDC();
    
    title->SetTextSize(0.32);
    stitle1->SetTextSize(0.35);
    stitle2->SetTextSize(0.35);
    stitle3->SetTextSize(0.35);
    stitle4->SetTextSize(0.35);
    stitle5->SetTextSize(0.35);
    stitle6->SetTextSize(0.35);
    
    canvas->cd();
    c1->Draw();
    stitlepad1->Draw();
    c2->Draw();
    stitlepad2->Draw();
    c3->Draw();
    stitlepad3->Draw();
    c4->Draw();
    stitlepad4->Draw();
    c5->Draw();
    stitlepad5->Draw();
    c6->Draw();
    stitlepad6->Draw();
    titlepad->Draw();

    //x plots
    c1->Divide(4,1);
    c1->cd(1);
    hCartesianAveDiffPhiR[0]->Draw("colz");
    c1->cd(2);
    hCartesianAveDiffPhiR[1]->Draw("colz");
    c1->cd(3);
    hCartesianShiftDifferencePhiR[0]->Draw();
    //c1->cd(4)->Clear();  
    c1->cd(4);
    //hCMmodelSliceRvTrue->Draw("colz");
    hSamplePerBinRZ->Draw("colz");
    
    //y plots
    c2->Divide(4,1);
    c2->cd(1);
    hCartesianAveDiffPhiR[2]->Draw("colz");
    c2->cd(2);
    hCartesianAveDiffPhiR[3]->Draw("colz");
    c2->cd(3);
    hCartesianShiftDifferencePhiR[1]->Draw();
    //c2->cd(4)->Clear();
    c2->cd(4);
    //hStripesPerBin->Draw("colz");
    hSamplePerBinPhiR->Draw("colz");
    
    //r cart
    c3->Divide(4,1);
    c3->cd(1);
    hCylindricalAveDiffPhiR[0]->Draw("colz");
    c3->cd(2);
    hCylindricalAveDiffPhiR[1]->Draw("colz");
    c3->cd(3);
    hCylindricalShiftDifferencePhiR[0]->Draw();
    c3->cd(4);
    hRShiftTrue->Draw();
    
    //phi cart
    c4->Divide(4,1);
    c4->cd(1);
    hCylindricalAveDiffPhiR[2]->Draw("colz");
    c4->cd(2);
    hCylindricalAveDiffPhiR[3]->Draw("colz");
    c4->cd(3);
    hCylindricalShiftDifferencePhiR[1]->Draw();
    c4->cd(4);
    hPhiShiftTrue->Draw();

    //r to true comparison
    c5->Divide(4,1);
    c5->cd(1);
    hCompareRTrue_PhiR->Draw("colz");
    c5->cd(2);
    hRDiffvR_PhiR->Draw("colz");
    c5->cd(3);
    hRDiffvZ_PhiR->Draw("colz");
    c5->cd(4);
    hRDiffvPhi_PhiR->Draw("colz");

    //phi to true comparison
    c6->Divide(4,1);
    c6->cd(1);
    hComparePhiTrue_PhiR->Draw("colz");
    c6->cd(2);
    hPhiDiffvR_PhiR->Draw("colz");
    c6->cd(3);
    hPhiDiffvZ_PhiR->Draw("colz");
    c6->cd(4);
    hPhiDiffvPhi_PhiR->Draw("colz");

    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.01,0.4,Form("Event %d; %s", ifile, sourcefilename.Data())); 
    title->Draw();
    
    stitlepad1->cd();
    stitlepad1->Clear();
    stitle1->DrawLatex(0.45,0.2,"X Model"); 
    stitle1->Draw();
     
    stitlepad2->cd();
    stitlepad2->Clear();
    stitle2->DrawLatex(0.45,0.2,"Y Model"); 
    stitle2->Draw();

    stitlepad3->cd();
    stitlepad3->Clear();
    stitle3->DrawLatex(0.45,0.2,"R Model"); 
    stitle3->Draw();

    stitlepad4->cd();
    stitlepad4->Clear();
    stitle4->DrawLatex(0.45,0.2,"Phi Model"); 
    stitle4->Draw();

    stitlepad5->cd();
    stitlepad5->Clear();
    stitle5->DrawLatex(0.4,0.2,"Comparing R Model to True"); 
    stitle5->Draw();

    stitlepad6->cd();
    stitlepad6->Clear();
    stitle6->DrawLatex(0.4,0.2,"Comparing Phi Model to True"); 
    stitle6->Draw();

    if(ifile == 0){ 
      //if(ifile == 1){
      canvas->Print("CMDistortionAnalysisPhiR.pdf(","pdf");
    } else if((ifile == 1) || (ifile == nEvents - 1)){
      canvas->Print("CMDistortionAnalysisPhiR.pdf","pdf");
    }
  }

  TCanvas *summary = new TCanvas("summary","ShiftPlotsSummary",2000,3000);

  TPad *sumtitlepad = new TPad("sumtitlepad","",0.0,0.96,1.0,1.0);
  TPad *sumplots = new TPad("sumplotspad","",0.0,0.0,1.0,0.96);

  TLatex *sumtitle = new TLatex(0.0,0.0,"");

  sumtitle->SetNDC();
  sumtitle->SetTextSize(0.4);

  summary->cd();
  sumplots->Draw();
  sumtitlepad->Draw();

  sumplots->Divide(4,6);
  sumplots->cd(1);
  hDifferenceMeanR->Draw();
  sumplots->cd(2);
  hDifferenceStdDevR->Draw();
  sumplots->cd(3);
  hTrueMeanR->Draw();
  sumplots->cd(4);
  hTrueStdDevR->Draw();
  sumplots->cd(5);
  hDifferenceMeanPhi->Draw();
  sumplots->cd(6);
  hDifferenceStdDevPhi->Draw();
  sumplots->cd(7);
  hTrueMeanPhi->Draw();
  sumplots->cd(8);
  hTrueStdDevPhi->Draw();
  sumplots->cd(9);
  sumplots->cd(10)->Clear();
  sumplots->cd(11)->Clear();
  sumplots->cd(12)->Clear();
  sumplots->cd(13)->Clear();
  sumplots->cd(14)->Clear();
  sumplots->cd(15)->Clear();
  sumplots->cd(16)->Clear();
  sumplots->cd(17)->Clear();
  sumplots->cd(18)->Clear();
  sumplots->cd(19)->Clear();
  sumplots->cd(20)->Clear();
  sumplots->cd(21)->Clear();
  sumplots->cd(22)->Clear();
  sumplots->cd(23)->Clear();
  sumplots->cd(24)->Clear();

  sumtitlepad->cd();
  sumtitlepad->Clear();
  sumtitle->DrawLatex(0.4,0.4,"Summary of Events"); 
  summary->Print("CMDistortionAnalysisPhiR.pdf)","pdf");

  return 0;
}
