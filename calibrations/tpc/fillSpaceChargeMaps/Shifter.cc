#include "Shifter.h"
#include <TVector3.h>
#include <TFile.h>
#include <TH3F.h>
#include <TString.h>


Shifter::Shifter(TString truthfilename, TString correctionfilename){
  //load a 'truth' distortion map and, optionally, a map of a measured correction to those distortions
  //this code is currently set up to load a particular correction map that doesn't have distortions
  // in X,Y, and Z components, but rather only in R, R*Phi, and Z components.
  
  //single event distortion file
  hasTruth=false;//assume the file doesn't load correctly, until we prove otherwise.
  forward=NULL;
  hX=NULL;hY=NULL;hZ=NULL;
  if (truthfilename!=""){
    forward=TFile::Open(truthfilename,"READ"); 
    if (forward!=NULL){
      hX=(TH3F*)forward->Get("hIntDistortionX");
      hY=(TH3F*)forward->Get("hIntDistortionY");
      hZ=(TH3F*)forward->Get("hIntDistortionZ");

      //not strictly needed, but handy:
      hR=(TH3F*)forward->Get("hIntDistortionR");
      hPhi=(TH3F*)forward->Get("hIntDistortionP");

    }
  }
  if (hX!=NULL && hY!=NULL && hZ!=NULL){
    hasTruth=true;
  }

   //single event distortion file
  hasCorrection=false;//assume the file doesn't load correctly, until we prove otherwise.
  average=NULL;
  hXave=NULL;hYave=NULL;hZave=NULL;
  if (correctionfilename!=""){
    //average=TFile::Open(correctionfilename,"READ");
     average=TFile::Open("/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps_rec/Distortions_full_realistic_micromegas_all-coarse.root","READ");
    if (average!=NULL){
      hXave=(TH3F*)average->Get("hIntDistortionX");
      hYave=(TH3F*)average->Get("hIntDistortionY");
      hZave=(TH3F*)average->Get("hIntDistortionZ");
     
  
  hRave=(TH3F*)average->Get("hIntDistortionR");
  hPhiave=(TH3F*)average->Get("hIntDistortionP");
    }
  }
  if (hXave!=NULL && hYave!=NULL && hZave!=NULL){
    hasCorrection=true;
  }

}

TVector3 Shifter::ShiftForward(TVector3 position){
double x, y, z, xshift, yshift, zshift;
  //const double mm = 1.0;
  //const double cm = 10.0;
  TVector3 shiftposition;

  x= position.X();
  y= position.Y();
  z= position.Z();

  double r=position.Perp();
  double phi=position.Phi();
  if(position.Phi() < 0.0){
    phi = position.Phi() + 2.0*TMath::Pi(); 
  }

  //distort coordinate of stripe
  xshift=0;
  yshift=0;
  zshift=0;
  if (hasTruth){
    xshift=hX->Interpolate(phi,r,z);
    yshift=hY->Interpolate(phi,r,z);
    zshift=hZ->Interpolate(phi,r,z);
  }
  
  //remove average distortion
  if (hasCorrection){
    double raveshift=hRave->Interpolate(phi,r,z);
    double paveshift=hPhiave->Interpolate(phi,r,z); //hugo confirms the units are cm
    double cosphi=cos(phi);
    double sinphi=sin(phi);
    xshift-=raveshift*cosphi-paveshift*sinphi;
    yshift-=raveshift*sinphi+paveshift*cosphi;
    
    zshift-=hZave->Interpolate(phi,r,z);
  }

  TVector3 forwardshift(x+xshift,y+yshift,z+zshift);

  return forwardshift;
}

TVector3 Shifter::ShiftBack(TVector3 forwardshift){
double x, y, z;
 // const double mm = 1.0;
  //const double cm = 10.0;
  TVector3 shiftposition;

  x= forwardshift.X();
  y= forwardshift.Y();
  z= forwardshift.Z();

  double rforward=forwardshift.Perp();
  double phiforward=forwardshift.Phi();
  if(forwardshift.Phi() < 0.0){
    phiforward += 2.0*TMath::Pi();
  }
  
  double xshiftback=-1*hXBack->Interpolate(phiforward,rforward,z);
  double yshiftback=-1*hYBack->Interpolate(phiforward,rforward,z);
  double zshiftback=-1*hZBack->Interpolate(phiforward,rforward,z);
    
  shiftposition.SetXYZ(x+xshiftback,y+yshiftback,z+zshiftback);

  return shiftposition;
}

TVector3 Shifter::Shift(TVector3 position){
  
  return ShiftBack(ShiftForward(position));
}
