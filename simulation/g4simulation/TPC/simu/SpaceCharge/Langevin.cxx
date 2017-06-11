#include "Langevin.h"
#include <string>

#include "TH3F.h"
#include "TFile.h"
#include "TMath.h"

//=====================
Langevin::Langevin() {
  fDebug = 0;
  fEr = NULL;
  fEp = NULL;
  fEz = NULL;
  fDeltaR = NULL;
  fRDeltaPHI = NULL;
  fMirrorZ=false;
  fInnerRadius=1.;
  fOutterRadius=2.;
  fHalfLength=1.;
  fNRadialSteps=1;
  fNAzimuthalSteps=1;
  fNLongitudinalSteps=1;
  fFileNameRoot="rho";
}
//=====================
Langevin::~Langevin() {
}
//=====================
void Langevin::Make() {
  InitMaps();

  //------------
  //STEER
  if(fEr==NULL) {
    printf("Input fields not found. A B O R T I N G\n");
    return;
  }
  if(fDebug>0) printf("Langevin solutions are being computed... \n");

  const float prec=1e-8;
  const float kE0 = 400; //V/cm
  //const float kB0 = 5; //kGauss ALICE / STAR
  const float kB0 = 15; //kGauss PHENIX
  //Langevin gas parameters
  // P9 gas
  //const float kT1 = 1.34; //NIM Phys Res A235,296 (1985)
  //const float kT2 = 1.11; //NIM Phys Res A235,296 (1985)
  // P10 gas
  //const float kT1 = 1.36; //measured by STAR
  //const float kT2 = 1.11; //measured by STAR
  // 85.7%Ne / 9.5%CO2 / 4.5%N2
  const float kT1 = 1.0; // measured by ALICE
  const float kT2 = 1.0; // measured by ALICE
  //const float kV0 = 4.0; // [cm/us] @ E=400V/cm (need fine tunning)
  const float kV0 = 6.0; // [cm/us] @ E=400V/cm (need fine tunning)
  float wt = -10*kB0*kV0/kE0; //[kGauss] * [cm/us] / [V/cm] //0.32 != -10*5*4/400
  float c0 = 1.0/(1.0+kT2*kT2*wt*wt);
  float c1 = kT1*wt/(1.0+kT1*kT1*wt*wt);
  float c2 = kT2*kT2*wt*wt*c0;
  if(fDebug>1) printf("wt = %f \n",wt);
  if(fDebug>1) printf("c0 = %f \n",c0);
  if(fDebug>1) printf("c1 = %f \n",c1);
  if(fDebug>1) printf("c2 = %f \n",c2);
  float dz = fEz->GetZaxis()->GetBinWidth(1);
  for(int r=0; r!=fNRadialSteps; ++r)
    for(int p=0; p!=fNAzimuthalSteps; ++p)
      for(int z=0; z!=fNLongitudinalSteps; ++z) {
	float xcurrz = fEz->GetZaxis()->GetBinCenter(z+1);
	if(fMirrorZ && xcurrz>0) continue;
	float intDR = 0, intRDPHI = 0;
	// integrating over drift ditance
	int minZ = 0;
	int maxZ = z+1;
	if( xcurrz>0 ) {
	  minZ=z;
	  maxZ=fNLongitudinalSteps;
	}
	for(int zprime=minZ; zprime!=maxZ; ++zprime) {
	  float dR, RdPHI;
	  float dEz = fEz->GetBinContent( r+1, p+1, z+1 ) + kE0;
	  dR = c0*fEr->GetBinContent( r+1, p+1, zprime+1 )/dEz;
	  dR += c1*fEp->GetBinContent( r+1, p+1, zprime+1 )/dEz;
	  RdPHI = -c1*fEr->GetBinContent( r+1, p+1, zprime+1 )/dEz;
	  RdPHI += c0*fEp->GetBinContent( r+1, p+1, zprime+1 )/dEz;
	  if(!TMath::AreEqualAbs( dR, 0, prec )) intDR += dR*dz;
	  if(!TMath::AreEqualAbs( RdPHI, 0, prec )) intRDPHI += RdPHI*dz;
	}
	fDeltaR->SetBinContent( r+1, p+1, z+1, intDR );
	fRDeltaPHI->SetBinContent( r+1, p+1, z+1, intRDPHI );
	if(fMirrorZ) {
	  fDeltaR->SetBinContent( r+1, p+1, fNLongitudinalSteps-z, intDR );
	  fRDeltaPHI->SetBinContent( r+1, p+1, fNLongitudinalSteps-z, intRDPHI );
	}
	if(fDebug>2) printf("@{Ir,Ip,Iz}={%d,%d,%d}, deltaR=%f\n",r,p,z,intDR);
	if(fDebug>2) printf("@{Ir,Ip,Iz}={%d,%d,%d}, RdeltaPHI=%f\n",r,p,z,intRDPHI);
      }
  if(fDebug>0) printf("[DONE]\n");
  //------------

  SaveMaps();
}
//=====================
void Langevin::ReadFile() {
  const char *inputfile = Form("%s_1.root",fFileNameRoot.data());
  if(fDebug>0) printf("Langevin is reading the input file %s... ",inputfile);
  TFile *ifile = new TFile(inputfile);
  fEr = (TH3F*) ifile->Get("Er");
  fEp = (TH3F*) ifile->Get("Ep");
  fEz = (TH3F*) ifile->Get("Ez");
  if(fDebug>0) printf("[DONE]\n");
}
//=====================
void Langevin::InitMaps() {
  ReadFile();
  if(fDebug>0) printf("Langevin is initializing distortion maps... ");
  fDeltaR = new TH3F("mapDeltaR","#delta_{R} [cm];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
		     fNRadialSteps,fInnerRadius,fOutterRadius,
		     fNAzimuthalSteps,0,TMath::TwoPi(),
		     fNLongitudinalSteps,-fHalfLength,+fHalfLength);
  fRDeltaPHI = new TH3F("mapRDeltaPHI","R#delta_{#varphi} [cm];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
		     fNRadialSteps,fInnerRadius,fOutterRadius,
		     fNAzimuthalSteps,0,TMath::TwoPi(),
		     fNLongitudinalSteps,-fHalfLength,+fHalfLength);
  if(fDebug>0) printf("[DONE]\n");
}
//=====================
void Langevin::SaveMaps() {
  const char *outputfile= Form("%s_2.root",fFileNameRoot.data());
  if(fDebug>0) printf("Langevin saving distortion maps... ");
  TFile *ofile = new TFile(outputfile,"RECREATE");
  ofile->WriteObject(fDeltaR,"mapDeltaR");
  ofile->WriteObject(fRDeltaPHI,"mapRDeltaPHI");
  ofile->Close();
  if(fDebug>0) printf("[DONE]\n");
}
