#include "QPileUp.h"
#include <string>

#include "TH3F.h"
#include "TFile.h"
#include "TMath.h"

//=====================
QPileUp::QPileUp() {
  fDebug = 0;
  fRho = NULL;
  fInnerRadius=1.;
  fOutterRadius=2.;
  fHalfLength=1.;
  fNRadialSteps=1;
  fNAzimuthalSteps=1;
  fNLongitudinalSteps=1;
  fFileNameRoot="rho";
}
//=====================
QPileUp::~QPileUp() {
}
//=====================
void QPileUp::Make() {
  InitMaps();
  SaveMaps();
}
//=====================
void QPileUp::InitMaps() {
  if(fDebug>0) printf("QPileUp is initializing rho map ... ");
  if(fDebug>1) {
    printf("\nTPC radial axis from %f to %f cm. spanned in %d steps\n",fInnerRadius,fOutterRadius,fNRadialSteps);
    printf("TPC longitudinal axis from %f to %f cm. spanned in %d steps\n",-fHalfLength,fHalfLength,fNLongitudinalSteps);
    printf("TPC azimuthal axis spanned in %d steps\n",fNAzimuthalSteps);
  }
  fRho = new TH3F("rho","ChargeDensity [fC/cm^3];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
		  fNRadialSteps,fInnerRadius,fOutterRadius,
		  fNAzimuthalSteps,0,TMath::TwoPi(),
		  fNLongitudinalSteps,-fHalfLength,+fHalfLength);
  if(fDebug>0) printf("[DONE]\n");
}
//=====================
void QPileUp::SaveMaps() {
  const char *outputfile= Form("%s_0.root",fFileNameRoot.data());
  if(fDebug>0) printf("QPileUp saving maps into %s ... ",outputfile);
  TFile *ofile = new TFile(outputfile,"RECREATE");
  ofile->WriteObject(fRho,"rho");
  ofile->Close();
  if(fDebug>0) printf("[DONE]\n");
}
