#include "FieldMaps.h"
#include <string>

#include "TH3F.h"
#include "TFile.h"
#include "TMath.h"

//=====================
FieldMaps::FieldMaps() {
  fDebug = 0;
  fRadialBin = -1;
  fEr = NULL;
  fEp = NULL;
  fEz = NULL;
  fRho = NULL;
  fInnerRadius=1.;
  fOutterRadius=2.;
  fHalfLength=1.;
  fNRadialSteps=1;
  fNAzimuthalSteps=1;
  fNLongitudinalSteps=1;
  fFileNameRoot="rho";
  fLSNameRoot="none";
  fMirrorZ=false;
}
//=====================
FieldMaps::~FieldMaps() {
}
//=====================
void FieldMaps::Make(int n) {
  fRadialBin = n;
  InitMaps();
  ComputeE();
  SaveMaps();
}
//=====================
void FieldMaps::InitMaps() {
  if(fDebug>0) printf("FieldMaps is initializing field maps... ");
  fEr = new TH3F("Er","Er [V/cm];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
		 fNRadialSteps,fInnerRadius,fOutterRadius,
		 fNAzimuthalSteps,0,TMath::TwoPi(),
		 fNLongitudinalSteps,-fHalfLength,+fHalfLength);
  fEp = new TH3F("Ep","Ephi [V/cm];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
		 fNRadialSteps,fInnerRadius,fOutterRadius,
		 fNAzimuthalSteps,0,TMath::TwoPi(),
		 fNLongitudinalSteps,-fHalfLength,+fHalfLength);
  fEz = new TH3F("Ez","Ez [V/cm];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
		 fNRadialSteps,fInnerRadius,fOutterRadius,
		 fNAzimuthalSteps,0,TMath::TwoPi(),
		 fNLongitudinalSteps,-fHalfLength,+fHalfLength);
  if(fDebug>0) printf("[DONE]\n");
  ReadFile();
}
//=====================
void FieldMaps::ReadFile() {
  const char *inputfile= Form("%s_0.root",fFileNameRoot.data());
  if(fDebug>0) printf("FieldMaps is reading from %s... ",inputfile);
  TFile *ifile = new TFile(inputfile);
  if(ifile->IsZombie()) ifile=NULL;
  if(!ifile) exit(1);
  fRho = (TH3F*) ifile->Get("rho");
  if(fDebug>0) printf("[DONE]\n");
  // PLUG Gr Gz HERE
}
//=====================
void FieldMaps::SaveMaps() {
  TString filename;
  if(fRadialBin<0) filename = Form("%s_1.root",fFileNameRoot.data());
  else filename = Form("%s_1_%d.root",fFileNameRoot.data(),fRadialBin);
  if(fDebug>1) printf("FieldMaps saving field maps into %s... ",filename.Data());
  TFile *ofile = new TFile(filename.Data(),"RECREATE");
  ofile->WriteObject(fEr,"Er");
  ofile->WriteObject(fEp,"Ep");
  ofile->WriteObject(fEz,"Ez");
  ofile->Close();
  if(fDebug>1) printf("[DONE]\n");
}
//=====================
float FieldMaps::ReadCharge(float rprime, float phiprime, float zprime, float dr, float dphi,float dz) {
  int brmin = fRho->GetXaxis()->FindBin(rprime - dr/2 + 1e-8);
  int brmax = fRho->GetXaxis()->FindBin(rprime + dr/2 - 1e-8);
  int bpmin = fRho->GetYaxis()->FindBin(phiprime - dphi/2 + 1e-8);
  int bpmax = fRho->GetYaxis()->FindBin(phiprime + dphi/2 - 1e-8);
  int bzmin = fRho->GetZaxis()->FindBin(zprime - dz/2 + 1e-8);
  int bzmax = fRho->GetZaxis()->FindBin(zprime + dz/2 - 1e-8);
  if(brmin<1||bpmin<1||bzmin<1) return 0;
  if(brmax<brmin) brmax = brmin +1;
  if(bpmax<bpmin) bpmax = bpmin +1;
  if(bzmax<bzmin) bzmax = bzmin +1;

  float ddr = fRho->GetXaxis()->GetBinWidth(1);
  float ddp = fRho->GetYaxis()->GetBinWidth(1);
  float ddz = fRho->GetZaxis()->GetBinWidth(1);
  float totalQ = 0;
  for(int br=brmin; br!=brmax+1; ++br)
    for(int bp=bpmin; bp!=bpmax+1; ++bp)
      for(int bz=bzmin; bz!=bzmax+1; ++bz) {
	float rho = fRho->GetBinContent(br,bp,bz);
	float r = fRho->GetXaxis()->GetBinCenter( br );
	float dv = r*ddp*ddr*ddz;
	totalQ += rho*dv;
      }
  //return rprime*fRho->Integral( brmin,brmax , bpmin,bpmax , bzmin,bzmax );
  return totalQ;
}
