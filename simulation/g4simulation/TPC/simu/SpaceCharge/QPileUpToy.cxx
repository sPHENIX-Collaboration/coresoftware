#include "QPileUpToy.h"
#include <string>

#include "TH3F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TMath.h"

//=====================
QPileUpToy::QPileUpToy(float gf, float mu, float rt, float eps, float rad) {
  fGasFactor = gf;
  fMultiplicity = mu;
  fDAQRate = rt;
  fEPS = eps;
  fRad = rad;
}
//=====================
void QPileUpToy::Make() {
  InitMaps();

  //------------
  //STEER
  if(fDebug>0) printf("QPileUp is being computed from TOY MODEL ... ");
  float e0 = 8.854187817e-3; //[C]/[Vm]*1e+9
  double a=fMultiplicity*fDAQRate*e0*fGasFactor; // [fC]/[m]*1e+9;
  float b=100.0/fHalfLength; //[1/m]
  float c=2.0/3.0*fEPS;
  float d=-fRad;
  float f2 = TMath::Power(fOutterRadius/100,-1)/(-1) - TMath::Power(fInnerRadius/100,-1)/(-1);
  float fn = TMath::Power(fOutterRadius/100,d+1)/(d+1) - TMath::Power(fInnerRadius/100,d+1)/(d+1);
  float fd = f2/fn;
  if(fDebug>1) printf("\n2PI*rho(r,z) = a fd (1 - b|z| + c) r^d\n");
  if(fDebug>1) printf("a = %f\n",a);
  if(fDebug>1) printf("b = %f\n",b);
  if(fDebug>1) printf("c = %f\n",c);
  if(fDebug>1) printf("fd = %f\n",fd);
  if(fDebug>1) printf("d = %f\n",d);
  for(int r=0; r!=fNRadialSteps; ++r) {
    float dr = fRho->GetXaxis()->GetBinCenter( r+1 )/100.0; //[m]
    for(int z=0; z!=fNLongitudinalSteps; ++z) {
      float dz = fRho->GetZaxis()->GetBinCenter( z+1 )/100.0; //[m]
      for(int p=0; p!=fNAzimuthalSteps; ++p) {
	float dp = fRho->GetYaxis()->GetBinCenter( p+1 );
	float dRho = a*(1-b*TMath::Abs(dz)+c)*fd*TMath::Power(dr,d); //[fC]/[cm^3]
	fRho->SetBinContent(r+1,p+1,z+1,dRho);
	if(fDebug>2) printf("@{Ir,Ip,Iz}={%d (%f),%d (%f),%d (%f)}, rho %f\n",r,dr,p,dp,z,dz,dRho);
      }
    }
  }
  if(fDebug>0) printf("[DONE]\n");
  //------------

  SaveMaps();
}
