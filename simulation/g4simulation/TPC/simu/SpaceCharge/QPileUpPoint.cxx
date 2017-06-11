#include "QPileUpPoint.h"
#include "Constants.h"
#include <string>

#include "TH3F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TMath.h"

using namespace Constants;

//=====================
QPileUpPoint::QPileUpPoint(float r, float p, float z, float q) {
  QPileUp();
  fRper = r;
  fPper = p;
  fZper = z;
  fQabs = q;
}
//=====================
void QPileUpPoint::Make() {
  InitMaps();

  //------------
  //STEER
  if(fDebug>0) printf("QPileUp is being computed with Point-like source... \n");
  int R = fRho->GetXaxis()->GetNbins()*fRper;
  int P = fRho->GetYaxis()->GetNbins()*fPper;
  int Z = fRho->GetZaxis()->GetNbins()*fZper;
  if(fDebug>1) printf("binR %d\n",R);
  if(fDebug>1) printf("binP %d\n",P);
  if(fDebug>1) printf("binZ %d\n",Z);
  for(int r=0; r!=kNRadialSteps; ++r)
    for(int p=0; p!=kNAzimuthalSteps; ++p)
      for(int z=0; z!=kNLongitudinalSteps; ++z) {
	float dRho = 0;
	if(R==r&&P==p&&Z==z) dRho = fQabs;
	fRho->SetBinContent(r+1,p+1,z+1,dRho);
	if(fDebug>2) printf("@{Ir,Ip,Iz}={%d,%d,%d}, rho %f\n",r,p,z,dRho);
      }
  if(fDebug>0) printf("[DONE]\n");
  //------------

  SaveMaps();
}
