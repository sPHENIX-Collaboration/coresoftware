#include "FieldMapsLaplace.h"
#include "LaplaceSolution.h"
#include <string>
#include <iostream>
#include "TH3F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TMath.h"


//=====================
void FieldMapsLaplace::ComputeE() {
  const float prec=1e-8;

  //------------
  //STEER
  if(fDebug>0) printf("FieldMaps is being computed from Laplace solutions... \n");
  LaplaceSolution *laplace_green_solution;
  if(0) {
    laplace_green_solution = new LaplaceSolution(fLSNameRoot);
  } else {
    laplace_green_solution = new LaplaceSolution(fInnerRadius/100/*[m]*/,fOutterRadius/100/*[m]*/,fHalfLength/100/*[m]*/);
  }
  int brprimeINI=0;
  int brprimeEND=fNRadialSteps;
  if(fRadialBin>-0.5) {
    brprimeINI = fRadialBin;
    brprimeEND = fRadialBin+1;
  }
  // constructing volume element
  float dr = (fOutterRadius-fInnerRadius)/fNRadialSteps;
  float dphi = TMath::TwoPi()/fNAzimuthalSteps;
  float dz = 2.0*fHalfLength/fNLongitudinalSteps;
  float dV = dr*dphi*dz;

  float e0 = 8.854187817e+1; //[fC]/[V.cm]
  // looping over volume
  for(int br=0; br!=fNRadialSteps; ++br) {
    float r = fEr->GetXaxis()->GetBinCenter( br+1 ); //[cm]
    if(fDebug>0) printf("FieldMaps: %.2f -> 1\n",float(br)/fNRadialSteps);
    for(int bp=0; bp!=fNAzimuthalSteps; ++bp) {
      float phi = fEr->GetYaxis()->GetBinCenter( bp+1 ); //[rad]
      for(int bz=0; bz!=fNLongitudinalSteps; ++bz) {
	float z = fEr->GetZaxis()->GetBinCenter( bz+1 ); //[cm]
	//std::cout << " z :" << z << std::endl;
	if(fDebug>1) printf("FieldMaps: Integral( dvprime Q(rprime) G_E(rprime,[%f,%f,%f]) )",r,phi,z);
	float intEr = 0, intEp = 0, intEz = 0;
	if(fMirrorZ) if(z>0) continue;
	//std::cout << " z :" << z << std::endl;
	// looping over volume prime
	for(int bzprime=0; bzprime!=fNLongitudinalSteps; ++bzprime) {
	  float zprime = fEr->GetZaxis()->GetBinCenter( bzprime+1 ); //[cm]
	  //std::cout << " z*z' :" << z*zprime << std::endl;
	  if(z*zprime < 0) continue; // only integrates same side along Z
	  //std::cout << " z :"<< z << " z' :" << zprime << " z*z' :" << z*zprime << std::endl;
	  for(int brprime=brprimeINI; brprime!=brprimeEND; ++brprime) {
	    float rprime = fEr->GetXaxis()->GetBinCenter( brprime+1 ); //[cm]
	    for(int bphiprime=0; bphiprime!=fNAzimuthalSteps; ++bphiprime) {
	      float phiprime = fEr->GetYaxis()->GetBinCenter( bphiprime+1 );
	      //float phiprime = fEr->GetYaxis()->GetBinCenter(1);
	      float charge = ReadCharge(rprime,phiprime,zprime,dr,dphi,dz); // fC
	      if(TMath::AreEqualAbs( charge, 0, prec )) continue;
	      float GR = laplace_green_solution->Er(r/100,phi,TMath::Abs(z)/100, rprime/100, phiprime, TMath::Abs(zprime)/100); // GRad
	      float GP = 0;//laplace_green_solution->Ephi(r/100,phi,TMath::Abs(z)/100, rprime/100, phi, TMath::Abs(zprime)/100); // GPhi
	      float GZ = laplace_green_solution->Ez(r/100,phi,TMath::Abs(z)/100, rprime/100, phiprime, TMath::Abs(zprime)/100); // GPhi
	      if(!TMath::AreEqualAbs( GR, 0, prec )) intEr += GR*charge; // [fC/m^2]	    
	      //std::cout << " z :"<< z << " z' :" << zprime << " GR :" << GR << " q :" << charge <<" intEr :" << intEr<< std::endl;
	      if(!TMath::AreEqualAbs( GP, 0, prec )) intEp += GP*charge; // [fC/m^2]
	      if(!TMath::AreEqualAbs( GZ, 0, prec )) intEz += GZ*charge; // [fC/m^2]
	    }//loop over phi'
	  }// loop over r'
	}//loop over z'
	if(fDebug>1) printf(" => Er=%f ( similar for Ep=%f, Ez=%f )\n",intEr,intEp,intEz);
	float toVperCM = 0.155*1e-4/e0;
	fEr->Fill(r,phi,z,intEr*toVperCM); //[V/cm]
	fEp->Fill(r,phi,z,intEp*toVperCM); //[V/cm]
	fEz->Fill(r,phi,z,intEz*toVperCM); //[V/cm]
      } // loop over z
    } // loop over phi
  } //loop over r
  if(fDebug>0) printf("[DONE]\n");
  //------------
}
