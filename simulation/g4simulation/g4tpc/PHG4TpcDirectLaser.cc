#include "PHG4TpcDirectLaser.h"

#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TH2F.h"


//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
//R__LOAD_LIBRARY(libphg4hit.so)

// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

using namespace std;

//PHG4TpcDirectLaser::

PHG4TpcDirectLaser::PHG4TpcDirectLaser()
{
  begin_CM = 20. *cm; // outer radius of IFC
  end_CM = 78. * cm; // inner radius of OFC
  halfwidth_CM=0.5*cm;

  for (int i=0;i<4*2;i++){
    //PHG4Hits.push_back(new PHG4Hitv1());//instantiate it
    AppendLaserTrack(TMath::Pi()/180.*10.,TMath::Pi()/180.*90,i); //add random starting laser tracks...
  }
  
  
 return;
}


void PHG4TpcDirectLaser::SetPhiStepping(int n, float min,float max){
  if (n<0 || max<min){
    printf("SetPhiStepping values weird.  not setting.\n");
    return;
  }
  nPhiSteps=n;
  minPhi=min;
  maxPhi=max;
  nTotalSteps=nThetaSteps*nPhiSteps;
  return;
}
void PHG4TpcDirectLaser::SetThetaStepping(int n, float min,float max){
  if (n<0 || max<min){
    printf("SetThetaStepping values weird.  not setting.\n");
    return;
  }
  nThetaSteps=n;
  minTheta=min;
  maxTheta=max;
  nTotalSteps=nThetaSteps*nPhiSteps;

  return;
}
void PHG4TpcDirectLaser::AimToThetaPhi(float theta, float phi){
  ClearHits();
  for (int i=0;i<4*2;i++){
    AppendLaserTrack(theta,phi,i);
  }
return;
}
void PHG4TpcDirectLaser::AimToPatternStep(int n){
  n=n%nTotalSteps;//trim against overflows
  currentPatternStep=n;
  int phiStep=n%nThetaSteps;
  int thetaStep=n/nPhiSteps;
  AimToThetaPhi((maxTheta-minTheta)/(nThetaSteps*1.)*thetaStep,(maxPhi-minPhi)/(nPhiSteps*1.)*phiStep);
  return;
}

void PHG4TpcDirectLaser::ClearHits(){
  for (int i =0; i< (int)(PHG4Hits.size());i++)
   {
     delete (PHG4Hits[i]);
   } 
   PHG4Hits.clear();
   return;
}

TVector3 PHG4TpcDirectLaser::GetCmStrike(TVector3 start, TVector3 direction){
  TVector3 ret;
  float end=halfwidth_CM;
  if (start.Z()<0) end=-halfwidth_CM;
  float dist=end-start.Z();
  float direction_scale=dist/direction.Z();
  ret=start+direction*direction_scale;
  
  
  return ret;
}
TVector3 PHG4TpcDirectLaser::GetFieldcageStrike(TVector3 start, TVector3 direction){
  TVector3 ofc_strike=GetCylinderStrike(start,direction,end_CM);
  TVector3 ifc_strike=GetCylinderStrike(start,direction,begin_CM);
  TVector3 no_strike(999,999,999);

  //measure which one occurs 'first' along the track, by dividing by the trajectory z.

  float ifc_dist=(ifc_strike.Z()-start.Z())/direction.Z();
  float ofc_dist=(ofc_strike.Z()-start.Z())/direction.Z();

  if (ifc_dist<0){
    if (ofc_dist>0) return ofc_strike;
    return no_strike;
  }
  //now ifc strike is guaranteed to be positive or zero
  if (ofc_dist<0){
    return ifc_strike;
  }
  //now they're both positive
  if (ifc_dist<ofc_dist){
    return ifc_strike;
  } else return ofc_strike;
  
}

TVector3  PHG4TpcDirectLaser::GetCylinderStrike(TVector3 s, TVector3 v, float radius){


  float R2=radius*radius;
		
  //Generalized Parameters for collision with cylinder of radius R:
  //from quadratic formula solutions of when a vector intersects a circle:
  float a = v.x()*v.x()+v.y()*v.y();
  float b = 2*(v.x()*s.x()+v.y()*s.y());
  float c = s.x()*s.x()+s.y()*s.y()-R2;

  float rootterm=b*b-4*a*c;

  //if a==0 then we are parallel and will have no solutions.
  //if the rootterm is negative, we will have no real roots -- we are outside the cylinder and pointing skew to the cylinder such that we never cross.
  float t1=-1,t2=-1; //this is the distance, in units of v, we must travel to find a collision
  if (rootterm >= 0 && a > 0) {
    //Find the (up to) two points where we collide with the cylinder:
    float sqrtterm=sqrt(rootterm);
     t1 = (-b+sqrtterm)/(2*a);
     t2 = (-b-sqrtterm)/(2*a);
  }

  //if either of the t's are nonzero, we have a collision.  the collision closest to the start (hence with the smallest t that is greater than zero) is the one that happens.
  float min_t=t1;
  if (t2>t1 && t2>0) min_t=t2;
  TVector3 ret=s+v*min_t;
  return ret;
}
  

void PHG4TpcDirectLaser::AppendLaserTrack(float theta, float phi, int laser)
{ //this function generates a series of PHG4 hits from the specified laser, in the specified direction.
  PHG4Hitv1 *hit;
  float direction=1;
  if (laser>3) direction=-1;
  TVector3 pos(60.*cm,0.*cm,105.5*cm*direction);
  TVector3 dir(0.,0.,-1.*direction);



  //adjust direction:
  dir.RotateY(theta*direction);
  dir.RotateZ(phi*direction);

  //rotate to the correct laser spot:
  pos.RotateZ(TMath::TwoPi()/4.*laser);
  dir.RotateZ(TMath::TwoPi()/4.*laser);

  //find collision point
  TVector3 cm_strike=GetCmStrike(pos,dir);
  TVector3 fc_strike=GetFieldcageStrike(pos,dir);
  TVector3 strike=cm_strike;
  char strikeChar='c';
  if( fc_strike.Z()!=999){
    if ((fc_strike.Z()-pos.Z())/dir.Z()<(cm_strike.Z()-pos.Z())/dir.Z()){
      strike=fc_strike;
    strikeChar='f';
    }
  }

  //find length
  TVector3 delta=(strike-pos);
  float fullLength=delta.Mag();
  int nHitSteps=fullLength/maxHitLength+1;

  TVector3 start=pos;
  TVector3 end=start;
  TVector3 step=dir*(maxHitLength/(dir.Mag()));
  float stepLength=0;
  for (int i=0;i<nHitSteps;i++){
    start=end;//new starting point is the previous ending point.
    if (i+1==nHitSteps){
      //last step is the remainder size
      end=strike;
      delta=start-end;
      stepLength=delta.Mag();
    } else {
      //all other steps are uniform length
      end=start+step;
      stepLength=step.Mag();
    }
    
  //now compute the laser hit:
  //printf("PHG4TpcDirectLaser::Adding New Hit(%1.2f,%1.2f,%d): (%1.2f,%1.2f,%1.2f) to (%1.2f,%1.2f,%1.2f)\n",theta,phi,laser,start.X()/cm,start.Y()/cm,start.Z()/cm,end.X()/cm,end.Y()/cm,end.Z()/cm);
      if (i+1==nHitSteps){
	//printf("PHG4TpcDirectLaser::Strike %c\n",strikeChar);
      }

  //from phg4tpcsteppingaction.cc
  hit = new PHG4Hitv1();
  hit->set_layer(99); // dummy number
  //here we set the entrance values in cm
  hit->set_x(0, start.X() / cm);
  hit->set_y(0, start.Y() / cm);
  hit->set_z(0, start.Z() / cm);
  //and the exist values
  hit->set_x(1, end.X() / cm);
  hit->set_y(1, end.Y() / cm);
  hit->set_z(1, end.Z() / cm);


  // momentum
  hit->set_px(0, 700.0); // GeV
  hit->set_py(0, 700.0);
  hit->set_pz(0, 700.0);
  
  // time in ns
  hit->set_t(0, 0.0); // nanosecond
  //set and save the track ID
  hit->set_trkid(-1); // dummy number

  
  hit->set_px(1, 700.0); // dummy large number, in GeV
  hit->set_py(1, 700.0);
  hit->set_pz(1, 700.0);
  
  hit->set_t(1, 0.0); // dummy number, nanosecond

  //calculate the total energy deposited

  //should calc this stuff in advance!
  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm

  //double Ne_NTotal = 43;    // Number/cm
  //double CF4_NTotal = 100;  // Number/cm
  //double Tpc_NTot = 0.90 * Ne_NTotal + 0.10 * CF4_NTotal;

  double Tpc_NTot = nElectrons;
  double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;

  //double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  //double Tpc_ElectronsPerGeV = Tpc_NTot / Tpc_dEdx*1e6; //electrons per gev.

  double electrons_per_gev = Tpc_dEdx*1e6 / Tpc_NTot; // GeV dep per electron

  float electrons_per_length=300./cm;//hardcoded
  float totalE=electrons_per_length*stepLength/electrons_per_gev;//rcc dummy hardcoded 300 electrons per cm!

  hit->set_eion(totalE);
  hit->set_edep(totalE);

  PHG4Hits.push_back(hit);
  }
  
  return;
}
