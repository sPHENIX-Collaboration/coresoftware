#include "BEmcRecFEMC.h"
#include "BEmcCluster.h"

#include <cmath>

void BEmcRecFEMC::Tower2Global(float E, float xC, float yC,
				 float& xA, float& yA, float& zA ) 
// With shower depth correction reflected in zA
{
  const float DZ = 18; // in cm, tower half length
  const float D = 13.3; // in cm, shower depth at 1 GeV relative to tower face; obtained from GEANT
  const float X0 = 2.2; // in cm; obtained from GEANT (should be ~ rad length)

  xA = yA = zA = 0;

  int ix = xC+0.5; // tower #
  if( ix<0 || ix >= fNx ) {
    printf("Error in BEmcRecFEMC::SectorToGlobal: wrong input x: %d\n",ix);
    return;
  }

  int iy = yC+0.5; // tower #
  if( iy<0 || iy >= fNy ) {
    printf("Error in BEmcRecFEMC::SectorToGlobal: wrong input y: %d\n",iy);
    return;
  }

  // Next tower in x
  TowerGeom geomx;
  int idx = 0;
  if( ix<fNx/2 ) {
    idx += 1;
    while(!GetTowerGeometry(ix+idx,iy,geomx) && idx< fNx/2) idx += 1;
  }
  else {
    idx -= 1;
    while(!GetTowerGeometry(ix+idx,iy,geomx) && idx>-fNx/2) idx -= 1;
  }
  if( idx>=fNx/2 || idx<=-fNx/2 ) {
    printf("Error in BEmcRecFEMC::Tower2Global: Error in geometry extraction for x= %f (ix,iy)=(%d,%d)\n",xC,ix,iy);
    return;
  }

  // Next tower in y
  TowerGeom geomy;
  int idy = 0;
  if( iy<fNy/2 ) {
    idy += 1;
    while(!GetTowerGeometry(ix,iy+idy,geomy) && idy< fNy/2) idy += 1;
  }
  else {
    idy -= 1;
    while(!GetTowerGeometry(ix,iy+idy,geomy) && idy>-fNy/2) idy -= 1;
  }
  if( idy>=fNy/2 || idy<=-fNy/2 ) {
    printf("Error in BEmcRecFEMC::Tower2Global: Error in geometry extraction for y= %f (ix,iy)=(%d,%d)\n",yC,ix,iy);
    return;
  }

  float dx, dy;
  float Zcenter;

  TowerGeom geom0;
  if( GetTowerGeometry(ix,iy,geom0) ) { // Found; this is for normal case 
    dx = (geomx.Xcenter - geom0.Xcenter)/float(idx);
    dy = (geomy.Ycenter - geom0.Ycenter)/float(idy);
    xA = geom0.Xcenter + (xC-ix)*dx;
    yA = geom0.Ycenter + (yC-iy)*dy;
    Zcenter = geom0.Zcenter;
  }
  else { // Not found; weird case (cluster center of gravity outside the EMCal)
    dx = (geomx.Xcenter - geomy.Xcenter)/float(idx);
    dy = (geomy.Ycenter - geomx.Ycenter)/float(idy);
    xA = geomx.Xcenter + (xC-ix)*dx - dx*idx;
    yA = geomx.Ycenter + (yC-iy)*dy; // Here I take Ycenter from geomx!
    Zcenter = geomx.Zcenter;
    //    printf("BEmcRecFEMC::Tower2Global: CG outside EMCal: input=(%f,%f), tower ind=(%d,%d); dd=(%d,%d); global=(%f,%f)\n",xC,yC,ix,iy,idx,idy,xA,yA);
  }

  float logE = log(0.1);
  if( E>0.1 ) logE = log(E);
  float ZcenterV = Zcenter - fVz;
  float cosT = fabs(ZcenterV)/sqrt(xA*xA+yA*yA+ZcenterV*ZcenterV);
  zA = (Zcenter-DZ) + (D+X0*logE)*cosT;

  //  zA = Zcenter; //!!!!!
}

void BEmcRecFEMC::CorrectEnergy(float Energy, float x, float y, 
				   float* Ecorr)
{
  *Ecorr = Energy;
}

float BEmcRecFEMC::GetImpactAngle(float e, float x, float y)
// Get impact angle, (x,y) - position in Sector frame (cm)
{
  /*
  float xVert, yVert, zVert;
  float vx, vy, vz;

  GlobalToSector( fVx, fVy, fVz, &xVert, &yVert, &zVert );
  vz = -zVert;
  vy = y - yVert;
  vx = x - xVert;
  // From this point X coord in sector frame is Z coord in Global Coord System !!!
  *sinT = sqrt((vx*vx+vy*vy)/(vx*vx+vy*vy+vz*vz));
  */
  return 0;
}

void BEmcRecFEMC::CorrectECore(float Ecore, float x, float y, float* Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers, 
  // long energy leakage and angle dependance
  //
  // (x,y) - shower CG in tower units (not projected anywhere!)
  
  //  *Ecorr = Ecore;
  float ec, ec2, corr;
  const float par1 = 0.938;
  const float par2 = 0.50;
  const float par3 = 0.067;

  *Ecorr = Ecore;
  if( Ecore < 0.01 ) return;

  float xA, yA, zA;
  Tower2Global(Ecore/0.91, x, y, xA, yA, zA); // 0.91 - average correction
  float tanT = sqrt(xA*xA+yA*yA)/fabs(zA-fVz);
  corr = par1 * ( 1 - tanT*tanT*tanT*(par2 + par3*log(Ecore)) );
  ec = Ecore/corr;

  //  CorrectEnergy( ec, x, y, &ec2);
  ec2 = ec; // !!!!! CorrectEnergy must be implemented
  *Ecorr = ec2;
}

void BEmcRecFEMC::CorrectPosition(float Energy, float x, float y,
				     float* pxc, float* pyc )
{
  // Corrects the Shower Center of Gravity for the systematic error due to 
  // the limited tower size
  //
  // Everything here is in tower units. 
  // (x,y) - CG position, (*pxc,*pyc) - corrected position

  float xZero, yZero, bx, by;
  float t, x0, y0;
  int ix0, iy0;

  *pxc = x;
  *pyc = y;
  //  return;

  if( Energy<0.01 ) return;

  float xA, yA, zA;
  Tower2Global(Energy, x, y, xA, yA, zA);
  zA -= fVz;
  float sinTx = xA/sqrt(xA*xA+zA*zA);
  float sinTy = yA/sqrt(yA*yA+zA*zA);
  float sin2Tx = sinTx*sinTx;
  float sin2Ty = sinTy*sinTy;

  if( sinTx > 0 ) xZero= -0.417 *sinTx - 1.500 * sin2Tx ;
  else            xZero= -0.417 *sinTx + 1.500 * sin2Tx ;

  if( sinTy > 0 ) yZero = -0.417 * sinTy - 1.500 * sin2Ty;
  else            yZero = -0.417 * sinTy + 1.500 * sin2Ty;

  t = 0.98 + 0.98*sqrt(Energy);
  t *= 0.7; // Temp correction
  bx = 0.17-0.009*log(Energy) + t*sin2Tx;
  by = 0.17-0.009*log(Energy) + t*sin2Ty;

  //  xZero *= 0.5;
  //  yZero *= 0.5;
  
  x0 = x + xZero;
  ix0 = EmcCluster::lowint(x0 + 0.5);

  if( EmcCluster::ABS(x0-ix0) <= 0.5 ) {
    x0 = (ix0-xZero)+bx*asinh( 2.*(x0-ix0)*sinh(0.5/bx) );
    *pxc = x0;
  }
  else {
    *pxc =  x;
    printf("????? Something wrong in BEmcRecFEMC::CorrectPosition: x=%f  dx=%f\n", x, x0-ix0);
  }

  y0 = y + yZero;
  iy0 = EmcCluster::lowint(y0 + 0.5);

  if( EmcCluster::ABS(y0-iy0) <= 0.5 ) {
    y0 = (iy0-yZero)+by*asinh( 2.*(y0-iy0)*sinh(0.5/by) );
    *pyc = y0;
  }
  else {
    *pyc = y;
    printf("????? Something wrong in BEmcRecFEMC::CorrectPosition: y=%f  dy=%f\n", y, y0-iy0);
  }
  
}
