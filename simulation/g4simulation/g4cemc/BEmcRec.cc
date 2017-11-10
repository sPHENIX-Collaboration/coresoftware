// Name: BEmcRec.h
// Author: A. Bazilevsky, Apr 2012
// Modified from EmcSectorRec.cxx and EmcScSectorRec.cxx

// BEmcRec -- base class for barrel EMcal for sPHENIX

// ///////////////////////////////////////////////////////////////////////////

#include "BEmcRec.h"
#include "BEmcCluster.h"
#include <TMath.h>

#include <cmath>

using namespace std;

// Define and initialize static members

// Minimal shower energy when splitting peakarea onto showers, used in Gamma()
float const BEmcRec::fgMinShowerEnergy=0.1;

// Max number of clusters in sector, used in Find_Clusters()
int const BEmcRec::fgMaxLen=10000;

// Default level, now for the Conf Level: 1% for GEANT, 2%-5% for TestBeam

float BEmcRec::fgChi2Level[50]={
    6.634899, 4.605171, 3.780564, 3.318915, 3.017103,    
    2.801872, 2.639259, 2.511249, 2.407341, 2.320967,    
    2.247720, 2.184744, 2.129863, 2.081515, 2.038526,    
    1.999994, 1.965214, 1.933627, 1.904781, 1.878311,    
    1.853912, 1.831334, 1.810365, 1.790825, 1.772564, 
    1.755449, 1.739367, 1.724222, 1.709926, 1.696406,    
    1.683593, 1.671430, 1.659864, 1.648850, 1.638344,    
    1.628311, 1.618716, 1.609528, 1.600721, 1.592268,    
    1.584148, 1.576338, 1.568822, 1.561579, 1.554596,    
    1.547856, 1.541346, 1.535055, 1.528968, 1.523077 };   

// For the Conf Level: 1% for GEANT, 2%-5% for TestBeam
float BEmcRec::fgChi2Level1[50]={
    6.634899, 4.605171, 3.780564, 3.318915, 3.017103,    
    2.801872, 2.639259, 2.511249, 2.407341, 2.320967,    
    2.247720, 2.184744, 2.129863, 2.081515, 2.038526,    
    1.999994, 1.965214, 1.933627, 1.904781, 1.878311,    
    1.853912, 1.831334, 1.810365, 1.790825, 1.772564, 
    1.755449, 1.739367, 1.724222, 1.709926, 1.696406,    
    1.683593, 1.671430, 1.659864, 1.648850, 1.638344,    
    1.628311, 1.618716, 1.609528, 1.600721, 1.592268,    
    1.584148, 1.576338, 1.568822, 1.561579, 1.554596,    
    1.547856, 1.541346, 1.535055, 1.528968, 1.523077 };   

// For the Conf Level: 2% for GEANT, 4%-7% for TestBeam
float BEmcRec::fgChi2Level2[50]={
    5.411895, 3.912024, 3.278443, 2.916812, 2.677547,    
    2.505458, 2.374582, 2.271008, 2.186567, 2.116065,    
    2.056169, 2.004491, 1.959343, 1.919481, 1.883964,    
    1.852072, 1.823237, 1.797008, 1.773021, 1.750981,    
    1.730640, 1.711795, 1.694274, 1.677931, 1.662643,    
    1.648301, 1.634814, 1.622101, 1.610093, 1.598727,    
    1.587948, 1.577709, 1.567968, 1.558684, 1.549824,    
    1.541357, 1.533256, 1.525494, 1.518051, 1.510903,    
    1.504033, 1.497424, 1.491059, 1.484924, 1.479006,    
    1.473292, 1.467771, 1.462433, 1.457267, 1.452265 };

// ///////////////////////////////////////////////////////////////////////////
// BEmcRec member functions

BEmcRec::BEmcRec()
{
  fModules=new vector<EmcModule>;
  fClusters=new vector<EmcCluster>;
  SetPeakThreshold(0.08);
  SetChi2Limit(2);
}

// ///////////////////////////////////////////////////////////////////////////

BEmcRec::~BEmcRec()
{

  if(fModules){

    fModules->clear();
    delete fModules;

  }

  if(fClusters){

    fClusters->clear();
    delete fClusters;

  }
}

// ///////////////////////////////////////////////////////////////////////////

void  BEmcRec::SetGeometry( int nx, int ny, float txsz, float tysz )
{
  fNx = nx;
  fNy = ny;
  fModSizex = txsz;
  fModSizey = tysz;
}

/*
void  BEmcRec::SetGeometry(SecGeom const &geom, PHMatrix * rm, PHVector * tr )
{
  // Sets the sector geometry.
  // Should be called before Find_Clusters(...)
  // Note: static data members are initialized => has to be called after
  // a new BEmcRec has been created.

  emcrm = *rm; 
  emctr = *tr; 
  PHFrame local;
  PHFrame global=PHGeometry::MatrixAndVector2frames(local,emcrm,emctr);
  PHGeometry::frames2MatrixAndVector(global,local,invemcrm,invemctr);

  // The number of towers in X and Y dir.
  fNx = geom.nx;
  fNy = geom.ny;

  // Tower size in X and Y dir.
  fModSizex = geom.Tower_xSize;
  fModSizey = geom.Tower_ySize;

}
*/

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetModules(vector<EmcModule> const *modules)
{

  *fModules=*modules;

#if 0
  // Make sure that each hit(fired module) knows what sector it belogs to
  vector<EmcModule>::iterator listIter;
  for(listIter=fModules->begin(); listIter!=fModules->end(); listIter++)
    listIter->fOwner=this;
#endif // #if 0

}

// ///////////////////////////////////////////////////////////////////////////

int BEmcRec::FindClusters()
{
  // Cluster search algorithm based on Lednev's one developed for GAMS.
  // Returns -1 if fgMaxLen parameter is too small (just increase it)
  
  int nhit, nCl;
  int LenCl[fgMaxLen];
  int next, ib, ie, iab, iae, last, LastCl, leng, ich;
  int ia = 0;
  
  EmcModule *vv;
  EmcModule *vhit, *vt;
  EmcCluster Clt(this);
  vector<EmcModule>::iterator ph;
  vector<EmcModule> hl;
  
  (*fClusters).erase(  (*fClusters).begin(),  (*fClusters).end() );
  nhit = (*fModules).size();
  
  if( nhit <= 0 ) return 0;
  if( nhit == 1 ) {
    Clt.ReInitialize( (*fModules) );
    fClusters->push_back(Clt);
    return 1;
  }

  vt = new EmcModule[nhit];
  vhit = new EmcModule[nhit];
  
  ph = (*fModules).begin();
  vv = vhit;
  while( ph != (*fModules).end() ) *vv++ = *ph++;
  
  qsort( vhit, nhit, sizeof(EmcModule), HitNCompare );
  
  nCl=0;
  next = 0;
  for( ich=1; ich<nhit+1; ich++ ){
    if( ich<nhit ) ia=vhit[ich].ich;

    // New subcluster
    //
    if( (ia-vhit[ich-1].ich > 1)  // the beginning of new subcluster
	|| (ich >= nhit)          // just finish defining last sub-cluster
	||(ia-ia/fNx*fNx == 0) ){ // new raw -> new subcluster

      ib=next;
      ie=ich-1;
      next=ich;
      if( nCl >= fgMaxLen ) {
	return -1;
      }
      nCl++;
      LenCl[nCl-1]=next-ib;
      if( nCl > 1 ) {

	// Job to glue the subclusters with common edge
	//
	iab=vhit[ib].ich;       // The last subcl begin
	iae=vhit[ie].ich;       // The last subcl end
	last=ib-1;              // The prelast subcl end
	LastCl=nCl-2;
	for( int iCl=LastCl; iCl>=0; iCl-- ){
	  leng=LenCl[iCl];

	  if( iab-vhit[last].ich > fNx ) goto new_ich;
	  for( int ichc=last; ichc >= last-leng+1; ichc-- ) {

	    //	    if( iab-vhit[ichc].ich >  fNx ) goto new_icl;

	    //	    if( iae-vhit[ichc].ich >= fNx
	    if( (vhit[ichc].ich+fNx <= iae && vhit[ichc].ich+fNx >= iab)
		|| ((iae%fNx == fNx-1) && (iae-vhit[ichc].ich == fNx-1))  // Only for CYLinder geom !!!!
		) {

	      // Swap iCl-cluster towers (of length "leng") with whatever was between it and the last subcluster (of length "ib-1-last") - to make it adjecent to the last subcluster
	      CopyVector( &vhit[last+1-leng], vt, leng );
	      CopyVector( &vhit[last+1], &vhit[last+1-leng], ib-1-last );
	      CopyVector( vt, &vhit[ib-leng], leng );

	      // Now the number of clusters is reduced by 1 and the length of the last one increased by iCl-cluster length "leng"
	      for( int i=iCl; i<nCl-2; i++ ) LenCl[i]=LenCl[i+1];
	      ib -= leng;
	      LenCl[nCl-2]=LenCl[nCl-1]+leng;
	      nCl--;
	      goto new_icl;

	    }
	  } // for( int ichc=last

	new_icl:           last=last-leng;
	} // for( int iCl=LastCl

      } // if( nCl > 1

    } // if( (ia-vhit

  new_ich:   continue;
  } //  for( ich=1

  if( nCl > 0 ) {
    ib=0;
    for( int iCl=0; iCl<nCl; iCl++ ) { 
      leng=LenCl[iCl];
      hl.erase( hl.begin(), hl.end() );
      for( ich=0; ich<leng; ich++ ) hl.push_back(vhit[ib+ich]);
      Clt.ReInitialize(hl);
      ib += LenCl[iCl];
      fClusters->push_back(Clt);
    }
  }
  delete [] vhit;
  delete [] vt;

  return nCl;
  
}

// ///////////////////////////////////////////////////////////////////////////
void BEmcRec::GetImpactAngle(float x, float y, float *sinT )
  // Get impact angle, (x,y) - position in Sector frame (cm)
{
  float xVert, yVert, zVert;
  float vx, vy, vz;

  GlobalToSector( fVx, fVy, fVz, &xVert, &yVert, &zVert );
  vz = -zVert;
  vy = y - yVert;
  vx = x - xVert;
  // From this point X coord in sector frame is Z coord in Global Coord System !!!
  *sinT = sqrt((vx*vx+vy*vy)/(vx*vx+vy*vy+vz*vz));
}


void BEmcRec::GlobalToSector(float xgl, float ygl, float zgl, float* px,
				 float* py, float* pz)
{

  *px = xgl + fModSizex*(fNx-1)/2.;
  *py = ygl + fModSizey*(fNy-1)/2.;
  *pz = zgl;
  /*
  PHPoint phnxHit(xgl, ygl, zgl);
  PHPoint emcHit  = PHGeometry::transformPoint(invemcrm, invemctr, phnxHit);
  *px =  emcHit.getX();
  *py =  emcHit.getY();
  *pz =  emcHit.getZ();
  */
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SectorToGlobal(float xsec, float ysec, float zsec,
				 float* px, float* py, float* pz ) 
{
  *px = xsec - fModSizex*(fNx-1)/2.;
  *py = ysec - fModSizey*(fNy-1)/2.;
  *pz = 0.;
  /*
  PHPoint emcHit(xsec, ysec, zsec);
  PHPoint phnxHit  = PHGeometry::transformPoint(emcrm, emctr, emcHit);
  *px =  phnxHit.getX();
  *py =  phnxHit.getY();
  *pz =  phnxHit.getZ();
  */  
}


// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SectorToGlobalErr( float dxsec, float dysec, float dzsec,
				     float* pdx, float* pdy, float* pdz ) 
{
  *pdx = 0.;
  *pdy = 0.;
  *pdz = 0.;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::Gamma(int nh, EmcModule* phit, float* pchi, float* pchi0,
			float* pe1, float* px1, float* py1, float* pe2,
			float* px2, float* py2, int &ndf)
{
  // Tests for 1 or 2 photon hypothesis by minimizing the chi2.
  // If the energy of one of two showers is less then fgMinShowerEnergy 
  // they are merged
  //
  // Returns two shower parameters (pe1,px1,py1,pe2,px2,py2). 
  // If one shower hypothesis accepted -> *pe2=0
  // *pchi  - chi2 after splitting
  // *pchi0 - chi2 before splitting
  // ndf contains number of degrees of freedom (not equal to number of cluster
  // modules for PbGl).

  float e1, x1, y1, e2, x2, y2;
  float chi, chi0, chi00, chisq0, chisave;
  float chir, chil, chiu, chid;
  int dof;
  float x0, y0, d2, xm2;
  float stepx, stepy, parx, pary;
  const float dxy=0.06;
  const float stepmin=0.01;
  const float zTG=100;
  const float xmcut=0.0015; // (GeV), for overlapping showers separation

  *pe1=0;
  *px1=0;
  *py1=0;
  *pe2=0;
  *px2=0;
  *py2=0;
  if( nh <= 0 ) return;
  Mom1(nh,phit,&e1,&x1,&y1);
  *pe1=e1;     
  if( e1 <= 0 ) return;
  
  SetProfileParameters(0, e1,x1,y1);

  chisave = *pchi;
  chi = *pchi;
  // ClusterChisq parameter list changed MV 28.01.00
  chi0 = ClusterChisq(nh, phit, e1, x1, y1, ndf);

  chisq0 = chi0;
  dof = ndf; // nh->ndf MV 28.01.00

  // ndf=0 means the cluster's chi2 cannot be found; in this case chi0=0.
  if( dof < 1 ) dof=1;
  chi = chisq0/dof;
  x0 = x1;
  y0 = y1;
  for(;;){

    chir = ClusterChisq(nh, phit, e1, x0+dxy, y0, ndf);
    chil = ClusterChisq(nh, phit, e1, x0-dxy, y0, ndf);
    chiu = ClusterChisq(nh, phit, e1, x0, y0+dxy, ndf);
    chid = ClusterChisq(nh, phit, e1, x0, y0-dxy, ndf);
    
    if( (chi0 > chir) || (chi0 > chil) ) {
      stepx = dxy;
      if( chir > chil ) stepx = -stepx;
    }
    else {
      stepx = 0;
      parx = chir+chil-2*chi0;
      if( parx > 0 ) stepx = -dxy*(chir-chil)/2/parx;
    }
    
    if( (chi0 > chiu) || (chi0 > chid) ) {
      stepy = dxy;
      if( chiu > chid ) stepy = -stepy;
    }
    else {
      stepy = 0;
      pary = chiu+chid-2*chi0;
      if( pary > 0 ) stepy = -dxy*(chiu-chid)/2/pary;
    }
    if( (EmcCluster::ABS(stepx) < stepmin) && (EmcCluster::ABS(stepy) < stepmin) ) break;
    chi00 = ClusterChisq(nh, phit, e1, x0+stepx, y0+stepy, ndf);

    if( chi00 >= chi0 ) break;
    chi0 = chi00;
    x0 += stepx;
    y0 += stepy;
  }
  if( chi0 < chisq0 ) {
    x1 = x0;
    y1 = y0;
    chi = chi0/dof;
  }
  
  *pchi0 = chi;
  *pchi = chi;
  *px1 = x1;
  *py1 = y1;

  if( e1 <= fgMinShowerEnergy ) return;
  
  if( chi > chisave ) {
    TwoGamma(nh,phit,&chi,&e1,&x1,&y1,&e2,&x2,&y2);
    if( e2 > 0 ) {
      d2 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))/zTG/zTG;
      xm2 = e1*e2*d2;
      if( xm2 > 0 ) xm2 = sqrt(xm2);
      if( xm2 > xmcut && e1 > fgMinShowerEnergy && e2 > fgMinShowerEnergy) {
	*pe1 = e1;
	*px1 = x1;
	*py1 = y1;
	*pe2 = e2;
	*px2 = x2;
	*py2 = y2;
	*pchi = chi;
      }
    }	
  }
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::Mom1(int nh, EmcModule* phit, float* pe, float* px,
		       float* py)
{
  // First momentum calculation

  float a, xw, yw, e;
  int ix, iy;
  EmcModule* p;
  
  *pe=0;
  *px=0;
  *py=0;
  if( nh <= 0 ) return;
  p=phit;
  xw=0;
  yw=0;
  e=0;
  for( int i=0; i<nh; i++ ) {
    a = p->amp;
    iy = p->ich / fNx;
    ix = p->ich - iy*fNx;
    e += a;
    xw += ix*a;
    yw += iy*a;
    p++;
  }
  *pe = e;
  if( e <= 0 ) return;
  *px = xw/e;
  *py = yw/e;

}

// ///////////////////////////////////////////////////////////////////////////

EmcModule BEmcRec::ShiftX(int ish, EmcModule ehit)
{
  EmcModule hh = ehit;
  int iy = hh.ich / fNx;
  int ix = hh.ich % fNx + ish;
  while(ix<0   ) ix+=fNx;
  while(ix>=fNx) ix-=fNx;
  int ich = iy*fNx + ix;
  hh.ich = ich;
  return hh;
}

int BEmcRec::ShiftX(int ishift, int nh, EmcModule* phit0, EmcModule* phit1)
// Used to avoid edge effect in numbering in X-direction (from fNx-1 to 0)
// (for CYL geometry)
//
// Shift cluster by "ishift" towers in X
// If ishift=0, let it decide how much to shift
{
  int ishift_def = fNx/2;
  int ix, iy, ich, ish;

  bool bshift = false;
  if( ishift==0 ) {
    for( int i=0; i<nh; i++ ) {
      phit1[i] = phit0[i];
      ix = phit0[i].ich % fNx;
      if( ix==0 ) bshift = true;
    }
    if( !bshift ) return 0;
  }

  if( ishift != 0 ) ish = ishift;
  else              ish = ishift_def;

  int ixmin =  999;
  int ixmax = -999;
  for( int i=0; i<nh; i++ ) {
    phit1[i] = phit0[i];
    iy = phit0[i].ich / fNx;
    ix = phit0[i].ich % fNx + ish;

    while(ix<0   ) ix+=fNx;
    while(ix>=fNx) ix-=fNx;

    if( ixmin>ix ) ixmin = ix;
    if( ixmax<ix ) ixmax = ix;

    ich = iy*fNx + ix;
    phit1[i].ich = ich;
  }

  if( ishift==0 && ixmax-ixmin > fNx/2 ) printf("!!! Warning: Too long cluster (%d towers): reconstruction may be wrong !!!\n",ixmax-ixmin+1);
  return ish;
}


void BEmcRec::Momenta(int nh, EmcModule* phit, float* pe, float* px,
			  float* py, float* pxx, float* pyy, float* pyx )
{
  // First and second momenta calculation
  
  float a, x, y, e, xx, yy, yx;
  int ix, iy, i;
  EmcModule* p;
  
  *pe=0;
  *px=0;
  *py=0;
  *pxx=0;
  *pyy=0;
  *pyx=0;
  if( nh <= 0 ) return;
  
  //  p=phit;
  EmcModule* phit1 = new EmcModule[nh];
  int ish = ShiftX(0, nh, phit, phit1);
  p = phit1;

  x=0;
  y=0;
  e=0;
  xx=0;
  yy=0;
  yx=0;
  for( i=0; i<nh; i++ ) {
    a = p->amp;
    iy = p->ich / fNx;
    ix = p->ich - iy*fNx;
    e += a;
    x += ix*a;
    y += iy*a;
    xx += a*ix*ix;
    yy += a*iy*iy;
    yx += a*ix*iy;
    p++;
  }
  *pe = e;

  if( e>0 ) {
    x /= e;
    y /= e;
    xx = xx/e - x*x;
    yy = yy/e - y*y;
    yx = yx/e - y*x;
    
    x -= float(ish);
    while(x<-0.5    ) x += float(fNx);
    while(x>=fNx-0.5) x -= float(fNx);
    
    *px = x;
    *py = y;
    *pxx = xx;
    *pyy = yy;
    *pyx = yx;
  }

  delete [] phit1;
}

// ///////////////////////////////////////////////////////////////////////////
// Static functions

int BEmcRec::HitNCompare(const void* h1, const void* h2)
{

  return ( ((EmcModule*)h1)->ich - ((EmcModule*)h2)->ich );

}

// ///////////////////////////////////////////////////////////////////////////

int BEmcRec::HitACompare(const void* h1, const void* h2)
{

  float amp1 = ((EmcModule*)h1)->amp;
  float amp2 = ((EmcModule*)h2)->amp;
  return (amp1<amp2) ? 1: (amp1>amp2) ? -1 : 0;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::ZeroVector(int* v, int N)
{

  int* p = v;
  for(int i=0; i<N; i++) *p++ = 0;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::ZeroVector(float* v, int N)
{

  float* p = v;
  for(int i=0; i<N; i++) *p++ = 0;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::ZeroVector(EmcModule* v, int N)
{
  //  memset(v, 0, N*sizeof(EmcModule));
  for(int i=0; i<N; i++){ v[i].ich=0; v[i].amp=0; v[i].tof=0; }
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::ResizeVector(int* Vector, int OldSize, int NewSize)
{

  int* vsave;

  if( OldSize <= 0 ) { Vector = new int[NewSize]; return; }
  vsave = new int[OldSize];
  CopyVector( Vector, vsave, OldSize );
  delete [] Vector;
  Vector = new int[NewSize];
  if( NewSize > OldSize ) CopyVector( vsave, Vector, OldSize );
  else 			CopyVector( vsave, Vector, NewSize );
  delete [] vsave;
  return;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CopyVector(int* from, int* to, int N)
{

  if( N <= 0 ) return;
  for( int i=0; i<N; i++ ) to[i]=from[i];

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CopyVector(EmcModule* from, EmcModule* to, int N)
{

  if( N <= 0 ) return;
  for( int i=0; i<N; i++ ) to[i]=from[i];

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::c3to5(float e0, float x0, float y0, float eps,
		      float dx, float dy,
		      float* pe1, float* px1, float* py1,
		      float* pe2, float* px2, float* py2)
{
  // 3 to 5-dim space conversion subroutine.
  // eps=(e1-e2)/e0,  (e0=e1+e2), x0*e0=x1*e1+x2*e2, dx=x1-x2

  *pe1 = e0*(1+eps)/2;
  *pe2 = e0 - *pe1;
  *px1 = x0 + dx*(1-eps)/2;
  *py1 = y0 + dy*(1-eps)/2;
  *px2 = x0 - dx*(1+eps)/2;
  *py2 = y0 - dy*(1+eps)/2;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetChi2Limit(int limit)
{
  // Sets the limit for PeakArea splitting onto 2 EMShowers:
  // limit=0 -> For the Conf Level: 0%
  // limit=1 -> For the Conf Level: 1% for GEANT, 2%-5% for TestBeam
  // limit=2 -> For the Conf Level: 2% for GEANT, 4%-7% for TestBeam
  
  int i;
  
  switch ( limit ) {

  case 0:
    for( i=0; i<50; i++ ) fgChi2Level[i]=9999.;
    break;
  case 1:
    for( i=0; i<50; i++ ) fgChi2Level[i]=fgChi2Level1[i];
    break;
  case 2:
    for( i=0; i<50; i++ ) fgChi2Level[i]=fgChi2Level2[i];
    break;
  default:
    for( i=0; i<50; i++ ) fgChi2Level[i]=fgChi2Level1[i];
    break;

  }
}



/////////////////////////////////////////////////////////////////////
// From EmcScSectorRec
/////////////////////////////////////////////////////////////////////



// Define and initialize static members

// Parameters for sigma in Hi2 calculations (p.36-37 v.3)
float BEmcRec::fgEpar00 = 0.005;
float BEmcRec::fgEpar0 = 0.0014;
float BEmcRec::fgEpar1 = 0.03;
float BEmcRec::fgEpar2 = -0.03;
// This is for PPLO mode !!!
float BEmcRec::fgEpar3 = 0.;
float BEmcRec::fgEpar4 = 4.0;

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CorrectEnergy(float Energy, float x, float y, 
				   float* Ecorr)
{
  // Corrects the EM Shower Energy for attenuation in fibers and 
  // long energy leakage
  //
  // (x,y) - impact position (cm) in Sector frame

  float sinT;
  float att, leak, corr;
  const float leakPar = 0.0033; // parameter from fit
  const float attPar = 120; // Attenuation in module (cm)
  const float X0 = 2; // radiation length (cm)

  *Ecorr = Energy;
  if( Energy < 0.01 ) return;

  GetImpactAngle(x, y, &sinT); // sinT not used so far
  leak = 2-sqrt(1+leakPar*log(1+Energy)*log(1+Energy));
  att = exp(log(Energy)*X0/attPar);
  corr = leak*att;
  *Ecorr = Energy/corr;
}

/////////////////////////////////////////////////////////////////

void BEmcRec::CorrectECore(float Ecore, float x, float y, float* Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers, 
  // long energy leakage and angle dependance
  //
  // (x,y) - impact position (cm) in Sector frame

  float ec, ec2, corr;
  float sinT;
  const float par1 = 0.918;
  const float par2 = 1.35;
  const float par3 = 0.003;

  *Ecorr = Ecore;
  if( Ecore < 0.01 ) return;

  GetImpactAngle(x, y, &sinT );
  corr = par1 * ( 1 - par2*sinT*sinT*sinT*sinT*(1 - par3*log(Ecore)) );
  ec = Ecore/corr;

  CorrectEnergy( ec, x, y, &ec2);
  *Ecorr = ec2;
}

/////////////////////////////////////////////////////////////////////

void BEmcRec::CorrectPosition(float Energy, float x, float y,
				     float* pxc, float* pyc, bool callSetPar)
{
  // Corrects the Shower Center of Gravity for the systematic error due to 
  // the limited tower size and angle shift
  //
  // Everything here is in cell units. 
  // (x,y) - CG position, (*pxc,*pyc) - corrected position

  float xShift, yShift, xZero, yZero, bx, by;
  float t, x0, y0;
  int ix0, iy0;
  int signx, signy;

  const float Xrad = 0.3; // !!!!! Need to put correct value
  const float Remc = 90.; // EMCal inner radius. !!!!! Should be obtained from geometry container

  SetProfileParameters( 0, Energy, x, y );
  if( fSinTx >= 0 ) signx =  1;
  else 	   signx = -1;
  if( fSinTy >= 0 ) signy =  1;
  else 	   signy = -1;
  t = 5.0+1.0*log(Energy); // In Rad Length units
  t *= ( Xrad/Remc/GetModSizex() ); // !!!!!
  xShift = t*fSinTx;
  yShift = t*fSinTy;
  xZero=xShift-(0.417*EmcCluster::ABS(fSinTx)+1.500*fSinTx*fSinTx)*signx;
  yZero=yShift-(0.417*EmcCluster::ABS(fSinTy)+1.500*fSinTy*fSinTy)*signy;
  xZero = xShift; // ...Somehow this works better !!!!!
  yZero = yShift; // ...Somehow this works better !!!!!
  t = 0.98 + 0.98*sqrt(Energy); // !!!!! Still from PHENIX
  bx = 0.15 + t*fSinTx*fSinTx;
  by = 0.15 + t*fSinTy*fSinTy;

  x0 = x;
  x0 = x0 - xShift + xZero;
  ix0 = EmcCluster::lowint(x0 + 0.5);
  if( EmcCluster::ABS(x0-ix0) <= 0.5 ) {
    x0 = (ix0-xZero)+bx*asinh( 2.*(x0-ix0)*sinh(0.5/bx) );
    *pxc = x0;
  }
  else {
    *pxc =  x - xShift;
    printf("????? Something wrong in CorrectPosition: x=%f  dx=%f\n", x, x0-ix0);
  }

  y0 = y;
  y0 = y0 - yShift + yZero;
  iy0 = EmcCluster::lowint(y0 + 0.5);
  if( EmcCluster::ABS(y0-iy0) <= 0.5 ) {
    y0 = (iy0-yZero)+by*asinh( 2.*(y0-iy0)*sinh(0.5/by) );
    *pyc = y0;
  }
  else {
    *pyc = y - yShift;
    printf("????? Something wrong in CorrectPosition: y=%f  dy=%f\n", y, y0-iy0);
  }
  
}


// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CalculateErrors( float e, float x, float y, float* pde,
				   float* pdx, float* pdy, float* pdz)
{
  // Returns the errors for the reconstructed energy and position 
  // (in the hypothesis of EM shower)
  // Should be called only just after CorrectPosition !!!

  float de, dy, dz, dxg, dyg, dzg;
  static float ae = 0.076, be = 0.022;  	// de/E = a/sqrt(E)&b
  static float a = 0.57, b = 0.155, d = 1.6;	// dx = a/sqrt(E)+b (cm)
  static float dx = 0.1;  // (cm)
  
  de = sqrt( ae*ae*e + be*be*e*e );
  dz = a/sqrt(e) + b;
  dy = dz;
  dz = sqrt( dz*dz + d*d*fSinTx*fSinTx );
  dy = sqrt( dy*dy + d*d*fSinTy*fSinTy );
  
  SectorToGlobalErr( dx, dy, dz, &dxg, &dyg, &dzg );
  
  *pde = de;
  *pdx = dxg;
  *pdy = dyg;
  *pdz = dzg;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetProfileParameters(int sec, float Energy, float x,
					float y )
{
  // Axis Z here means X in two dim Sector coord system !!! 
  // So one has to supply (x,y) coordinates in cell units (x instead of z)
  // If sec < 0 this routine changes only Energy dependent parameters - 
  // the angle dependent ones are set in the previous call

  static float sin2ax, sin2ay, sin2a, lgE;
  //  float vx, vy, vz;
  //  float xVert, yVert, zVert;
  int sign;

  if( sec >= 0 ) {
    // Non-ortogonality in phi in a sector (of 8 rows)
    int ix = x+0.5;
    float dx = x+0.5 - ix; // from 0 to 1
    int ix8 = ix%8;
    dx += ix8;
    fSinTx = (dx-4)*fModSizex;
    fSinTy = 0;
    sin2a = fSinTx*fSinTx + fSinTy*fSinTy;
    fSin4T = sin2a*sin2a;
    sin2ax = fSinTx*fSinTx;
    sin2ay = fSinTy*fSinTy;
  }
  
  if( Energy <= 1.e-10 ) lgE=0;
  else lgE=log(Energy);
  
  fPpar1=0.59-(1.45+0.13*lgE)*sin2a;
  fPpar2=0.265+(0.80+0.32*lgE)*sin2a;
  fPpar3=0.25+(0.45-0.036*lgE)*sin2a;
  fPpar4=0.42;
  
  if( fSinTx > 0 ) sign = 1;
  else sign = -1;
  fPshiftx = (1.05+0.12*lgE) * sin2ax * sign;
  fPshiftx = 0; // !!!!! Untill tuned ... may not be necessary
  
  if( fSinTy > 0 ) sign = 1;
  else sign = -1;
  fPshifty = (1.05+0.12*lgE) * sin2ay * sign;
  fPshifty = 0; // !!!!! Untill tuned ... may not be necessary
}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::PredictEnergy(float xc, float yc, float en)
{
  // Calculates the energy deposited in the tower, the distance between 
  // its center and shower Center of Gravity being (xc,yc)
  // en - shower energy
  // If en<0 -> no Shower Profile parameters change is needed

  float dx, dy, r1, r2, r3, e;
  
  if( en > 0 ) SetProfileParameters(-1,en,xc,yc);
  dx=fabs(xc-fPshiftx);
  dy=EmcCluster::ABS(yc-fPshifty);
  e=0;
  r2=dx*dx+dy*dy;
  r1=sqrt(r2);
  r3=r2*r1;
  e=fPpar1*exp(-r3/fPpar2)+fPpar3*exp(-r1/fPpar4);
  
  return e;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::TwoGamma(int nh, EmcModule* phit, float* pchi, float* pe1,
			   float* px1, float* py1, float* pe2, float* px2,
			   float* py2)
{

  float e0, x0, y0, xx, yy, yx;
  float dxy, rsg2, rsq;
  float dxc, dyc, r, epsc;
  int ix, iy, ixy, in, iter, dof;
  float step, cosi, chisq2, u;
  float e1c, x1c, y1c, e2c, x2c, y2c;
  float eps0 = 0.0;
  float eps1, eps2, chisqc, ex;
  float dx1, dy1, dx2, dy2, a0, d;
  float dchi, dchi0, dd, dchida, a1, a2;
  float gr = 0.0;
  float grec, grxc, gryc, grc, gx1, gx2, gy1, gy2;
  float gre = 0.0;
  float grx = 0.0;
  float gry = 0.0;
  float scal;
  float dx0 = 0.0;
  float dy0 = 0.0;
  
  const float epsmax=0.9999;
  const float stpmin=0.025;
  const float delch=2;
  
  Momenta(nh,phit,&e0,&x0,&y0,&xx,&yy,&yx);
  *pe2 = 0;
  *px2 = 0;
  *py2 = 0;
  if( nh <= 0 ) return;
  //  choosing of the starting point
  dxy = xx-yy;
  rsg2 = dxy*dxy + 4*yx*yx;
  if( rsg2 < 1e-20 ) rsg2 = 1e-20;
  rsq = sqrt(rsg2);
  dxc = -sqrt((rsq+dxy)*2);
  dyc =  sqrt((rsq-dxy)*2);
  if( yx >= 0 ) dyc = -dyc;
  r = sqrt(dxc*dxc + dyc*dyc);
  epsc = 0;
  for( in=0; in<nh; in++ ) {
    ixy = phit[in].ich;
    iy = ixy/fNx;
    ix = ixy - iy*fNx;
    u = (ix-x0)*dxc/r + (iy-y0)*dyc/r;
    epsc -= phit[in].amp * u * EmcCluster::ABS(u);
  }
  epsc /= (e0*rsq);
  if( epsc >  0.8 ) epsc = 0.8;
  if( epsc < -0.8 ) epsc =-0.8;
  dxc /= sqrt(1-epsc*epsc);
  dyc /= sqrt(1-epsc*epsc);
  //  Start of iterations
  step = 0.1;
  cosi = 0;
  chisq2 = 1.e35;
  for( iter=0; iter<100; iter++)
    {
      c3to5(e0,x0,y0,epsc,dxc,dyc,&e1c,&x1c,&y1c,&e2c,&x2c,&y2c);
      eps1 = (1+epsc)/2;
      eps2 = (1-epsc)/2;
      chisqc = 0;
      for( in=0; in<nh; in++ ) {
	ex = phit[in].amp;
	ixy = phit[in].ich;
	iy = ixy/fNx;
	ix = ixy - iy*fNx;
	dx1 = x1c - ix;
	dy1 = y1c - iy;
	dx2 = x2c - ix;
	dy2 = y2c - iy;
	a0 = e1c*PredictEnergy(dx1, dy1, e1c) + e2c*PredictEnergy(dx2, dy2, e2c);
	d = fgEpar00*fgEpar00 + e0*( fgEpar1*a0/e0 + fgEpar2*a0*a0/e0/e0 +fgEpar3*a0*a0*a0/e0/e0/e0 ) + e0*sqrt(e0)*fgEpar4*a0/e0*(1-a0/e0)*fSin4T + e0*e0*fgEpar0*fgEpar0;
	chisqc += (a0-ex)*(a0-ex)/d;
      }
      if( chisqc >= chisq2 ) {
	if( iter > 0 ) {
	  dchi = chisqc-chisq2;
	  dchi0 = gr*step;
	  step /= (2*sqrt(1+dchi/dchi0));
	}
	step /= 2;
      }
      else {
	// Calculation of gradient
	grec = 0;
	grxc = 0;
	gryc = 0;
	for( in=0; in<nh; in++ ) {
	  ex = phit[in].amp;
	  ixy = phit[in].ich;
	  iy = ixy/fNx;
	  ix = ixy - iy*fNx;
	  dx1 = x1c - ix;
	  dy1 = y1c - iy;
	  dx2 = x2c - ix;
	  dy2 = y2c - iy;
	  a1 = e1c*PredictEnergy(dx1,dy1,e1c);
	  a2 = e2c*PredictEnergy(dx2,dy2,e2c);
	  a0 = a1 + a2;
	  d = fgEpar00*fgEpar00 + e0*( fgEpar1*a0/e0 + fgEpar2*a0*a0/e0/e0 +fgEpar3*a0*a0*a0/e0/e0/e0 ) + e0*sqrt(e0)*fgEpar4*a0/e0*(1-a0/e0)*fSin4T + e0*e0*fgEpar0*fgEpar0;
	  dd = (a0-ex)/d;
	  dchida = dd*( 2 - dd*(fgEpar1 + 2*fgEpar2*a0/e0 + 3*fgEpar3*a0*a0/e0/e0 + e0*sqrt(e0)*fgEpar4*fSin4T*(1-2*a0/e0) + 2*fgEpar0*fgEpar0*a0) );
	  gx1 = ( e1c*PredictEnergy(x1c+0.05-ix,dy1,e1c) - a1 )*20;
	  gx2 = ( e2c*PredictEnergy(x2c+0.05-ix,dy2,e2c) - a2 )*20;
	  gy1 = ( e1c*PredictEnergy(dx1, y1c+0.05-iy,e1c) - a1 )*20;
	  gy2 = ( e2c*PredictEnergy(dx2, y2c+0.05-iy,e2c) - a2 )*20;
	  grec += (dchida*((a1/e1c-a2/e2c)*e0 - (gx1+gx2)*dxc -(gy1+gy2)*dyc)/2);
	  grxc += (dchida*(gx1*eps2-gx2*eps1));
	  gryc += (dchida*(gy1*eps2-gy2*eps1));
	}
	grc = sqrt(grec*grec + grxc*grxc + gryc*gryc);
	if( grc < 1e-10 ) grc = 1e-10;
	if( iter > 0 ) {
	  cosi = (gre*grec + grx*grxc + gry*gryc ) / (gr*grc);
	  scal = EmcCluster::ABS(gr/grc - cosi);
	  if( scal < 0.1 ) scal = 0.1;
	  step /= scal;
	}
	chisq2 = chisqc;
	eps0 = epsc;
	dx0 = dxc;
	dy0 = dyc;
	gre = grec;
	grx = grxc;
	gry = gryc;
	gr = grc;
      }
      epsc = eps0 - step*gre/gr;
      while( EmcCluster::ABS(epsc) >= epsmax ) {
	step /= 2;
	epsc = eps0 - step*gre/gr;
      }
      dxc = dx0 - step*grx/gr;
      dyc = dy0 - step*gry/gr;
      if( step*gr < stpmin ) break;
    }
  if( (*pchi)*nh-chisq2 < delch ) return;
  dof = nh;
  if( dof < 1 ) dof = 1;
  *pchi = chisq2/dof;
  c3to5(e0,x0,y0,eps0,dx0,dy0,pe1,px1,py1,pe2,px2,py2);

}

// ///////////////////////////////////////////////////////////////////////////

int BEmcRec::GetTowerID( int iy, int iz, int nn, int* iyy, int* izz, float* ee )
{
  while(iz<0   ) iz+=fNx;
  while(iz>=fNx) iz-=fNx;

  for( int i=0; i<nn; i++ ) {
    if( iy==iyy[i] && iz==izz[i]  ) return i;
  }
  return -1;
}

float BEmcRec::GetProb(vector<EmcModule> HitList, float &chi2, int &ndf)
{
  const float thresh = 0.01;
  const int DXY=3; // 2 is for 5x5 matrix; 3 for 7x7 matrix
  const int Nmax = 1000;
  float ee[Nmax];
  int iyy[Nmax];
  int izz[Nmax];

  int ich;
  vector<EmcModule>::iterator ph = HitList.begin();

  chi2 = 0;
  ndf = 0;

  int nn = 0;

  while( ph != HitList.end() ) {
    ee[nn] = ph->amp;
    if( ee[nn]>thresh ) {
      ich = ph->ich;
      izz[nn] = ich%fNx;
      iyy[nn] = ich/fNx;
      nn++;
      if( nn>=Nmax ) {
	printf("BEmcRec::GetProb: Cluster size is too big. Skipping the rest of the towers\n");
	break;
      }
    } // if( ee[nn]
    ++ph;
  } // while( ph 

  if( nn<=0 ) return -1;

  int iy0=-1, iz0=-1;
  float emax=0;

  for( int i=0; i<nn; i++ ) {
    if( ee[i]>emax ) {
      emax = ee[i];
      iy0 = iyy[i];
      iz0 = izz[i];
    }
  }

  if( emax<=0 ) return -1;

  int id;
  float etot=0;
  float sz=0;
  float sy=0;

  for( int idz=-DXY; idz<=DXY; idz++ ) {
    for( int idy=-DXY; idy<=DXY; idy++ ) {
      id = GetTowerID(iy0+idy, iz0+idz, nn, iyy, izz, ee);
      if( id>=0 ) {
	etot += ee[id];
	sz += ee[id]*(iz0+idz);
	sy += ee[id]*(iy0+idy);
      }
    }
  }
  float zcg = sz/etot; // Here cg allowed to be out of range
  float ycg = sy/etot;
  int iz0cg = int(zcg+0.5);
  int iy0cg = int(ycg+0.5);
  float ddz = fabs(zcg-iz0cg);
  float ddy = fabs(ycg-iy0cg);

  int isz = 1;
  if( zcg-iz0cg < 0 ) isz = -1;
  int isy = 1;
  if( ycg-iy0cg < 0 ) isy = -1;

  // 4 central towers: 43  
  //                   12
  // Tower 1 - central one  
  float e1, e2, e3, e4;
  e1 = e2 = e3 = e4 = 0;
  id = GetTowerID(iy0cg,    iz0cg,     nn, iyy, izz, ee);
  if( id>=0 ) e1 = ee[id];
  id = GetTowerID(iy0cg,    iz0cg+isz, nn, iyy, izz, ee);
  if( id>=0 ) e2 = ee[id];
  id = GetTowerID(iy0cg+isy,iz0cg+isz, nn, iyy, izz, ee);
  if( id>=0 ) e3 = ee[id];
  id = GetTowerID(iy0cg+isy,iz0cg,     nn, iyy, izz, ee);
  if( id>=0 ) e4 = ee[id];

  float e1t = (e1+e2+e3+e4)/etot;
  float e2t = (e1+e2-e3-e4)/etot;
  float e3t = (e1-e2-e3+e4)/etot;
  float e4t = (e3)/etot;
  //  float e5t = (e2+e4)/etot;

  float rr = sqrt((0.5-ddz)*(0.5-ddz)+(0.5-ddy)*(0.5-ddy));

  float c1, c2, c11;

  float logE = log(etot);

  // e1 energy is the most effective for PID if properly tuned !
  // Discrimination power is very sensitive to paramter c1: the bigger it is 
  // the better discrimination;
  c1 = 0.95;
  c2 = 0.0066364*logE+0.00466667;
  if( c2<0 ) c2=0;
  float e1p = c1 - c2*rr*rr;
  c1 = 0.034 - 0.01523*logE + 0.0029*logE*logE;
  float err1 = c1;

  // For e2
  c1 = 0.00844086 + 0.00645359*logE - 0.00119381*logE*logE;
  if( etot>15 ) c1 = 0.00844086 + 0.00645359*log(15.) - 0.00119381*log(15.)*log(15.); // Const at etot>15GeV
  if( c1<0 ) c1 = 0;
  c2 = 3.9; // Fixed
  float e2p = sqrt(c1+0.25*c2)-sqrt(c1+c2*ddy*ddy); // =0 at ddy=0.5

  c1 = 0.0212333 + 0.0420473/etot;
  c2 = 0.090; // Fixed
  float err2 = c1+c2*ddy;
  if( ddy>0.3 ) err2 = c1+c2*0.3; // Const at ddy>0.3

  // For e3
  c1 = 0.0107857 + 0.0056801*logE - 0.000892016*logE*logE;
  if( etot>15 ) c1 = 0.0107857 +  0.0056801*log(15.) - 0.000892016*log(15.)*log(15.); // Const at etot>15GeV
  if( c1<0 ) c1 = 0;
  c2 = 3.9; // Fixed
  float e3p = sqrt(c1+0.25*c2)-sqrt(c1+c2*ddz*ddz); // =0 at ddz=0.5

  //  c1 = 0.0200 + 0.042/etot;
  c1 = 0.0167 + 0.058/etot;
  c2 = 0.090; // Fixed
  float err3 = c1+c2*ddz;
  if( ddz>0.3 ) err3 = c1+c2*0.3; // Const at ddz>0.3

  // For e4
  float e4p = 0.25-0.668*rr+0.460*rr*rr;
  c11 = 0.171958 + 0.0142421*logE - 0.00214827*logE*logE;
  //  c11 = 0.171085 + 0.0156215*logE - -0.0025809*logE*logE;
  float err4 = 0.102-1.43*c11*rr+c11*rr*rr; // Min is set to x=1.43/2.
  err4 *= 1.1;

  chi2 = 0.;
  chi2 += (e1p-e1t)*(e1p-e1t)/err1/err1;
  chi2 += (e2p-e2t)*(e2p-e2t)/err2/err2;
  chi2 += (e3p-e3t)*(e3p-e3t)/err3/err3;
  chi2 += (e4p-e4t)*(e4p-e4t)/err4/err4;
  ndf = 4;

  //  chi2 /= 1.1;
  float prob = TMath::Prob(chi2, ndf);

  return prob;
}

float BEmcRec::ClusterChisq(int nh, EmcModule* phit, float e, float x,
				float y, int &ndf)
{

  float chi=0;
  int ixy, ix, iy;
  float et, a, d;
  
  for( int in=0; in<nh; in++ ) {
    ixy = phit[in].ich;
    iy = ixy/fNx;
    ix = ixy - iy*fNx;
    et = phit[in].amp;
    a = PredictEnergy(x-ix, y-iy, -1);
    d = fgEpar00*fgEpar00 + e*(fgEpar1*a + fgEpar2*a*a + fgEpar3*a*a*a) + 
      e*sqrt(e)*fgEpar4*a*(1-a)*fSin4T + e*e*fgEpar0*fgEpar0;
    a *= e;
    chi += (et-a)*(et-a)/d;
  }

  ndf=nh; // change needed for PbGl MV 28.01.00
  return chi;

}

// ///////////////////////////////////////////////////////////////////////////


float BEmcRec::Chi2Limit(int ND)
{
  //  Here the reverse Chi2Correct function is used
  
  float rn, a, b, chi2;
  
  if( ND < 1 ) return 9999.;  // Should we put 0. here?
  
  chi2 = fgChi2Level[EmcCluster::min(ND,50)-1];
  if( chi2 > 100 ) return 9999.; // Why should chi2 ever be >100?
  
  rn = ND;
  b = 0.072*sqrt(rn);
  a = 6.21/(sqrt(rn)+4.7);
  
  return chi2*a/(1.-chi2*b);

}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::Chi2Correct(float Chi2, int ND)
{
  // Chi2 - is reduced Chi2: Chi2/ND !!
  // MV 02.22.2000: Actually the above is not true. The supplied value of Chi2
  // has been already divided by ND. So Chi2 here is only corrected.

  float rn, a, b, c;
  
  if( ND < 1 ) return 9999.; // Should we put 0. here?
  
  rn = ND;
  b = 0.072*sqrt(rn);
  a = 6.21/(sqrt(rn)+4.7);
  c = a + b*Chi2;
  if( c < 1 ) c=1;
  
  return Chi2/c;

}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetTowerThreshold(float Thresh)
{
  fgTowerThresh = Thresh;
  fgEpar0 = Thresh*0.07;
  fgEpar00 = EmcCluster::max( (double)Thresh/3, 0.005 );
}

// **********************************************************************

void BEmcRec::getTowerPos(int ix, int iy, float &x, float & y){
  x = 2.859+5.562*ix+int(ix/12)*0.256;
  y = 2.859+5.562*iy+int(iy/12)*0.156;
}

// **********************************************************************


/// Converts coordinates in units of towers into cm's (Local coord. system)
void   BEmcRec::TowersToSector(float xT, float yT, float & xS, float & yS){
  int   x  = int(xT);
  float dx = xT - x;
  int   y  = int(yT);
  float dy = yT - y;
  xS = fModSizex*(x+0.5) + int(xT/12)*0.256 + fModSizex*dx;
  yS = fModSizey*(y+0.5) + int(yT/12)*0.156 + fModSizey*dy;
}

// **********************************************************************
/// Returns  coordinates of the tower centers in cm's (Local coord. system)
void   BEmcRec::TowersToSector(int xT,   int yT,   float & xS, float & yS){
    xS = fModSizex*(xT+0.5) + int(xT/12)*0.256;
    yS = fModSizey*(yT+0.5) + int(yT/12)*0.156;
}

// **********************************************************************
/// Converts Local Sector coordinates in cm into integer tower numbers
void   BEmcRec::SectorToTowers(float xS, float yS, int & xT,   int & yT){
  // PbSc
  xT = int(((xS-int(xS/67.0)*67.0)-0.078)/fModSizex) + 12*int(xS/67.0);
  yT = int(((yS-int(yS/66.9)*66.9)-0.078)/fModSizey) + 12*int(yS/66.9);

}

// ///////////////////////////////////////////////////////////////////////////

/* Future improvements:

1. FindClusters(): to ensure that all EmcModules are above energy threshold 
set by SetThreshold routine (or default one)

*/

// ///////////////////////////////////////////////////////////////////////////
// EOF
