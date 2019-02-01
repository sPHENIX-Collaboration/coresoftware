// Name: EmcCluster.cxx
// Author: A. Bazilevsky (RIKEN-BNL)
// Major modifications by M. Volkov (RRC KI) Jan 27 2000

#include "BEmcCluster.h"
#include "BEmcRec.h"

#include <TMath.h>

#include <cstring>

using namespace std;

// Define and initialize static members

// Max number of peaks in cluster; used in EmcCluster::GetPeaks(...)
int const EmcCluster::fgMaxNofPeaks=1000;

// Used in EmcCluster::GetPeaks(...), it is the number of iterations
// when fitting photons to peak areas
int const EmcCluster::fgPeakIter=6;

// Emin cuts cluster energy: Ecl >= Epk1+Epk2 +...+Epkn !!!!!!!!!!!!!
float const EmcCluster::fgEmin=0.002;

// chisq=3 devides about 1% of single showers (now it isn't used)
float const EmcCluster::fgChisq=3.;

// define meaningless values for (x,y)
float const EmcCluster::fgXABSURD=-999999.;
float const EmcCluster::fgYABSURD=-999999.;

EmcModule::EmcModule() 
  : ich(0),
    softKey(0),
    amp(0),
    tof(0),
    deadmap(0),
    warnmap(0),
    adc(0),
    tac(0)
{
}

//_____________________________________________________________________________
EmcModule::EmcModule(int ich_, int softkey_, float amp_, float tof_,
		     int deadmap_, int warnmap_, float adc_, float tac_) 
  : ich(ich_),
    softKey(softkey_),
    amp(amp_),
    tof(tof_),
    deadmap(deadmap_),
    warnmap(warnmap_),
    adc(adc_),
    tac(tac_)		     	            
{  
}

// ///////////////////////////////////////////////////////////////////////////
// EmcCluster member functions

void EmcCluster::GetCorrPos(float* px, float* py)
// Returns the cluster corrected position in tower units
// Corrected for S-oscilations, not for shower depth
// Shower depth (z-coord) is defined in fOwner->Tower2Global()
{
  
  float e, x, y, xx, yy, xy;
  
  e = GetTotalEnergy();
  GetMoments(&x, &y, &xx, &xy, &yy);
  fOwner->CorrectPosition(e, x, y, px, py);
}

// ///////////////////////////////////////////////////////////////////////////

void EmcCluster::GetGlobalPos( float& xA, float& yA, float& zA )
// Returns the cluster position in PHENIX global coord system
{

   float xc, yc;
   float e = GetTotalEnergy();
   GetCorrPos( &xc, &yc );
   fOwner->Tower2Global( e, xc, yc, xA, yA, zA );
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTowerEnergy( int ich )
// Returns the energy of the ich-tower (0 if ich not found in the fHitList)
{
     vector<EmcModule>::iterator ph;
     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) {
       if( (*ph).ich == ich ) return (*ph).amp;
	ph++;
     }
     return 0;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTowerEnergy( int ix, int iy )
// Returns the energy of the tower ix,iy (0 if tower not found in the fHitList)
{
     int ich, ixl, iyl;
     vector<EmcModule>::iterator ph;

     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) {
       ich = (*ph).ich;
       iyl = ich/fOwner->GetNx();
       ixl = ich % fOwner->GetNx();
       if( ixl == ix && iyl == iy ) return (*ph).amp;
	ph++;
     }
     return 0;
}

// ///////////////////////////////////////////////////////////////////////////
float EmcCluster::GetTowerToF( int ich )
// Returns the ToF of the ich-tower (0 if ich not found in the fHitList)
{
     vector<EmcModule>::iterator ph;
     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) {
       if( (*ph).ich == ich ) return (*ph).tof;
	ph++;
     }
     return 0;
}

// ///////////////////////////////////////////////////////////////////////////

int EmcCluster::GetTowerDeadMap( int ich )
// Returns the Dead Map of the ich-tower (0 if ich not found in the fHitList)
{
     vector<EmcModule>::iterator ph;
     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) {
       if( (*ph).ich == ich ) return (*ph).deadmap;
       ph++;
     }
     return 0;
}

// ///////////////////////////////////////////////////////////////////////////

int EmcCluster::GetTowerWarnMap( int ich )
  // Returns the Warning Map of the ich-tower (0 if ich not found in the fHitList)
{
     vector<EmcModule>::iterator ph;
     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) {
       if( (*ph).ich == ich ) return (*ph).warnmap;
       ph++;
     }
     return 0;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTowerADC( int ich )
  // Returns ADC of the ich-tower (0 if ich not found in the fHitList)
{
     vector<EmcModule>::iterator ph;
     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) {
       if( (*ph).ich == ich ) return (*ph).adc;
       ph++;
     }
     return 0.;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTowerTAC( int ich )
  // Returns ADC of the ich-tower (0 if ich not found in the fHitList)
{
     vector<EmcModule>::iterator ph;
     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) {
       if( (*ph).ich == ich ) return (*ph).tac;
       ph++;
     }
     return 0.;
}

// ///////////////////////////////////////////////////////////////////////////

// Returns Number of Dead Channels in 3x3 Modules around MaxTower
// see emc-calib/Calib/emcQAs.C for deadmap defn.
int EmcCluster::GetNDead()
{
  EmcModule hmax = GetMaxTower();
  int dead = hmax.deadmap;
  int ndead = 0;

  dead = dead>>4;		// start at bit 4
  for ( int iy=0; iy<3; iy++ )
    {
      for ( int iz=0; iz<3; iz++ )
        {
          ndead += dead&1;
          dead = dead>>1;	// shift to next tower
        }
      dead = dead>>2;		// shift up to next row
    }
  return ndead;
}

// ///////////////////////////////////////////////////////////////////////////

int EmcCluster::GetDeadMap()
{
  // MV 2001/12/06 Returns the deadmap of the dominant tower
  
  EmcModule hmax = GetMaxTower();
  return hmax.deadmap;

}

// ///////////////////////////////////////////////////////////////////////////

int EmcCluster::GetWarnMap()
{
  // MV 2001/12/06 Returns the warnmap of the dominant tower
  
  EmcModule hmax = GetMaxTower();
  return hmax.warnmap;

}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTotalEnergy()
// Returns the cluster total energy
{
     vector<EmcModule>::iterator ph;
     float et=0;
     if( fHitList.empty() ) return 0;
     ph = fHitList.begin();
     while( ph != fHitList.end() ) { 
	et += (*ph).amp;
	ph++;
     }
     return et;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetECoreCorrected()
// Returns the energy in core towers around the cluster Center of Gravity
// Corrected for energy leak sidewise from core towers
{
  float x, y, xx, yy, xy;
  float ecore, ecorecorr;
  ecore = GetECore();
  GetMoments(&x, &y, &xx, &xy, &yy);
  fOwner->CorrectECore(ecore, x, y, &ecorecorr);
  return ecorecorr;
}

float EmcCluster::GetECore()
// Returns the energy in core towers around the cluster Center of Gravity
{
    vector<EmcModule>::iterator ph;
    float xcg, ycg, xx, xy, yy, dx, dy;
    float energy, es, et;
    int ix, iy, ixy;

    GetMoments( &xcg, &ycg, &xx, &xy, &yy ); // Now it is in cell units
    //    xcg /= fOwner->GetModSizex();
    //    ycg /= fOwner->GetModSizey();
    energy = GetTotalEnergy();
    fOwner->SetProfileParameters(0,energy,xcg,ycg);

    es=0;
    if( fHitList.empty() ) return 0;
    ph = fHitList.begin();
    while( ph != fHitList.end() ) { 
      ixy = (*ph).ich;
      iy = ixy/fOwner->GetNx();
      ix = ixy - iy*fOwner->GetNx();
      //      dx = xcg - ix;
      dx = fOwner->fTowerDist(float(ix),xcg);
      dy = ycg - iy;
      et = fOwner->PredictEnergy(dx, dy, -1);
      if( et > 0.01 ) es += (*ph).amp;
      ph++;
    }
    return es;

}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetE4()
// Returns the energy in 2x2 towers around the cluster Center of Gravity
{
    float xcg, ycg, xx, xy, yy;
    float e1, e2, e3, e4;
    int ix0, iy0, isx, isy;

    GetMoments( &xcg, &ycg, &xx, &xy, &yy );
    xcg /= fOwner->GetModSizex();
    ycg /= fOwner->GetModSizey();
    ix0 = int(xcg+0.5);
    iy0 = int(ycg+0.5);

    isx = 1;
    if( xcg-ix0 < 0 ) isx = -1;
    isy = 1;
    if( ycg-iy0 < 0 ) isy = -1;

    e1 = GetTowerEnergy(ix0,     iy0    );
    e2 = GetTowerEnergy(ix0+isx, iy0    );
    e3 = GetTowerEnergy(ix0+isx, iy0+isy);
    e4 = GetTowerEnergy(ix0    , iy0+isy);

    return (e1 + e2 + e3 + e4);
}
// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetE9()
// Returns the energy in 3x3 towers around the cluster Center of Gravity
{
     float xcg, ycg, xx, xy, yy;
     int ich, ix0, iy0, nhit;
     float e;

     e=0;

     nhit = fHitList.size();

     if( nhit <= 0 ) return e;

     GetMoments( &xcg, &ycg, &xx, &xy, &yy );
     //     xcg /= fOwner->GetModSizex();
     //     ycg /= fOwner->GetModSizey();
     ix0 = int(xcg+0.5);
     iy0 = int(ycg+0.5);
     ich = iy0*fOwner->GetNx() + ix0;

     return GetE9(ich);

}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetE9( int ich )
// Returns the energy in 3x3 towers around the tower ich
{
  int ix, iy, ixy, idx, idy;
  vector<EmcModule>::iterator ph;
  
  int iy0 = ich/fOwner->GetNx();
  int ix0 = ich - iy0*fOwner->GetNx();

  float es=0;
  if( fHitList.empty() ) return 0;
  ph = fHitList.begin();
  while( ph != fHitList.end() ) { 
    ixy = (*ph).ich;
    iy = ixy/fOwner->GetNx();
    ix = ixy - iy*fOwner->GetNx();
    idx = fOwner->iTowerDist(ix0,ix);
    idy = iy-iy0;
    if( abs(idx)<=1 && abs(idy)<=1 ) es += (*ph).amp;
    ph++;
  }
  return es;
}

// ///////////////////////////////////////////////////////////////////////////

EmcModule EmcCluster::GetImpactTower()
// Returns the EmcModule corresponding to the reconstructed impact tower
{
     float x, y;
     int ix, iy, ich;
     EmcModule ht;

     GetCorrPos( &x, &y );
     ix = lowint(x/fOwner->GetModSizex() + 0.5);
     iy = lowint(y/fOwner->GetModSizey() + 0.5);
     if( ix < 0 || ix > fOwner->GetNx()-1 || iy < 0 || iy > fOwner->GetNy()-1 ) {
	printf("????? EmcClusterChi2: Something wrong in GetImpactTower: (x,y)=(%f,%f)  (ix,iy)=(%d,%d) \n", x, y, ix, iy);
	//!!!!!	memset(&ht, 0, sizeof(EmcModule)); // MV 2002/03/12 bugfix
	ht.ich = -1;
	ht.amp = 0;
	ht.tof = 0;
	return ht;
     }
     else {
	ich = iy*fOwner->GetNx() + ix;
	ht.ich = ich;
	ht.amp = GetTowerEnergy( ich );
	ht.tof = GetTowerToF( ich );
	ht.deadmap = GetTowerDeadMap( ich );
	ht.warnmap = GetTowerWarnMap( ich ); // MV 2002/02/18 bugfix
	ht.adc=GetTowerADC(ich); // MV 2002/03/12 bugfix
	ht.tac=GetTowerTAC(ich); // MV 2002/03/12 bugfix
	return ht;
     }
}

// ///////////////////////////////////////////////////////////////////////////

EmcModule EmcCluster::GetMaxTower()
// Returns the EmcModule with the maximum energy
{
     vector<EmcModule>::iterator ph;
     float emax=0;
     EmcModule ht;

     //!!!!!     memset(&ht, 0, sizeof(EmcModule)); // MV 2002/03/12 bugfix
     ht.ich = -1;
     ht.amp = 0;
     ht.tof = 0;
     if( fHitList.empty() ) return ht;

     ph = fHitList.begin();
     while( ph != fHitList.end() ) { 
	if( (*ph).amp > emax ) { emax = (*ph).amp; ht = *ph; }
	ph++;
     }
     return ht;
}

// ///////////////////////////////////////////////////////////////////////////

void EmcCluster::GetHits(EmcModule* phit, int n)
// Returns n EmcModules (sorted) with the maximum energy
{
     int nhit;
     EmcModule *hlist, *vv;
     vector<EmcModule>::iterator ph;
     vector<EmcModule> hl;

     fOwner->ZeroVector(phit, n);
     nhit = fHitList.size();

     if( nhit <= 0 ) return;

     hlist = new EmcModule[nhit];

     ph = fHitList.begin();
     vv = hlist;
     while( ph != fHitList.end() ) *vv++ = *ph++;

     qsort( hlist, nhit, sizeof(EmcModule), fOwner->HitACompare );
     for( int i=0; i<min(nhit,n); i++ ) phit[i]=hlist[i];
     delete [] hlist;
}

// ///////////////////////////////////////////////////////////////////////////

void EmcCluster::GetMoments( float* px, float* py, float* pxx, float* pxy, float* pyy )
//  Returns cluster 1-st (px,py) and 2-d momenta (pxx,pxy,pyy) in cell unit
{
     vector<EmcModule>::iterator ph;
     float e, x, y, xx, yy, xy;
     int nhit;
     EmcModule *phit, *p;

     *px = fgXABSURD; *py = fgYABSURD;
     *pxx = 0; *pxy = 0; *pyy = 0;
     nhit=fHitList.size();
     if( nhit <= 0 ) return;

     phit = new EmcModule[nhit];
     ph = fHitList.begin();
     p=phit;
     while( ph != fHitList.end() ) { 
	p->ich = (*ph).ich;
	p->amp = (*ph).amp;
	ph++;
	p++;
     }
     fOwner->Momenta( nhit, phit, &e, &x, &y, &xx, &yy, &xy );
     /*
     *px = x*fOwner->GetModSizex();
     *py = y*fOwner->GetModSizey();
     *pxx = xx*fOwner->GetModSizex()*fOwner->GetModSizex();
     *pxy = xy*fOwner->GetModSizex()*fOwner->GetModSizey();
     *pyy = yy*fOwner->GetModSizey()*fOwner->GetModSizey();
     */
     *px = x;
     *py = y;
     *pxx = xx;
     *pxy = xy;
     *pyy = yy;

     delete [] phit;
}

// ///////////////////////////////////////////////////////////////////////////

void EmcCluster::GetErrors( float* pde, float* pdx, float* pdy, float* pdz)
{
  //  Returns the errors for the reconstructed energy and position

    float e, x, y;

    e = GetTotalEnergy();
    GetCorrPos( &x, &y );
    fOwner->CalculateErrors( e, x, y, pde, pdx, pdy, pdz);

}

float EmcCluster::GetProb(float &chi2, int &ndf)
{
  return fOwner->GetProb(fHitList, chi2, ndf);
}

// ///////////////////////////////////////////////////////////////////////////

void EmcCluster::GetChar( float* pe, 
			  float* pxcg, float* pycg, 
			  float* pxc, float* pyc, 
			  float* pxg, float* pyg, float* pzg,
			  float* pxx, float* pxy, float* pyy,
			  float* pde, float* pdx, float* pdy, float* pdz  )
{
  // This method replaces methods GetTotalEnergy, GetMoments, 
  // GetCorrPos, GetGlobalPos, GetErrors

  *pe = GetTotalEnergy();
  GetMoments( pxcg, pycg, pxx, pxy, pyy );
  fOwner->CorrectPosition( *pe, *pxcg, *pycg, pxc, pyc );
  fOwner->SectorToGlobal( *pxc, *pyc, 0, pxg, pyg, pzg );
  fOwner->CalculateErrors( *pe, *pxc, *pyc, pde, pdx, pdy, pdz);

}

// ///////////////////////////////////////////////////////////////////////////

int EmcCluster::GetPeaks( vector<EmcPeakarea> *PkList, vector<EmcModule> *ppeaks )
{
  // Splits the cluster onto peakarea's
  // The number of peakarea's is equal to the number of Local Maxima 
  // in the cluster. Local Maxima can have the energy not less then 
  // defined in fgMinPeakEnergy
  //
  // Output: PkList - vector of peakarea's
  //         ppeaks - vector of peak EmcModules (one for each peakarea)
  //
  // Returns: >= 0 Number of Peaks;
  //	      -1 The number of Peaks is greater then fgMaxNofPeaks;
  //		 (just increase parameter fgMaxNofPeaks)
  
  int npk, ipk, nhit;
  int ixypk, ixpk, iypk, in, nh, ic;
  int ixy, ix, iy, nn;
  int idx, idy;
  int ig, ng, igmpk1[fgMaxNofPeaks], igmpk2[fgMaxNofPeaks];
  int PeakCh[fgMaxNofPeaks];
  float epk[fgMaxNofPeaks*2], xpk[fgMaxNofPeaks*2], ypk[fgMaxNofPeaks*2];
  float ratio, eg, dx, dy, a;
  float *Energy[fgMaxNofPeaks], *totEnergy, *tmpEnergy;
  EmcModule *phit, *hlist, *vv;
  EmcPeakarea peak(fOwner);
  vector<EmcModule>::iterator ph;
  vector<EmcModule> hl;
  
  PkList->clear();
  ppeaks->clear();

  nhit = fHitList.size();
  
  if( nhit <= 0 ) return 0;

  hlist  = new EmcModule[nhit];
  
  ph = fHitList.begin();
  vv = hlist;
  while( ph != fHitList.end() ) *vv++ = *ph++;

  // sort by linear channel number
  qsort( hlist, nhit, sizeof(EmcModule), fOwner->HitNCompare );

  //
  //  Find peak (maximum) position (towers with local maximum amp)
  //
  npk=0;
  for( ic=0; ic<nhit; ic++ ) {
    float amp = hlist[ic].amp;
    if( amp > fOwner->GetPeakThreshold() ) {
      int ich = hlist[ic].ich;
      //old      int ichmax = ich + fOwner->GetNx() + 1;
      //old      int ichmin = ich - fOwner->GetNx() - 1;
      // Look into three raws only
      int ichmax = (ich/fOwner->GetNx()+2)*fOwner->GetNx() - 1;
      int ichmin = (ich/fOwner->GetNx()-1)*fOwner->GetNx();
      int ixc = ich - ich/fOwner->GetNx()*fOwner->GetNx();
      int in = ic + 1;
      // check right, and 3 towers above if there is another max
      while( (in < nhit) && (hlist[in].ich <= ichmax) ) {
	int ix = hlist[in].ich - hlist[in].ich/fOwner->GetNx()*fOwner->GetNx();
	//old	if( (abs(ixc-ix) <= 1) && (hlist[in].amp >= amp) ) goto new_ic;
	if( (abs(fOwner->iTowerDist(ix,ixc)) <= 1) && (hlist[in].amp >= amp) ) goto new_ic;
	in++;
      }
      in = ic - 1;
      // check left, and 3 towers below if there is another max
      while( (in >= 0) && (hlist[in].ich >= ichmin) ) {
	int ix = hlist[in].ich - hlist[in].ich/fOwner->GetNx()*fOwner->GetNx();
	//old	if( (abs(ixc-ix) <= 1) && (hlist[in].amp > amp) ) goto new_ic;
	if( (abs(fOwner->iTowerDist(ix,ixc)) <= 1) && (hlist[in].amp > amp) ) goto new_ic;
	in--;
      }
      if( npk >= fgMaxNofPeaks ) { 
	delete [] hlist; 
	printf("!!! Error in EmcCluster::GetPeaks(): too many peaks in a cluster (>%d). May need tower energy threshold increase for clustering.\n",fgMaxNofPeaks);
	return -1; 
      }

      // ic is a maximum in a 3x3 tower group
      PeakCh[npk]=ic;
      npk++;
    }
  new_ic: continue;
  }
  /*
  for( ipk=0; ipk<npk; ipk++ ) {
    ic = PeakCh[ipk];
    printf("  %d: E=%f\n", ipk, hlist[ic].amp);
  }
  */  

  // there was only one peak
  if( npk <= 1 ) {
    hl.clear();
    for( int ich=0; ich<nhit; ich++ ) hl.push_back(hlist[ich]);
    peak.ReInitialize(hl);
    PkList->push_back(peak);

    if( npk == 1 ) ppeaks->push_back(hlist[PeakCh[0]]);
    else           ppeaks->push_back(GetMaxTower());

    delete [] hlist;
    return 1;
  }
 
  // there were more than one peak

  for( ipk=0; ipk<npk; ipk++ ) { Energy[ipk] = new float[nhit]; }
  totEnergy = new float[nhit];
  tmpEnergy = new float[nhit];
  for ( int i = 0; i < nhit; ++i ) 
    {
      totEnergy[i]=0.0;
      tmpEnergy[i]=0.0;
    }
  //
  // Divide energy in towers among photons positioned in every peak
  //
  ratio = 1.;
  for( int iter=0; iter<fgPeakIter; iter++ ) {
    fOwner->ZeroVector( tmpEnergy, nhit );
    for( ipk=0; ipk<npk; ipk++ ) {
      
      ic = PeakCh[ipk];
      if( iter > 0 ) ratio = Energy[ipk][ic]/totEnergy[ic];
      eg = hlist[ic].amp * ratio;
      ixypk = hlist[ic].ich;
      iypk = ixypk/fOwner->GetNx();
      ixpk = ixypk - iypk*fOwner->GetNx();
      epk[ipk] = eg;
      xpk[ipk] = 0.;	// CoG from the peak tower
      ypk[ipk] = 0.;	// CoG from the peak tower
      
      int ichmax = (iypk+2)*fOwner->GetNx() - 1;
      int ichmin = (iypk-1)*fOwner->GetNx();

      // add up energies of tower to the right and up
      if( ic < nhit-1 ){
	for( in=ic+1; in<nhit; in++ ) {
	  ixy = hlist[in].ich;
	  iy = ixy/fOwner->GetNx();
	  ix = ixy - iy*fOwner->GetNx();

	  //old	  if( (ixy-ixypk) > fOwner->GetNx()+1 ) break;
	  if( ixy > ichmax ) break;
	  //old	  if( abs(ix-ixpk) <= 1 ) {
	  idx = fOwner->iTowerDist(ixpk,ix);
	  idy = iy - iypk;
	  if( abs(idx) <= 1 ) {
	    if( iter > 0 ) ratio = Energy[ipk][in]/totEnergy[in];
	    eg = hlist[in].amp * ratio;
	    epk[ipk] += eg;
	    xpk[ipk] += eg*idx;
	    ypk[ipk] += eg*idy;
	  }
	}
      } // if ic
      
      // add up energies of tower to the left and down
      if( ic >= 1 ){
	for( in=ic-1; in >= 0; in-- ) {
	  ixy = hlist[in].ich;
	  iy = ixy/fOwner->GetNx();
	  ix = ixy - iy*fOwner->GetNx();

	  //old	  if( (ixypk-ixy) > fOwner->GetNx()+1 ) break;
	  if( ixy < ichmin ) break;
	  //old	  if( abs(ix-ixpk) <= 1 ) {
	  idx = fOwner->iTowerDist(ixpk,ix);
	  idy = iy - iypk;
	  if( abs(idx) <= 1 ) {
	    if( iter > 0 ) ratio = Energy[ipk][in]/totEnergy[in];
	    eg = hlist[in].amp * ratio;
	    epk[ipk] += eg;
	    xpk[ipk] += eg*idx;
	    ypk[ipk] += eg*idy;
	  }
	}
      } // if ic
      
      xpk[ipk] = xpk[ipk]/epk[ipk] + ixpk;
      ypk[ipk] = ypk[ipk]/epk[ipk] + iypk;
      fOwner->SetProfileParameters(0, epk[ipk], xpk[ipk], ypk[ipk]);

      for( in=0; in<nhit; in++ ) {
	ixy = hlist[in].ich;
	iy = ixy/fOwner->GetNx();
	ix = ixy - iy*fOwner->GetNx();
	dx = fOwner->fTowerDist(float(ix),xpk[ipk]);
	dy = ypk[ipk]-iy;
	a = 0;
	
        // predict energy within 2.5 cell square around local peak
	if( ABS(dx) < 2.5 && ABS(dy) < 2.5 )
	  a = epk[ipk]*fOwner->PredictEnergy(dx, dy, -1);
	
	Energy[ipk][in] = a;
	tmpEnergy[in] += a;
      }
      
    } // for ipk
    float flim = 0.000001;
    for ( in = 0; in < nhit; in++ )
      {
	if (tmpEnergy[in] > flim)
	  {
	    totEnergy[in] = tmpEnergy[in];
	  }
	else
	  {
	    totEnergy[in] = flim;
	  }
      }
  } // for iter
  
  phit = new EmcModule[nhit];
  
  ng=0;
  for( ipk=0; ipk<npk; ipk++ ) {
    nh=0;

    // fill phit, the tower distribution for a peak
    for( in=0; in<nhit; in++ ) {
      if( tmpEnergy[in] > 0 ) {
	ixy = hlist[in].ich;
	a = hlist[in].amp * Energy[ipk][in]/tmpEnergy[in];
	if( a > fgEmin ) {
	  phit[nh].ich=ixy;
	  phit[nh].amp=a;
	  phit[nh].tof= hlist[in].tof;  // Not necessary here
	  phit[nh].deadmap= hlist[in].deadmap;  // Not necessary here
	  phit[nh].warnmap= hlist[in].warnmap;  // MV 2002/02/18 bugfix

	  // MV 2002/03/12 bugfix: let adc, tdc=0 for split clusters
	  //	  phit[nh].adc=0.; // MV 2002/03/12 bugfix
	  //	  phit[nh].tac=0.; // MV 2002/03/12 bugfix
	  phit[nh].adc=hlist[in].adc; // GD 04/12/2006
	  phit[nh].tac=hlist[in].tac; // GD 04/12/2006

	  nh++;
	}
      }
    }

    // if number hits > 0
    if( nh>0 ) {
      // !!!!! Exclude Gamma() for now until it is proved working and useful
      /*
      chi=fOwner->Chi2Limit(nh)*10;
      int ndf; // just a plug for changed Gamma parameter list MV 28.01.00

      fOwner->Gamma(nh, phit,&chi, &chi0, &epk[ng], &xpk[ng], &ypk[ng],
		    &epk[ng+1], &xpk[ng+1], &ypk[ng+1], ndf);

      igmpk1[ipk]=ng;
      igmpk2[ipk]=ng;
      if( epk[ng+1]>0 ) { ng++; igmpk2[ipk]=ng;}
      */
      igmpk1[ipk]=ng;
      igmpk2[ipk]=ng;
      
      ng++;
    }
    else {
      igmpk1[ipk]=-1;
      igmpk2[ipk]=-1;
    }
  } // for( ipk=
  
  fOwner->ZeroVector( tmpEnergy, nhit );
  for( ipk=0; ipk<npk; ipk++ ) {
    ig=igmpk1[ipk];
    if( ig >= 0 ) {
      //      printf("  %d: X=%f Y=%f\n",ipk,xpk[ig], ypk[ig]);
      fOwner->SetProfileParameters(0, epk[ig], xpk[ig], ypk[ig]);
      for( in=0; in<nhit; in++ ) {
	Energy[ipk][in]=0;
	for(ig=igmpk1[ipk]; ig<=igmpk2[ipk]; ig++){
	  ixy = hlist[in].ich;
	  iy = ixy/fOwner->GetNx();
	  ix = ixy - iy*fOwner->GetNx();
	  dx = fOwner->fTowerDist(float(ix),xpk[ig]);
	  dy = ypk[ig]-iy;
	  a = epk[ig]*fOwner->PredictEnergy(dx, dy, epk[ig]);
	  Energy[ipk][in] += a;
	  tmpEnergy[in] += a;
	}
      } // for( in
    } // if( ig >= 0
  } // for( ipk

  nn=0;
  for( ipk=0; ipk<npk; ipk++ ) {
    nh=0;
    for( in=0; in<nhit; in++ ) {
      if( tmpEnergy[in] > 0 ) {
	ixy = hlist[in].ich;
	a = hlist[in].amp * Energy[ipk][in]/tmpEnergy[in];
	if( a > fgEmin ) {
	  phit[nh].ich=ixy;
	  phit[nh].amp=a;
	  phit[nh].tof=hlist[in].tof;
	  phit[nh].deadmap=hlist[in].deadmap;
	  phit[nh].warnmap=hlist[in].warnmap; // MV 2002/02/18 bugfix

	  // MV 2002/03/12 bugfix: let adc, tdc=0 for split clusters
	  //	  phit[nh].adc=0.; // MV 2002/03/12 bugfix
	  //	  phit[nh].tac=0.; // MV 2002/03/12 bugfix
	  phit[nh].adc=hlist[in].adc; // GD 04/12/2006
	  phit[nh].tac=hlist[in].tac; // GD 04/12/2006

	  nh++;
	}
      }
    } // for( in
    if( nh>0 ) {
      //      *ip++ = hlist[PeakCh[ipk]];
      ppeaks->push_back(hlist[PeakCh[ipk]]);
      hl.clear();
      for( in=0; in<nh; in++ ) hl.push_back(phit[in]);
      peak.ReInitialize(hl);
      PkList->push_back(peak);
      nn++;
    }
  } // for( ipk
  
  delete [] phit;
  delete [] hlist;
  for( ipk=0; ipk<npk; ipk++ ) delete [] Energy[ipk];
  delete [] totEnergy;
  delete [] tmpEnergy;

  return nn;
}

// ///////////////////////////////////////////////////////////////////////////
// EmcPeakarea member functions

void EmcPeakarea::GetChar( float* pe, float* pec, 
			   float* pecore, float* pecorec, 
			   float* pxcg, float* pycg, 		// center of gravity
			   float* pxcgmin, float* pycgmin, 	// 
			   float* pxc, float* pyc, 		// Local (Sector) coords
			   float* pxg, float* pyg, float* pzg,	// Global coords
			   float* pxx, float* pxy, float* pyy,	// moments
			   float* pchi,
			   float* pde, float* pdx, float* pdy, float* pdz  )
// This method replaces "cluster" methods GetTotalEnergy, GetMoments,
// GetCorrPos, GetGlobalPos, GetErrors and "EmcPeakarea" methods GetCGmin, GetChi2
{
    float chi, chi0;
    float e1, x1, y1, e2, x2, y2;
    int nh;
    EmcModule *phit, *vv;
    vector<EmcModule>::iterator ph;
    vector<EmcModule> hl;
    float tmplvalue; // temporary lvalue

    *pe = 0;
    hl = fHitList;
    nh = hl.size();
    if( nh <= 0 ) return;

    phit = new EmcModule[nh];

    ph = hl.begin();
    vv = phit;
    while( ph != hl.end() ) *vv++ = *ph++;

    chi = fgChisq*1000;
    int ndf; // Gamma parameter list changed MV 28.01.00
    fOwner->Gamma(nh, phit, &chi, &chi0, &e1, &x1, &y1, &e2, &x2, &y2, ndf);
    fNdf=ndf;
    *pchi=chi0;

    // Calculate CL
    tmplvalue=fOwner->Chi2Correct(chi0, ndf)*ndf; // nh->ndf MV 28.01.00
    if(tmplvalue>0.) fCL=TMath::Prob(tmplvalue, ndf);
    else fCL=1.; // nothing to say about this peak area

    // Shower Center of Gravity after shower profile fit
    *pxcgmin = x1*fOwner->GetModSizex();
    *pycgmin = y1*fOwner->GetModSizey();

    *pe = GetTotalEnergy();
    *pecore = GetECore();
    GetMoments( pxcg, pycg, pxx, pxy, pyy );
    fOwner->CorrectPosition(*pe, *pxcgmin, *pycgmin, pxc, pyc);
    fOwner->SectorToGlobal( *pxc, *pyc, 0, pxg, pyg, pzg );
    fOwner->CalculateErrors( *pe, *pxc, *pyc, pde, pdx, pdy, pdz);
    //    fOwner->CorrectEnergy( *pe, *pxcg, *pycg, pec ); // MM 02.11.2000
    fOwner->CorrectECore( *pecore, *pxc, *pyc, pecorec );

    delete [] phit;

}

// ///////////////////////////////////////////////////////////////////////////

float EmcPeakarea::GetChi2()
// Returns Hi2 after its minimization fluctuating CG position 
// (i.e. after shower profile fit)
{
    float chi, chi0;
    float e1, x1, y1, e2, x2, y2;
    int nh;
    EmcModule *phit, *vv;
    vector<EmcModule>::iterator ph;
    vector<EmcModule> hl;

    hl = fHitList;
    nh = hl.size();
    if( nh <= 0 ) return 0;

    phit = new EmcModule[nh];

    ph = hl.begin();
    vv = phit;
    while( ph != hl.end() ) *vv++ = *ph++;

    chi = fgChisq*1000;
    int ndf; // Gamma parameter list changed MV 28.01.00
    fOwner->Gamma(nh, phit, &chi, &chi0, &e1, &x1, &y1, &e2, &x2, &y2, ndf);
    fNdf=ndf;
    delete [] phit;
    return chi0;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcPeakarea::GetCLNew()
// Conf. Level based on new chi2 calculation
{
    float xcg, ycg, xx, xy, yy;
    float e1, e2, e3, e4;
    float e1m, e2m, e3m, e4m;
    float e1p, e2p, e3p, e4p;
    float s1, s2, s3, s4;
    float dx, dy, dt, et, etot, sc;
    float chi2, pr;
    int ix0, iy0, isx, isy, ndf;

    etot = GetTotalEnergy();
    GetMoments( &xcg, &ycg, &xx, &xy, &yy );
    xcg /= fOwner->GetModSizex();
    ycg /= fOwner->GetModSizey();
    ix0 = int(xcg+0.5);
    iy0 = int(ycg+0.5);
    dx = ABS(xcg-ix0);
    dy = ABS(ycg-iy0);

    isx = 1;
    if( xcg-ix0 < 0 ) isx = -1;
    isy = 1;
    if( ycg-iy0 < 0 ) isy = -1;

    e1 = GetTowerEnergy(ix0,     iy0    );
    e2 = GetTowerEnergy(ix0+isx, iy0    );
    e3 = GetTowerEnergy(ix0+isx, iy0+isy);
    e4 = GetTowerEnergy(ix0    , iy0+isy);

    if( dy > dx ) {
      et = e2;
      e2 = e4;
      e4 = et;
      dt = dx;
      dx = dy;
      dy = dt;
    }

    e1m = e1 + e2 + e3 + e4;
    e2m = e1 + e2 - e3 - e4;
    e3m = e1 - e2 - e3 + e4;
    e4m = e4 - e3;

    e1p = 0.932;
    e2p = 0.835 - 2*dy*dy/(dy+0.099);
    e3p = 0.835 - 2*dx*dx/(dx+0.099);
    e4p = 0.02;

    sc = sqrt(0.1*0.1/etot + 0.03*0.03)/0.04;
    s1 = sc*0.02;
    sc = sqrt(0.1*0.1/etot + 0.02*0.02)/0.04;
    s2 = sc*(0.056-0.026*e2p);
    s3 = sc*(0.056-0.026*e2p);
    s4 = sc*0.03;

    chi2 = 0.;
    chi2 += (e1p*etot-e1m)*(e1p*etot-e1m)/s1/s1/etot/etot;
    chi2 += (e2p*etot-e2m)*(e2p*etot-e2m)/s2/s2/etot/etot;
    chi2 += (e3p*etot-e3m)*(e3p*etot-e3m)/s3/s3/etot/etot;
    chi2 += (e4p*etot-e4m)*(e4p*etot-e4m)/s4/s4/etot/etot;
    chi2 /= 0.7;

    ndf = 4;
    pr = TMath::Prob(chi2, ndf);
    return pr;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcPeakarea::GetChi2New()
// Conf. Level based on new chi2 calculation
{
    float xcg, ycg, xx, xy, yy;
    float e1, e2, e3, e4;
    float e1m, e2m, e3m, e4m;
    float e1p, e2p, e3p, e4p;
    float s1, s2, s3, s4;
    float dx, dy, dt, et, etot, sc;
    float chi2;
    int ix0, iy0, isx, isy, ndf;

    etot = GetTotalEnergy();
    GetMoments( &xcg, &ycg, &xx, &xy, &yy );
    xcg /= fOwner->GetModSizex();
    ycg /= fOwner->GetModSizey();
    ix0 = int(xcg+0.5);
    iy0 = int(ycg+0.5);
    dx = ABS(xcg-ix0);
    dy = ABS(ycg-iy0);

    isx = 1;
    if( xcg-ix0 < 0 ) isx = -1;
    isy = 1;
    if( ycg-iy0 < 0 ) isy = -1;

    e1 = GetTowerEnergy(ix0,     iy0    );
    e2 = GetTowerEnergy(ix0+isx, iy0    );
    e3 = GetTowerEnergy(ix0+isx, iy0+isy);
    e4 = GetTowerEnergy(ix0    , iy0+isy);

    if( dy > dx ) {
      et = e2;
      e2 = e4;
      e4 = et;
      dt = dx;
      dx = dy;
      dy = dt;
    }

    e1m = e1 + e2 + e3 + e4;
    e2m = e1 + e2 - e3 - e4;
    e3m = e1 - e2 - e3 + e4;
    e4m = e4 - e3;

    e1p = 0.932;
    e2p = 0.835 - 2*dy*dy/(dy+0.099);
    e3p = 0.835 - 2*dx*dx/(dx+0.099);
    e4p = 0.02;

    sc = sqrt(0.1*0.1/etot + 0.03*0.03)/0.04;
    s1 = sc*0.02;
    sc = sqrt(0.1*0.1/etot + 0.02*0.02)/0.04;
    s2 = sc*(0.056-0.026*e2p);
    s3 = sc*(0.056-0.026*e2p);
    s4 = sc*0.03;

    chi2 = 0.;
    chi2 += (e1p*etot-e1m)*(e1p*etot-e1m)/s1/s1/etot/etot;
    chi2 += (e2p*etot-e2m)*(e2p*etot-e2m)/s2/s2/etot/etot;
    chi2 += (e3p*etot-e3m)*(e3p*etot-e3m)/s3/s3/etot/etot;
    chi2 += (e4p*etot-e4m)*(e4p*etot-e4m)/s4/s4/etot/etot;
    chi2 /= 0.7;

    ndf = 4;
    return chi2/ndf;
}

// ///////////////////////////////////////////////////////////////////////////

void EmcPeakarea::GetCGmin( float* px, float* py )
{
  // Gets CG coordinates corresponding to min Hi2 (after shower shape fit)

  float chi, chi0;
  float e1, x1, y1, e2, x2, y2;
  int nh;
  EmcModule *phit, *vv;
  vector<EmcModule>::iterator ph;
  vector<EmcModule> hl;
  
  *px = fgXABSURD;
  *py = fgYABSURD;
  hl = fHitList;
  nh = hl.size();
  if( nh <= 0 ) return;
  
  phit = new EmcModule[nh];
  
  ph = hl.begin();
  vv = phit;
  while( ph != hl.end() ) *vv++ = *ph++;
  
  chi = fgChisq*1000;
  int ndf; // Gamma parameter list changed MV 28.01.00
  fOwner->Gamma(nh, phit, &chi, &chi0, &e1, &x1, &y1, &e2, &x2, &y2, ndf);
  fNdf=ndf;
  *px = x1*fOwner->GetModSizex();
  *py = y1*fOwner->GetModSizey();
  
  delete [] phit;

}

