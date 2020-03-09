// Name: BEmcRec.h
// Author: A. Bazilevsky, Apr 2012
// Modified from EmcSectorRec.cxx and EmcScSectorRec.cxx

// BEmcRec -- base class for sPHENIX EMCal

// ///////////////////////////////////////////////////////////////////////////

#include "BEmcRec.h"
#include "BEmcCluster.h"

#include <TMath.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>

using namespace std;

// Define and initialize static members

// Minimal shower energy when splitting peakarea onto showers, used in Gamma()
float const BEmcRec::fgMinShowerEnergy = 0.1;

// Max number of clusters in sector, used in Find_Clusters()
int const BEmcRec::fgMaxLen = 10000;

// Default level, now for the Conf Level: 1% for GEANT, 2%-5% for TestBeam

float BEmcRec::fgChi2Level[50] = {
    6.634899, 4.605171, 3.780564, 3.318915, 3.017103,
    2.801872, 2.639259, 2.511249, 2.407341, 2.320967,
    2.247720, 2.184744, 2.129863, 2.081515, 2.038526,
    1.999994, 1.965214, 1.933627, 1.904781, 1.878311,
    1.853912, 1.831334, 1.810365, 1.790825, 1.772564,
    1.755449, 1.739367, 1.724222, 1.709926, 1.696406,
    1.683593, 1.671430, 1.659864, 1.648850, 1.638344,
    1.628311, 1.618716, 1.609528, 1.600721, 1.592268,
    1.584148, 1.576338, 1.568822, 1.561579, 1.554596,
    1.547856, 1.541346, 1.535055, 1.528968, 1.523077};

// For the Conf Level: 1% for GEANT, 2%-5% for TestBeam
float BEmcRec::fgChi2Level1[50] = {
    6.634899, 4.605171, 3.780564, 3.318915, 3.017103,
    2.801872, 2.639259, 2.511249, 2.407341, 2.320967,
    2.247720, 2.184744, 2.129863, 2.081515, 2.038526,
    1.999994, 1.965214, 1.933627, 1.904781, 1.878311,
    1.853912, 1.831334, 1.810365, 1.790825, 1.772564,
    1.755449, 1.739367, 1.724222, 1.709926, 1.696406,
    1.683593, 1.671430, 1.659864, 1.648850, 1.638344,
    1.628311, 1.618716, 1.609528, 1.600721, 1.592268,
    1.584148, 1.576338, 1.568822, 1.561579, 1.554596,
    1.547856, 1.541346, 1.535055, 1.528968, 1.523077};

// For the Conf Level: 2% for GEANT, 4%-7% for TestBeam
float BEmcRec::fgChi2Level2[50] = {
    5.411895, 3.912024, 3.278443, 2.916812, 2.677547,
    2.505458, 2.374582, 2.271008, 2.186567, 2.116065,
    2.056169, 2.004491, 1.959343, 1.919481, 1.883964,
    1.852072, 1.823237, 1.797008, 1.773021, 1.750981,
    1.730640, 1.711795, 1.694274, 1.677931, 1.662643,
    1.648301, 1.634814, 1.622101, 1.610093, 1.598727,
    1.587948, 1.577709, 1.567968, 1.558684, 1.549824,
    1.541357, 1.533256, 1.525494, 1.518051, 1.510903,
    1.504033, 1.497424, 1.491059, 1.484924, 1.479006,
    1.473292, 1.467771, 1.462433, 1.457267, 1.452265};

// ///////////////////////////////////////////////////////////////////////////
// BEmcRec member functions

BEmcRec::BEmcRec()
  : bCYL(true)
  , fNx(-99999)
  , fNy(-99999)
  , fModSizex(NAN)
  , fModSizey(NAN)
  , fVx(NAN)
  , fVy(NAN)
  , fVz(NAN)
  , fgTowerThresh(NAN)
  , fgMinPeakEnergy(NAN)
  , fSin4T(NAN)
  , fSinTx(NAN)
  , fSinTy(NAN)
  , fPpar1(NAN)
  , fPpar2(NAN)
  , fPpar3(NAN)
  , fPpar4(NAN)
  , fPshiftx(NAN)
  , fPshifty(NAN)
    //  , _emcprof(nullptr)
{
  fModules = new vector<EmcModule>;
  fClusters = new vector<EmcCluster>;
  SetPeakThreshold(0.08);
  SetChi2Limit(2);
}

// ///////////////////////////////////////////////////////////////////////////

BEmcRec::~BEmcRec()
{
  //  if (_emcprof) delete _emcprof;

  if (fModules)
  {
    fModules->clear();
    delete fModules;
  }

  if (fClusters)
  {
    fClusters->clear();
    delete fClusters;
  }
}

// ///////////////////////////////////////////////////////////////////////////
void BEmcRec::LoadProfile(const char *fname) 
{
  printf("Warning from BEmcRec::LoadProfile(): No acton defined for shower profile evaluation; should be defined in a detector specific module BEmcRec{Name}\n");
}


void BEmcRec::SetGeometry(int nx, int ny, float txsz, float tysz)
{
  fNx = nx;
  fNy = ny;
  fModSizex = txsz;
  fModSizey = tysz;
}

void BEmcRec::PrintTowerGeometry(const char* fname)
{
  FILE* pf = fopen(fname, "w");
  if (!pf)
  {
    printf("Error in BEmcRec::PrintTowerGeometry(): Failed to open file %s\n", fname);
    return;
  }

  //  printf("Info: Print from BEmcRec::PrintTowerGeometry():\n");
  //  printf("      Number of bins: %d %d\n",fNx,fNy);
  fprintf(pf, "Number of bins:\n%d %d\n", fNx, fNy);
  fprintf(pf, "ix iy x y z dx0 dy0 dz0 dx1 dy1 dz1\n");
  int ich;
  TowerGeom geom;
  std::map<int, TowerGeom>::iterator it;
  for (int iy = 0; iy < fNy; iy++)
  {
    for (int ix = 0; ix < fNx; ix++)
    {
      ich = iy * fNx + ix;
      it = fTowerGeom.find(ich);
      if (it != fTowerGeom.end())
      {
        geom = it->second;
        //	printf("       %d %d: %f %f %f\n",ix,iy,geom.Xcenter,geom.Ycenter,geom.Zcenter);
        fprintf(pf, "%d %d %f %f %f %f %f %f %f %f %f\n", ix, iy, geom.Xcenter, geom.Ycenter, geom.Zcenter, geom.dX[0],  geom.dY[0],  geom.dZ[0], geom.dX[1], geom.dY[1],  geom.dZ[1]);
      }
    }
  }

  fclose(pf);
}

bool BEmcRec::GetTowerGeometry(int ix, int iy, TowerGeom& geom)
{
  if (ix < 0 || ix >= fNx || iy < 0 || iy >= fNy) return false;

  int ich = iy * fNx + ix;
  std::map<int, TowerGeom>::iterator it = fTowerGeom.find(ich);
  if (it == fTowerGeom.end()) return false;

  geom = it->second;
  return true;
}

bool BEmcRec::SetTowerGeometry(int ix, int iy, float xx, float yy, float zz)
{
  if (ix < 0 || ix >= fNx || iy < 0 || iy >= fNy) return false;

  TowerGeom geom;
  geom.Xcenter = xx;
  geom.Ycenter = yy;
  geom.Zcenter = zz;
  geom.dX[0] = geom.dX[1] = 0; // These should be calculated by CompleteTowerGeometry()
  geom.dY[0] = geom.dY[1] = 0;
  geom.dZ[0] = geom.dZ[1] = 0;

  int ich = iy * fNx + ix;
  fTowerGeom[ich] = geom;
  return true;
}

bool BEmcRec::CompleteTowerGeometry()
// Calculates tower front size from coordinates of tower center coordinates
{
  if( fTowerGeom.empty() || fNx <= 0 ) {
    printf("Error in BEmcRec::CalculateTowerSize(): Tower geometry not well setup (NX=%d)\n",fNx);
    return false;
  }

  std::map<int, TowerGeom>::iterator it;

  for( it = fTowerGeom.begin(); it != fTowerGeom.end(); it++ ){

    int ich = it->first;
    TowerGeom geom0 = it->second;
    int ix = ich%fNx;
    int iy = ich/fNx;

    // Next tower in x
    TowerGeom geomx;
    int idx = 0;
    if (ix < fNx / 2)
      {
	idx += 1;
	while (!GetTowerGeometry(ix + idx, iy, geomx) && idx < fNx / 2) idx += 1;
      }
    else
      {
	idx -= 1;
	while (!GetTowerGeometry(ix + idx, iy, geomx) && idx > -fNx / 2) idx -= 1;
      }
    if (idx >= fNx / 2 || idx <= -fNx / 2)
      {
	printf("Error in BEmcRec::CalculateTowerSize(): Error when locating neighbour for (ix,iy)=(%d,%d)\n", ix, iy);
	return false;
      }
    
    // Next tower in y
    TowerGeom geomy;
    int idy = 0;
    if (iy < fNy / 2)
      {
	idy += 1;
	while (!GetTowerGeometry(ix, iy + idy, geomy) && idy < fNy / 2) idy += 1;
      }
    else
      {
	idy -= 1;
	while (!GetTowerGeometry(ix, iy + idy, geomy) && idy > -fNy / 2) idy -= 1;
      }
    if (idy >= fNy / 2 || idy <= -fNy / 2)
      {
	printf("Error in BEmcRec::CalculateTowerSize(): Error when locating neighbour for (ix,iy)=(%d,%d)\n", ix, iy);
	return false;
      }

    geom0.dX[0] = (geomx.Xcenter - geom0.Xcenter) / float(idx);
    geom0.dY[0] = (geomx.Ycenter - geom0.Ycenter) / float(idx);
    geom0.dZ[0] = (geomx.Zcenter - geom0.Zcenter) / float(idx);
    geom0.dX[1] = (geomy.Xcenter - geom0.Xcenter) / float(idy);
    geom0.dY[1] = (geomy.Ycenter - geom0.Ycenter) / float(idy);
    geom0.dZ[1] = (geomy.Zcenter - geom0.Zcenter) / float(idy);

    it->second = geom0;

  } // it = fTowerGeom.begin()
  
  return true;
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

void BEmcRec::Tower2Global(float E, float xC, float yC,
                               float& xA, float& yA, float& zA)
// xC and yC are local position in tower units
// For CYL geometry (xC,yC) is actually (phiC,zC)
{

  xA = 0;
  yA = 0;
  zA = 0;

  int ix = xC + 0.5;  // tower #
  if (ix < 0 || ix >= fNx)
  {
    printf("Error in BEmcRec::SectorToGlobal: wrong input x: %d\n", ix);
    return;
  }

  int iy = yC + 0.5;  // tower #
  if (iy < 0 || iy >= fNy)
  {
    printf("Error in BEmcRec::SectorToGlobal: wrong input y: %d\n", iy);
    return;
  }

  // Get tower where the shower is positioned
  TowerGeom geom0;

  if (!GetTowerGeometry(ix, iy, geom0)) { 
    // Weird case: cluster center of gravity outside the EMCal, take geometry from the neighbouring tower
    int idx[4] = {1,0,-1, 0};
    int idy[4] = {0,1, 0,-1};
    int ii = 0;
    while( ii<4 && !GetTowerGeometry(ix+idx[ii], iy+idy[ii], geom0) ) ii++;
    if( ii >= 4 ) {
      printf("Error in BEmcRec::SectorToGlobal: can not identify neighbour for tower (%d,%d)\n", ix,iy);
      return;
    }
    float Xc = geom0.Xcenter - idx[ii]*geom0.dX[0] - idy[ii]*geom0.dX[1];
    float Yc = geom0.Ycenter - idx[ii]*geom0.dY[0] - idy[ii]*geom0.dY[1];
    float Zc = geom0.Zcenter - idx[ii]*geom0.dZ[0] - idy[ii]*geom0.dZ[1];
    geom0.Xcenter = Xc;
    geom0.Ycenter = Yc;
    geom0.Zcenter = Zc;
  }

  float xt = geom0.Xcenter + (xC - ix) * geom0.dX[0] + (yC - iy) * geom0.dX[1];
  float yt = geom0.Ycenter + (xC - ix) * geom0.dY[0] + (yC - iy) * geom0.dY[1];
  float zt = geom0.Zcenter + (xC - ix) * geom0.dZ[0] + (yC - iy) * geom0.dZ[1];

  CorrectShowerDepth(E, xt, yt, zt, xA, yA, zA);

  //  rA = sqrt(xA*xA+yA*yA);
  //  phiA = atan2(yA, xA);
  
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetModules(vector<EmcModule> const* modules)
{
  *fModules = *modules;

#if 0
  // Make sure that each hit(fired module) knows what sector it belogs to
  vector<EmcModule>::iterator listIter;
  for(listIter=fModules->begin(); listIter!=fModules->end(); listIter++)
    listIter->fOwner=this;
#endif  // #if 0
}

// ///////////////////////////////////////////////////////////////////////////

int BEmcRec::iTowerDist(int ix1, int ix2)
// Distrance in tower units
{
  int idist = ix2 - ix1;
  if (bCYL)
  {
    int idistr = fNx - abs(idist);  // Always >0
    if (idistr < abs(idist))
    {  // Then count in opposite direction
      if (idist < 0)
        idist = idistr;
      else
        idist = -idistr;
    }
  }
  //  printf("Dist %d %d: %d\n",ix1,ix2,idist);
  return idist;
}

float BEmcRec::fTowerDist(float x1, float x2)
{
  float dist = x2 - x1;
  if (bCYL)
  {
    float distr = fNx - fabs(dist);  // Always >0
    if (distr < abs(dist))
    {  // Then count in opposite direction
      if (dist < 0)
        dist = distr;
      else
        dist = -distr;
    }
  }
  return dist;
}

// ///////////////////////////////////

int BEmcRec::FindClusters()
{
  // Cluster search algorithm based on Lednev's one developed for GAMS.
  // Returns -1 if fgMaxLen parameter is too small (just increase it)

  int nhit, nCl;
  int LenCl[fgMaxLen];
  int next, ib, ie, iab, iae, last, LastCl, leng, ich;
  int ia = 0;

  EmcModule* vv;
  EmcModule *vhit, *vt;
  EmcCluster Clt(this);
  vector<EmcModule>::iterator ph;
  vector<EmcModule> hl;

  (*fClusters).clear();
  nhit = (*fModules).size();

  if (nhit <= 0) return 0;
  if (nhit == 1)
  {
    Clt.ReInitialize((*fModules));
    fClusters->push_back(Clt);
    return 1;
  }

  vt = new EmcModule[nhit];
  vhit = new EmcModule[nhit];

  ph = (*fModules).begin();
  vv = vhit;
  while (ph != (*fModules).end()) *vv++ = *ph++;

  qsort(vhit, nhit, sizeof(EmcModule), HitNCompare);

  nCl = 0;
  next = 0;
  for (ich = 1; ich < nhit + 1; ich++)
  {
    if (ich < nhit) ia = vhit[ich].ich;

    // New subcluster
    //
    if ((ia - vhit[ich - 1].ich > 1)  // the beginning of new subcluster
        || (ich >= nhit)              // just finish defining last sub-cluster
        || (ia - ia / fNx * fNx == 0))
    {  // new raw -> new subcluster

      ib = next;
      ie = ich - 1;
      next = ich;
      if (nCl >= fgMaxLen)
      {
        delete[] vhit;
        delete[] vt;
        return -1;
      }
      nCl++;
      LenCl[nCl - 1] = next - ib;
      if (nCl > 1)
      {
        // Job to glue the subclusters with common edge
        //
        iab = vhit[ib].ich;  // The last subcl begin
        iae = vhit[ie].ich;  // The last subcl end
        last = ib - 1;       // The prelast subcl end
        LastCl = nCl - 2;
        for (int iCl = LastCl; iCl >= 0; iCl--)
        {
          leng = LenCl[iCl];

          if (iab - vhit[last].ich > fNx) goto new_ich;
          for (int ichc = last; ichc >= last - leng + 1; ichc--)
          {
            //	    if( iab-vhit[ichc].ich >  fNx ) goto new_icl; // From PHENIX version !!! This may be not right for complicated clusters, where tower ordering is not conserved

            //	    if( iae-vhit[ichc].ich >= fNx // From PHENIX version
            if ((vhit[ichc].ich + fNx <= iae && vhit[ichc].ich + fNx >= iab) || (bCYL && (iae % fNx == fNx - 1) && (iae - vhit[ichc].ich == fNx - 1))  // Only for CYLinder geom !!!!
            )
            {
              // Swap iCl-cluster towers (of length "leng") with whatever was between it and the last subcluster (of length "ib-1-last") - to make it adjecent to the last subcluster
              CopyVector(&vhit[last + 1 - leng], vt, leng);
              CopyVector(&vhit[last + 1], &vhit[last + 1 - leng], ib - 1 - last);
              CopyVector(vt, &vhit[ib - leng], leng);

              // Now the number of clusters is reduced by 1 and the length of the last one increased by iCl-cluster length "leng"
              for (int i = iCl; i < nCl - 2; i++) LenCl[i] = LenCl[i + 1];
              ib -= leng;
              LenCl[nCl - 2] = LenCl[nCl - 1] + leng;
              nCl--;
              goto new_icl;
            }
          }  // for( int ichc=last

        new_icl:
          last = last - leng;
        }  // for( int iCl=LastCl

      }  // if( nCl > 1

    }  // if( (ia-vhit

  new_ich:
    continue;
  }  //  for( ich=1

  if (nCl > 0)
  {
    ib = 0;
    for (int iCl = 0; iCl < nCl; iCl++)
    {
      leng = LenCl[iCl];
      hl.clear();
      for (ich = 0; ich < leng; ich++) hl.push_back(vhit[ib + ich]);
      Clt.ReInitialize(hl);
      ib += LenCl[iCl];
      fClusters->push_back(Clt);
    }
  }
  delete[] vhit;
  delete[] vt;

  return nCl;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::GlobalToSector(float xgl, float ygl, float zgl, float* px,
                             float* py, float* pz)
{
  *px = xgl + fModSizex * (fNx - 1) / 2.;
  *py = ygl + fModSizey * (fNy - 1) / 2.;
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

void BEmcRec::SectorToGlobal(float xC, float yC, float zC,
                             float* px, float* py, float* pz)
{
  *px = *py = *pz;

  int ich;
  std::map<int, TowerGeom>::iterator it;

  int ix = xC + 0.5;  // tower #
  if (ix < 0 || ix >= fNx)
  {
    printf("Error in BEmcRec::SectorToGlobal: wrong input x: %d\n", ix);
    return;
  }

  int iy = yC + 0.5;  // tower #
  if (iy < 0 || iy >= fNy)
  {
    printf("Error in BEmcRec::SectorToGlobal: wrong input y: %d\n", iy);
    return;
  }

  ich = iy * fNx + ix;
  it = fTowerGeom.find(ich);
  if (it == fTowerGeom.end())
  {
    printf("Error in BEmcRec::SectorToGlobal: wrong input (x,y): %f %f\n", xC, yC);
    return;
  }
  TowerGeom geom0 = it->second;

  // Next tower in x
  ich = iy * fNx + (ix + 1);
  it = fTowerGeom.find(ich);
  if (it == fTowerGeom.end())
  {
    ich = iy * fNx + (ix - 1);
    it = fTowerGeom.find(ich);
    if (it == fTowerGeom.end())
    {
      printf("Error in BEmcRec::SectorToGlobal: Error in geometery extraction for x= %f\n", xC);
      return;
    }
  }
  TowerGeom geomx = it->second;

  // Next tower in y
  ich = (iy + 1) * fNx + ix;
  it = fTowerGeom.find(ich);
  if (it == fTowerGeom.end())
  {
    ich = (iy - 1) * fNx + ix;
    it = fTowerGeom.find(ich);
    if (it == fTowerGeom.end())
    {
      printf("Error in BEmcRec::SectorToGlobal: Error in geometery extraction for y= %f\n", yC);
      return;
    }
  }
  TowerGeom geomy = it->second;

  float dx = fabs(geom0.Xcenter - geomx.Xcenter);
  float dy = fabs(geom0.Ycenter - geomy.Ycenter);
  *px = geom0.Xcenter + (xC - ix) * dx;
  *py = geom0.Ycenter + (yC - iy) * dy;
  *pz = geom0.Zcenter;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SectorToGlobalErr(float dxsec, float dysec, float dzsec,
                                float* pdx, float* pdy, float* pdz)
{
  *pdx = 0.;
  *pdy = 0.;
  *pdz = 0.;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::Gamma(int nh, EmcModule* phit0, float* pchi, float* pchi0,
                    //void BEmcRec::Gamma(int nh, EmcModule* phit, float* pchi, float* pchi0,
                    float* pe1, float* px1, float* py1, float* pe2,
                    float* px2, float* py2, int& ndf)
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
  float chi, chi0, chisq0, chisave;
  int dof;
  float x0, y0;
  float stepx, stepy, parx, pary;
  const float dxy = 0.06;
  const float stepmin = 0.01;
  const float zTG = 100;
  const float xmcut = 0.0015;  // (GeV), for overlapping showers separation

  *pe1 = 0;
  *px1 = 0;
  *py1 = 0;
  *pe2 = 0;
  *px2 = 0;
  *py2 = 0;
  if (nh <= 0) return;

  EmcModule* phit = new EmcModule[nh];
  int ish = ShiftX(0, nh, phit0, phit);
  if (ish < -fNx)
  {
    delete[] phit;
    return;
  }

  //  int ish = 0;
  float xx, yy, xy;  // Not used anywhere
  Momenta(nh, phit, &e1, &x1, &y1, &xx, &yy, &xy);
  *pe1 = e1;

  //!!!!! Begin: Exclude chi2 minimizitation and 2-phot splitting
  // Until this is proved working and useful (and ClusterChisq() properly tuned) !!!!!
  /*
  *px1=x1-ish;
  *py1=y1;
  delete [] phit; 
  return;
  */
  //!!!!! End: Exclude chi2 minimizitation and 2-phot splitting

  if (e1 <= 0)
  {
    delete[] phit;
    return;
  }
  //  if( e1 <= 0 ) return;

  SetProfileParameters(0, e1, x1, y1);

  chisave = *pchi;
  chi = *pchi;
  // ClusterChisq parameter list changed MV 28.01.00
  chi0 = ClusterChisq(nh, phit, e1, x1, y1, ndf);

  chisq0 = chi0;
  dof = ndf;  // nh->ndf MV 28.01.00

  // ndf=0 means the cluster's chi2 cannot be found; in this case chi0=0.
  if (dof < 1) dof = 1;
  chi = chisq0 / dof;
  x0 = x1;
  y0 = y1;
  for (;;)
  {
    double chir = ClusterChisq(nh, phit, e1, x0 + dxy, y0, ndf);
    double chil = ClusterChisq(nh, phit, e1, x0 - dxy, y0, ndf);
    double chiu = ClusterChisq(nh, phit, e1, x0, y0 + dxy, ndf);
    double chid = ClusterChisq(nh, phit, e1, x0, y0 - dxy, ndf);

    if ((chi0 > chir) || (chi0 > chil))
    {
      stepx = dxy;
      if (chir > chil) stepx = -stepx;
    }
    else
    {
      stepx = 0;
      parx = chir + chil - 2 * chi0;
      if (parx > 0) stepx = -dxy * (chir - chil) / 2 / parx;
    }

    if ((chi0 > chiu) || (chi0 > chid))
    {
      stepy = dxy;
      if (chiu > chid) stepy = -stepy;
    }
    else
    {
      stepy = 0;
      pary = chiu + chid - 2 * chi0;
      if (pary > 0) stepy = -dxy * (chiu - chid) / 2 / pary;
    }
    if ((EmcCluster::ABS(stepx) < stepmin) && (EmcCluster::ABS(stepy) < stepmin)) break;
    double chi00 = ClusterChisq(nh, phit, e1, x0 + stepx, y0 + stepy, ndf);

    if (chi00 >= chi0) break;
    chi0 = chi00;
    x0 += stepx;
    y0 += stepy;
  }
  if (chi0 < chisq0)
  {
    x1 = x0;
    y1 = y0;
    chi = chi0 / dof;
  }

  *pchi0 = chi;
  *pchi = chi;
  *py1 = y1;
  x1 -= float(ish);
  while (x1 < -0.5) x1 += float(fNx);
  while (x1 >= fNx - 0.5) x1 -= float(fNx);
  *px1 = x1;

  if (e1 <= fgMinShowerEnergy)
  {
    delete[] phit;
    return;
  }
  //  if( e1 <= fgMinShowerEnergy ) return;

  if (chi > chisave)
  {
    TwoGamma(nh, phit, &chi, &e1, &x1, &y1, &e2, &x2, &y2);
    //    printf("Chi=%f E=%f %f X=%f %f Y=%f %f\n",chi,e1,e2,x1,x2,y1,y2);
    if (e2 > 0)
    {
      double d2 = ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / zTG / zTG;
      double xm2 = e1 * e2 * d2;
      if (xm2 > 0) xm2 = sqrt(xm2);
      if (xm2 > xmcut && e1 > fgMinShowerEnergy && e2 > fgMinShowerEnergy)
      {
        x1 -= float(ish);
        while (x1 < -0.5) x1 += float(fNx);
        while (x1 >= fNx - 0.5) x1 -= float(fNx);
        x2 -= float(ish);
        while (x2 < -0.5) x2 += float(fNx);
        while (x2 >= fNx - 0.5) x2 -= float(fNx);

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
  delete[] phit;
}

// ///////////////////////////////////////////////////////////////////////////

EmcModule BEmcRec::ShiftX(int ish, EmcModule& ehit)
{
  EmcModule hh = ehit;
  int iy = hh.ich / fNx;
  int ix = hh.ich % fNx + ish;
  while (ix < 0) ix += fNx;
  while (ix >= fNx) ix -= fNx;
  int ich = iy * fNx + ix;
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
  int ishift_def = fNx / 2;
  int ix, ish;

  if (ishift == 0)
  {
    bool bshift = false;
    for (int i = 0; i < nh; i++)
    {
      phit1[i] = phit0[i];
      ix = phit0[i].ich % fNx;
      if (ix == 0) bshift = true;
    }
    if (!bshift) return 0;
  }

  if (ishift != 0)
    ish = ishift;
  else
    ish = ishift_def;

  int ixmin = 999;
  int ixmax = -999;
  for (int i = 0; i < nh; i++)
  {
    phit1[i] = phit0[i];
    int iy = phit0[i].ich / fNx;
    ix = phit0[i].ich % fNx + ish;

    while (ix < 0) ix += fNx;
    while (ix >= fNx) ix -= fNx;

    if (ixmin > ix) ixmin = ix;
    if (ixmax < ix) ixmax = ix;

    int ich = iy * fNx + ix;
    phit1[i].ich = ich;
  }

  if (ishift == 0 && ixmax - ixmin > fNx / 2)
  {
    printf("!!! Error BEmcRec::ShiftX(): Too long cluster (>%d towers in phi): reconstruction may be wrong. May need tower energy threshold increase for clustering.\n", fNx / 2);
    return -999;
  }

  return ish;
}

void BEmcRec::Momenta(int nh, EmcModule* phit, float* pe, float* px,
                      float* py, float* pxx, float* pyy, float* pyx)
{
  // First and second momenta calculation

  float a, x, y, e, xx, yy, yx;
  EmcModule* p;

  *pe = 0;
  *px = 0;
  *py = 0;
  *pxx = 0;
  *pyy = 0;
  *pyx = 0;
  if (nh <= 0) return;

  // Find max energy tower
  //
  p = phit;
  float emax = 0;
  int ichmax = 0;
  for (int i = 0; i < nh; i++)
  {
    a = p->amp;
    if (a > emax)
    {
      emax = a;
      ichmax = p->ich;
    }
    p++;
  }
  if (emax <= 0) return;

  int iymax = ichmax / fNx;
  int ixmax = ichmax - iymax * fNx;

  // Calculate CG relative to max energy tower

  x = 0;
  y = 0;
  e = 0;
  xx = 0;
  yy = 0;
  yx = 0;
  p = phit;
  for (int i = 0; i < nh; i++)
  {
    a = p->amp;
    int iy = p->ich / fNx;
    int ix = p->ich - iy * fNx;
    int idx = iTowerDist(ixmax, ix);
    int idy = iy - iymax;
    e += a;
    x += idx * a;
    y += idy * a;
    xx += a * idx * idx;
    yy += a * idy * idy;
    yx += a * idx * idy;
    p++;
  }

  *pe = e;

  if (e > 0)
  {
    x /= e;
    y /= e;
    xx = xx / e - x * x;
    yy = yy / e - y * y;
    yx = yx / e - y * x;

    x += ixmax;
    y += iymax;

    while (x < -0.5) x += float(fNx);
    while (x >= fNx - 0.5) x -= float(fNx);

    *px = x;
    *py = y;
    *pxx = xx;
    *pyy = yy;
    *pyx = yx;
  }
}

// ///////////////////////////////////////////////////////////////////////////
// Static functions

int BEmcRec::HitNCompare(const void* h1, const void* h2)
{
  return (static_cast<const EmcModule*>(h1)->ich - static_cast<const EmcModule*>(h2)->ich);
}

// ///////////////////////////////////////////////////////////////////////////

int BEmcRec::HitACompare(const void* h1, const void* h2)
{
  float amp1 = static_cast<const EmcModule*>(h1)->amp;
  float amp2 = static_cast<const EmcModule*>(h2)->amp;
  return (amp1 < amp2) ? 1 : (amp1 > amp2) ? -1 : 0;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::ZeroVector(int* v, int N)
{
  int* p = v;
  for (int i = 0; i < N; i++) *p++ = 0;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::ZeroVector(float* v, int N)
{
  float* p = v;
  for (int i = 0; i < N; i++) *p++ = 0;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::ZeroVector(EmcModule* v, int N)
{
  //  memset(v, 0, N*sizeof(EmcModule));
  for (int i = 0; i < N; i++)
  {
    v[i].ich = 0;
    v[i].amp = 0;
    v[i].tof = 0;
  }
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CopyVector(int* from, int* to, int N)
{
  if (N <= 0) return;
  for (int i = 0; i < N; i++) to[i] = from[i];
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CopyVector(EmcModule* from, EmcModule* to, int N)
{
  if (N <= 0) return;
  for (int i = 0; i < N; i++) to[i] = from[i];
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::c3to5(float e0, float x0, float y0, float eps,
                    float dx, float dy,
                    float* pe1, float* px1, float* py1,
                    float* pe2, float* px2, float* py2)
{
  // 3 to 5-dim space conversion subroutine.
  // eps=(e1-e2)/e0,  (e0=e1+e2), x0*e0=x1*e1+x2*e2, dx=x1-x2

  *pe1 = e0 * (1 + eps) / 2;
  *pe2 = e0 - *pe1;
  *px1 = x0 + dx * (1 - eps) / 2;
  *py1 = y0 + dy * (1 - eps) / 2;
  *px2 = x0 - dx * (1 + eps) / 2;
  *py2 = y0 - dy * (1 + eps) / 2;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetChi2Limit(int limit)
{
  // Sets the limit for PeakArea splitting onto 2 EMShowers:
  // limit=0 -> For the Conf Level: 0%
  // limit=1 -> For the Conf Level: 1% for GEANT, 2%-5% for TestBeam
  // limit=2 -> For the Conf Level: 2% for GEANT, 4%-7% for TestBeam

  int i;

  switch (limit)
  {
  case 0:
    for (i = 0; i < 50; i++) fgChi2Level[i] = 9999.;
    break;
  case 1:
    for (i = 0; i < 50; i++) fgChi2Level[i] = fgChi2Level1[i];
    break;
  case 2:
    for (i = 0; i < 50; i++) fgChi2Level[i] = fgChi2Level2[i];
    break;
  default:
    for (i = 0; i < 50; i++) fgChi2Level[i] = fgChi2Level1[i];
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

/////////////////////////////////////////////////////////////////
/*
void BEmcRec::CorrectPosition(float Energy, float x, float y,
				     float* pxc, float* pyc)
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
*/
// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CalculateErrors(float e, float x, float y, float* pde,
                              float* pdx, float* pdy, float* pdz)
{
  // Returns the errors for the reconstructed energy and position
  // (in the hypothesis of EM shower)
  // Should be called only just after CorrectPosition !!!

  float de, dy, dz, dxg, dyg, dzg;
  static float ae = 0.076, be = 0.022;        // de/E = a/sqrt(E)&b
  static float a = 0.57, b = 0.155, d = 1.6;  // dx = a/sqrt(E)+b (cm)
  static float dx = 0.1;                      // (cm)

  de = sqrt(ae * ae * e + be * be * e * e);
  dz = a / sqrt(e) + b;
  dy = dz;
  dz = sqrt(dz * dz + d * d * fSinTx * fSinTx);
  dy = sqrt(dy * dy + d * d * fSinTy * fSinTy);

  SectorToGlobalErr(dx, dy, dz, &dxg, &dyg, &dzg);

  *pde = de;
  *pdx = dxg;
  *pdy = dyg;
  *pdz = dzg;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetProfileParameters(int sec, float Energy, float x,
                                   float y)
{
  // Axis Z here means X in two dim Sector coord system !!!
  // So one has to supply (x,y) coordinates in cell units (x instead of z)
  // If sec < 0 this routine changes only Energy dependent parameters -
  // the angle dependent ones are set in the previous call

  static float sin2ax, sin2ay, lgE;
  //  float vx, vy, vz;
  //  float xVert, yVert, zVert;
  int sign;

  if (sec >= 0)
  {
    // Non-ortogonality in phi in a sector (of 8 rows)
    int ix = x + 0.5;
    float dx = x + 0.5 - ix;  // from 0 to 1
    int ix8 = ix % 8;
    dx += ix8;
    fSinTx = (dx - 4) * fModSizex;
    fSinTx = 0;
    fSinTy = 0;
    static float sin2a = fSinTx * fSinTx + fSinTy * fSinTy;
    fSin4T = sin2a * sin2a;
    sin2ax = fSinTx * fSinTx;
    sin2ay = fSinTy * fSinTy;
  }

  if (Energy <= 1.e-10)
    lgE = 0;
  else
    lgE = log(Energy);
  /*  
  fPpar1=0.59-(1.45+0.13*lgE)*sin2a;
  fPpar2=0.265+(0.80+0.32*lgE)*sin2a;
  fPpar3=0.25+(0.45-0.036*lgE)*sin2a;
  fPpar4=0.42;
  */
  fPpar1 = 0.549;
  fPpar2 = 0.304;
  fPpar3 = 0.389;
  fPpar4 = 0.326;

  if (fSinTx > 0)
    sign = 1;
  else
    sign = -1;
  fPshiftx = (1.05 + 0.12 * lgE) * sin2ax * sign;
  fPshiftx = 0;  // !!!!! Untill tuned ... may not be necessary

  if (fSinTy > 0)
    sign = 1;
  else
    sign = -1;
  fPshifty = (1.05 + 0.12 * lgE) * sin2ay * sign;
  fPshifty = 0;  // !!!!! Untill tuned ... may not be necessary
}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::PredictEnergy(float xc, float yc, float en)
{
  // Calculates the energy deposited in the tower, the distance between
  // its center and shower Center of Gravity being (xc,yc)
  // en - shower energy
  // If en<0 -> no Shower Profile parameters change is needed

  float dx, dy, r1, r2, r3;

  if (en > 0) SetProfileParameters(-1, en, xc, yc);
  dx = fabs(xc - fPshiftx);
  dy = EmcCluster::ABS(yc - fPshifty);
  r2 = dx * dx + dy * dy;
  r1 = sqrt(r2);
  r3 = r2 * r1;
  double e = fPpar1 * exp(-r3 / fPpar2) + fPpar3 * exp(-r1 / fPpar4);

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
  float step, cosi;
  float e1c, x1c, y1c, e2c, x2c, y2c;
  float eps0 = 0.0;
  float ex;
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

  const float epsmax = 0.9999;
  const float stpmin = 0.025;
  const float delch = 2;

  Momenta(nh, phit, &e0, &x0, &y0, &xx, &yy, &yx);
  *pe2 = 0;
  *px2 = 0;
  *py2 = 0;
  if (nh <= 0) return;
  //  choosing of the starting point
  dxy = xx - yy;
  rsg2 = dxy * dxy + 4 * yx * yx;
  if (rsg2 < 1e-20) rsg2 = 1e-20;
  rsq = sqrt(rsg2);
  dxc = 0;
  if (rsq + dxy > 0) dxc = -sqrt((rsq + dxy) * 2);  // To avoid the possibility of small negative number due to precision
  dyc = 0;
  if (rsq - dxy > 0) dyc = sqrt((rsq - dxy) * 2);  // To avoid the possibility of small negative number due to precision
  if (yx >= 0) dyc = -dyc;
  r = sqrt(dxc * dxc + dyc * dyc);
  epsc = 0;
  for (in = 0; in < nh; in++)
  {
    ixy = phit[in].ich;
    iy = ixy / fNx;
    ix = ixy - iy * fNx;
    double u = (ix - x0) * dxc / r + (iy - y0) * dyc / r;
    epsc -= phit[in].amp * u * EmcCluster::ABS(u);
  }
  epsc /= (e0 * rsq);
  if (epsc > 0.8) epsc = 0.8;
  if (epsc < -0.8) epsc = -0.8;
  dxc /= sqrt(1 - epsc * epsc);
  dyc /= sqrt(1 - epsc * epsc);
  //  Start of iterations
  step = 0.1;
  cosi = 0;
  double chisq2 = 1.e35;
  for (iter = 0; iter < 100; iter++)
  {
    c3to5(e0, x0, y0, epsc, dxc, dyc, &e1c, &x1c, &y1c, &e2c, &x2c, &y2c);
    double eps1 = (1 + epsc) / 2;
    double eps2 = (1 - epsc) / 2;
    double chisqc = 0;
    for (in = 0; in < nh; in++)
    {
      ex = phit[in].amp;
      ixy = phit[in].ich;
      iy = ixy / fNx;
      ix = ixy - iy * fNx;
      dx1 = x1c - ix;
      dy1 = y1c - iy;
      dx2 = x2c - ix;
      dy2 = y2c - iy;
      a0 = e1c * PredictEnergy(dx1, dy1, e1c) + e2c * PredictEnergy(dx2, dy2, e2c);
      d = fgEpar00 * fgEpar00 + e0 * (fgEpar1 * a0 / e0 + fgEpar2 * a0 * a0 / e0 / e0 + fgEpar3 * a0 * a0 * a0 / e0 / e0 / e0) + e0 * sqrt(e0) * fgEpar4 * a0 / e0 * (1 - a0 / e0) * fSin4T + e0 * e0 * fgEpar0 * fgEpar0;
      chisqc += (a0 - ex) * (a0 - ex) / d;
    }
    if (chisqc >= chisq2)
    {
      if (iter > 0)
      {
        dchi = chisqc - chisq2;
        dchi0 = gr * step;
        step /= (2 * sqrt(1 + dchi / dchi0));
      }
      step /= 2;
    }
    else
    {
      // Calculation of gradient
      grec = 0;
      grxc = 0;
      gryc = 0;
      for (in = 0; in < nh; in++)
      {
        ex = phit[in].amp;
        ixy = phit[in].ich;
        iy = ixy / fNx;
        ix = ixy - iy * fNx;
        dx1 = x1c - ix;
        dy1 = y1c - iy;
        dx2 = x2c - ix;
        dy2 = y2c - iy;
        a1 = e1c * PredictEnergy(dx1, dy1, e1c);
        a2 = e2c * PredictEnergy(dx2, dy2, e2c);
        a0 = a1 + a2;
        d = fgEpar00 * fgEpar00 + e0 * (fgEpar1 * a0 / e0 + fgEpar2 * a0 * a0 / e0 / e0 + fgEpar3 * a0 * a0 * a0 / e0 / e0 / e0) + e0 * sqrt(e0) * fgEpar4 * a0 / e0 * (1 - a0 / e0) * fSin4T + e0 * e0 * fgEpar0 * fgEpar0;
        dd = (a0 - ex) / d;
        dchida = dd * (2 - dd * (fgEpar1 + 2 * fgEpar2 * a0 / e0 + 3 * fgEpar3 * a0 * a0 / e0 / e0 + e0 * sqrt(e0) * fgEpar4 * fSin4T * (1 - 2 * a0 / e0) + 2 * fgEpar0 * fgEpar0 * a0));
        gx1 = (e1c * PredictEnergy(x1c + 0.05 - ix, dy1, e1c) - a1) * 20;
        gx2 = (e2c * PredictEnergy(x2c + 0.05 - ix, dy2, e2c) - a2) * 20;
        gy1 = (e1c * PredictEnergy(dx1, y1c + 0.05 - iy, e1c) - a1) * 20;
        gy2 = (e2c * PredictEnergy(dx2, y2c + 0.05 - iy, e2c) - a2) * 20;
        grec += (dchida * ((a1 / e1c - a2 / e2c) * e0 - (gx1 + gx2) * dxc - (gy1 + gy2) * dyc) / 2);
        grxc += (dchida * (gx1 * eps2 - gx2 * eps1));
        gryc += (dchida * (gy1 * eps2 - gy2 * eps1));
      }
      grc = sqrt(grec * grec + grxc * grxc + gryc * gryc);
      if (grc < 1e-10) grc = 1e-10;
      if (iter > 0)
      {
        cosi = (gre * grec + grx * grxc + gry * gryc) / (gr * grc);
        scal = EmcCluster::ABS(gr / grc - cosi);
        if (scal < 0.1) scal = 0.1;
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
    epsc = eps0 - step * gre / gr;
    while (EmcCluster::ABS(epsc) >= epsmax)
    {
      step /= 2;
      epsc = eps0 - step * gre / gr;
    }
    dxc = dx0 - step * grx / gr;
    dyc = dy0 - step * gry / gr;
    if (step * gr < stpmin) break;
  }
  if ((*pchi) * nh - chisq2 < delch) return;
  dof = nh;
  if (dof < 1) dof = 1;
  *pchi = chisq2 / dof;
  c3to5(e0, x0, y0, eps0, dx0, dy0, pe1, px1, py1, pe2, px2, py2);
}

// ///////////////////////////////////////////////////////////////////////////

int BEmcRec::GetTowerID(int iy, int iz, int nn, int* iyy, int* izz, float* ee)
{
  while (iz < 0) iz += fNx;
  while (iz >= fNx) iz -= fNx;

  for (int i = 0; i < nn; i++)
  {
    if (iy == iyy[i] && iz == izz[i]) return i;
  }
  return -1;
}

float BEmcRec::GetProb(vector<EmcModule> HitList, float et, float xg, float yg, float zg, float& chi2, int& ndf)
// Do nothing; should be defined in a detector specific module BEmcRec<Name>
{
  chi2 = 0;
  ndf = 0;
  return -1;
}

float BEmcRec::ClusterChisq(int nh, EmcModule* phit, float e, float x,
                            float y, int& ndf)
{
  double chi = 0;

  for (int in = 0; in < nh; in++)
  {
    int ixy = phit[in].ich;
    int iy = ixy / fNx;
    int ix = ixy - iy * fNx;
    double et = phit[in].amp;
    double a = PredictEnergy(x - ix, y - iy, -1);
    double d = fgEpar00 * fgEpar00 + e * (fgEpar1 * a + fgEpar2 * a * a + fgEpar3 * a * a * a) +
               e * sqrt(e) * fgEpar4 * a * (1 - a) * fSin4T + e * e * fgEpar0 * fgEpar0;
    a *= e;
    chi += (et - a) * (et - a) / d;
  }

  ndf = nh;  // change needed for PbGl MV 28.01.00
  return chi;
}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::Chi2Limit(int ND)
{
  //  Here the reverse Chi2Correct function is used

  float rn, a, b, chi2;

  if (ND < 1) return 9999.;  // Should we put 0. here?

  chi2 = fgChi2Level[EmcCluster::min(ND, 50) - 1];
  if (chi2 > 100) return 9999.;  // Why should chi2 ever be >100?

  rn = ND;
  b = 0.072 * sqrt(rn);
  a = 6.21 / (sqrt(rn) + 4.7);

  return chi2 * a / (1. - chi2 * b);
}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::Chi2Correct(float Chi2, int ND)
{
  // Chi2 - is reduced Chi2: Chi2/ND !!
  // MV 02.22.2000: Actually the above is not true. The supplied value of Chi2
  // has been already divided by ND. So Chi2 here is only corrected.

  float rn, a, b, c;

  if (ND < 1) return 9999.;  // Should we put 0. here?

  rn = ND;
  b = 0.072 * sqrt(rn);
  a = 6.21 / (sqrt(rn) + 4.7);
  c = a + b * Chi2;
  if (c < 1) c = 1;

  return Chi2 / c;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::SetTowerThreshold(float Thresh)
{
  fgTowerThresh = Thresh;
  fgEpar0 = Thresh * 0.07;
  fgEpar00 = EmcCluster::max((double) Thresh / 3, 0.005);
}

// **********************************************************************

void BEmcRec::getTowerPos(int ix, int iy, float& x, float& y)
{
  x = 2.859 + 5.562 * ix + int(ix / 12) * 0.256;
  y = 2.859 + 5.562 * iy + int(iy / 12) * 0.156;
}

// **********************************************************************

/// Converts coordinates in units of towers into cm's (Local coord. system)
void BEmcRec::TowersToSector(float xT, float yT, float& xS, float& yS)
{
  int x = int(xT);
  float dx = xT - x;
  int y = int(yT);
  float dy = yT - y;
  xS = fModSizex * (x + 0.5) + int(xT / 12) * 0.256 + fModSizex * dx;
  yS = fModSizey * (y + 0.5) + int(yT / 12) * 0.156 + fModSizey * dy;
}

// **********************************************************************
/// Returns  coordinates of the tower centers in cm's (Local coord. system)
void BEmcRec::TowersToSector(int xT, int yT, float& xS, float& yS)
{
  xS = fModSizex * (xT + 0.5) + int(xT / 12) * 0.256;
  yS = fModSizey * (yT + 0.5) + int(yT / 12) * 0.156;
}

// **********************************************************************
/// Converts Local Sector coordinates in cm into integer tower numbers
void BEmcRec::SectorToTowers(float xS, float yS, int& xT, int& yT)
{
  // PbSc
  xT = int(((xS - int(xS / 67.0) * 67.0) - 0.078) / fModSizex) + 12 * int(xS / 67.0);
  yT = int(((yS - int(yS / 66.9) * 66.9) - 0.078) / fModSizey) + 12 * int(yS / 66.9);
}

// ///////////////////////////////////////////////////////////////////////////

/* Future improvements:

1. FindClusters(): to ensure that all EmcModules are above energy threshold 
set by SetThreshold routine (or default one)

*/

// ///////////////////////////////////////////////////////////////////////////
// EOF
