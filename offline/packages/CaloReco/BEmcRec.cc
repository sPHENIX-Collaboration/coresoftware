// Name: BEmcRec.h
// Author: A. Bazilevsky, Apr 2012
// Modified from EmcSectorRec.cxx and EmcScSectorRec.cxx

// BEmcRec -- base class for sPHENIX EMCal

// ///////////////////////////////////////////////////////////////////////////

#include "BEmcRec.h"
#include "BEmcCluster.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>

using namespace std;

// Define and initialize static members

// Max number of clusters, used in FindClusters(), automatically extended when needed
int const BEmcRec::fgMaxLen = 1000;

// ///////////////////////////////////////////////////////////////////////////
// BEmcRec member functions

BEmcRec::BEmcRec()
  : bCYL(true)
  , fNx(-1)
  , fNy(-1)
  , fVx(0)
  , fVy(0)
  , fVz(0)
  , fgTowerThresh(0.01)
  , fgMinPeakEnergy(0.08)
  , m_ThisName("NOTSET")
//  , _emcprof(nullptr)
{
  fTowerGeom.clear();
  fModules = new vector<EmcModule>;
  fClusters = new vector<EmcCluster>;
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

void BEmcRec::LoadProfile(const string& fname)
{
  cout << "Warning from BEmcRec::LoadProfile(): No acton defined for shower profile evaluation; should be defined in a detector specific module " << Name() << endl;
}

void BEmcRec::PrintTowerGeometry(const string& fname)
{
  ofstream outfile(fname);
  if (!outfile.is_open())
  {
    cout << "Error in BEmcRec::PrintTowerGeometry(): Failed to open file "
         << fname << endl;
    return;
  }
  outfile << "Number of bins:" << endl;
  outfile << fNx << " " << fNy << endl;
  outfile << "ix iy x y z dx0 dy0 dz0 dx1 dy1 dz1" << endl;
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
        outfile << ix << " " << iy << " " << geom.Xcenter << " "
                << geom.Ycenter << " " << geom.Zcenter << " " << geom.dX[0] << " "
                << geom.dY[0] << " " << geom.dZ[0] << " " << geom.dX[1] << " "
                << geom.dY[1] << " " << geom.dZ[1] << endl;
      }
    }
  }
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
  geom.dX[0] = geom.dX[1] = 0;  // These should be calculated by CompleteTowerGeometry()
  geom.dY[0] = geom.dY[1] = 0;
  geom.dZ[0] = geom.dZ[1] = 0;

  int ich = iy * fNx + ix;
  fTowerGeom[ich] = geom;
  return true;
}

bool BEmcRec::CompleteTowerGeometry()
// Calculates tower front size from coordinates of tower center coordinates
{
  if (fTowerGeom.empty() || fNx <= 0)
  {
    cout << "Error in BEmcRec::CalculateTowerSize(): Tower geometry not well setup (NX = "
         << fNx << ")" << endl;
    return false;
  }

  std::map<int, TowerGeom>::iterator it;

  for (it = fTowerGeom.begin(); it != fTowerGeom.end(); it++)
  {
    int ich = it->first;
    TowerGeom geom0 = it->second;
    int ix = ich % fNx;
    int iy = ich / fNx;

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
      cout << "Error in BEmcRec::CompleteTowerGeometry(): Error when locating neighbour for (ix,iy)=("
           << ix << "," << iy << ")" << endl;
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
      cout << "Error in BEmcRec::CompleteTowerGeometry(): Error when locating neighbour for (ix,iy)=("
           << ix << "," << iy << ")" << endl;
      return false;
    }

    geom0.dX[0] = (geomx.Xcenter - geom0.Xcenter) / float(idx);
    geom0.dY[0] = (geomx.Ycenter - geom0.Ycenter) / float(idx);
    geom0.dZ[0] = (geomx.Zcenter - geom0.Zcenter) / float(idx);
    geom0.dX[1] = (geomy.Xcenter - geom0.Xcenter) / float(idy);
    geom0.dY[1] = (geomy.Ycenter - geom0.Ycenter) / float(idy);
    geom0.dZ[1] = (geomy.Zcenter - geom0.Zcenter) / float(idy);

    it->second = geom0;

  }  // it = fTowerGeom.begin()

  return true;
}

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
    cout << "Error in BEmcRec::Tower2Global: wrong input x: " << ix << endl;
    return;
  }

  int iy = yC + 0.5;  // tower #
  if (iy < 0 || iy >= fNy)
  {
    cout << "Error in BEmcRec::Tower2Global: wrong input y: " << iy << endl;
    return;
  }

  // Get tower where the shower is positioned
  TowerGeom geom0;

  if (!GetTowerGeometry(ix, iy, geom0))
  {
    // Weird case: cluster center of gravity outside the EMCal, take geometry from the neighbouring tower
    int idx[4] = {1, 0, -1, 0};
    int idy[4] = {0, 1, 0, -1};
    int ii = 0;
    while (ii < 4 && !GetTowerGeometry(ix + idx[ii], iy + idy[ii], geom0)) ii++;
    if (ii >= 4)
    {
      cout << "Error in BEmcRec::Tower2Global: can not identify neighbour for tower ("
           << ix << "," << iy << ")" << endl;
      return;
    }
    float Xc = geom0.Xcenter - idx[ii] * geom0.dX[0] - idy[ii] * geom0.dX[1];
    float Yc = geom0.Ycenter - idx[ii] * geom0.dY[0] - idy[ii] * geom0.dY[1];
    float Zc = geom0.Zcenter - idx[ii] * geom0.dZ[0] - idy[ii] * geom0.dZ[1];
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
  //  cout << "Dist " << ix1 << " " << ix2 << ": " << idist << endl;
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
// Cluster search algorithm based on Lednev's one developed for GAMS.
// Returns number of clusters found
{
  int nhit, nCl;
  //  int LenCl[fgMaxLen];
  int* LenCl;
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

  int MaxLen = fgMaxLen;
  LenCl = new int[MaxLen];
  ZeroVector(LenCl, MaxLen);

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
      if (nCl >= MaxLen)
      {
        //        delete[] vhit;
        //        delete[] vt;
        //        return -1;
        int* LenCltmp = new int[MaxLen];
        CopyVector(LenCl, LenCltmp, MaxLen);
        delete[] LenCl;
        LenCl = new int[MaxLen * 2];
        ZeroVector(LenCl, MaxLen * 2);
        CopyVector(LenCltmp, LenCl, MaxLen);
        delete[] LenCltmp;
        MaxLen *= 2;
        //	cout << "Extend array size to " << MaxLen << endl;
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
  delete[] LenCl;
  delete[] vhit;
  delete[] vt;

  return nCl;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::Momenta(vector<EmcModule>* phit, float& pe, float& px,
                      float& py, float& pxx, float& pyy, float& pyx)
{
  // First and second momenta calculation

  float a, x, y, e, xx, yy, yx;
  vector<EmcModule>::iterator ph;

  pe = 0;
  px = 0;
  py = 0;
  pxx = 0;
  pyy = 0;
  pyx = 0;
  if (phit->empty()) return;

  // Find max energy tower
  //
  ph = phit->begin();
  float emax = 0;
  int ichmax = 0;
  while (ph != phit->end())
  {
    a = ph->amp;
    if (a > emax)
    {
      emax = a;
      ichmax = ph->ich;
    }
    ++ph;
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
  ph = phit->begin();
  while (ph != phit->end())
  {
    a = ph->amp;
    int iy = ph->ich / fNx;
    int ix = ph->ich - iy * fNx;
    int idx = iTowerDist(ixmax, ix);
    int idy = iy - iymax;
    e += a;
    x += idx * a;
    y += idy * a;
    xx += a * idx * idx;
    yy += a * idy * idy;
    yx += a * idx * idy;
    ph++;
  }
  pe = e;

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

    px = x;
    py = y;
    pxx = xx;
    pyy = yy;
    pyx = yx;
  }
}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::PredictEnergy(float xc, float yc, float en)
{
  // Calculates the energy deposited in the tower, the distance between
  // its center and shower Center of Gravity being (xc,yc)
  // en - shower energy

  float dx, dy, r1, r2, r3;
  float fPpar1, fPpar2, fPpar3, fPpar4;

  float fPshiftx = 0;  // !!!!! Untill tuned ... may not be necessary
  float fPshifty = 0;  // !!!!! Untill tuned ... may not be necessary

  /*
  float lgE;
  if (en <= 1.e-10)
    lgE = 0;
  else
    lgE = log(en);
  fPpar1=0.59-(1.45+0.13*lgE)*sin2a;
  fPpar2=0.265+(0.80+0.32*lgE)*sin2a;
  fPpar3=0.25+(0.45-0.036*lgE)*sin2a;
  fPpar4=0.42;
  */
  fPpar1 = 0.549;
  fPpar2 = 0.304;
  fPpar3 = 0.389;
  fPpar4 = 0.326;
  /*
  fPpar1 = 0.486;
  fPpar2 = 0.302;
  fPpar3 = 0.354;
  fPpar4 = 0.407;
  */
  /*
  fPpar1 = 0.343;
  fPpar2 = 0.509;
  fPpar3 = 0.199;
  fPpar4 = 0.548;
  */

  //  if (en > 0) SetProfileParameters(-1, en, xc, yc);

  dx = fabs(xc - fPshiftx);
  dy = fabs(yc - fPshifty);
  r2 = dx * dx + dy * dy;
  r1 = sqrt(r2);
  r3 = r2 * r1;
  double e = fPpar1 * exp(-r3 / fPpar2) + fPpar3 * exp(-r1 / fPpar4);

  return e;
}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::GetProb(vector<EmcModule> HitList, float et, float xg, float yg, float zg, float& chi2, int& ndf)
// Do nothing; should be defined in a detector specific module BEmcRec{Name}
{
  chi2 = 0;
  ndf = 0;
  return -1;
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

/* Future improvements:

1. FindClusters(): to ensure that all EmcModules are above energy threshold 
set by SetThreshold routine (or default one)

*/

// ///////////////////////////////////////////////////////////////////////////
// EOF
