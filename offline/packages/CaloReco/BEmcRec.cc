// Name: BEmcRec.h
// Author: A. Bazilevsky, Apr 2012
// Modified from EmcSectorRec.cxx and EmcScSectorRec.cxx

// BEmcRec -- base class for sPHENIX EMCal

// ///////////////////////////////////////////////////////////////////////////

#include "BEmcRec.h"
#include "BEmcCluster.h"
#include "BEmcProfile.h"

#include <TMath.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>

// ///////////////////////////////////////////////////////////////////////////
// BEmcRec member functions

BEmcRec::BEmcRec()
{
  fTowerGeom.clear();
  fModules = new std::vector<EmcModule>;
  fClusters = new std::vector<EmcCluster>;
}

// ///////////////////////////////////////////////////////////////////////////

BEmcRec::~BEmcRec()
{
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

  delete _emcprof;
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::LoadProfile(const std::string& /*fname*/)
{
  std::cout << "Warning from BEmcRec::LoadProfile(): No acton defined for shower profile evaluation; should be defined in a detector specific module " << Name() << std::endl;
}

void BEmcRec::PrintTowerGeometry(const std::string& fname)
{
  std::ofstream outfile(fname);
  if (!outfile.is_open())
  {
    std::cout << "Error in BEmcRec::PrintTowerGeometry(): Failed to open file "
              << fname << std::endl;
    return;
  }
  outfile << "Number of bins:" << std::endl;
  outfile << fNx << " " << fNy << std::endl;
  outfile << "ix iy x y z dx0 dy0 dz0 dx1 dy1 dz1" << std::endl;
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
                << geom.dY[1] << " " << geom.dZ[1] << std::endl;
        //	std::cout << "Z0: " << geom.dZ[0] << " || Z1: " << geom.dZ[1] << std::endl;
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
    std::cout << "Error in BEmcRec::CalculateTowerSize(): Tower geometry not well setup (NX = "
              << fNx << ")" << std::endl;
    return false;
  }

  const int nb = 8;
  int idx[nb] = {0, 1, 0, -1, -1, 1, 1, -1};
  int idy[nb] = {-1, 0, 1, 0, -1, -1, 1, 1};

  std::map<int, TowerGeom>::iterator it;

  for (it = fTowerGeom.begin(); it != fTowerGeom.end(); ++it)
  {
    int ich = it->first;
    TowerGeom geom0 = it->second;
    int ix = ich % fNx;
    int iy = ich / fNx;

    TowerGeom geomx;
    int inx = 0;

    while (inx < nb && (idx[inx] == 0 || !GetTowerGeometry(ix + idx[inx], iy + idy[inx], geomx))) inx++;
    if (inx >= nb)
    {
      std::cout << "Error in BEmcRec::CompleteTowerGeometry(): Error when locating neighbour for (ix,iy)=("
                << ix << "," << iy << ")" << std::endl;
      return false;
    }

    TowerGeom geomy;
    int iny = 0;

    while (iny < nb && (idy[iny] == 0 || !GetTowerGeometry(ix + idx[iny], iy + idy[iny], geomy))) iny++;
    if (iny >= nb)
    {
      std::cout << "Error in BEmcRec::CompleteTowerGeometry(): Error when locating neighbour for (ix,iy)=("
                << ix << "," << iy << ")" << std::endl;
      return false;
    }

    geom0.dX[0] = (geomx.Xcenter - geom0.Xcenter) / float(idx[inx]);
    geom0.dY[0] = (geomx.Ycenter - geom0.Ycenter) / float(idx[inx]);
    geom0.dZ[0] = (geomx.Zcenter - geom0.Zcenter) / float(idx[inx]);
    geom0.dX[1] = (geomy.Xcenter - geom0.Xcenter) / float(idy[iny]);
    geom0.dY[1] = (geomy.Ycenter - geom0.Ycenter) / float(idy[iny]);
    geom0.dZ[1] = (geomy.Zcenter - geom0.Zcenter) / float(idy[iny]);

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
    std::cout << m_ThisName << " Error in BEmcRec::Tower2Global: wrong input x: " << ix << std::endl;
    return;
  }

  int iy = yC + 0.5;  // tower #
  if (iy < 0 || iy >= fNy)
  {
    std::cout << "Error in BEmcRec::Tower2Global: wrong input y: " << iy << std::endl;
    return;
  }

  // Get tower where the shower is positioned
  TowerGeom geom0;

  if (!GetTowerGeometry(ix, iy, geom0))
  {
    // Weird case: cluster center of gravity outside the EMCal, take geometry from the neighbouring tower
    const int idx[4] = {1, 0, -1, 0};
    const int idy[4] = {0, 1, 0, -1};
    int ii = 0;
    while (ii < 4 && !GetTowerGeometry(ix + idx[ii], iy + idy[ii], geom0)) ii++;
    if (ii >= 4)
    {
      std::cout << "Error in BEmcRec::Tower2Global: can not identify neighbour for tower ("
                << ix << "," << iy << ")" << std::endl;
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
  //  std::cout << "Dist " << ix1 << " " << ix2 << ": " << idist << std::endl;
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
  std::vector<EmcModule>::iterator ph;
  std::vector<EmcModule> hl;

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
        //	std::cout << "Extend array size to " << MaxLen << std::endl;
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

void BEmcRec::Momenta(std::vector<EmcModule>* phit, float& pe, float& px,
                      float& py, float& pxx, float& pyy, float& pyx,
                      float thresh)
{
  // First and second momenta calculation

  float a, x, y, e, xx, yy, yx;
  std::vector<EmcModule>::iterator ph;

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
    if (a > thresh)
    {
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
    }
    ++ph;
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

float BEmcRec::PredictEnergy(float en, float xcg, float ycg, int ix, int iy)
{
  if (_emcprof != nullptr && bProfileProb) return PredictEnergyProb(en, xcg, ycg, ix, iy);

  float dx = fabs(fTowerDist(float(ix), xcg));
  float dy = ycg - iy;
  return PredictEnergyParam(en, dx, dy);
}

float BEmcRec::PredictEnergyParam(float /*en*/, float xc, float yc)
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

float BEmcRec::PredictEnergyProb(float en, float xcg, float ycg, int ix, int iy)
// Predict tower energy from profiles used in GetProb()
// This is expected to be used in BEmcCluster::GetSubClusters
{
  if (_emcprof == nullptr) return -1;

  while (xcg < -0.5) xcg += float(fNx);
  while (xcg >= fNx - 0.5) xcg -= float(fNx);

  int ixcg = int(xcg + 0.5);
  int iycg = int(ycg + 0.5);
  float ddx = fabs(xcg - ixcg);
  float ddy = fabs(ycg - iycg);

  float xg, yg, zg;
  Tower2Global(en, xcg, ycg, xg, yg, zg);

  float theta, phi;
  GetImpactThetaPhi(xg, yg, zg, theta, phi);

  int isx = 1;
  if (xcg - ixcg < 0) isx = -1;
  int isy = 1;
  if (ycg - iycg < 0) isy = -1;

  int idx = iTowerDist(ixcg, ix) * isx;
  int idy = (iy - iycg) * isy;

  int id = -1;
  if (idx == 0 && idy == 0)
    id = 0;
  else if (idx == 1 && idy == 0)
    id = 1;
  else if (idx == 1 && idy == 1)
    id = 2;
  else if (idx == 0 && idy == 1)
    id = 3;

  if (id < 0)
  {
    float dx = fabs(fTowerDist(xcg, float(ix)));
    float dy = fabs(iy - ycg);
    float rr = sqrt(dx * dx + dy * dy);
    //    return PredictEnergyParam(en, dx, dy);
    return _emcprof->PredictEnergyR(en, theta, phi, rr);
  }

  float ep[4], err[4];
  for (int ip = 0; ip < 4; ip++)
  {
    _emcprof->PredictEnergy(ip, en, theta, phi, ddx, ddy, ep[ip], err[ip]);
  }

  float eout;

  if (id == 0)
    eout = (ep[1] + ep[2]) / 2. + ep[3];
  else if (id == 1)
    eout = (ep[0] - ep[2]) / 2. - ep[3];
  else if (id == 3)
    eout = (ep[0] - ep[1]) / 2. - ep[3];
  else
    eout = ep[3];

  //  if( eout<0 ) printf("id=%d eout=%f: ep= %f %f %f %f Input: E=%f xcg=%f ycg=%f\n",id,eout,ep[0],ep[1],ep[2],ep[3],en,xcg,ycg);
  if (eout < 0) eout = 1e-6;

  return eout;
}

// ///////////////////////////////////////////////////////////////////////////

float BEmcRec::GetTowerEnergy(int iy, int iz, std::vector<EmcModule>* plist)
{
  int nn = plist->size();
  if (nn <= 0) return 0;

  for (int i = 0; i < nn; i++)
  {
    int ich = (*plist)[i].ich;
    int iyt = ich / fNx;
    int izt = ich % fNx;
    if (iy == iyt && iz == izt)
    {
      return (*plist)[i].amp;
    }
  }
  return 0;
}

// !!!!! Change here to a ponter to HitList
float BEmcRec::GetProb(std::vector<EmcModule> HitList, float en, float xg, float yg, float zg, float& chi2, int& ndf)
// Do nothing; should be defined in a detector specific module BEmcRec{Name}
{
  //  float enoise = 0.01;  // 10 MeV per tower
  float enoise = GetProbNoiseParam();
  //  float thresh = 0.01;
  float thresh = GetTowerThreshold();

  chi2 = 0;
  ndf = 0;
  if (_emcprof == nullptr) return -1;

  if (!(_emcprof->IsLoaded()))
  {
    return -1;
  }

  int nn = HitList.size();
  if (nn <= 0) return -1;

  float theta, phi;
  GetImpactThetaPhi(xg, yg, zg, theta, phi);

  // z coordinate below means x coordinate

  float etot;
  float zcg, ycg;
  float zz, yy, yz;
  Momenta(&HitList, etot, zcg, ycg, zz, yy, yz, thresh);

  int iz0cg = int(zcg + 0.5);
  int iy0cg = int(ycg + 0.5);
  float ddz = fabs(zcg - iz0cg);
  float ddy = fabs(ycg - iy0cg);

  int isz = 1;
  if (zcg - iz0cg < 0) isz = -1;
  int isy = 1;
  if (ycg - iy0cg < 0) isy = -1;

  // 4 central towers: 43
  //                   12
  // Tower 1 - central one
  float e1, e2, e3, e4;
  e1 = GetTowerEnergy(iy0cg, iz0cg, &HitList);
  e2 = GetTowerEnergy(iy0cg, iz0cg + isz, &HitList);
  e3 = GetTowerEnergy(iy0cg + isy, iz0cg + isz, &HitList);
  e4 = GetTowerEnergy(iy0cg + isy, iz0cg, &HitList);
  if (e1 < thresh) e1 = 0;
  if (e2 < thresh) e2 = 0;
  if (e3 < thresh) e3 = 0;
  if (e4 < thresh) e4 = 0;

  float e1t = (e1 + e2 + e3 + e4) / etot;
  float e2t = (e1 + e2 - e3 - e4) / etot;
  float e3t = (e1 - e2 - e3 + e4) / etot;
  float e4t = (e3) / etot;
  //  float rr = sqrt((0.5-ddz)*(0.5-ddz)+(0.5-ddy)*(0.5-ddy));

  // Predicted values
  const int NP = 4;  // From BEmcProfile
  float ep[NP];
  float err[NP];
  for (int ip = 0; ip < NP; ip++)
  {
    _emcprof->PredictEnergy(ip, en, theta, phi, ddz, ddy, ep[ip], err[ip]);
    if (ep[ip] < 0)
    {
      return -1;
    }
    if (ip < 3)
    {
      err[ip] = sqrt(err[ip] * err[ip] + 4 * enoise * enoise / etot / etot);
    }
    else
    {
      err[ip] = sqrt(err[ip] * err[ip] + 1 * enoise * enoise / etot / etot);
    }
  }

  chi2 = 0.;
  chi2 += (ep[0] - e1t) * (ep[0] - e1t) / err[0] / err[0];
  chi2 += (ep[1] - e2t) * (ep[1] - e2t) / err[1] / err[1];
  chi2 += (ep[2] - e3t) * (ep[2] - e3t) / err[2] / err[2];
  chi2 += (ep[3] - e4t) * (ep[3] - e4t) / err[3] / err[3];
  ndf = 4;

  chi2 /= 1.5;

  float prob = TMath::Prob(chi2, ndf);

  return prob;
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

void BEmcRec::CopyVector(const int* from, int* to, int N)
{
  if (N <= 0) return;
  for (int i = 0; i < N; i++) to[i] = from[i];
}

// ///////////////////////////////////////////////////////////////////////////

void BEmcRec::CopyVector(const EmcModule* from, EmcModule* to, int N)
{
  if (N <= 0)
  {
    return;
  }
  for (int i = 0; i < N; i++)
  {
    to[i] = from[i];
  }
}

// ///////////////////////////////////////////////////////////////////////////

/* Future improvements:

1. FindClusters(): to ensure that all EmcModules are above energy threshold 
set by SetThreshold routine (or default one)

*/

// ///////////////////////////////////////////////////////////////////////////
// EOF
