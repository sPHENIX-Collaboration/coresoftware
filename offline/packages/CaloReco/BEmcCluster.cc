// Name: EmcCluster.cxx
// Author: A. Bazilevsky (RIKEN-BNL)
// Major modifications by M. Volkov (RRC KI) Jan 27 2000

#include "BEmcCluster.h"
#include "BEmcRec.h"

#include <iostream>

// Define and initialize static members

// Max number of peaks in cluster; used in EmcCluster::GetSubClusters(...)
int const EmcCluster::fgMaxNofPeaks = 1000;

// Used in EmcCluster::GetSubClusters(...), it is the number of iterations
// when fitting photons to peak areas
int const EmcCluster::fgPeakIter = 6;

// Emin cuts cluster energy: Ecl >= Epk1+Epk2 +...+Epkn !!!!!!!!!!!!!
float const EmcCluster::fgEmin = 0.002;

EmcModule::EmcModule()
  : ich(0)
  , amp(0)
  , tof(0)
{
}

//_____________________________________________________________________________
EmcModule::EmcModule(int ich_, float amp_, float tof_)
  : ich(ich_)
  , amp(amp_)
  , tof(tof_)
{
}

// ///////////////////////////////////////////////////////////////////////////
// EmcCluster member functions

void EmcCluster::GetCorrPos(float& xc, float& yc)
// Returns the cluster corrected position in tower units
// Corrected for S-oscilations, not for shower depth
// Shower depth (z-coord) is defined in fOwner->Tower2Global()
{
  float e, x, y, xx, yy, xy;

  fOwner->Momenta(&fHitList, e, x, y, xx, yy, xy);
  fOwner->CorrectPosition(e, x, y, xc, yc);
}

// ///////////////////////////////////////////////////////////////////////////

void EmcCluster::GetGlobalPos(float& xA, float& yA, float& zA)
// Returns the cluster position in PHENIX global coord system
{
  float xc, yc;
  float e = GetTotalEnergy();
  GetCorrPos(xc, yc);
  fOwner->Tower2Global(e, xc, yc, xA, yA, zA);
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTowerEnergy(int ich)
// Returns the energy of the ich-tower (0 if ich not found in the fHitList)
{
  std::vector<EmcModule>::iterator ph;
  if (fHitList.empty()) return 0;
  ph = fHitList.begin();
  while (ph != fHitList.end())
  {
    if ((*ph).ich == ich) return (*ph).amp;
    ++ph;
  }
  return 0;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTowerEnergy(int ix, int iy)
// Returns the energy of the tower ix,iy (0 if tower not found in the fHitList)
{
  std::vector<EmcModule>::iterator ph;

  if (fHitList.empty()) return 0;
  ph = fHitList.begin();
  int ich = iy * fOwner->GetNx() + ix;
  return GetTowerEnergy(ich);
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTowerToF(int ich)
// Returns the ToF of the ich-tower (0 if ich not found in the fHitList)
{
  std::vector<EmcModule>::iterator ph;
  if (fHitList.empty()) return 0;
  ph = fHitList.begin();
  while (ph != fHitList.end())
  {
    if ((*ph).ich == ich) return (*ph).tof;
    ++ph;
  }
  return 0;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetTotalEnergy()
// Returns the cluster total energy
{
  std::vector<EmcModule>::iterator ph;
  float et = 0;
  if (fHitList.empty()) return 0;
  ph = fHitList.begin();
  while (ph != fHitList.end())
  {
    et += (*ph).amp;
    ++ph;
  }
  return et;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetECoreCorrected()
// Returns the energy in core towers around the cluster Center of Gravity
// Corrected for energy leak sidewise from core towers
{
  float e, x, y, xx, yy, xy;
  float ecore, ecorecorr;
  ecore = GetECore();
  fOwner->Momenta(&fHitList, e, x, y, xx, yy, xy);
  fOwner->CorrectECore(ecore, x, y, ecorecorr);
  return ecorecorr;
}

float EmcCluster::GetECore()
// Returns the energy in core towers around the cluster Center of Gravity
{
  const float thresh = 0.01;

  std::vector<EmcModule>::iterator ph;
  float xcg, ycg, xx, xy, yy;
  float energy, es;

  fOwner->Momenta(&fHitList, energy, xcg, ycg, xx, yy, xy);
  //  fOwner->SetProfileParameters(0, energy, xcg, ycg);

  es = 0;
  if (fHitList.empty()) return 0;
  ph = fHitList.begin();
  while (ph != fHitList.end())
  {
    int ixy = (*ph).ich;
    int iy = ixy / fOwner->GetNx();
    int ix = ixy - iy * fOwner->GetNx();
    //      dx = xcg - ix;
    //    float dx = fOwner->fTowerDist(float(ix), xcg);
    //    float dy = ycg - iy;
    //    float et = fOwner->PredictEnergy(dx, dy, energy);
    float et = fOwner->PredictEnergy(energy, xcg, ycg, ix, iy);
    if (et > thresh) es += (*ph).amp;
    ++ph;
  }
  return es;
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetE4()
// Returns the energy in 2x2 towers around the cluster Center of Gravity
{
  float et, xcg, ycg, xx, xy, yy;
  float e1, e2, e3, e4;
  int ix0, iy0, isx, isy;

  fOwner->Momenta(&fHitList, et, xcg, ycg, xx, yy, xy);
  ix0 = int(xcg + 0.5);
  iy0 = int(ycg + 0.5);

  isx = 1;
  if (xcg - ix0 < 0) isx = -1;
  isy = 1;
  if (ycg - iy0 < 0) isy = -1;

  e1 = GetTowerEnergy(ix0, iy0);
  e2 = GetTowerEnergy(ix0 + isx, iy0);
  e3 = GetTowerEnergy(ix0 + isx, iy0 + isy);
  e4 = GetTowerEnergy(ix0, iy0 + isy);

  return (e1 + e2 + e3 + e4);
}
// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetE9()
// Returns the energy in 3x3 towers around the cluster Center of Gravity
{
  float et, xcg, ycg, xx, xy, yy;
  int ich, ix0, iy0, nhit;

  nhit = fHitList.size();

  if (nhit <= 0) return 0;

  fOwner->Momenta(&fHitList, et, xcg, ycg, xx, yy, xy);
  ix0 = int(xcg + 0.5);
  iy0 = int(ycg + 0.5);
  ich = iy0 * fOwner->GetNx() + ix0;

  return GetE9(ich);
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetE9(int ich)
// Returns the energy in 3x3 towers around the tower ich
{
  std::vector<EmcModule>::iterator ph;

  int iy0 = ich / fOwner->GetNx();
  int ix0 = ich - iy0 * fOwner->GetNx();

  float es = 0;
  if (fHitList.empty()) return 0;
  ph = fHitList.begin();
  while (ph != fHitList.end())
  {
    int ixy = (*ph).ich;
    int iy = ixy / fOwner->GetNx();
    int ix = ixy - iy * fOwner->GetNx();
    int idx = fOwner->iTowerDist(ix0, ix);
    int idy = iy - iy0;
    if (abs(idx) <= 1 && abs(idy) <= 1) es += (*ph).amp;
    ++ph;
  }
  return es;
}

// ///////////////////////////////////////////////////////////////////////////

EmcModule EmcCluster::GetMaxTower()
// Returns the EmcModule with the maximum energy
{
  std::vector<EmcModule>::iterator ph;
  float emax = 0;
  EmcModule ht;

  ht.ich = -1;
  ht.amp = 0;
  ht.tof = 0;
  if (fHitList.empty()) return ht;

  ph = fHitList.begin();
  while (ph != fHitList.end())
  {
    if ((*ph).amp > emax)
    {
      emax = (*ph).amp;
      ht = *ph;
    }
    ++ph;
  }
  return ht;
}

// ///////////////////////////////////////////////////////////////////////////

void EmcCluster::GetMoments(float& x, float& y, float& xx, float& xy, float& yy)
//  Returns cluster 1-st (px,py) and 2-d momenta (pxx,pxy,pyy) in cell unit
{
  float e;
  fOwner->Momenta(&fHitList, e, x, y, xx, yy, xy);
}

// ///////////////////////////////////////////////////////////////////////////

float EmcCluster::GetProb(float& chi2, int& ndf)
{
  float e, xg, yg, zg;
  e = GetTotalEnergy();
  GetGlobalPos(xg, yg, zg);
  return fOwner->GetProb(fHitList, e, xg, yg, zg, chi2, ndf);
}

// ///////////////////////////////////////////////////////////////////////////

int EmcCluster::GetSubClusters(std::vector<EmcCluster>& PkList, std::vector<EmcModule>& ppeaks)
{
  // Splits the cluster onto subclusters
  // The number of subclusters is equal to the number of Local Maxima in a cluster.
  // Local Maxima can have the energy not less then
  // defined in fgMinPeakEnergy
  //
  // Output: PkList - vector of subclusters
  //         ppeaks - vector of peak EmcModules (one for each subcluster)
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
  float epk[fgMaxNofPeaks * 2], xpk[fgMaxNofPeaks * 2], ypk[fgMaxNofPeaks * 2];
  float ratio, eg, dx, dy, a;
  float *Energy[fgMaxNofPeaks], *totEnergy, *tmpEnergy;
  EmcModule *phit, *hlist, *vv;
  EmcCluster peak(fOwner);
  std::vector<EmcModule>::iterator ph;
  std::vector<EmcModule> hl;

  PkList.clear();
  ppeaks.clear();

  nhit = fHitList.size();

  if (nhit <= 0) return 0;

  hlist = new EmcModule[nhit];

  ph = fHitList.begin();
  vv = hlist;
  while (ph != fHitList.end()) *vv++ = *ph++;

  // sort by linear channel number
  qsort(hlist, nhit, sizeof(EmcModule), fOwner->HitNCompare);

  //
  //  Find peak (maximum) position (towers with local maximum amp)
  //
  npk = 0;
  for (ic = 0; ic < nhit; ic++)
  {
    float amp = hlist[ic].amp;
    if (amp > fOwner->GetPeakThreshold())
    {
      int ich = hlist[ic].ich;
      //old      int ichmax = ich + fOwner->GetNx() + 1;
      //old      int ichmin = ich - fOwner->GetNx() - 1;
      // Look into three raws only
      int ichmax = (ich / fOwner->GetNx() + 2) * fOwner->GetNx() - 1;
      int ichmin = (ich / fOwner->GetNx() - 1) * fOwner->GetNx();
      int ixc = ich - ich / fOwner->GetNx() * fOwner->GetNx();
      int inA = ic + 1;
      // check right, and 3 towers above if there is another max
      while ((inA < nhit) && (hlist[inA].ich <= ichmax))
      {
        int ixA = hlist[inA].ich - hlist[inA].ich / fOwner->GetNx() * fOwner->GetNx();
        //old	if( (abs(ixc-ixA) <= 1) && (hlist[inA].amp >= amp) ) goto new_ic;
        if ((abs(fOwner->iTowerDist(ixA, ixc)) <= 1) && (hlist[inA].amp >= amp)) goto new_ic;
        inA++;
      }
      inA = ic - 1;
      // check left, and 3 towers below if there is another max
      while ((inA >= 0) && (hlist[inA].ich >= ichmin))
      {
        int ixA = hlist[inA].ich - hlist[inA].ich / fOwner->GetNx() * fOwner->GetNx();
        //old	if( (abs(ixc-ixA) <= 1) && (hlist[inA].amp > amp) ) goto new_ic;
        if ((abs(fOwner->iTowerDist(ixA, ixc)) <= 1) && (hlist[inA].amp > amp)) goto new_ic;
        inA--;
      }
      if (npk >= fgMaxNofPeaks)
      {
        delete[] hlist;
        std::cout << "!!! Error in EmcCluster::GetSubClusters(): too many peaks in a cluster (>"
                  << fgMaxNofPeaks
                  << "). May need tower energy threshold increase for clustering." << std::endl;
        return -1;
      }

      // ic is a maximum in a 3x3 tower group
      PeakCh[npk] = ic;
      npk++;
    }
  new_ic:
    continue;
  }
  /*
  for( ipk=0; ipk<npk; ipk++ ) {
    ic = PeakCh[ipk];
    std::cout << "  " << ipk << ": E = " << hlist[ic].amp << std::endl;
  }
  */

  // there was only one peak
  if (npk <= 1)
  {
    hl.clear();
    for (int ich = 0; ich < nhit; ich++) hl.push_back(hlist[ich]);
    peak.ReInitialize(hl);
    PkList.push_back(peak);

    if (npk == 1)
      ppeaks.push_back(hlist[PeakCh[0]]);
    else
      ppeaks.push_back(GetMaxTower());

    delete[] hlist;
    return 1;
  }

  // there were more than one peak

  for (ipk = 0; ipk < npk; ipk++)
  {
    Energy[ipk] = new float[nhit];
  }
  totEnergy = new float[nhit];
  tmpEnergy = new float[nhit];
  for (int i = 0; i < nhit; ++i)
  {
    totEnergy[i] = 0.0;
    tmpEnergy[i] = 0.0;
  }
  //
  // Divide energy in towers among photons positioned in every peak
  //
  ratio = 1.;
  for (int iter = 0; iter < fgPeakIter; iter++)
  {
    fOwner->ZeroVector(tmpEnergy, nhit);
    for (ipk = 0; ipk < npk; ipk++)
    {
      ic = PeakCh[ipk];
      if (iter > 0) ratio = Energy[ipk][ic] / totEnergy[ic];
      eg = hlist[ic].amp * ratio;
      ixypk = hlist[ic].ich;
      iypk = ixypk / fOwner->GetNx();
      ixpk = ixypk - iypk * fOwner->GetNx();
      epk[ipk] = eg;
      xpk[ipk] = 0.;  // CoG from the peak tower
      ypk[ipk] = 0.;  // CoG from the peak tower

      int ichmax = (iypk + 2) * fOwner->GetNx() - 1;
      int ichmin = (iypk - 1) * fOwner->GetNx();

      // add up energies of tower to the right and up
      if (ic < nhit - 1)
      {
        for (in = ic + 1; in < nhit; in++)
        {
          ixy = hlist[in].ich;
          iy = ixy / fOwner->GetNx();
          ix = ixy - iy * fOwner->GetNx();

          //old	  if( (ixy-ixypk) > fOwner->GetNx()+1 ) break;
          if (ixy > ichmax) break;
          //old	  if( abs(ix-ixpk) <= 1 ) {
          idx = fOwner->iTowerDist(ixpk, ix);
          idy = iy - iypk;
          if (abs(idx) <= 1)
          {
            if (iter > 0) ratio = Energy[ipk][in] / totEnergy[in];
            eg = hlist[in].amp * ratio;
            epk[ipk] += eg;
            xpk[ipk] += eg * idx;
            ypk[ipk] += eg * idy;
          }
        }
      }  // if ic

      // add up energies of tower to the left and down
      if (ic >= 1)
      {
        for (in = ic - 1; in >= 0; in--)
        {
          ixy = hlist[in].ich;
          iy = ixy / fOwner->GetNx();
          ix = ixy - iy * fOwner->GetNx();

          //old	  if( (ixypk-ixy) > fOwner->GetNx()+1 ) break;
          if (ixy < ichmin) break;
          //old	  if( abs(ix-ixpk) <= 1 ) {
          idx = fOwner->iTowerDist(ixpk, ix);
          idy = iy - iypk;
          if (abs(idx) <= 1)
          {
            if (iter > 0) ratio = Energy[ipk][in] / totEnergy[in];
            eg = hlist[in].amp * ratio;
            epk[ipk] += eg;
            xpk[ipk] += eg * idx;
            ypk[ipk] += eg * idy;
          }
        }
      }  // if ic

      xpk[ipk] = xpk[ipk] / epk[ipk] + ixpk;
      ypk[ipk] = ypk[ipk] / epk[ipk] + iypk;
      //      fOwner->SetProfileParameters(0, epk[ipk], xpk[ipk], ypk[ipk]);

      for (in = 0; in < nhit; in++)
      {
        ixy = hlist[in].ich;
        iy = ixy / fOwner->GetNx();
        ix = ixy - iy * fOwner->GetNx();
        dx = fOwner->fTowerDist(float(ix), xpk[ipk]);
        dy = ypk[ipk] - iy;
        a = 0;

        // predict energy within 2.5 cell square around local peak
        if (ABS(dx) < 2.5 && ABS(dy) < 2.5)
          //          a = epk[ipk] * fOwner->PredictEnergy(dx, dy, epk[ipk]);
          a = epk[ipk] * fOwner->PredictEnergy(epk[ipk], xpk[ipk], ypk[ipk], ix, iy);

        Energy[ipk][in] = a;
        tmpEnergy[in] += a;
      }

    }  // for ipk
    float flim = 0.000001;
    for (in = 0; in < nhit; in++)
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
  }  // for iter

  phit = new EmcModule[nhit];

  ng = 0;
  for (ipk = 0; ipk < npk; ipk++)
  {
    nh = 0;

    // fill phit, the tower distribution for a peak
    for (in = 0; in < nhit; in++)
    {
      if (tmpEnergy[in] > 0)
      {
        ixy = hlist[in].ich;
        a = hlist[in].amp * Energy[ipk][in] / tmpEnergy[in];
        if (a > fgEmin)
        {
          phit[nh].ich = ixy;
          phit[nh].amp = a;
          phit[nh].tof = hlist[in].tof;  // Not necessary here

          nh++;
        }
      }
    }

    // if number hits > 0
    if (nh > 0)
    {
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
      igmpk1[ipk] = ng;
      igmpk2[ipk] = ng;

      ng++;
    }
    else
    {
      igmpk1[ipk] = -1;
      igmpk2[ipk] = -1;
    }
  }  // for( ipk=

  fOwner->ZeroVector(tmpEnergy, nhit);
  for (ipk = 0; ipk < npk; ipk++)
  {
    ig = igmpk1[ipk];
    if (ig >= 0)
    {
      //      std::cout << "  " << ipk << ": X = " << xpk[ig]
      //           << " Y = " << ypk[ig] << std::endl;
      //      fOwner->SetProfileParameters(0, epk[ig], xpk[ig], ypk[ig]);
      for (in = 0; in < nhit; in++)
      {
        Energy[ipk][in] = 0;
        for (ig = igmpk1[ipk]; ig <= igmpk2[ipk]; ig++)
        {
          ixy = hlist[in].ich;
          iy = ixy / fOwner->GetNx();
          ix = ixy - iy * fOwner->GetNx();
          dx = fOwner->fTowerDist(float(ix), xpk[ig]);
          dy = ypk[ig] - iy;
          //          a = epk[ig] * fOwner->PredictEnergy(dx, dy, epk[ig]);
          a = epk[ig] * fOwner->PredictEnergy(epk[ig], xpk[ig], ypk[ig], ix, iy);
          Energy[ipk][in] += a;
          tmpEnergy[in] += a;
        }
      }  // for( in
    }    // if( ig >= 0
  }      // for( ipk

  nn = 0;
  for (ipk = 0; ipk < npk; ipk++)
  {
    nh = 0;
    for (in = 0; in < nhit; in++)
    {
      if (tmpEnergy[in] > 0)
      {
        ixy = hlist[in].ich;
        a = hlist[in].amp * Energy[ipk][in] / tmpEnergy[in];
        if (a > fgEmin)
        {
          phit[nh].ich = ixy;
          phit[nh].amp = a;
          phit[nh].tof = hlist[in].tof;

          nh++;
        }
      }
    }  // for( in
    if (nh > 0)
    {
      //      *ip++ = hlist[PeakCh[ipk]];
      ppeaks.push_back(hlist[PeakCh[ipk]]);
      hl.clear();
      for (in = 0; in < nh; in++) hl.push_back(phit[in]);
      peak.ReInitialize(hl);
      PkList.push_back(peak);
      nn++;
    }
  }  // for( ipk

  delete[] phit;
  delete[] hlist;
  for (ipk = 0; ipk < npk; ipk++) delete[] Energy[ipk];
  delete[] totEnergy;
  delete[] tmpEnergy;

  return nn;
}
