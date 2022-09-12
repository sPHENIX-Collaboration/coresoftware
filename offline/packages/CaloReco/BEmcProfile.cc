#include "BEmcProfile.h"
#include "BEmcCluster.h"

#include <TFile.h>
#include <TH1.h>  // for TH1F
#include <TMath.h>
#include <TROOT.h>

#include <boost/format.hpp>

#include <cmath>    // for sqrt, log, pow, fabs
#include <cstdlib>  // for abs
#include <iostream>
#include <memory>  // for allocator_traits<>::value_type

const int NP = 4;  // Number of profiles in a bin

BEmcProfile::BEmcProfile(const std::string& fname)
  : bloaded(false)
  , thresh(0.01)
  , nen(0)
  , nth(0)
  , energy_array(nullptr)
  , theta_array(nullptr)
  , hmean(nullptr)
  , hsigma(nullptr)
  , m_Verbosity(0)
{
  TFile* f = new TFile(fname.c_str());
  if (f->IsZombie())
  {
    std::cout << "BEmcProfile: Error when opening profile data file " << fname << std::endl;
    return;
  }

  gROOT->cd();

  TH1F* hen = static_cast<TH1F*>(f->Get("hen"));
  if (!hen)
  {
    std::cout << "BEmcProfile: Error when loading profile data: hen" << std::endl;
    f->Close();
    return;
  }

  TH1F* hth = static_cast<TH1F*>(f->Get("hth"));
  if (!hth)
  {
    std::cout << "BEmcProfile: Error when loading profile data: hth" << std::endl;
    f->Close();
    return;
  }

  nen = hen->GetNbinsX();
  nth = hth->GetNbinsX();

  energy_array = new float[nen];
  for (int i = 0; i < nen; i++) energy_array[i] = hen->GetBinContent(i + 1);
  theta_array = new float[nth];
  for (int i = 0; i < nth; i++) theta_array[i] = hth->GetBinContent(i + 1);

  if (Verbosity())
  {
    std::cout << "BEmcProfile: Loading profile data from file " << fname << std::endl;

    std::cout << " " << nen << " bins in energy: ";
    for (int i = 0; i < nen; i++)
    {
      std::cout << boost::format("%4.1f, ") % energy_array[i];
    }
    std::cout << std::endl;
    std::cout << "  " << nth << " bins in theta:  ";
    for (int i = 0; i < nth; i++)
    {
      std::cout << boost::format("%4.2f, ") % theta_array[i];
    }
    std::cout << std::endl;
  }
  hmean = new TH1F*[nen * nth * NP];
  hsigma = new TH1F*[nen * nth * NP];
  hr4 = new TH1F*[nen * nth];

  TH1F* hh;
  int ii = 0;
  int ii2 = 0;

  std::string hname;
  for (int it = 0; it < nth; it++)
  {
    for (int ie = 0; ie < nen; ie++)
    {
      for (int ip = 0; ip < NP; ip++)
      {
        hname = boost::str(boost::format("hmean%d_en%d_t%d") % (ip + 1) % ie % it);
        hh = static_cast<TH1F*>(f->Get(hname.c_str()));
        if (!hh)
        {
          std::cout << "BEmcProfile: Could not load histogram " << hname
                    << ", Error when loading profile data for hmean it = "
                    << it << ", ie = " << ie << "ip = " << ip << std::endl;
          f->Close();
          return;
        }
        hmean[ii] = static_cast<TH1F*>(hh->Clone());

        hname = boost::str(boost::format("hsigma%d_en%d_t%d") % (ip + 1) % ie % it);
        hh = static_cast<TH1F*>(f->Get(hname.c_str()));
        if (!hh)
        {
          std::cout << "BEmcProfile: Could not load histogram " << hname
                    << ", Error when loading profile data for hsigma it = "
                    << it << ", ie = " << ie << ", ip = " << ip << std::endl;
          f->Close();
          return;
        }
        hsigma[ii] = static_cast<TH1F*>(hh->Clone());

        ii++;
      }

      hname = boost::str(boost::format("hr4_en%d_t%d") % ie % it);
      hh = static_cast<TH1F*>(f->Get(hname.c_str()));

      if (!hh)
      {
        std::cout << "BEmcProfile: Error when loading profile data for hr4 it = "
                  << it << ", ie = " << ie << std::endl;
        f->Close();
        return;
      }
      hr4[ii2] = static_cast<TH1F*>(hh->Clone());
      ii2++;
    }
  }

  f->Close();
  bloaded = true;
}

BEmcProfile::~BEmcProfile()
{
  if (bloaded)
  {
    delete[] energy_array;
    delete[] theta_array;
    for (int i = 0; i < nth * nen * NP; i++)
    {
      delete hmean[i];
      delete hsigma[i];
    }
    delete[] hmean;
    delete[] hsigma;

    for (int i = 0; i < nth * nen; i++)
    {
      delete hr4[i];
    }
    delete[] hr4;
  }
}

float BEmcProfile::GetProb(std::vector<EmcModule>* plist, int NX, float en, float theta, float phi)
{
  float enoise = 0.01;  // 10 MeV per tower

  if (!bloaded)
  {
    //    std::cout << "Error in BEmcProfile::GetProb: profiles not loaded" << std::endl;
    return -1;
  }

  int nn = plist->size();
  if (nn <= 0) return -1;

  // z coordinate below means x coordinate

  float ee;
  int ich;  // iy, iz;

  int iy0 = -1, iz0 = -1;
  float emax = 0;
  for (int i = 0; i < nn; i++)
  {
    ee = (*plist)[i].amp;
    if (ee > emax)
    {
      emax = ee;
      ich = (*plist)[i].ich;
      iy0 = ich / NX;
      iz0 = ich % NX;
    }
  }
  if (emax <= 0) return -1;

  float etot = 0;
  float sz = 0;
  float sy = 0;
  for (int i = 0; i < nn; i++)
  {
    ee = (*plist)[i].amp;
    ich = (*plist)[i].ich;
    int iy = ich / NX;
    int iz = ich % NX;
    if (ee > thresh && abs(iz - iz0) <= 3 && abs(iy - iy0) <= 3)
    {
      etot += ee;
      sz += iz * ee;
      sy += iy * ee;
    }
  }
  float zcg = sz / etot;
  float ycg = sy / etot;
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
  e1 = GetTowerEnergy(iy0cg, iz0cg, plist, NX);
  e2 = GetTowerEnergy(iy0cg, iz0cg + isz, plist, NX);
  e3 = GetTowerEnergy(iy0cg + isy, iz0cg + isz, plist, NX);
  e4 = GetTowerEnergy(iy0cg + isy, iz0cg, plist, NX);
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
  float ep[NP];
  float err[NP];
  for (int ip = 0; ip < NP; ip++)
  {
    PredictEnergy(ip, en, theta, phi, ddz, ddy, ep[ip], err[ip]);
    if (ep[ip] < 0) return -1;
    if (ip < 3)
      err[ip] = sqrt(err[ip] * err[ip] + 4 * enoise * enoise / etot / etot);
    else
      err[ip] = sqrt(err[ip] * err[ip] + 1 * enoise * enoise / etot / etot);
  }

  float chi2 = 0.;
  chi2 += (ep[0] - e1t) * (ep[0] - e1t) / err[0] / err[0];
  chi2 += (ep[1] - e2t) * (ep[1] - e2t) / err[1] / err[1];
  chi2 += (ep[2] - e3t) * (ep[2] - e3t) / err[2] / err[2];
  chi2 += (ep[3] - e4t) * (ep[3] - e4t) / err[3] / err[3];
  int ndf = 4;

  float prob = TMath::Prob(chi2, ndf);

  return prob;
}

/*
float BEmcProfile::GetProbTest(std::vector<EmcModule>* plist, int NX, float en, float theta, float& test_rr, float& test_et, float& test_ep, float& test_err)
// This is only for test purposes. Can be removed if not needed
{
  int nn = plist->size();
  if( nn <= 0 ) return -1;

  // z coordinate below means x coordinate

  float ee;
  int ich, iy, iz;

  int iy0=-1, iz0=-1;
  float emax=0;
  for( int i=0; i<nn; i++ ) {
    ee = (*plist)[i].amp;
    if( ee>emax ) {
      emax = ee;
      ich = (*plist)[i].ich;
      iy0 = ich/NX;
      iz0 = ich%NX;
    }
  }
  if( emax<=0 ) return -1;

  float etot=0;
  float sz=0;
  float sy=0;
  for( int i=0; i<nn; i++ ) {
    ee = (*plist)[i].amp;
    ich = (*plist)[i].ich;
    iy = ich/NX;
    iz = ich%NX;
    if( ee>thresh && abs(iz-iz0)<=3 && abs(iy-iy0)<=3 ) {
      etot += ee;
      sz += iz*ee;
      sy += iy*ee;
    }
  }
  float zcg = sz/etot;
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
  e1 = GetTowerEnergy(iy0cg,    iz0cg,     plist, NX);
  e2 = GetTowerEnergy(iy0cg,    iz0cg+isz, plist, NX);
  e3 = GetTowerEnergy(iy0cg+isy,iz0cg+isz, plist, NX);
  e4 = GetTowerEnergy(iy0cg+isy,iz0cg,     plist, NX);
  if( e1<thresh ) e1 = 0;
  if( e2<thresh ) e2 = 0;
  if( e3<thresh ) e3 = 0;
  if( e4<thresh ) e4 = 0;

  float e1t = (e1+e2+e3+e4)/etot;
  float e2t = (e1+e2-e3-e4)/etot;
  float e3t = (e1-e2-e3+e4)/etot;
  float e4t = (e3)/etot;
  float rr = sqrt((0.5-ddz)*(0.5-ddz)+(0.5-ddy)*(0.5-ddy));

  // Predicted values
  float ep[NP];
  float err[NP];
  for( int ip=0; ip<NP; ip++ ) 
    PredictEnergy(ip, en, theta, ddz, ddy, ep[ip], err[ip]);

  float chi2 = 0.;
  chi2 += (ep[0]-e1t)*(ep[0]-e1t)/err[0]/err[0];
  chi2 += (ep[1]-e2t)*(ep[1]-e2t)/err[1]/err[1];
  chi2 += (ep[2]-e3t)*(ep[2]-e3t)/err[2]/err[2];
  chi2 += (ep[3]-e4t)*(ep[3]-e4t)/err[3]/err[3];
  int ndf = 4;

  float prob = TMath::Prob(chi2, ndf);

  test_rr = rr;
  test_et = e1t;
  test_ep = ep[0];
  test_err = err[0];

  return prob;
}
*/

void BEmcProfile::PredictEnergy(int ip, float energy, float theta, float /*phi*/, float ddz, float ddy, float& ep, float& err)
// ip changes from 0 to NP-1, meaning the profile index 1,2,..,NP
{
  ep = err = -1;

  if (!bloaded)
  {
    //std::cout << "Error in BEmcProfile::PredictEnergy: profiles not loaded" << std::endl;
    return;
  }

  if (ip < 0 || ip >= NP)
  {
    std::cout << "Error in BEmcProfile::PredictEnergy: profile index = "
              << ip << " but should be from 0 to " << NP - 1 << std::endl;
    return;
  }

  if (energy <= 0)
  {
    std::cout << "Error in BEmcProfile::PredictEnergy: energy = "
              << energy << " but should be >0" << std::endl;
    return;
  }
  if (theta < 0)
  {
    std::cout << "Error in BEmcProfile::PredictEnergy: theta = "
              << theta << " but should be >=0" << std::endl;
    return;
  }

  if (ddz < 0 || ddz > 0.5)
  {
    std::cout << "Error in BEmcProfile::PredictEnergy: ddz = "
              << ddz << " but should be from 0 to 0.5" << std::endl;
  }

  if (ddy < 0 || ddy > 0.5)
  {
    std::cout << "Error in BEmcProfile::PredictEnergy: ddy = "
              << ddy << " but should be from 0 to 0.5" << std::endl;
  }

  // Safety margin (slightly away from bin edge)
  //  if (ddz < 0.021) ddz = 0.021;
  //  if (ddz > 0.489) ddz = 0.489;
  //  if (ddy < 0.021) ddy = 0.021;
  //  if (ddy > 0.489) ddy = 0.489;

  // Energy bin
  int ie2 = 0;
  while (ie2 < nen && energy > energy_array[ie2]) ie2++;
  if (ie2 == 0)
    ie2 = 1;
  else if (ie2 >= nen)
    ie2 = nen - 1;
  int ie1 = ie2 - 1;
  //  int ie1 = ie2-2; // For a test()

  // Theta bin
  int it2 = 0;
  while (it2 < nth && theta > theta_array[it2]) it2++;
  if (it2 == 0)
    it2 = 1;
  else if (it2 >= nth)
    it2 = nth - 1;
  int it1 = it2 - 1;
  //  int it1 = it2-2; // For a test()

  //  std::cout << "Energy bin= " << ie1 << " " << ie2 << " ("
  //       << energy << ")  Theta bin= " << it1 << " " << it2
  //       << " (" << theta << ")" << std::endl;

  float rr = sqrt((0.5 - ddz) * (0.5 - ddz) + (0.5 - ddy) * (0.5 - ddy));

  float xx = rr;
  if (ip == 1)
    xx = ddy;
  else if (ip == 2)
    xx = ddz;

  float en1 = energy_array[ie1];
  float en2 = energy_array[ie2];
  float th1 = theta_array[it1];
  float th2 = theta_array[it2];

  // 1st index - ie, second index - it
  int ii11 = ip + ie1 * NP + it1 * nen * NP;
  int ii21 = ip + ie2 * NP + it1 * nen * NP;
  int ii12 = ip + ie1 * NP + it2 * nen * NP;
  int ii22 = ip + ie2 * NP + it2 * nen * NP;

  int ibin = hmean[ii11]->FindBin(xx);

  // Log (1/sqrt) energy dependence of mean (sigma)
  //
  float pr11 = hmean[ii11]->GetBinContent(ibin);
  float pr21 = hmean[ii21]->GetBinContent(ibin);
  float prt1 = pr11 + (pr21 - pr11) / (log(en2) - log(en1)) * (log(energy) - log(en1));
  if (prt1 < 0) prt1 = 0;

  float er11 = hsigma[ii11]->GetBinContent(ibin);
  float er21 = hsigma[ii21]->GetBinContent(ibin);
  float ert1 = er11 + (er21 - er11) / (1. / sqrt(en2) - 1. / sqrt(en1)) * (1. / sqrt(energy) - 1. / sqrt(en1));
  if (ert1 < 0) ert1 = 0;

  float pr12 = hmean[ii12]->GetBinContent(ibin);
  float pr22 = hmean[ii22]->GetBinContent(ibin);
  float prt2 = pr12 + (pr22 - pr12) / (log(en2) - log(en1)) * (log(energy) - log(en1));
  if (prt2 < 0) prt2 = 0;

  float er12 = hsigma[ii12]->GetBinContent(ibin);
  float er22 = hsigma[ii22]->GetBinContent(ibin);
  float ert2 = er12 + (er22 - er12) / (1. / sqrt(en2) - 1. / sqrt(en1)) * (1. / sqrt(energy) - 1. / sqrt(en1));
  if (ert2 < 0) ert2 = 0;

  // Quadratic theta dependence of mean and sigma
  //
  float pr = prt1 + (prt2 - prt1) / (pow(th2, 2) - pow(th1, 2)) * (pow(theta, 2) - pow(th1, 2));
  if (pr < 0) pr = 0;
  float er = ert1 + (ert2 - ert1) / (pow(th2, 2) - pow(th1, 2)) * (pow(theta, 2) - pow(th1, 2));
  if (er < 0) er = 0;

  // Additional error due to binning in xx
  //
  int ibin1 = ibin;
  if (ibin > 1) ibin1 = ibin - 1;
  int ibin2 = ibin;
  if (ibin < hmean[ii11]->GetNbinsX())
  {
    if (hmean[ii11]->GetBinContent(ibin + 1) > 0)
    {
      ibin2 = ibin + 1;
    }
  }
  float dd = (hmean[ii11]->GetBinContent(ibin2) -
              hmean[ii11]->GetBinContent(ibin1)) /
             2.;
  //  if( fabs(dd)>er )
  // {
  //   std::cout << "ie = " << ie1 << ", it = " << it1 << ", bin = "
  // 	 << ibin << ": " << er << " " << dd << std::endl;
  // }
  er = sqrt(er * er + dd * dd);

  ep = pr;
  err = er;
}

float BEmcProfile::PredictEnergyR(float energy, float theta, float /*phi*/, float rr)
{
  if (!bloaded)
  {
    //    std::cout << "Error in BEmcProfile::PredictEnergyR: profiles not loaded" << std::endl;
    return 0;
  }

  if (energy <= 0)
  {
    std::cout << "Error in BEmcProfile::PredictEnergyR: energy = " << energy
              << " but should be >0" << std::endl;
    return 0;
  }
  if (theta < 0)
  {
    std::cout << "Error in BEmcProfile::PredictEnergyR: theta = " << theta
              << " but should be >=0" << std::endl;
    return 0;
  }

  if (rr < 1)
  {
    std::cout << "Error in BEmcProfile::PredictEnergyR: rr = " << rr
              << " but should be >1" << std::endl;
  }

  // Energy bin
  int ie2 = 0;
  while (ie2 < nen && energy > energy_array[ie2]) ie2++;
  if (ie2 == 0)
    ie2 = 1;
  else if (ie2 >= nen)
    ie2 = nen - 1;
  int ie1 = ie2 - 1;
  //  int ie1 = ie2-2; // For a test()

  // Theta bin
  int it2 = 0;
  while (it2 < nth && theta > theta_array[it2]) it2++;
  if (it2 == 0)
    it2 = 1;
  else if (it2 >= nth)
    it2 = nth - 1;
  int it1 = it2 - 1;
  //  int it1 = it2-2; // For a test()

  //  printf("Energy bin= %d %d (%f)  Theta bin= %d %d (%f)\n",ie1,ie2,energy,it1,it2,theta);

  float en1 = energy_array[ie1];
  float en2 = energy_array[ie2];
  float th1 = theta_array[it1];
  float th2 = theta_array[it2];

  // 1st index - ie, second index - it
  int ii11 = ie1 + it1 * nen;
  int ii21 = ie2 + it1 * nen;
  int ii12 = ie1 + it2 * nen;
  int ii22 = ie2 + it2 * nen;

  int ibin = hr4[ii11]->FindBin(rr);

  // Log (1/sqrt) energy dependence of mean (sigma)
  //
  float pr11 = hr4[ii11]->GetBinContent(ibin);
  float pr21 = hr4[ii21]->GetBinContent(ibin);
  float prt1 = pr11 + (pr21 - pr11) / (log(en2) - log(en1)) * (log(energy) - log(en1));
  if (prt1 < 0) prt1 = 0;

  float pr12 = hr4[ii12]->GetBinContent(ibin);
  float pr22 = hr4[ii22]->GetBinContent(ibin);
  float prt2 = pr12 + (pr22 - pr12) / (log(en2) - log(en1)) * (log(energy) - log(en1));
  if (prt2 < 0) prt2 = 0;

  // Quadratic theta dependence of mean and sigma
  //
  float pr = prt1 + (prt2 - prt1) / (pow(th2, 2) - pow(th1, 2)) * (pow(theta, 2) - pow(th1, 2));
  if (pr < 0) pr = 0;

  return pr;
}

float BEmcProfile::GetTowerEnergy(int iy, int iz, std::vector<EmcModule>* plist, int NX)
{
  int nn = plist->size();
  if (nn <= 0) return 0;

  for (int i = 0; i < nn; i++)
  {
    int ich = (*plist)[i].ich;
    int iyt = ich / NX;
    int izt = ich % NX;
    if (iy == iyt && iz == izt)
    {
      return (*plist)[i].amp;
    }
  }
  return 0;
}
