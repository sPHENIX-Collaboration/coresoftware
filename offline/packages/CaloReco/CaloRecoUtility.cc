#include "CaloRecoUtility.h"

#include "BEmcCluster.h"
#include "BEmcRec.h"  // for BEmcRec
#include "BEmcRecCEMC.h"

#include <calobase/RawCluster.h>
#include <calobase/RawTowerDefs.h>  // for decode_index1, decode_index2

#include <ffamodules/CDBInterface.h>

#include <phool/phool.h>  // for PHWHERE

#include <cmath>    // for atan2, log, cos, fabs, sin, sqrt
#include <cstdlib>  // for exit, getenv
#include <iostream>
#include <map>  // for _Rb_tree_const_iterator, operat...
#include <string>
#include <utility>  // for pair
#include <vector>   // for vector

void CaloRecoUtility::ShowerDepthCorrZVertex(RawCluster* clus, float vz)
{
  //

  float xA = clus->get_x();
  float yA = clus->get_y();
  float zA = clus->get_z();

  float zC = zA;

  //      bemcC.CorrectShowerDepth(clus->get_energy(),xt,yt,zt,xC,yC, zC);
  // Correction in z
  // Just tuned for sim data ... don't fully understand why it works like that

  float logE = log(0.1);
  if (clus->get_energy() > 0.1)
  {
    logE = std::log(clus->get_energy());
  }

  float rA = std::sqrt(xA * xA + yA * yA);
  //  float theta_twr = GetTowerTheta(xA,yA,zA);
  float theta_twr;
  if (std::fabs(zA) <= 15)
  {
    theta_twr = 0;
  }
  else if (zA > 15)
  {
    theta_twr = std::atan2(zA - 15, rA);
  }
  else
  {
    theta_twr = std::atan2(zA + 15, rA);
  }

  //  fVz = 0;

  float theta_tr = std::atan2(zA - vz, rA);
  float L = -1.3 + 0.7 * logE;  // Shower CG in long. direction
  float dz = L * std::sin(theta_tr - theta_twr) / std::cos(theta_twr);

  dz -= vz * 0.10;

  zC = zA - dz;

  clus->set_z(zC);
}

void CaloRecoUtility::ProbCorrsZVertex(RawCluster* clus, float vz)
{
  if (!_profLoaded)
  {
    LoadProfile();
  }

  float vvv[3];
  vvv[0] = 0;
  vvv[1] = 0;
  vvv[2] = vz;

  _bemc->SetVertex(vvv);
  // construct hitlist, which is just the towers hit in the EmcModule data structure format
  //   (for this cluster)

  EmcModule vhit;
  std::vector<EmcModule> HitList;
  HitList.erase(HitList.begin(), HitList.end());
  int ich;

  RawCluster::TowerConstRange towers = clus->get_towers();
  RawCluster::TowerConstIterator toweriter;

  for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
  {
    int ix = RawTowerDefs::decode_index2(toweriter->first);  // index2 is phi in CYL
    int iy = RawTowerDefs::decode_index1(toweriter->first);  // index1 is eta in CYL
    /*
    ix -= BINX0;
    iy -= BINY0;

    if (ix >= 0 && ix < NBINX && iy >= 0 && iy < NBINY)
      {

    */
    //      ich = iy * NBINX+ ix;
    ich = iy * 256 + ix;
    // add key field to vhit
    vhit.ich = ich;
    vhit.amp = toweriter->second;
    // vhit.amp = tower->get_energy() * fEnergyNorm;  // !!! Global Calibration
    vhit.tof = 0.01;  // not used

    HitList.push_back(vhit);

    //} // if check NBINS
  }

  float chi2 = 0;
  int ndf = 0;
  float prob = _bemc->GetProb(HitList, clus->get_energy(), clus->get_x(), clus->get_y(), clus->get_z(), chi2, ndf);

  clus->set_prob(prob);
  if (ndf > 0)
  {
    clus->set_chi2(chi2 / ndf);
  }
  else
  {
    clus->set_chi2(0);
  }
}

void CaloRecoUtility::LoadProfile()
{
  const char* calibroot = getenv("CALIBRATIONROOT");

  if (calibroot == nullptr)
  {
    std::cout << PHWHERE << "CALIBRATIONROOT env var not set" << std::endl;
    exit(1);
  }
  std::string fname_emc_prof = calibroot;
  fname_emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";

  std::cout << "CaloRecoUtility:::loading emc_prof from " << fname_emc_prof << std::endl;

  std::string url = CDBInterface::instance()->getUrl("EMCPROFILE", fname_emc_prof);
  _bemc->LoadProfile(url);

  _profLoaded = true;
}

CaloRecoUtility::CaloRecoUtility()
  : _profLoaded(false)
{
  _bemc = new BEmcRecCEMC();

  _bemc->SetDim(256, 96);

  _bemc->SetTowerThreshold(0.030);

  float fProbNoiseParam = 0.04;
  _bemc->SetProbNoiseParam(fProbNoiseParam);
}

// this two stupid functions are  only here because of a
// cppcheck warning that should have been suppressed
//  "recommended to have a copy constructor/op= because of
//  dynamic alloc resource"
CaloRecoUtility::CaloRecoUtility(CaloRecoUtility& cru)
{
  _profLoaded = false;

  if (cru._bemc == nullptr)
  {
    _bemc = nullptr;
    return;
  }

  _bemc = new BEmcRecCEMC();

  _bemc->SetDim(256, 96);

  _bemc->SetTowerThreshold(0.030);

  float fProbNoiseParam = 0.04;
  _bemc->SetProbNoiseParam(fProbNoiseParam);
}

CaloRecoUtility& CaloRecoUtility::operator=(const CaloRecoUtility& cru)
{
  if (this == &cru)
  {
    return *this;
  }

  _profLoaded = false;
  if (cru._bemc == nullptr)
  {
    _bemc = nullptr;
    return *this;
  }

  _bemc = new BEmcRecCEMC();

  _bemc->SetDim(256, 96);

  _bemc->SetTowerThreshold(0.030);

  float fProbNoiseParam = 0.04;
  _bemc->SetProbNoiseParam(fProbNoiseParam);

  return *this;
}

CaloRecoUtility::~CaloRecoUtility()
{
  delete _bemc;
}
