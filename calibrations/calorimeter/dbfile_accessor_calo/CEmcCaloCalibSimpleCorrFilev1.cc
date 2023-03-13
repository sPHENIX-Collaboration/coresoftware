#include "CEmcCaloCalibSimpleCorrFilev1.h"

#include <calobase/RawTowerDefs.h>  // for keytype

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>  // for exit
#include <fstream>
#include <iostream>
#include <map>  // for _Rb_tree_iterator

class TBranch;

//CEmcCaloCalibSimpleCorrFilev1::CEmcCaloCalibSimpleCorrFilev1(RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE)

void CEmcCaloCalibSimpleCorrFilev1::Open(const std::string &CalibrationFileName)
{
  //  _corrs[0] = 9.92939;

  if (!CalibrationFileName.empty())
  {
    TFile mef(CalibrationFileName.c_str());
    TTree *corrtree = (TTree *) mef.Get(m_CalTreeName.c_str());
    if (!corrtree)
    {
      std::cout << "CEmcCCSCorrFilev1.cc:  getting simple emcal db calibrations from file "
                << CalibrationFileName
                << "; treename: " << m_CalTreeName
                << std::endl;
      return;
    }

    TBranch *b_corr;   //!
    TBranch *b_towid;  //!
    float corr;
    int towid;

    corrtree->SetBranchAddress("corr", &corr, &b_corr);
    corrtree->SetBranchAddress("towid", &towid, &b_towid);

    /*
      int etabin = -1;
      int phibin = -1;
      double recal = 1.;
      */

    // see elsewhere in this file comments about this specific key
    // lines with key should be updated to use more standard key
    // if this is not a standard tower id convention
    //int inkey = etabin*100+phibin;

    long int nentries = corrtree->GetEntriesFast();

    for (long int jentry = 0; jentry < nentries; jentry++)
    {
      //	  corrtree->LoadTree(jentry);
      corrtree->GetEntry(jentry);

      if (!std::isfinite(corr))
      {
        std::cout << "Calibration constant at towid " << towid
                  << " in " << CalibrationFileName
                  << " is not finite: " << corr << std::endl;
        gSystem->Exit(1);
        exit(1);
      }

      // _corrs stlmap
      _corrs[towid] = corr;

      // at() does a bounds check
      //	      m_RecalArray.at(etabin).at(phibin) = recal;
      // aboveline  is Chris's pull request updated version of Ejiro original code
    }
  }
}

void CEmcCaloCalibSimpleCorrFilev1::View()
{
  std::cout << "corrs View " << _corrs[0] << std::endl;
}

float CEmcCaloCalibSimpleCorrFilev1::getCorr(const unsigned int ieta, const unsigned int iphi)
{
  RawTowerDefs::keytype key = TowKey(ieta, iphi);

  if (_corrs.find(key) != _corrs.end())
  {
    return _corrs[key];
  }
  else
  {
    std::cout << "calibrations/hcalCCSCFv1:: "
              << "corr not found for key " << key
              << " ieta " << ieta << "  " << iphi
              << ", returning -999" << std::endl;
    return -999;
  }
}

CaloCalibSimpleCorrFile::ConstIterator CEmcCaloCalibSimpleCorrFilev1::AddCorr(const unsigned int ieta, const unsigned int iphi, float corr)
{
  //this is just test code it is probably not correct and needs reimplemented
  // the expected change is:
  // for key it should use one of the standard id's for each hcal channel
  // if this is not

  //  int key = ieta*100+iphi;
  RawTowerDefs::keytype key = TowKey(ieta, iphi);
  _corrs[key] = corr;
  return _corrs.find(key);
}

void CEmcCaloCalibSimpleCorrFilev1::ViewReadable()
{
  std::cout << "corrs View Readable " << _corrs[0] << std::endl;
}
