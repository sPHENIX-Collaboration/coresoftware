#include "HcalCaloCalibSimpleCorrFilev1.h"

#include <TSystem.h>

#include <cmath>
#include <cstdlib>  // for exit
#include <fstream>
#include <iostream>
#include <map>  // for _Rb_tree_iterator

void HcalCaloCalibSimpleCorrFilev1::Open(const std::string &CalibrationFileName)
{
  //  _corrs[0] = 9.92939;

  if (!CalibrationFileName.empty())
  {
    std::ifstream calibrate_tower;
    calibrate_tower.open(CalibrationFileName, std::ifstream::in);
    if (calibrate_tower.is_open())
    {
      int etabin = -1;
      int phibin = -1;
      double recal = 1.;
      calibrate_tower >> etabin >> phibin >> recal;

      std::cout << "HcalCCSCorrFilev1.cc:  getting simple hcal db calibrations from file "
                << CalibrationFileName
                << " :  " << etabin << " " << phibin << " " << recal << " ... "
                << std::endl;

      while (!calibrate_tower.eof())
      {
        if (!std::isfinite(recal))
        {
          std::cout << "Calibration constant at etabin " << etabin
                    << ", phibin " << phibin << " in " << CalibrationFileName
                    << " is not finite: " << recal << std::endl;
          gSystem->Exit(1);
          exit(1);
        }

        // see elsewhere in this file comments about this specific key
        // lines with key should be updated to use more standard key
        // if this is not a standard tower id convention
        int inkey = etabin * 100 + phibin;

        // _corrs stlmap
        _corrs[inkey] = recal;

        // at() does a bounds check
        //	      m_RecalArray.at(etabin).at(phibin) = recal;
        // aboveline  is Chris's pull request updated version of Ejiro original code
        calibrate_tower >> etabin >> phibin >> recal;
      }
      calibrate_tower.close();
    }
  }
}

void HcalCaloCalibSimpleCorrFilev1::View()
{
  std::cout << "corrs View " << _corrs[0] << std::endl;
}

float HcalCaloCalibSimpleCorrFilev1::getCorr(const unsigned int ieta, const unsigned int iphi)
{
  int key = ieta * 100 + iphi;

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

CaloCalibSimpleCorrFile::ConstIterator HcalCaloCalibSimpleCorrFilev1::AddCorr(const unsigned int ieta, const unsigned int iphi, float corr)
{
  //this is just test code it is probably not correct and needs reimplemented
  // the expected change is:
  // for key it should use one of the standard id's for each hcal channel
  // if this is not

  int key = ieta * 100 + iphi;
  _corrs[key] = corr;
  return _corrs.find(key);
}

void HcalCaloCalibSimpleCorrFilev1::ViewReadable()
{
  std::cout << "corrs View Readable " << _corrs[0] << std::endl;
}
