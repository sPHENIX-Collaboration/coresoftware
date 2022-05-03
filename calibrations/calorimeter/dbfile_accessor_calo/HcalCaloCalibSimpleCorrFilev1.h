#ifndef CALOCALIB_HCALSIMPLECORRFILEv1_H
#define CALOCALIB_HCALSIMPLECORRFILEv1_H

#include "CaloCalibSimpleCorrFile.h"

#include <calobase/RawTowerDefs.h>

#include <cstddef>
#include <iostream>
#include <map>
#include <utility>
#include <string>

class HcalCaloCalibSimpleCorrFilev1 : public CaloCalibSimpleCorrFile
{
 public:
  HcalCaloCalibSimpleCorrFilev1(){}
  ~HcalCaloCalibSimpleCorrFilev1() override {}

  /*
  void Reset();
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;
  */

  void Open(const char * ) override;
  void View() override;
  void ViewReadable() override;
  
  float getCorr(const unsigned int ieta, const unsigned int iphi) override;


  ConstIterator AddCorr(const unsigned int ieta, const unsigned int iphi, float corr) override;
  void set_CalibrationFileName(const char * inFileName) {m_CalibrationFileName = inFileName;}

  

 protected:
  //  std::array< std::array<double,64>, 24> m_RecalArray;
  // use towerid  map instead of 2d array

  std::string m_CalibrationFileName;
//    ClassDefOverride(HcalCaloCalibSimpleCorrFilev1, 2);

};

#endif
