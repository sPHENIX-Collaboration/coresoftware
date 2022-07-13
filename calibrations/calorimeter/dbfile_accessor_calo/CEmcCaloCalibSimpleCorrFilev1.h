#ifndef CALOCALIB_CEmcSIMPLECORRFILEv1_H
#define CALOCALIB_CEmcSIMPLECORRFILEv1_H

#include "CaloCalibSimpleCorrFile.h"

#include <calobase/RawTowerDefs.h>

#include <string>

class CEmcCaloCalibSimpleCorrFilev1 : public CaloCalibSimpleCorrFile
{
 public:
  // CEmcCaloCalibSimpleCorrFilev1() : m_CalTreeName("emc_corr_tree") {}

  CEmcCaloCalibSimpleCorrFilev1(RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE)
    : CaloCalibSimpleCorrFile(caloid)
    , m_CalTreeName("emc_corr_tree")
  {
  }
  ~CEmcCaloCalibSimpleCorrFilev1() override {}

  /*
  void Reset();
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;
  */

  void Open(const std::string &) override;
  void View() override;
  void ViewReadable() override;

  float getCorr(const unsigned int ieta, const unsigned int iphi) override;

  RawTowerDefs::keytype TowKey(const unsigned int ieta, const unsigned int iphi)
  {
    return ieta * 1000 + iphi;
  }

  ConstIterator AddCorr(const unsigned int ieta, const unsigned int iphi, float corr) override;

  void set_CalTreeName(const char *inTreename) { m_CalTreeName = inTreename; }

 protected:
  //  std::array< std::array<double,64>, 24> m_RecalArray;
  // use towerid  map instead of 2d array

  //  TNtuple * _readNtup;

  std::string m_CalTreeName;
  //    ClassDefOverride(CEmcCaloCalibSimpleCorrFilev1, 2);
};

#endif
