#ifndef CALOBASE_RAWTOWERDEADMAPV1_H
#define CALOBASE_RAWTOWERDEADMAPV1_H

#include "RawTowerDeadMap.h"
#include "RawTowerDefs.h"

#include <iostream>

class RawTowerDeadMapv1 : public RawTowerDeadMap
{
 public:
  RawTowerDeadMapv1(RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE)
    : _caloid(caloid)
  {
  }
  ~RawTowerDeadMapv1() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;

  void setCalorimeterID(RawTowerDefs::CalorimeterId caloid) override { _caloid = caloid; }
  RawTowerDefs::CalorimeterId getCalorimeterID() override { return _caloid; }
  void addDeadTower(const unsigned int ieta, const unsigned int iphi) override;
  void addDeadTower(RawTowerDefs::keytype key) override;

  bool isDeadTower(RawTowerDefs::keytype key) override;
  bool isDeadTower(const unsigned int ieta, const unsigned int iphi) override;
  //! return all towers
  const Map &getDeadTowers(void) const override;
  Map &getDeadTowers(void) override;

  unsigned int size() const override { return m_DeadTowers.size(); }

 private:
  RawTowerDefs::CalorimeterId _caloid;
  Map m_DeadTowers;

  ClassDefOverride(RawTowerDeadMapv1, 1)
};

#endif
