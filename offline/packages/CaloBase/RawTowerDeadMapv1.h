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
  virtual ~RawTowerDeadMapv1() {}

  virtual void Reset() override;
  virtual int isValid() const override;
  virtual void identify(std::ostream &os = std::cout) const override;

  virtual void setCalorimeterID(RawTowerDefs::CalorimeterId caloid) override { _caloid = caloid; }
  virtual RawTowerDefs::CalorimeterId getCalorimeterID() override { return _caloid; }
  virtual void addDeadTower(const unsigned int ieta, const unsigned int iphi) override;
  virtual void addDeadTower(RawTowerDefs::keytype key) override;

  virtual bool isDeadTower(RawTowerDefs::keytype key) override;
  virtual bool isDeadTower(const unsigned int ieta, const unsigned int iphi) override;
  //! return all towers
  virtual const Map &getDeadTowers(void) const override;
  virtual Map &getDeadTowers(void) override;

  virtual unsigned int size() const override { return m_DeadTowers.size(); }

 private:
  RawTowerDefs::CalorimeterId _caloid;
  Map m_DeadTowers;

  ClassDefOverride(RawTowerDeadMapv1, 1)
};

#endif
