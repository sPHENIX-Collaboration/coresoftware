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
  virtual void Reset();
  virtual int isValid() const;

  virtual void identify(std::ostream &os = std::cout) const;

  virtual void setCalorimeterID(RawTowerDefs::CalorimeterId caloid) { _caloid = caloid; }
  virtual RawTowerDefs::CalorimeterId getCalorimeterID() { return _caloid; }
  virtual void addDeadTower(const unsigned int ieta, const unsigned int iphi);
  virtual void addDeadTower(RawTowerDefs::keytype key);

  virtual bool isDeadTower(RawTowerDefs::keytype key);
  virtual bool isDeadTower(const unsigned int ieta, const unsigned int iphi);
  //! return all towers
  virtual const Map &getDeadTowers(void) const;
  virtual Map &getDeadTowers(void);

  virtual unsigned int size() const { return m_DeadTowers.size(); }

 private:
  RawTowerDefs::CalorimeterId _caloid;
  Map m_DeadTowers;

  ClassDef(RawTowerDeadMapv1, 1)
};

#endif
