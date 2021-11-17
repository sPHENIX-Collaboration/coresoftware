#ifndef CALOBASE_RAWTOWERDEADMAP_H
#define CALOBASE_RAWTOWERDEADMAP_H

#include "RawTowerDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <set>

class RawTowerDeadMap : public PHObject
{
 public:
  typedef std::set<RawTowerDefs::keytype> Map;

  ~RawTowerDeadMap() override {}

  void Reset() override;
  int isValid() const override;

  void identify(std::ostream &os = std::cout) const override;

  virtual void setCalorimeterID(RawTowerDefs::CalorimeterId /*caloid*/) {}
  virtual RawTowerDefs::CalorimeterId getCalorimeterID() { return RawTowerDefs::NONE; }
  virtual void addDeadTower(const unsigned int ieta, const unsigned int iphi);
  virtual void addDeadTower(RawTowerDefs::keytype key);

  virtual bool isDeadTower(RawTowerDefs::keytype key);
  virtual bool isDeadTower(const unsigned int ieta, const unsigned int iphi);
  //! return all towers
  virtual const Map &getDeadTowers(void) const;
  virtual Map &getDeadTowers(void);

  virtual unsigned int size() const { return 0; }

 protected:
  RawTowerDeadMap(RawTowerDefs::CalorimeterId /*caloid*/ = RawTowerDefs::NONE)
  {
  }

 private:
  ClassDefOverride(RawTowerDeadMap, 1)
};

#endif
