#ifndef CALOBASE_RAWTOWERCONTAINER_H
#define CALOBASE_RAWTOWERCONTAINER_H

#include "RawTowerDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>

class RawTower;

class RawTowerContainer : public PHObject
{
 public:
  typedef std::map<RawTowerDefs::keytype, RawTower *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  RawTowerContainer(RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE)
    : _caloid(caloid)
  {
  }

  ~RawTowerContainer() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;

  void setCalorimeterID(RawTowerDefs::CalorimeterId caloid) { _caloid = caloid; }
  RawTowerDefs::CalorimeterId getCalorimeterID() { return _caloid; }

  ConstIterator AddTower(const unsigned int ieta, const unsigned int iphi, RawTower *twr);
  ConstIterator AddTower(RawTowerDefs::keytype key, RawTower *twr);

  RawTower *getTower(RawTowerDefs::keytype key);
  const RawTower *getTower(RawTowerDefs::keytype key) const;

  RawTower *getTower(const unsigned int ieta, const unsigned int iphi);
  const RawTower *getTower(const unsigned int ieta, const unsigned int iphi) const;

  RawTower *getTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il);
  const RawTower *getTower(const unsigned int ieta, const unsigned int iphi, const unsigned int il) const;

  //! return all towers
  ConstRange getTowers(void) const;
  Range getTowers(void);

  unsigned int size() const { return _towers.size(); }
  void compress(const double emin);
  double getTotalEdep() const;

 protected:
  RawTowerDefs::CalorimeterId _caloid;
  Map _towers;

  ClassDefOverride(RawTowerContainer, 1)
};

#endif
