#ifndef TOWERINFOCONTAINER_H
#define TOWERINFOCONTAINER_H

#include "TowerInfo.h"

#include <phool/PHObject.h>

#include <climits>
#include <map>

class TowerInfoContainer : public PHObject
{
 public:
  typedef std::map<int, TowerInfo*> TowerMap;
  typedef TowerMap::const_iterator ConstIter;
  typedef TowerMap::iterator Iter;

  enum DETECTOR
  {
    EMCAL = 0,
    HCAL = 1,
    SEPD = 2
  };

  TowerInfoContainer() = default;
  ~TowerInfoContainer() override = default;
  typedef std::map<unsigned int, TowerInfo *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::pair<Iterator, Iterator> Range;

  virtual void Reset() override {}
  virtual void add(TowerInfo* /*ti*/, int /*pos*/) {}
  virtual TowerInfo* at(int /*pos*/) { return nullptr; }
  virtual unsigned int encode_key(unsigned int /*towerIndex*/) { return UINT_MAX; }

  virtual size_t size() { return 0; }

  virtual unsigned int getTowerPhiBin(unsigned int /*towerIndex*/) { return UINT_MAX; }
  virtual unsigned int getTowerEtaBin(unsigned int /*towerIndex*/) { return UINT_MAX; }

  virtual ConstIter begin() const;
  virtual ConstIter find(int key) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(int key);
  virtual Iter end();

 private:
  ClassDefOverride(TowerInfoContainer, 1);
};

#endif
