#ifndef TOWERINFOCONTAINER_H
#define TOWERINFOCONTAINER_H

#include "TowerInfo.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>

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

  virtual void Reset() override;
  virtual void add(TowerInfo* /*ti*/, int /*pos*/);
  virtual TowerInfo* at(int /*pos*/);
  virtual int encode_key(int /*towerIndex*/) { return 0; }
  virtual TowerMap getTowerMap();

  virtual size_t size() { return 0; }

  virtual int getTowerPhiBin(int /*towerIndex*/) { return -1; }
  virtual int getTowerEtaBin(int /*towerIndex*/) { return -1; }

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
