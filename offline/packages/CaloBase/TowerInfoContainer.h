#ifndef TOWERINFOCONTAINER_H
#define TOWERINFOCONTAINER_H

#include "TowerInfo.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>


class TowerInfoContainer : public PHObject
{
 public:
  TowerInfoContainer() = default;
  ~TowerInfoContainer() override = default;

  virtual void Reset() override;
  virtual void add(TowerInfo* /*ti*/, int /*pos*/);
  virtual TowerInfo* at(int /*pos*/);
  virtual size_t size() { return 0; }

 private:
  ClassDefOverride(TowerInfoContainer, 1);
};

#endif
