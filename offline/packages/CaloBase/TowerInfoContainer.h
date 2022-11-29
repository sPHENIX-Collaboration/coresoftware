#ifndef TOWERINFOCONTAINER_H
#define TOWERINFOCONTAINER_H

#include <cstddef>  // for size_t
#include <iostream>
#include <set>

#include <phool/PHObject.h>
#include <TClonesArray.h>

#include "TowerInfo.h"

class TowerInfoContainer : public PHObject
{
 public:
  TowerInfoContainer();
  ~TowerInfoContainer() override;

  void Reset() override;
  void add(TowerInfo *ci, int pos);
  TowerInfo* at(int pos);
  int size() { return _clones->GetEntries(); }

 private:
  TClonesArray *_clones;


  ClassDefOverride(TowerInfoContainer, 1);
};

#endif
