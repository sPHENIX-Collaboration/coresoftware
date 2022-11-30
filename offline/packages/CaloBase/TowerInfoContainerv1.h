#ifndef TOWERINFOCONTAINERV1_H
#define TOWERINFOCONTAINERV1_H

#include "TowerInfov1.h"
#include "TowerInfoContainer.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>


class TowerInfoContainerv1 : public TowerInfoContainer
{
 public:
  TowerInfoContainerv1();
  ~TowerInfoContainerv1() override;

  void Reset() override;
  void add(TowerInfov1 *ti, int pos);
  TowerInfov1* at(int pos) override;
  size_t size() override { return _clones->GetEntries(); }

 protected:
  TClonesArray *_clones;

  ClassDefOverride(TowerInfoContainerv1, 1);
};

#endif
