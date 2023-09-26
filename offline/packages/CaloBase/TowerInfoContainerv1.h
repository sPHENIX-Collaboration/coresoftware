#ifndef TOWERINFOCONTAINERV1_H
#define TOWERINFOCONTAINERV1_H

#include "TowerInfoContainer.h"
#include "TowerInfov1.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>

class TowerInfoContainerv1 : public TowerInfoContainer
{
 public:
  TowerInfoContainerv1(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerv1() {}

  ~TowerInfoContainerv1() override;

  void Reset() override;
  TowerInfov1 *get_tower_at_channel(int pos) override;
  TowerInfov1 *get_tower_at_key(int pos) override;

  size_t size() override { return _clones->GetEntries(); }

 protected:
  TClonesArray *_clones = nullptr;

 private:
  ClassDefOverride(TowerInfoContainerv1, 1);
};

#endif
