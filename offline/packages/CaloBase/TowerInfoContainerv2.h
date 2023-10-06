#ifndef TOWERINFOCONTAINERV2_H
#define TOWERINFOCONTAINERV2_H

#include "TowerInfoContainer.h"
#include "TowerInfov2.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>
// this is basically a copy of TowerInfoContainerv1.h, but with TowerInfov2...
class TowerInfoContainerv2 : public TowerInfoContainer
{
 public:
  TowerInfoContainerv2(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerv2() {}

  ~TowerInfoContainerv2() override;

  void Reset() override;
  TowerInfov2 *get_tower_at_channel(int pos) override;
  TowerInfov2 *get_tower_at_key(int pos) override;

  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  size_t size() override { return _clones->GetEntries(); }

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

 private:
  ClassDefOverride(TowerInfoContainerv2, 1);
};

#endif