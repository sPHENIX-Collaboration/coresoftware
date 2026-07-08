#ifndef TOWERINFOCONTAINERV5_H
#define TOWERINFOCONTAINERV5_H

#include "TowerInfoContainer.h"
#include "TowerInfov5.h"

#include <TClonesArray.h>

#include <cstddef>
#include <iostream>

class PHObject;

class TowerInfoContainerv5 : public TowerInfoContainer
{
 public:
  TowerInfoContainerv5(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerv5() = default;
  PHObject *CloneMe() const override { return new TowerInfoContainerv5(*this); }
  TowerInfoContainerv5(const TowerInfoContainerv5 &);

  ~TowerInfoContainerv5() override;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;
  TowerInfov5 *get_tower_at_channel(int pos) override;
  TowerInfov5 *get_tower_at_key(int pos) override;

  size_t size() const override { return _clones->GetEntries(); }
  DETECTOR get_detectorid() const override { return _detector; }

 protected:
  TClonesArray *_clones{nullptr};
  DETECTOR _detector{DETECTOR_INVALID};

 private:
  ClassDefOverride(TowerInfoContainerv5, 1);
};

#endif
