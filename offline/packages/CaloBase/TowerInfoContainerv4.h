#ifndef TOWERINFOCONTAINERV4_H
#define TOWERINFOCONTAINERV4_H

#include "TowerInfoContainer.h"
#include "TowerInfov4.h"

#include <TClonesArray.h>

#include <cstddef>
#include <iostream>

class PHObject;

// this is basically a copy of TowerInfoContainerv1.h, but with TowerInfov4...
class TowerInfoContainerv4 : public TowerInfoContainer
{
 public:
  TowerInfoContainerv4(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerv4() {}
  PHObject *CloneMe() const override { return new TowerInfoContainerv4(*this); }
  TowerInfoContainerv4(const TowerInfoContainerv4 &);

  ~TowerInfoContainerv4() override;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;
  TowerInfov4 *get_tower_at_channel(int pos) override;
  TowerInfov4 *get_tower_at_key(int pos) override;

  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  size_t size() const override { return _clones->GetEntries(); }
  DETECTOR get_detectorid() const override { return _detector; }

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

 private:
  ClassDefOverride(TowerInfoContainerv4, 1);
};

#endif
