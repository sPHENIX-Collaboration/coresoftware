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
  PHObject *CloneMe() const override { return new TowerInfoContainerv1(*this); }
  TowerInfoContainerv1(const TowerInfoContainerv1 &);

  ~TowerInfoContainerv1() override;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;
  TowerInfov1 *get_tower_at_channel(int pos) override;
  TowerInfov1 *get_tower_at_key(int pos) override;

  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  size_t size() const override { return _clones->GetEntries(); }

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

 private:
  ClassDefOverride(TowerInfoContainerv1, 1);
};

#endif
