#ifndef TOWERINFOCONTAINERV3_H
#define TOWERINFOCONTAINERV3_H

#include "TowerInfoContainer.h"
#include "TowerInfov3.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>

class TowerInfoContainerv3 : public TowerInfoContainer
{
 public:
  TowerInfoContainerv3(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerv3() {}
  PHObject *CloneMe() const override { return new TowerInfoContainerv3(*this); }
  TowerInfoContainerv3(const TowerInfoContainerv3 &);

  ~TowerInfoContainerv3() override;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;
  TowerInfov3 *get_tower_at_channel(int pos) override;
  TowerInfov3 *get_tower_at_key(int pos) override;

  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  size_t size() const override { return _clones->GetEntries(); }
  DETECTOR get_detectorid() const override {return _detector;}

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

 private:
  ClassDefOverride(TowerInfoContainerv3, 1);
};

#endif
