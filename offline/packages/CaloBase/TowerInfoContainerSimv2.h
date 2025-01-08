#ifndef TOWERINFOCONTAINERSIMV2_H
#define TOWERINFOCONTAINERSIMV2_H

#include "TowerInfoContainer.h"
#include "TowerInfoSimv2.h" 

#include <TClonesArray.h>

class PHObject;

class TowerInfoContainerSimv2 : public TowerInfoContainer
{
 public:
  TowerInfoContainerSimv2(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerSimv2() {}
  PHObject *CloneMe() const override { return new TowerInfoContainerSimv2(*this); }
  TowerInfoContainerSimv2(const TowerInfoContainerSimv2 &);

  ~TowerInfoContainerSimv2() override;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;
  TowerInfoSimv2 *get_tower_at_channel(int pos) override;
  TowerInfoSimv2 *get_tower_at_key(int pos) override;

  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  size_t size() const override { return _clones->GetEntries(); }
  DETECTOR get_detectorid() const override { return _detector; }

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

 private:
  ClassDefOverride(TowerInfoContainerSimv2, 1);
};

#endif
