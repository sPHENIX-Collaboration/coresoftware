#ifndef TOWERINFOCONTAINERSIMV1_H
#define TOWERINFOCONTAINERSIMV1_H

#include "TowerInfoContainer.h"
#include "TowerInfoSimv1.h" 

#include <TClonesArray.h>

class PHObject;

class TowerInfoContainerSimv1 : public TowerInfoContainer
{
 public:
  TowerInfoContainerSimv1(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerSimv1() {}
  PHObject *CloneMe() const override { return new TowerInfoContainerSimv1(*this); }
  TowerInfoContainerSimv1(const TowerInfoContainerSimv1 &);

  ~TowerInfoContainerSimv1() override;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;
  TowerInfoSimv1 *get_tower_at_channel(int pos) override;
  TowerInfoSimv1 *get_tower_at_key(int pos) override;

  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  size_t size() const override { return _clones->GetEntries(); }
  DETECTOR get_detectorid() const override { return _detector; }

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

 private:
  ClassDefOverride(TowerInfoContainerSimv1, 1);
};

#endif
