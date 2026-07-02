#ifndef TOWERINFOCONTAINERSIMV3_H
#define TOWERINFOCONTAINERSIMV3_H

#include "TowerInfoContainer.h"
#include "TowerInfoSimv3.h"

#include <TClonesArray.h>

class PHObject;

class TowerInfoContainerSimv3 : public TowerInfoContainer
{
 public:
  TowerInfoContainerSimv3(DETECTOR detec);

  // default constructor for ROOT IO
  TowerInfoContainerSimv3() = default;
  PHObject *CloneMe() const override { return new TowerInfoContainerSimv3(*this); }
  TowerInfoContainerSimv3(const TowerInfoContainerSimv3 &);
  TowerInfoContainerSimv3 &operator=(const TowerInfoContainerSimv3 &) = delete;

  ~TowerInfoContainerSimv3() override;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;
  TowerInfoSimv3 *get_tower_at_channel(int pos) override;
  TowerInfoSimv3 *get_tower_at_key(int pos) override;

  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  size_t size() const override { return _clones->GetEntries(); }
  DETECTOR get_detectorid() const override { return _detector; }

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

 private:
  ClassDefOverride(TowerInfoContainerSimv3, 1);
};

#endif
