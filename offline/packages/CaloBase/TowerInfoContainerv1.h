#ifndef TOWERINFOCONTAINERV1_H
#define TOWERINFOCONTAINERV1_H

#include "TowerInfoContainer.h"
#include "TowerInfov1.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>

class TowerInfoContainerv1 : public TowerInfoContainer
{
 public:
  TowerInfoContainerv1(DETECTOR detec = DETECTOR::EMCAL);
  ~TowerInfoContainerv1() override;

  void Reset() override;
  void add(TowerInfov1 *ti, int pos);
  TowerInfov1 *at(int pos) override;
  unsigned int encode_key(unsigned int towerIndex) override;
  TowerMap getTowerMap() override;

  size_t size() override { return _clones->GetEntries(); }

  unsigned int getTowerPhiBin(unsigned int towerIndex) override;
  unsigned int getTowerEtaBin(unsigned int towerIndex) override;

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(int key) const override { return _map.find(key); }
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(int key) override { return _map.find(key); }
  Iter end() override { return _map.end(); }

 protected:
  TClonesArray *_clones;
  DETECTOR _detector;
  TowerMap _map;

 private:
  using TowerInfoContainer::add;

  ClassDefOverride(TowerInfoContainerv1, 1);
};

#endif
