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
  typedef std::map<unsigned int, TowerInfo *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::pair<Iterator, Iterator> Range;

  void Reset() override;
  TowerInfov1 *at(int pos) override;
  unsigned int encode_key(unsigned int towerIndex) override;
  unsigned int decode_key(unsigned int tower_key) override;

  unsigned int encode_epd(unsigned int towerIndex) override;
  unsigned int encode_hcal(unsigned int towerIndex) override;
  unsigned int encode_emcal(unsigned int towerIndex) override;

  unsigned int decode_epd(unsigned int towerIndex) override;
  unsigned int decode_hcal(unsigned int towerIndex) override;
  unsigned int decode_emcal(unsigned int towerIndex) override;

  Range getTowers(void);

  size_t size() override { return _clones->GetEntries(); }

  unsigned int getTowerPhiBin(unsigned int towerIndex) override;
  unsigned int getTowerEtaBin(unsigned int towerIndex) override;

 protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

  //! static Tower index map, not saved on DST output and constructed on the fly
  Map _towers;  //!

 private:
  ClassDefOverride(TowerInfoContainerv1, 1);
};

#endif
