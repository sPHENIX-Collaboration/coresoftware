#ifndef TOWERINFOCONTAINER_H
#define TOWERINFOCONTAINER_H

#include "TowerInfo.h"

#include <phool/PHObject.h>

#include <climits>
#include <map>

class TowerInfoContainer : public PHObject
{
 public:
  typedef std::map<unsigned int, TowerInfo*> TowerMap;
  typedef TowerMap::const_iterator ConstIter;
  typedef TowerMap::iterator Iter;

  enum DETECTOR
  {
    EMCAL = 0,
    HCAL = 1,
    SEPD = 2,
    MBD = 3,
    ZDC = 4,
    DETECTOR_INVALID = 9999
  };

  TowerInfoContainer() = default;
  ~TowerInfoContainer() override = default;
  void identify(std::ostream& os = std::cout) const override;

  virtual void Reset() override {}
  virtual TowerInfo* get_tower_at_channel(int /*index*/) { return nullptr; }
  virtual TowerInfo* get_tower_at_key(int /*key*/) { return nullptr; }
  virtual size_t size() const { return 0; }

  virtual unsigned int encode_key(unsigned int /*towerIndex*/) { return UINT_MAX; }
  virtual unsigned int decode_key(unsigned int /*towerIndex*/) { return UINT_MAX; }

  virtual unsigned int encode_epd(unsigned int /*towerIndex*/);
  virtual unsigned int encode_hcal(unsigned int /*towerIndex*/);
  virtual unsigned int encode_emcal(unsigned int /*towerIndex*/);
  virtual unsigned int encode_mbd(unsigned int /*towerIndex*/);
  virtual unsigned int encode_zdc(unsigned int /*towerIndex*/);

  virtual unsigned int decode_epd(unsigned int /*towerIndex*/);
  virtual unsigned int decode_hcal(unsigned int /*towerIndex*/);
  virtual unsigned int decode_emcal(unsigned int /*towerIndex*/);
  virtual unsigned int decode_mbd(unsigned int /*towerIndex*/);
  virtual unsigned int decode_zdc(unsigned int /*towerIndex*/);

  virtual unsigned int getTowerPhiBin(unsigned int /*towerIndex*/);
  virtual unsigned int getTowerEtaBin(unsigned int /*towerIndex*/);

 private:
  ClassDefOverride(TowerInfoContainer, 1);
};

#endif
