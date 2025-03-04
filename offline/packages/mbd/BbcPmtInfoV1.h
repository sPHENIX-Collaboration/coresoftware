#ifndef MBD_BBCPMTINFOV1_H
#define MBD_BBCPMTINFOV1_H

#include <calobase/TowerInfo.h>

#include <phool/PHObject.h>

#include <limits>

class BbcPmtInfoV1 : public TowerInfo
{
 public:
  BbcPmtInfoV1() = default;
  // BbcPmtInfoV1(TowerInfo& tower);
  ~BbcPmtInfoV1() override = default;
  void Reset() override;

  //! Clear is used by TClonesArray to reset the tower to initial state without calling destructor/constructor
  void Clear(Option_t* = "") override;

  //! Prints out exact identity of object
  void identify(std::ostream& out = std::cout) const override;

  //! isValid returns non zero if object contains vailid data
  virtual int isValid() const override
  {
    if (bq == std::numeric_limits<Float_t>::quiet_NaN()) return 0;
    return 1;
  }

  // short get_time() override { return get_t(); } // deprecated for MBD since MBD returns time in ns
  Float_t get_energy() override { return get_q(); }

  // void set_time(short t) override;  // too simple to be used for MBD
  void set_energy(float energy) override { set_q(energy); }

  Short_t get_pmt() { return bpmt; }
  Float_t get_q() { return bq; }
  Float_t get_t() { return btq; }
  Float_t get_tt() { return btt; }
  Float_t get_tq() { return btq; }

  void set_pmt(const Short_t pmt) { bpmt = pmt; }
  void set_q(const Float_t q) { bq = q; }
  void set_tt(const Float_t t) { btt = t; }
  void set_tq(const Float_t t) { btq = t; }
  void set_pmt(const Short_t pmt, const Float_t q, const Float_t tt, const Float_t tq)
  {
    bpmt = pmt;
    bq = q;
    btt = tt;
    btq = tq;
  }

 private:
  Short_t bpmt{-1};  ///
  Float_t bq{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t btt{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t btq{std::numeric_limits<Float_t>::quiet_NaN()};

  ClassDefOverride(BbcPmtInfoV1, 1);
};

#endif  // __BBCPMTINFOV1_H__
