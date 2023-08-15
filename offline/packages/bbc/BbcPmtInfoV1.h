#ifndef __BBCPMTINFOV1_H__
#define __BBCPMTINFOV1_H__

#include <calobase/TowerInfo.h>

#include <phool/PHObject.h>

#include <cmath>

class BbcPmtInfoV1 : public TowerInfo
{
 public:
  BbcPmtInfoV1() {}
  //BbcPmtInfoV1(TowerInfo& tower);
  ~BbcPmtInfoV1() override = default;
  void Reset() override;

  //! Clear is used by TClonesArray to reset the tower to initial state without calling destructor/constructor
  void Clear(Option_t* = "") override;

  //! Identify object
  void identify(std::ostream& os = std::cout) const override;

  //! isValid returns non zero if object contains vailid data
  virtual int isValid() const override { if ( bq == NAN ) return 0; return 1; }

  /*
  void set_time(short t) override;
  short get_time() override { return _time; }
  void set_energy(float energy) override { _energy = energy; }
  float get_energy() override { return _energy; }
  */

  Short_t get_pmt() { return bpmt; }
  Float_t get_q()   { return bq; }
  Float_t get_t()   { return btq; }
  Float_t get_tt()  { return btt; }
  Float_t get_tq()  { return btq; }

  void set_pmt(const Short_t pmt) { bpmt = pmt; }
  void set_q(const Float_t q) { bq = q; }
  void set_tt(const Float_t t) { btt = t; }
  void set_tq(const Float_t t) { btq = t; }
  void set_pmt(const Short_t pmt, const Float_t q, const Float_t tt, const Float_t tq) {
    bpmt = pmt;
    bq = q;
    btt = tt;
    btq = tq;
  }

 private:
  Short_t bpmt {-1};
  Float_t bq   {NAN};
  Float_t btt  {NAN};
  Float_t btq  {NAN};

  ClassDefOverride(BbcPmtInfoV1, 1);
};

#endif  // __BBCPMTINFOV1_H__
