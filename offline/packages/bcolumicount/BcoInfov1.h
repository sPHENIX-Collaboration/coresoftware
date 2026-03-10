// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BCOLUMICOUNT_BCOINFOV1_H
#define BCOLUMICOUNT_BCOINFOV1_H

#include "BcoInfo.h"

#include <array>
#include <iostream>

class BcoInfov1 : public BcoInfo
{
 public:
  /// ctor
  BcoInfov1() = default;

  /// dtor
  ~BcoInfov1() override = default;

  ///  Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& out = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  uint64_t get_previous_bco() const override { return bco[0]; }
  uint64_t get_current_bco() const override { return bco[1]; }
  uint64_t get_future_bco() const override { return bco[2]; }

  void set_previous_bco(uint64_t val) override { bco[0] = val; }
  void set_current_bco(uint64_t val) override { bco[1] = val; }
  void set_future_bco(uint64_t val) override { bco[2] = val; }

  int get_previous_evtno() const override { return evtno[0]; }
  int get_current_evtno() const override { return evtno[1]; }
  int get_future_evtno() const override { return evtno[2]; }

  void set_previous_evtno(int val) override { evtno[0] = val; }
  void set_current_evtno(int val) override { evtno[1] = val; }
  void set_future_evtno(int val) override { evtno[2] = val; }

 private:
  std::array<uint64_t, 3> bco{0};
  std::array<int, 3> evtno{0};

  ClassDefOverride(BcoInfov1, 1)
};

#endif
