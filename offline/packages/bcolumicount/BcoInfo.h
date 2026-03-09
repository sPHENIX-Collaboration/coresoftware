// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BCOLLUMICOUNT_BCOINFO_H
#define BCOLLUMICOUNT_BCOINFO_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>

///
class BcoInfo : public PHObject
{
 public:
  /// ctor - daughter class copy ctor needs this
  BcoInfo() = default;
  /// dtor
  ~BcoInfo() override = default;
  /// Clear Sync
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

    uint64_t get_previous_bco() const {return 0;}
  uint64_t get_current_bco() const {return 0;}
  uint64_t get_future_bco() const {return 0;}

  void set_previous_bco(uint64_t /*val*/) {return;}
  void set_current_bco(uint64_t /*val*/) {return;}
  void set_future_bco(uint64_t /*val*/) {return;}


 private:

  ClassDefOverride(BcoInfo, 1)
};

#endif
