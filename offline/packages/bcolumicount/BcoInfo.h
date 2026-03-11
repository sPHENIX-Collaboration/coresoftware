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

  virtual uint64_t get_previous_bco() const { return 0; }
  virtual uint64_t get_current_bco() const { return 0; }
  virtual uint64_t get_future_bco() const { return 0; }

  virtual void set_previous_bco(uint64_t /*val*/) { return; }
  virtual void set_current_bco(uint64_t /*val*/) { return; }
  virtual void set_future_bco(uint64_t /*val*/) { return; }

  virtual int get_previous_evtno() const { return 0; }
  virtual int get_current_evtno() const { return 0; }
  virtual int get_future_evtno() const { return 0; }

  virtual void set_previous_evtno(int /*val*/) { return; }
  virtual void set_current_evtno(int /*val*/) { return; }
  virtual void set_future_evtno(int /*val*/) { return; }

 private:
  ClassDefOverride(BcoInfo, 1)
};

#endif
