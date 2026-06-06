// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BCOLLUMICOUNT_STREAMINGBCOINFO_H
#define BCOLLUMICOUNT_STREAMINGBCOINFO_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>
#include <utility>


///
class StreamingBcoInfo : public PHObject
{
 public:
  /// ctor - daughter class copy ctor needs this
  StreamingBcoInfo() = default;
  /// dtor
  ~StreamingBcoInfo() override = default;
  /// Clear Sync
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  //int isValid() const override;

  virtual uint64_t get_bco() const { return 0; }

  virtual void set_bco(uint64_t /*val*/) { return; }

  virtual int get_evtno() const { return 0; }

  virtual void set_evtno(int /*val*/) { return; }

  virtual bool get_usable_bco_tag() const { return 0; }

  virtual void set_usable_bco_tag(bool /*val*/) { return; }

  virtual std::pair<uint64_t, uint64_t> get_bco_streaming_window() const { return std::make_pair(0, 0); }

  virtual void set_bco_streaming_window(std::pair<uint64_t, uint64_t> /*val*/) { return; }


 private:
  ClassDefOverride(StreamingBcoInfo, 1)
};

#endif
