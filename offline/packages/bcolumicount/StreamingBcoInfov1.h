// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BCOLLUMICOUNT_STREAMINGBCOINFOV1_H
#define BCOLLUMICOUNT_STREAMINGBCOINFOV1_H

#include "StreamingBcoInfo.h"


#include <iostream>
#include <limits>
#include <utility>


///
class StreamingBcoInfov1 : public StreamingBcoInfo
{
 public:
  /// ctor - daughter class copy ctor needs this
  StreamingBcoInfov1() = default;
  /// dtor
  ~StreamingBcoInfov1() override = default;
  /// Clear Sync
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  //int isValid() const override;

  virtual uint64_t get_bco() const override { return m_bco; }
  virtual void set_bco(uint64_t val) override { m_bco = val; }

  virtual int get_evtno() const override { return m_evtno; }
  virtual void set_evtno(int val) override { m_evtno = val; }

  virtual bool get_usable_bco_tag() const override { return m_usable_bco_tag; }
  virtual void set_usable_bco_tag(bool val) override { m_usable_bco_tag = val; }

  virtual std::pair<uint64_t, uint64_t> get_bco_streaming_window() const override { return m_bco_streaming_window; }
  virtual void set_bco_streaming_window(std::pair<uint64_t, uint64_t> val) override { m_bco_streaming_window = val; }


 private:
  uint64_t m_bco{0};
  int m_evtno{0};
  bool m_usable_bco_tag{false};
  std::pair<uint64_t, uint64_t> m_bco_streaming_window;

  ClassDefOverride(StreamingBcoInfov1, 1)
};

#endif
