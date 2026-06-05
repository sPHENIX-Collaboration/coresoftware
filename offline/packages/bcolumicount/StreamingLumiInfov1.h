// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BCOLLUMICOUNT_STREAMINGLUMIINFOV1_H
#define BCOLLUMICOUNT_STREAMINGLUMIINFOV1_H

#include "StreamingLumiInfo.h"


#include <iostream>
#include <limits>
#include <utility>


///
class StreamingLumiInfov1 : public StreamingLumiInfo
{
 public:
  /// ctor - daughter class copy ctor needs this
  StreamingLumiInfov1() = default;
  /// dtor
  ~StreamingLumiInfov1() override = default;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  //int isValid() const override;

  virtual double get_lumi_raw() const override { return m_lumi_raw; }
  virtual void set_lumi_raw(double val) override { m_lumi_raw = val; }

  virtual double get_lumi_live() const override { return m_lumi_live; }
  virtual void set_lumi_live(double val) override { m_lumi_live = val; }

  virtual double get_lumi_scaled() const override { return m_lumi_scaled; }
  virtual void set_lumi_scaled(double val) override { m_lumi_scaled = val; }


 private:
  double m_lumi_raw{0.};
  double m_lumi_live{0.};
  double m_lumi_scaled{0.};



  ClassDefOverride(StreamingLumiInfov1, 1)
};

#endif
