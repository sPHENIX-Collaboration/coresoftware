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

  virtual const std::array<double, 120> get_bunchnumber_lumi() const override { return m_bunchnumber_lumi; }
  virtual void set_bunchnumber_lumi(const std::array<double, 120>& vals) override { m_bunchnumber_lumi = vals; }

  virtual double get_lumi() const override { return m_lumi; }
  virtual void set_lumi(double val) override { m_lumi = val; }


 private:
 // lumi values saved in pb^-1
  std::array<double, 120> m_bunchnumber_lumi{0.};

  double m_lumi{0.};



  ClassDefOverride(StreamingLumiInfov1, 1)
};

#endif
