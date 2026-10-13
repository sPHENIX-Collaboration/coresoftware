// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BCOLLUMICOUNT_STREAMINGLUMIINFO_H
#define BCOLLUMICOUNT_STREAMINGLUMIINFO_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>
#include <utility>
#include <array>


///
class StreamingLumiInfo : public PHObject
{
 public:
  /// ctor - daughter class copy ctor needs this
  StreamingLumiInfo() = default;
  /// dtor
  ~StreamingLumiInfo() override = default;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  //int isValid() const override;

  virtual const std::array<double, 120> get_bunchnumber_lumi() const { return std::array<double, 120>{}; }
  virtual void set_bunchnumber_lumi(const std::array<double, 120>& /*vals*/) { return; }

  virtual double get_lumi() const { return 0; }
  virtual void set_lumi(double /*val*/) { return; }

 private:
  ClassDefOverride(StreamingLumiInfo, 1)
};

#endif
