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

  virtual const std::array<double, 120> get_bunchnumber_lumi_raw() const { return std::array<double, 120>{}; }
  virtual void set_bunchnumber_lumi_raw(const std::array<double, 120>& /*vals*/) { return; }

  virtual const std::array<double, 120> get_bunchnumber_lumi_live() const { return std::array<double, 120>{}; }
  virtual void set_bunchnumber_lumi_live(const std::array<double, 120>& /*vals*/) { return; }

  virtual const std::array<double, 120> get_bunchnumber_lumi_scaled() const { return std::array<double, 120>{}; }
  virtual void set_bunchnumber_lumi_scaled(const std::array<double, 120>& /*vals*/) { return; }

  virtual double get_lumi_raw() const { return 0; }
  virtual void set_lumi_raw(double /*val*/) { return; }

  virtual double get_lumi_live() const { return 0; }
  virtual void set_lumi_live(double /*val*/) { return; }

  virtual double get_lumi_scaled() const { return 0; }
  virtual void set_lumi_scaled(double /*val*/) { return; }


 private:
  ClassDefOverride(StreamingLumiInfo, 1)
};

#endif
