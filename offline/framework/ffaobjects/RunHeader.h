// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_RUNHEADER_H
#define FFAOBJECTS_RUNHEADER_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <string>  // for string

///
class RunHeader : public PHObject
{
 public:
  /// dtor
  ~RunHeader() override {}

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  /// get Run Number
  virtual int get_RunNumber() const;
  /// set Run Number
  virtual void set_RunNumber(const int run);

  virtual void set_floatval(const std::string & /*name*/, const float /*fval*/) { return; }
  virtual float get_floatval(const std::string & /*name*/) const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void set_intval(const std::string & /*name*/, const int /*ival*/) { return; }
  virtual int get_intval(const std::string & /*name*/) const { return std::numeric_limits<int>::min(); }

  virtual void set_bor_timestamp(const int time ) {set_intval("BORTIME",time);}
  virtual int get_bor_timestamp() const {return get_intval("BORTIME");}

  virtual void set_eor_timestamp(const int time ) {set_intval("EORTIME",time);}
  virtual int get_eor_timestamp() const {return get_intval("EORTIME");}
  /// switches off the pesky virtual warning messages
  void NoWarning(const int i = 1);

 private:
  void warning(const std::string &func) const;

  ClassDefOverride(RunHeader, 1)
};

#endif
