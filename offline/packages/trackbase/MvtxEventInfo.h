// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTXEVENTINFO_H
#define MVTXEVENTINFO_H

     /********************************/
    /* Subsystem event header class */
   /*          Cameron Dean        */
  /*      MIT (ctdean@mit.edu)    */
 /*             29/09/2023       */
/********************************/

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <string>  // for string
#include <map>

///
class MvtxEventInfo : public PHObject
{
 public:
  MvtxEventInfo() = default;

  /// dtor
  virtual ~MvtxEventInfo() = default;

  PHObject *CloneMe() const override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  void set_floatval(const std::string & /*name*/, const float /*fval*/);
  float get_floatval(const std::string & /*name*/) const;

  void set_intval(const std::string & /*name*/, const int32_t /*ival*/);
  int get_intval(const std::string & /*name*/) const;

  void set_int64val(const std::string & /*name*/, const int64_t /*ival*/);
  int get_int64val(const std::string & /*name*/) const;

  void set_uintval(const std::string & /*name*/, const uint32_t /*ival*/);
  int get_uintval(const std::string & /*name*/) const;

  void set_uint64val(const std::string & /*name*/, const uint64_t /*ival*/);
  int get_uint64val(const std::string & /*name*/) const;

  void set_stringval(const std::string & /*name*/, const std::string & /*ival*/);
  std::string get_stringval(const std::string & /*name*/) const;

  /// switches off the pesky virtual warning messages
  void NoWarning(const int i = 1);

 protected:
  std::map<std::string, int32_t> m_IntEventProperties;
  std::map<std::string, int64_t> m_Int64EventProperties;
  std::map<std::string, uint32_t> m_UintEventProperties;
  std::map<std::string, uint64_t> m_Uint64EventProperties;
  std::map<std::string, float> m_FloatEventProperties;
  std::map<std::string, std::string> m_StringEventProperties;

 private:
  void warning(const std::string &func) const;

  //ClassDefOverride(MvtxEventInfo, 1)
};

#endif
