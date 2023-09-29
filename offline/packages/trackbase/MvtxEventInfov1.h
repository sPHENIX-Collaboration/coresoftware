// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTXEVENTINFOV1_H
#define MVTXEVENTINFOV1_H

     /***************************/
    /* MVTX event header class */
   /*       Cameron Dean      */
  /*   MIT (ctdean@mit.edu)  */
 /*          29/09/2023     */
/***************************/

//#include <phool/PHObject.h>
//
//#include <cmath>
#include <iostream>
//#include <string>  // for string
//#include <map>

#include <SubdetectorEventInfo.h>

///
class MvtxEventInfov1 : public SubdetectorEventInfo 
{
 public:
  MvtxEventInfov1() = default;

  /// dtor
  virtual ~MvtxEventInfov1() = default;

  PHObject *CloneMe() const override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  void set_number_LL1(const int /*ival*/);
  int get_number_LL1() const;

  void set_strobe_BCO(const uint64_t /*ival*/);
  uint64_t get_strobe_BCO(const uint32_t /*ival*/) const;
  uint32_t get_number_strobe_BCO() const;

  /// switches off the pesky virtual warning messages
  void NoWarning(const int i = 1);

 private:
  void warning(const std::string &func) const;

  int m_number_LL1 = -1;
  std::map<uint32_t, uint64_t> m_strobe_BCO;

  std::string m_number_LL1_name = "Number LL1";

  //ClassDefOverride(MvtxEventInfov1, 1)
};

#endif
