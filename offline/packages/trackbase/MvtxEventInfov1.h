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

#include <iostream>

#include <set>

#include "MvtxEventInfo.h"

typedef std::pair<uint64_t, uint64_t> strobe_L1_pair;

///
class MvtxEventInfov1 : public MvtxEventInfo 
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

  void set_number_HB(const int /*ival*/);
  int get_number_HB() const;

  void set_strobe_BCO_L1_BCO(const uint64_t strobe_BCO, const uint64_t L1_BCO);
 
  unsigned int get_number_strobes() const;
  unsigned int get_number_L1s() const;

  std::set<uint64_t> get_strobe_BCOs() const;
  std::set<uint64_t> get_L1_BCOs() const;

  std::set<uint64_t> get_strobe_BCO_from_L1_BCO(const uint64_t ival) const;
  std::set<uint64_t> get_L1_BCO_from_strobe_BCO(const uint64_t ival) const;

  /// switches off the pesky virtual warning messages
  void NoWarning(const int i = 1);

 private:
  void warning(const std::string &func) const;

  int m_number_HB = -1;
  std::set<strobe_L1_pair> m_strobe_BCO_L1_BCO;

  std::string m_number_L1_name = "Number L1";
  std::string m_number_HB_name = "Number HB";

  //ClassDefOverride(MvtxEventInfov1, 1)
};

#endif
