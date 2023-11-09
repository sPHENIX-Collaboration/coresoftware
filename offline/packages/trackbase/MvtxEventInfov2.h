// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTXEVENTINFOV2_H
#define MVTXEVENTINFOV2_H

     /***************************/
    /* MVTX event header class */
   /*       Cameron Dean      */
  /*   MIT (ctdean@mit.edu)  */
 /*          09/11/2023     */
/***************************/

#include <iostream>

#include <set>

#include "MvtxEventInfov1.h"

typedef std::pair<uint64_t, uint64_t> strobe_L1_pair;

///
class MvtxEventInfov2 : public MvtxEventInfov1 
{
 public:
  MvtxEventInfov2() = default;

  /// dtor
  virtual ~MvtxEventInfov2() = default;

  PHObject *CloneMe() const override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  void set_strobe_BCO(const uint64_t strobe_BCO) { m_strobe_BCO = strobe_BCO; }

  uint64_t get_strobe_BCO() {return  m_strobe_BCO; }

  /// switches off the pesky virtual warning messages
  void NoWarning(const int i = 1);

 private:
  void warning(const std::string &func) const;

  uint64_t m_strobe_BCO = 0;

  //ClassDefOverride(MvtxEventInfov2, 1)
};

#endif
