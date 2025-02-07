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

///
class MvtxEventInfov3 : public MvtxEventInfo
{
 public:
  //! ctor
  MvtxEventInfov3() = default;

  //! cp/mv ctor
  MvtxEventInfov3(const MvtxEventInfov3 &) = default;
  MvtxEventInfov3(MvtxEventInfov3 &&) = default;

  //! cp/mv asignment
  MvtxEventInfov3 &operator=(const MvtxEventInfov3 &) = default;
  MvtxEventInfov3 &operator=(MvtxEventInfov3 &&) = default;

  //! dtor
  ~MvtxEventInfov3() override = default;

  PHObject *CloneMe() const override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  unsigned int get_number_strobes() const override { return m_strobe_BCOs.size(); }
  unsigned int get_number_L1s() const override { return m_L1_BCOs.size(); }

  std::set<uint64_t> get_strobe_BCOs() const override { return m_strobe_BCOs; }
  std::set<uint64_t> get_L1_BCOs() const override { return m_L1_BCOs; }

  void add_strobe_BCO(const uint64_t &strb_val) override { m_strobe_BCOs.insert(strb_val); }
  void add_L1_BCO(const uint64_t &strb_val) override { m_L1_BCOs.insert(strb_val); }

 protected:
 private:
  void warning(const std::string &func) const;

  std::set<uint64_t> m_strobe_BCOs;
  std::set<uint64_t> m_L1_BCOs;

  ClassDefOverride(MvtxEventInfov3, 1)
};

#endif
