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

#include "MvtxEventInfo.h"

typedef std::pair<uint64_t, uint64_t> strobe_L1_pair;

///
class MvtxEventInfov2 : public MvtxEventInfo
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

  void set_number_HB(const int /*ival*/) override;
  int get_number_HB() const override;

  void set_strobe_BCO(const uint64_t strobe_BCO) override;

  void set_strobe_BCO_L1_BCO(const uint64_t strobe_BCO, const uint64_t L1_BCO) override;

  unsigned int get_number_strobes() const override;
  unsigned int get_number_L1s() const override;

  // std::set<uint64_t> get_strobe_BCOs() const;
  std::set<uint64_t> get_strobe_BCOs() const override;
  std::set<uint64_t> get_L1_BCOs() const override;

  std::set<uint64_t> get_strobe_BCO_from_L1_BCO(const uint64_t ival) const override;
  std::set<uint64_t> get_L1_BCO_from_strobe_BCO(const uint64_t ival) const override;

 protected:
  std::set<strobe_L1_pair> m_strobe_BCO_L1_BCO;
  std::set<uint64_t> m_strobe_BCOs;

  std::string m_number_L1_name = "Number L1";
  std::string m_number_HB_name = "Number HB";

 private:
  void warning(const std::string &func) const;

  int m_number_HB = -1;

  ClassDefOverride(MvtxEventInfov2, 1)
};

#endif
