#ifndef FUN4ALLRAW_MVTXRAWEVTHEADERV1_H
#define FUN4ALLRAW_MVTXRAWEVTHEADERV1_H

#include "MvtxRawEvtHeader.h"

#include <phool/PHObject.h>

#include <cstdint>
#include <set>

class MvtxRawEvtHeaderv1 : public MvtxRawEvtHeader
{
 public:
  MvtxRawEvtHeaderv1() = default;
  ~MvtxRawEvtHeaderv1() = default;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  void AddFeeId(const int& feeid) override { m_MvtxFeeIdSet.insert(feeid); };
  void AddL1Trg(const uint64_t& gtmL1_bco) override { m_MvtxL1TrgSet.insert(gtmL1_bco); };

  void AddFeeId(const std::set<uint16_t>& mvtxFeeIds) override;
  void AddL1Trg(const std::set<uint64_t>& mvtxL1TrgSet) override;

  std::set<uint16_t>& getMvtxFeeIdSet() override { return m_MvtxFeeIdSet; };
  std::set<uint64_t>& getMvtxLvL1BCO() override { return m_MvtxL1TrgSet; };

 private:
  std::set<uint16_t> m_MvtxFeeIdSet;
  std::set<uint64_t> m_MvtxL1TrgSet;

  ClassDefOverride(MvtxRawEvtHeaderv1, 1)
};

#endif
