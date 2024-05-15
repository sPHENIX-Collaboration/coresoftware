#ifndef FUN4ALLRAW_MVTXRAWEVTHEADERV2_H
#define FUN4ALLRAW_MVTXRAWEVTHEADERV2_H

#include "MvtxRawEvtHeaderv1.h"

#include <phool/PHObject.h>

class MvtxRawEvtHeaderv2 : public MvtxRawEvtHeaderv1
{
 public:
  MvtxRawEvtHeaderv2() = default;
  ~MvtxRawEvtHeaderv2() = default;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

//  void AddFeeId(const int& feeid) override { m_MvtxFeeIdSet.insert(feeid); };
  void AddDetField(const uint16_t& feeId,
                   const uint32_t& detField) override { m_MvtxDetField[feeId] = detField; }
  void AddL1Trg(const uint64_t& gtmL1_bco) override { m_MvtxL1TrgSet.insert(gtmL1_bco); };

//  void AddFeeId(const std::set<uint16_t>& mvtxFeeIds) override;
  void AddDetField(const std::map<uint16_t, uint32_t>& ) override;
  void AddL1Trg(const std::set<uint64_t>& mvtxL1TrgSet) override;

  std::map<uint16_t, uint32_t>& getMvtxDetField() override { return m_MvtxDetField; };
  std::set<uint64_t>& getMvtxLvL1BCO() override { return m_MvtxL1TrgSet; };

 private:
  std::map<uint16_t, uint32_t> m_MvtxDetField;
  std::set<uint64_t> m_MvtxL1TrgSet;

  ClassDefOverride(MvtxRawEvtHeaderv2, 1)
};

#endif
