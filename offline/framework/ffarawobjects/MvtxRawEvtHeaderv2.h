#ifndef FUN4ALLRAW_MVTXRAWEVTHEADERV2_H
#define FUN4ALLRAW_MVTXRAWEVTHEADERV2_H

#include "MvtxRawEvtHeader.h"

#include <phool/PHObject.h>

class MvtxFeeIdInfo;
class TCloneArray;

class MvtxRawEvtHeaderv2 : public MvtxRawEvtHeader
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

  MvtxFeeIdInfo *AddFeeIdInfo() override;
  MvtxFeeIdInfo *AddFeeIdInfo(const MvtxFeeIdInfo *feeIdInfo) override;

  void AddL1Trg(const uint64_t& gtmL1_bco) override { m_MvtxL1TrgSet.insert(gtmL1_bco); };
  void AddL1Trg(const std::set<uint64_t>& mvtxL1TrgSet) override;

  std::set<uint64_t>& getMvtxLvL1BCO() override { return m_MvtxL1TrgSet; };

 private:
  TCloneArray *MvtxFeeIdInfoTCArray = nullptr;

  std::set<uint64_t> m_MvtxL1TrgSet;

  ClassDefOverride(MvtxRawEvtHeaderv2, 1)
};

#endif
