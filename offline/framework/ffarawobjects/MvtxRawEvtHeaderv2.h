#ifndef FUN4ALLRAW_MVTXRAWEVTHEADERV2_H
#define FUN4ALLRAW_MVTXRAWEVTHEADERV2_H

#include "MvtxRawEvtHeader.h"

#include <phool/PHObject.h>

class MvtxFeeIdInfo;
class TClonesArray;

class MvtxRawEvtHeaderv2 : public MvtxRawEvtHeader
{
 public:
  //! ctor
  MvtxRawEvtHeaderv2();

  //! cp/mv ctor
  MvtxRawEvtHeaderv2(const MvtxRawEvtHeaderv2 &) = default;
  MvtxRawEvtHeaderv2(MvtxRawEvtHeaderv2 &&) = default;

  //! cp/mv assignment
  MvtxRawEvtHeaderv2 &operator=(const MvtxRawEvtHeaderv2 &) = default;
  MvtxRawEvtHeaderv2 &operator=(MvtxRawEvtHeaderv2 &&) = default;

  //! dtor
  ~MvtxRawEvtHeaderv2() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  MvtxFeeIdInfo *AddFeeIdInfo() override;
  MvtxFeeIdInfo *AddFeeIdInfo(MvtxFeeIdInfo *feeIdInfo) override;

  uint64_t get_nFeeIdInfo() override;
  MvtxFeeIdInfo *get_feeIdInfo(unsigned int index) override;

  void AddL1Trg(const uint64_t &gtmL1_bco) override { m_MvtxL1TrgSet.insert(gtmL1_bco); };
  void AddL1Trg(const std::set<uint64_t> &mvtxL1TrgSet) override;

  std::set<uint64_t> &getMvtxLvL1BCO() override { return m_MvtxL1TrgSet; };

 private:
  TClonesArray *m_MvtxFeeIdInfoTCArray = nullptr;

  std::set<uint64_t> m_MvtxL1TrgSet;

  ClassDefOverride(MvtxRawEvtHeaderv2, 1)
};

#endif
