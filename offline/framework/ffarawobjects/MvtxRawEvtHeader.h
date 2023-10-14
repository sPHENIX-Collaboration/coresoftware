#ifndef FUN4ALLRAW_MVTXRAWEVTHEADER_H
#define FUN4ALLRAW_MVTXRAWEVTHEADER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <set>

class  MvtxRawEvtHeader: public PHObject
{

public:
  MvtxRawEvtHeader() = default;
  virtual ~MvtxRawEvtHeader() = default;

  ///Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  void AddFeeId(const int& feeid)
  {
    m_MvtxFeeIdSet.insert(feeid);
  };
  void AddFeeId(const std::set<int>& mvtxFeeIds);

  std::set<int>& getMvtxFeeIdSet() { return m_MvtxFeeIdSet; };

private:
  std::set<int> m_MvtxFeeIdSet;

  ClassDefOverride(MvtxRawEvtHeader,1)
};

#endif
