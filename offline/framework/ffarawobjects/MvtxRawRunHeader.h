#ifndef FUN4ALLRAW_MVTXRAWRUNHEADER_H
#define FUN4ALLRAW_MVTXRAWRUNHEADER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <set>

class  MvtxRawRunHeader: public PHObject
{

public:
  MvtxRawRunHeader() = default;
  virtual ~MvtxRawRunHeader() = default;

  ///Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  void AddL1Trg(const uint64_t& gtmL1_bco)
  {
    m_MvtxL1TrgSet.insert(gtmL1_bco);
  };
  void AddL1Trg(const std::set<uint64_t>& mvtxL1TrgSet);

  std::set<uint64_t>& getMvtxLvL1BCO() { return m_MvtxL1TrgSet; };

private:
  std::set<uint64_t> m_MvtxL1TrgSet;

  ClassDefOverride(MvtxRawRunHeader,1)
};

#endif
