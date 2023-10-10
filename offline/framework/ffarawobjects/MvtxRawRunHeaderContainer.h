#ifndef FUN4ALLRAW_MVTXRAWRUNHEADERCONTAINER_H
#define FUN4ALLRAW_MVTXRAWRUNHEADERCONTAINER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <set>

class  MvtxRawRunHeaderContainer: public PHObject
{

public:
  MvtxRawRunHeaderContainer() = default;
  virtual ~MvtxRawRunHeaderContainer() = default;

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

private:
  std::set<uint64_t> m_MvtxL1TrgSet;

  ClassDefOverride(MvtxRawRunHeaderContainer,1)
};

#endif
