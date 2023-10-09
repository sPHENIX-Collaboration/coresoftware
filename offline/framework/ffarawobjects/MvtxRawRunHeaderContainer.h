#ifndef FUN4ALLRAW_MVTXRAWRUNHEADERCONTAINER_H
#define FUN4ALLRAW_MVTXRAWRUNHEADERCONTAINER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <set>

struct MvtxL1Trg
{
  uint64_t bco : 40; // 40 bits gl1/gtm bco
  uint16_t bc  : 12; // 12 bits mvtx internal bc counter

  MvtxL1Trg() = default;
  ~MvtxL1Trg() = default;

  MvtxL1Trg(uint64_t _bco, uint16_t _bc) : bco(_bco), bc(_bc) {};

  bool operator<(const MvtxL1Trg& other) const
  {
    return (bco == other.bco) ? (bc < other.bc) : (bco < other.bco);
  }

  bool operator==(const MvtxL1Trg& other) const
  {
    return (bco == other.bco) ? (bc == other.bc) : false;
  }

};

class  MvtxRawRunHeaderContainer: public PHObject
{
  using l1_type_set = std::set<MvtxL1Trg>;

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

  void AddL1Trg(const MvtxL1Trg& gtmL1Trg)
  {
    m_MvtxL1TrgSet.insert(gtmL1Trg);
  };
  void AddL1Trg(const l1_type_set& mvtxL1TrgSet);

private:
  l1_type_set m_MvtxL1TrgSet;

  ClassDefOverride(MvtxRawRunHeaderContainer,1)
};

#endif
