#ifndef INTT_MAP_H
#define INTT_MAP_H

#include <phool/PHObject.h>

class InttMap : public PHObject
{
 public:
  typedef int field_t;
  field_t static const Wildcard = 0xffff;

  struct Online_s
  {
    field_t lyr = 0;
    field_t ldr = 0;
    field_t arm = 0;
    field_t chp = 0;
    field_t chn = 0;
  };

  struct RawData_s
  {
    field_t pid = 0;
    field_t fee = 0;
    field_t chp = 0;
    field_t chn = 0;
  };

  struct Offline_s
  {
    field_t layer = 0;
    field_t ladder_phi = 0;
    field_t ladder_z = 0;
    field_t strip_phi = 0;
    field_t strip_z = 0;
  };

  struct OnlineComparator
  {
    bool operator()(struct Online_s const&, struct Online_s const&) const;
  };

  struct OnlineWildcardComparator
  {
    bool operator()(struct Online_s const&, struct Online_s const&) const;
  };

  struct RawDataComparator
  {
    bool operator()(struct RawData_s const&, struct RawData_s const&) const;
  };

  struct RawDataWildcardComparator
  {
    bool operator()(struct RawData_s const&, struct RawData_s const&) const;
  };

  struct OfflineComparator
  {
    bool operator()(struct Offline_s const&, struct Offline_s const&) const;
  };

  struct OfflineWildcardComparator
  {
    bool operator()(struct Offline_s const&, struct Offline_s const&) const;
  };

 private:
  ClassDefOverride(InttMap, 1);
};

#endif  // INTT_MAP_H

// |<-            South                        |                       North            ->|
// |         B      |            A             |             A           |       B        |
// +----------------+--------------------------+-------------------------+----------------+
// |         1      |            0             |             2           |       3        | ladder_z  (offline)
// |  0  1  2  3  4 |  0  1  2  3  4  5  6  7  |  0  1  2  3  4  5  6  7 |  0  1  2  3  4 | strip_z   (offline)
// +----------------+--------------------------+-------------------------+----------------+
// |  0  1  2  3  4 |  5  6  7  8  9 10 11 12  | 25 24 23 22 21 20 19 18 | 17 16 15 14 13 | chp       (online)
// | 13 14 15 16 17 | 18 19 20 21 22 23 24 25  | 12 11 10  9  8  7  6  5 |  4  3  2  1  0 | chp       (online)
