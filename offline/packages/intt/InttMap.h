#ifndef INTT_MAP_H
#define INTT_MAP_H

#include <iostream>
#include <phool/PHObject.h>

class InttMap : public PHObject
{
 public:
  typedef int field_t;
  field_t static const Wildcard = 0xffff;

  struct Online_s
  {
    field_t lyr = Wildcard;
    field_t ldr = Wildcard;
    field_t arm = Wildcard;
    field_t chp = Wildcard;
    field_t chn = Wildcard;

	Online_s& operator++();
	Online_s& operator--();
  };
  static bool IsValid(struct Online_s const&);
  static bool IsValidOrWildcard(struct Online_s const&);

  const static struct Online_s OnlineBegin;
  const static struct Online_s OnlineEnd;
  const static struct Online_s OnlineRBegin;
  const static struct Online_s OnlineREnd;

  struct RawData_s
  {
    field_t pid = Wildcard;
    field_t fee = Wildcard;
    field_t chp = Wildcard;
    field_t chn = Wildcard;

	RawData_s& operator++();
	RawData_s& operator--();
  };
  static bool IsValid(struct RawData_s const&);
  static bool IsValidOrWildcard(struct RawData_s const&);

  const static struct RawData_s RawDataBegin;
  const static struct RawData_s RawDataEnd;
  const static struct RawData_s RawDataRBegin;
  const static struct RawData_s RawDataREnd;

  struct Offline_s
  {
    field_t layer = Wildcard;
    field_t ladder_phi = Wildcard;
    field_t ladder_z = Wildcard;
    field_t strip_z = Wildcard;
    field_t strip_phi = Wildcard;

	Offline_s& operator++();
	Offline_s& operator--();
  };
  static bool IsValid(struct Offline_s const&);
  static bool IsValidOrWildcard(struct Offline_s const&);

  const static struct Offline_s OfflineBegin;
  const static struct Offline_s OfflineEnd;
  const static struct Offline_s OfflineRBegin;
  const static struct Offline_s OfflineREnd;

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

bool operator<(struct InttMap::Online_s const&, InttMap::Online_s const&);
bool operator==(struct InttMap::Online_s const&, InttMap::Online_s const&);
bool operator!=(struct InttMap::Online_s const&, InttMap::Online_s const&);
bool operator>(struct InttMap::Online_s const&, InttMap::Online_s const&);
bool operator<=(struct InttMap::Online_s const&, InttMap::Online_s const&);
bool operator>=(struct InttMap::Online_s const&, InttMap::Online_s const&);

bool operator<(struct InttMap::RawData_s const&, InttMap::RawData_s const&);
bool operator==(struct InttMap::RawData_s const&, InttMap::RawData_s const&);
bool operator!=(struct InttMap::RawData_s const&, InttMap::RawData_s const&);
bool operator>(struct InttMap::RawData_s const&, InttMap::RawData_s const&);
bool operator<=(struct InttMap::RawData_s const&, InttMap::RawData_s const&);
bool operator>=(struct InttMap::RawData_s const&, InttMap::RawData_s const&);

bool operator<(struct InttMap::Offline_s const&, InttMap::Offline_s const&);
bool operator==(struct InttMap::Offline_s const&, InttMap::Offline_s const&);
bool operator!=(struct InttMap::Offline_s const&, InttMap::Offline_s const&);
bool operator>(struct InttMap::Offline_s const&, InttMap::Offline_s const&);
bool operator<=(struct InttMap::Offline_s const&, InttMap::Offline_s const&);
bool operator>=(struct InttMap::Offline_s const&, InttMap::Offline_s const&);

std::ostream& operator<<(std::ostream&, struct InttMap::Online_s const&);
std::ostream& operator<<(std::ostream&, struct InttMap::RawData_s const&);
std::ostream& operator<<(std::ostream&, struct InttMap::Offline_s const&);

#endif  // INTT_MAP_H

// |<-            South                        |                       North            ->|
// |         B      |            A             |             A           |       B        |
// +----------------+--------------------------+-------------------------+----------------+
// |         1      |            0             |             2           |       3        | ladder_z  (offline)
// |  0  1  2  3  4 |  0  1  2  3  4  5  6  7  |  0  1  2  3  4  5  6  7 |  0  1  2  3  4 | strip_z   (offline)
// +----------------+--------------------------+-------------------------+----------------+
// |  0  1  2  3  4 |  5  6  7  8  9 10 11 12  | 25 24 23 22 21 20 19 18 | 17 16 15 14 13 | chp       (online)
// | 13 14 15 16 17 | 18 19 20 21 22 23 24 25  | 12 11 10  9  8  7  6  5 |  4  3  2  1  0 | chp       (online)
