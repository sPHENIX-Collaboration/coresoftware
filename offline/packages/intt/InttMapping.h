#ifndef INTT_MAPPING_H
#define INTT_MAPPING_H

#include <cstdint>
#include <iostream>
#include <map>

class InttRawHit;

namespace InttNameSpace
{
  /// RawData_s (kept as aggregate for convenience)
  struct RawData_s
  {
    int felix_server = 0;
    int felix_channel = 0;
    int chip = 0;
    int channel = 0;

	using value_type = RawData_s;
	using pointer = RawData_s*;
	using reference = RawData_s&;

	RawData_s& operator*() { return *this; }
	RawData_s& operator++();

    friend bool operator<(InttNameSpace::RawData_s const&, InttNameSpace::RawData_s const&);
    friend std::ostream& operator<<(std::ostream&, InttNameSpace::RawData_s const&);
  };
  bool operator<(RawData_s const&, RawData_s const&);
  std::ostream& operator<<(std::ostream&, InttNameSpace::RawData_s const&);
  class AllRawDataChannels /// For range-based for loops
  {
  public:
	  static const RawData_s begin() { return {.felix_server = 0, .felix_channel = 0, .chip = 0, .channel = 0}; }
	  static const RawData_s end() { return {.felix_server = 8, .felix_channel = 0, .chip = 0, .channel = 0}; }
  };

  /// Online_s (kept as aggregate for convenience)
  struct Online_s
  {
    int lyr = 0;
    int ldr = 0;
    int arm = 0;
    int chp = 0;
    int chn = 0;

	using value_type = Online_s;
	using pointer = Online_s*;
	using reference = Online_s&;

	Online_s& operator*() { return *this; }
	Online_s& operator++();

    friend bool operator<(InttNameSpace::Online_s const&, InttNameSpace::Online_s const&);
    friend std::ostream& operator<<(std::ostream&, InttNameSpace::Online_s const&);
  };
  bool operator<(Online_s const&, Online_s const&);
  std::ostream& operator<<(std::ostream&, InttNameSpace::Online_s const&);
  class AllOnlineChannels /// For range-based for loops
  {
  public:
	  static const Online_s begin() { return {.lyr = 0, .ldr = 0, .arm = 0, .chp = 0, .chn = 0}; }
	  static const Online_s end() { return {.lyr = 4, .ldr = 0, .arm = 0, .chp = 0, .chn = 0}; }
  };

  /// Offline_s (kept as aggregate for convenience)
  struct Offline_s
  {
    int layer = 0;
    int ladder_phi = 0;
    int ladder_z = 0;
    int strip_x = 0; // row, global phi
    int strip_y = 0; // col, global z

	using value_type = Offline_s;
	using pointer = Offline_s*;
	using reference = Offline_s&;

	Offline_s& operator*() { return *this; }
	Offline_s& operator++();

    friend bool operator<(InttNameSpace::Offline_s const&, InttNameSpace::Offline_s const&);
    friend std::ostream& operator<<(std::ostream&, InttNameSpace::Online_s const&);
  };
  bool operator<(Offline_s const&, Offline_s const&);
  std::ostream& operator<<(std::ostream&, InttNameSpace::Offline_s const&);
  class AllOfflineChannels /// For range-based for loops
  {
  public:
	  static const Offline_s begin() { return {.layer = 3, .ladder_phi = 0, .ladder_z = 0, .strip_x = 0, .strip_y = 0}; }
	  static const Offline_s end() { return {.layer = 7, .ladder_phi = 0, .ladder_z = 0, .strip_x = 0, .strip_y = 0}; }
  };

  /// Methods
  RawData_s RawFromHit(InttRawHit*);

  Online_s ToOnline(Offline_s const&);
  Offline_s ToOffline(Online_s const&);

  RawData_s ToRawData(Online_s const&);
  Online_s ToOnline(RawData_s const&);

  inline RawData_s ToRawData(Offline_s const& offline) { return ToRawData(ToOnline(offline)); }
  inline Offline_s ToOffline(RawData_s const& rawdata) { return ToOffline(ToOnline(rawdata)); }
};

/// Relational operators with macro helper
#define DEFINE_RELATIONAL_OPERATORS(ClassName) \
	inline bool operator>(InttNameSpace::ClassName const& lhs, InttNameSpace::ClassName const& rhs) { return rhs < lhs; } \
	inline bool operator<=(InttNameSpace::ClassName const& lhs, InttNameSpace::ClassName const& rhs) { return !(rhs < lhs); } \
	inline bool operator>=(InttNameSpace::ClassName const& lhs, InttNameSpace::ClassName const& rhs) { return !(lhs < rhs); } \
	inline bool operator!=(InttNameSpace::ClassName const& lhs, InttNameSpace::ClassName const& rhs) { return (lhs < rhs) || (rhs < lhs); } \
	inline bool operator==(InttNameSpace::ClassName const& lhs, InttNameSpace::ClassName const& rhs) { return !(lhs < rhs) && !(rhs < lhs); }
DEFINE_RELATIONAL_OPERATORS(RawData_s)
DEFINE_RELATIONAL_OPERATORS(Online_s)
DEFINE_RELATIONAL_OPERATORS(Offline_s)
#undef DEFINE_RELATIONAL_OPERATORS

#endif  // INTT_MAPPING_H

//      Ladder Structure                //
//======================================//
//      U14     U1      Type B  North   //
//      U15     U2        .       .     //
//      U16     U3        .       .     //
//      U17     U4        .       .     //
//      U18     U5      Type B    .     //
//------------------------------  .     //
//      U19     U6      Type A    .     //
//      U20     U7        .       .     //
//      U21     U8        .       .     //
//      U22     U9        .       .     //
//      U23     U10       .       .     //
//      U24     U11       .       .     //
//      U25     U12       .       .     //
//      U26     U13     Type A  North   //
//--------------------------------------//
//      U13     U26     Type A  South   //
//      U12     U25       .       .     //
//      U11     U24       .       .     //
//      U10     U23       .       .     //
//      U9      U22       .       .     //
//      U8      U21       .       .     //
//      U7      U20       .       .     //
//      U6      U19     Type A    .     //
//------------------------------  .     //
//      U5      U18     Type B    .     //
//      U4      U17       .       .     //
//      U3      U16       .       .     //
//      U2      U15       .       .     //
//      U1      U14     Type B  South   //
//======================================//

//|<-            South                        |                       North            ->|
//|         B      |            A             |             A           |       B        |
//+----------------+--------------------------+-------------------------+----------------+
//|         1      |            0             |             2           |       3        | ladder_z  (offline)
//|  0  1  2  3  4 |  0  1  2  3  4  5  6  7  |  0  1  2  3  4  5  6  7 |  0  1  2  3  4 | strip_y   (offline)
//+----------------+--------------------------+-------------------------+----------------+
//|  0  1  2  3  4 |  5  6  7  8  9 10 11 12  | 25 24 23 22 21 20 19 18 | 17 16 15 14 13 | chp       (online)
//| 13 14 15 16 17 | 18 19 20 21 22 23 24 25  | 12 11 10  9  8  7  6  5 |  4  3  2  1  0 | chp       (online)
