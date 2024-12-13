// @file RDH.h
// @brief Classes for interpretation of MVTX Raw Data Header (FELIX + RU Header)
// @sa <O2/DataFormats/Headers/include/Headers/RAWDataHeader.h>

#ifndef MVTXDECODER_RDH_H
#define MVTXDECODER_RDH_H

#include <cstdint>
#include <type_traits>

namespace mvtx {

struct RDHv8 {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
  union {
    // default value
    uint64_t word0 = 0x0;
    struct {
      uint64_t flxHdrCntr : 32;        ///
      uint64_t zero0 : 32;
    };
  };
  union {
    uint64_t word1 = 0x0;
    struct {
      uint64_t zero1 :  64;
    };
  };
  union{
    uint64_t word2 = 0x0;
    struct {
      uint64_t zero2 : 56;
      uint64_t flxId : 8;
    };
  };
  union {
    uint64_t word3 = 0xab01200000000000;
    struct {
      uint64_t zero3 : 8;
      uint64_t pageSize : 12;
      uint64_t zero4 : 12;
      uint64_t gbtLink : 5;
      uint64_t zero5: 3;
      uint64_t flxHdrSize: 8;
      uint64_t flxHdrVersion : 8;
      uint64_t flxHdrCode : 8;
    };
  };
  ///
  union{
    struct __attribute__((packed)) {
      uint64_t word4 = 0x000000ffff2008;
      uint16_t word5 = 0x0;
    };
    struct __attribute__((packed)) {
      uint64_t rdhVersion : 8;         /// bit  0 to  7: RDH version
      uint64_t rdhSize : 8;            /// bit  8 to 15: RDH size
      uint64_t feeId : 16;             /// bit 16 to 31: FEE identifier
      uint64_t sourceId: 8;            /// bit 32 to 39: Source Id
      uint64_t detectorField : 32;     /// bit 40 to 71: Detector Id
      uint64_t zero6 : 8;
    };
  };
  union {
    struct __attribute__((packed)) {
      uint64_t word6 = 0x0;
      uint16_t word7 = 0x0;
    };
    struct __attribute__((packed)) {
      uint64_t bc : 12;
      uint64_t zero7 : 20;
      uint64_t orbit : 40;
      uint64_t zero8 : 8;
    };
  };
  union {
    struct __attribute__((packed)) {
      uint64_t word8 = 0x0;
      uint32_t word9 = 0x0;
    };
    struct __attribute__((packed)) {
      uint64_t trgType : 32;
      uint64_t packetCounter : 16;
      uint64_t stopBit : 8;
      uint64_t priority : 8;
      uint64_t zero9 : 16;
      uint64_t rdhGBTcounter : 16;
    };
  };
}; // RDH
#pragma GCC diagnostic pop

struct RDHAny {

  uint64_t word0 = 0x0;
  uint64_t word1 = 0x0;
  uint64_t word2 = 0x0;
  uint64_t word3 = 0x0;
  uint64_t word4 = 0x0;
  uint64_t word5 = 0x0;
  uint64_t word6 = 0x0;
  uint64_t word7 = 0x0;

  RDHAny(int v = 0); // 0 for default version

  template <typename H>
  RDHAny(const H& rdh);

  template <typename H>
  RDHAny& operator=(const H& rdh);

  //------------------ service methods
  using RDHv8 = mvtx::RDHv8; // V8

  /// make sure we RDH is a legitimate RAWDataHeader
  template <typename rdh>
  static constexpr void sanityCheckStrict()
  {
    static_assert(std::is_same<rdh, RDHv8>::value,
                  "not an RDH");
  }

  /// make sure we RDH is a legitimate RAWDataHeader or generic RDHAny placeholder
  template <typename rdh>
  static constexpr void sanityCheckLoose()
  {
    static_assert(std::is_same<rdh, RDHv8>::value || std::is_same<rdh, RDHAny>::value,
                  "not an RDH or RDHAny");
  }

  template <typename H>
  static const void* voidify(const H& rdh)
  {
    sanityCheckLoose<H>();
    return reinterpret_cast<const void*>(&rdh);
  }

  template <typename H>
  static void* voidify(H& rdh)
  {
    sanityCheckLoose<H>();
    return reinterpret_cast<void*>(&rdh);
  }

  const void* voidify() const { return voidify(*this); }
  void* voidify() { return voidify(*this); }

  template <typename H>
  H* as_ptr()
  {
    sanityCheckLoose<H>();
    return reinterpret_cast<H*>(this);
  }

  template <typename H>
  H& as_ref()
  {
    sanityCheckLoose<H>();
    return reinterpret_cast<H&>(this);
  }

 protected:
  void copyFrom(const void* rdh);
};

using RDH = RDHv8;

struct RDHUtils {

// dereference SRC pointer as DST type reference
#define TOREF(DST, SRC) *reinterpret_cast<DST*>(SRC)
// dereference SRC pointer as DST type const reference
#define TOCREF(DST, SRC) *reinterpret_cast<const DST*>(SRC)

  using RDHDef = mvtx::RDH; // wathever is default
  using RDHAny = mvtx::RDHAny;

  ///_______________________________
  template <typename H>
  static constexpr int getVersion()
  {
    RDHAny::sanityCheckStrict<H>();
    if (std::is_same<H, RDHv8>::value) {
      return 8;
    }
    return -1; // dummy value as this method will be used on the CPU only
  }

  ///_______________________________
  template <typename H>
  static uint8_t getVersion(const H& rdh)
  {
    return rdh.rdhVersion;
  } // same for all
  static uint8_t getVersion(const RDHAny& rdh) { return getVersion(rdh.voidify()); }
  static uint8_t getVersion(const void* rdhP) { return getVersion(TOCREF(RDHDef, rdhP)); }

  ///_______________________________
  static void printRDH(const RDHv8& rdh);
  static void printRDH(const RDHAny& rdh) { printRDH(rdh.voidify()); }
  static void printRDH(const void* rdhP);

  ///_______________________________
  static bool checkRDH(const RDHv8& rdh, bool verbose = true, bool checkZeros = false);
  static bool checkRDH(const RDHAny& rdh, bool verbose = true, bool checkZeros = false)
  {
     return checkRDH(rdh.voidify(), verbose, checkZeros);
  }
  static bool checkRDH(const void* rdhP, bool verbose = true, bool checkZeros = false);

  ///_______________________________
  template <typename H>
  static void dumpRDH(const H& rdh)
  {
    dumpRDH(reinterpret_cast<const void*>(&rdh));
  }
  static void dumpRDH(const void* rdhP);
};

///_________________________________
/// create from arbitrary RDH version
template <typename H>
inline RDHAny::RDHAny(const H& rdh)
{
  sanityCheckLoose<H>();
  copyFrom(&rdh);
}

///_________________________________
/// copy from arbitrary RDH version
template <typename H>
inline RDHAny& RDHAny::operator=(const H& rdh)
{
  sanityCheckLoose<H>();
  if (this != voidify(rdh)) {
    copyFrom(&rdh);
  }
  return *this;
}

}  // namespace mvtx

#endif // MVTXDECODER_RDH_H
