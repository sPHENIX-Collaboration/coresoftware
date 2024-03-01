// @file mvtx_utils.h
// @brief Declarations of helper classes for the MVTX

#ifndef MVTXDECODER_UTILS_H
#define MVTXDECODER_UTILS_H

#include <cstdint>
#include <cassert>
#include <iostream>

namespace mvtx_utils
{

#define clean_errno() (errno == 0 ? "None" : strerror(errno)) << " "
#define log_error std::cerr << "[ERROR] (" << __FILE__ << ":" << __LINE__ << ":errno: " << clean_errno()

  constexpr uint8_t FLXWordLength = 32;

  struct RdhExt_t
  {
    // FLX header
    uint8_t flxId;    // [23]
    uint16_t pageSize; // [25]
    uint16_t gbtLink;
    uint16_t flxHdrSize;
    uint16_t flxHdrVersion;
    // RU header
    uint8_t rdhVersion;
    uint8_t rdhSize;
    uint16_t feeId;
    uint8_t sourceId;
    uint32_t detectorField;
    uint16_t bc;
    uint64_t orbit;
    uint32_t trgType;
    uint16_t packetCounter;
    uint8_t  stopBit;
    uint8_t  priority;
    uint16_t rdhGBTcounter; // 10 bits

    RdhExt_t() = default;
    ~RdhExt_t() = default;

    void decode(const uint8_t* rdh_ptr)
    {
     // FELIX header
      flxId         = *(reinterpret_cast<const uint8_t*>(rdh_ptr +  23) ) & 0xFF;
      pageSize      = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 25) ) & 0x7FF;
      gbtLink       = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 28) ) & 0x7FF;
      flxHdrSize    = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 29) ) & 0x7FF;
      flxHdrVersion = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 30) ) & 0xFFFF;
      // RU header
      rdhVersion    = *(reinterpret_cast<const uint8_t*>(rdh_ptr + 32) ) & 0xFF;
      rdhSize       = *(reinterpret_cast<const uint8_t*>(rdh_ptr + 33) ) & 0xFF;
      feeId         = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 34) ) & 0xFFFF;
      sourceId      = *(reinterpret_cast<const uint8_t*>(rdh_ptr + 36) ) & 0xFF;
      detectorField = *(reinterpret_cast<const uint32_t*>(rdh_ptr + 37) ) & 0xFFFFFFFF;
      bc            = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 42) ) & 0xFFF;
      orbit         = *(reinterpret_cast<const uint64_t*>(rdh_ptr + 46) ) & 0xFFFFFFFFFF;
      trgType       = *(reinterpret_cast<const uint32_t*>(rdh_ptr + 52) ) & 0xFFFFFFFF;
      packetCounter = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 56) ) & 0xFFFF;
      stopBit       = *(reinterpret_cast<const uint8_t*>(rdh_ptr + 58) ) & 0xFF;
      priority      = *(reinterpret_cast<const uint8_t*>(rdh_ptr + 59) ) & 0xFF;
      rdhGBTcounter = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 62) ) & 0xFFFF;
    }
  };

  template < typename A, typename B >
  bool comp(A a, B b)
  {
    return a.second < b.second;
  }

} //namespace mvtx_utils

#endif
