// @file mvtx_utils.h
// @brief Declarations of helper classes for the MVTX

#ifndef MVTXDECODER_UTILS_H
#define MVTXDECODER_UTILS_H

#include <cstdint>
#include <cassert>
#include <iostream>
#include <limits>
#include <array>
#include <memory>
#include <phool/recoConsts.h>

namespace mvtx_utils
{

#define clean_errno() (errno == 0 ? "None" : strerror(errno)) << " "
#define log_error std::cerr << "[ERROR] (" << __FILE__ << ":" << __LINE__ << ":errno: " << clean_errno()

  constexpr uint8_t FLXWordLength = 32;

  struct RdhExt_t
  {
    // FLX header
    uint8_t flxId = std::numeric_limits<uint8_t>::quiet_NaN();    // [23]
    uint16_t pageSize = std::numeric_limits<uint16_t>::quiet_NaN(); // [25]
    uint16_t gbtLink = std::numeric_limits<uint16_t>::quiet_NaN();
    uint8_t  flxHdrSize = std::numeric_limits<uint8_t>::quiet_NaN();
    uint16_t flxHdrVersion = std::numeric_limits<uint16_t>::quiet_NaN();
    // RU header
    uint8_t rdhVersion = std::numeric_limits<uint8_t>::quiet_NaN();
    uint8_t rdhSize = std::numeric_limits<uint8_t>::quiet_NaN();
    uint16_t feeId = std::numeric_limits<uint16_t>::quiet_NaN();
    uint8_t sourceId = std::numeric_limits<uint8_t>::quiet_NaN();
    uint32_t detectorField = std::numeric_limits<uint32_t>::quiet_NaN();
    uint16_t bc = std::numeric_limits<uint16_t>::quiet_NaN();
    uint64_t orbit  = std::numeric_limits<uint64_t>::quiet_NaN();
    uint32_t trgType  = std::numeric_limits<uint32_t>::quiet_NaN();
    uint16_t packetCounter = std::numeric_limits<uint16_t>::quiet_NaN();
    uint8_t  stopBit = std::numeric_limits<uint8_t>::quiet_NaN();
    uint8_t  priority = std::numeric_limits<uint8_t>::quiet_NaN();
    uint16_t rdhGBTcounter = std::numeric_limits<uint16_t>::quiet_NaN(); // 10 bits

    RdhExt_t() = default;
    ~RdhExt_t() = default;

    void decode(const uint8_t* rdh_ptr)
    {
     // FELIX header
      flxId         = *(reinterpret_cast<const uint8_t*>(rdh_ptr +  23) ) & 0xFF;
      pageSize      = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 25) ) & 0x7FF;
      gbtLink       = *(reinterpret_cast<const uint16_t*>(rdh_ptr + 28) ) & 0x7FF;
      flxHdrSize    = *(reinterpret_cast<const uint8_t*>(rdh_ptr  + 29) ) & 0xFF;
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

    bool checkRDH( const bool verbose )
    {
      // check if rdh conform with RDH8 fields
      bool ok = true;
      if ( flxHdrSize != 0x20 )
      {
        if ( verbose )
        {
          std::cout << "RDH FLX Header Size 0x20 is expected instead of " << int(flxHdrSize) << std::endl;
        }
        ok = false;
      }

      if ( flxHdrVersion != 0xAB01 )
      {
      if ( verbose )
        {
          std::cout << "RDH FLX Header version 0x01AB is expected instead of " << int(flxHdrVersion) << std::endl;
        }
        ok = false;
      }

      if ( rdhVersion != 0x08 )
      {
        if ( verbose )
        {
          std::cout << "RDH version 8 is expected instead of " << int(rdhVersion) << std::endl;
        }
        ok = false;
      }

      if ( rdhSize != 32 )
      {
      if ( verbose )
        {
          std::cout << "RDH with header size of 64 B is expected instead of " << int(rdhSize) << std::endl;
        }
        ok = false;
      }

      return ok;
    }
  };

  template < typename A, typename B >
  bool comp(A a, B b)
  {
    return a.second < b.second;
  }

float getStrobeLength(const int& runNumber);


} //namespace mvtx_utils

#endif
