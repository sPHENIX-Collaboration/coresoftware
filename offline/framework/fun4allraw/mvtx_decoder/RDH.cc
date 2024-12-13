#include "RDH.h"

#include <stdexcept>
#include <cstring>
#include <iostream>
#include <iomanip>

#if __cplusplus >= 202002L
#include <format>
#endif

using namespace mvtx;

//_________________________________________________
/// placeholder copied from specific version
RDHAny::RDHAny(int v)
{
  if (v == 0)
  {
    *this = RDH{};
  }
  else if (v == 8)
  {
    *this = RDHv8{};
  }
  else
  {
    throw std::runtime_error(std::string("unsupported RDH version ") + std::to_string(v));
  }
}

//_________________________________________________
void RDHAny::copyFrom(const void* rdh)
{
  std::memcpy(this, rdh, sizeof(RDHAny));
}

#if __cplusplus >= 202002L
//_________________________________________________
void RDHUtils::printRDH(const RDHv8& rdh)
{
 // std::bitset<32> trb(rdh.trgType);
  std::cout << std::format(
    "FlxHdrCntr:{:d} flxId:{:d} pageSize:0x{:03x} gbtLink:0x{:02x} ",
    (rdh.flxHdrCntr), (rdh.flxId), (rdh.pageSize), (rdh.gbtLink)).c_str();
  std::cout << std::format(
    "FlxHdrSize:0x{:02x} FlxHdrVersion:0x{:02x} FlxHdrCode:0x{:02x}",
    (rdh.flxHdrSize), (rdh.flxHdrVersion), (rdh.flxHdrCode)).c_str();
  std::cout << std::endl;
  std::cout << std::format(
    "RDHversion:0x{:02x} RDHsize:0x{:02x} FEEID:0x{:04x} SrcID:0x{:02x} ",
    (rdh.rdhVersion), (rdh.rdhSize), (rdh.feeId), (rdh.sourceId)).c_str();
  std::cout << std::endl;
  std::cout << std::format(
    "detField:0x{:08x} bc:0x{:03x} orbit:0x{:010x}",
    (rdh.detectorField), (rdh.bc), (rdh.orbit)).c_str();
  std::cout << std::endl;
  std::cout << std::format(
    "trgType:0x{:08x} packetCounter:0x{:04x} stopBit:0x{:01x} priority:0x{:01x} RDHGBTCntr:0x{:04X} ",
    (rdh.trgType), (rdh.packetCounter), (rdh.stopBit), (rdh.priority), (rdh.rdhGBTcounter)).c_str();
  std::cout << std::endl;
}
#else
void RDHUtils::printRDH(const RDHv8& /*rdh*/)
{
  std::cerr << "RDHUtils::printRDH(const RDHv8&) only implemented in c++20." << std::endl;
}
#endif

//_________________________________________________
void RDHUtils::printRDH(const void* rdhP)
{
  int version = getVersion(rdhP);
  if(version==8){
      printRDH(*reinterpret_cast<const RDHv8*>(rdhP));
  }
  else{
      std::cerr << "Unexpected RDH version " << version << " from";
      dumpRDH(rdhP);
      throw std::runtime_error("invalid RDH provided");
  }
  
}

//_________________________________________________
#if __cplusplus >= 202002L
void RDHUtils::dumpRDH(const void* rdhP)
{
  const uint32_t* w32 = reinterpret_cast<const uint32_t*>(rdhP);
  for (int i = 0; i < 4; i++)
  {
    int l = 4 * i;
    std::cout << std::format(
      "[rdh{:d}] 0x{:08x} 0x{:08x} 0x{:08x} 0x{:08x}",
      i, w32[l + 3], w32[l + 2], w32[l + 1], w32[l]).c_str();
    std::cout << std::endl;
  }
}
#else
void RDHUtils::dumpRDH(const void* /*rdhP*/)
{
  std::cerr << "void RDHUtils::dumpRDH(const void*) only implemented in c++20." << std::endl;
}
#endif

//_________________________________________________
bool RDHUtils::checkRDH(const void* rdhP, bool verbose, bool checkZeros)
{
  int version = getVersion(rdhP);
  bool ok = true;

  if(version==8){
      ok = checkRDH(*reinterpret_cast<const RDHv8*>(rdhP), verbose, checkZeros);
    }
  else {
      ok = false;
      if (verbose)
      {
        std::cerr << "WARNING: "
                  << "Unexpected RDH version " << version << " from" << std::endl;
      }
  }
  if (!ok && verbose) {
    dumpRDH(rdhP);
  }
  return ok;
}

//_____________________________________________________________________
bool RDHUtils::checkRDH(const RDHv8& rdh, bool verbose, bool checkZeros)
{
  // check if rdh conforms with RDH8 fields
  bool ok = true;
  if (rdh.flxHdrSize != 0x20)
  {
    if (verbose)
    {
      std::cerr << "WARNING: "
                << "RDH FLX Header Size 0x20 is expected instead of "
                << int(rdh.flxHdrSize) << std::endl;
    }
    ok = false;
  }
  if (rdh.flxHdrVersion != 0x01)
  {
    if (verbose)
    {
      std::cerr << "WARNING: "
                << "RDH FLX Header version 0x01 is expected instead of "
                << int(rdh.flxHdrVersion) << std::endl;
    }
    ok = false;
  }
  if (rdh.flxHdrCode != 0xAB)
  {
    if (verbose)
    {
      std::cerr << "WARNING: "
                << "RDH FLX Header Code 0xAB is expected instead of "
                << int(rdh.flxHdrVersion) << std::endl;
    }
    ok = false;
  }
  if (rdh.rdhVersion != 8)
  {
    if (verbose) {
      std::cerr << "WARNING: "
                << "RDH version 7 is expected instead of "
                << int(rdh.rdhVersion) << std::endl;
    }
    ok = false;
  }
  if (rdh.rdhSize != 32)
  {
    if (verbose) {
      std::cerr << "WARNING: "
                << "RDH with header size of 32 B is expected instead of "
                << int(rdh.rdhSize) << std::endl;
    }
    ok = false;
  }
  if ((! rdh.packetCounter) && (rdh.stopBit))
  {
    if (verbose) {
      std::cerr << "WARNING: "
                << "RDH stopBit is not expected in packetCounter 0 "
                << int(rdh.packetCounter) << int(rdh.stopBit) <<  std::endl;
    }
    ok = false;
  }
  if (checkZeros && (rdh.zero0 || rdh.zero1 || rdh.zero2 || rdh.zero3 || rdh.zero4 || rdh.zero5))
  {
    if (verbose)
    {
      std::cerr << "WARNING: "
                << "Some reserved fields of RDH v7 are not empty"
                << std::endl;
    }
    ok = false;
  }
  if (!ok && verbose)
  {
    dumpRDH(rdh);
  }
  return ok;
}
