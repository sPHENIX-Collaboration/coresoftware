#ifndef __MVTX_POOL_H__
#define __MVTX_POOL_H__

#include "mvtx_decoder/PayLoadCont.h"
#include "mvtx_decoder/GBTLink.h"

#include <set>
#include <unordered_map>
#include <vector>

#include <cstdint>

class Packet;

class  mvtx_pool {

 public:
  mvtx_pool() = default;
  virtual ~mvtx_pool();

  virtual int addPacket(Packet *);

  int iValue(const int ,const char * what);
  int iValue(const int, const int, const char* what);
  int iValue(const int, const int, const int, const char* what);

  long long int lValue(const int, const int, const char* what);

 protected:
  int mvtx_decode();
  bool m_is_decoded = false;

    struct dumpEntry
  {
    int entry = -1;
  };

  void loadInput(Packet *);
  void setupLinks();

  mvtx::PayLoadCont mBuffer = {};
  std::unordered_map<uint16_t, dumpEntry> mFeeId2LinkID; // link fee_id to GBTLinks
  std::vector<mvtx::GBTLink> mGBTLinks;
  std::set<uint16_t> feeid_set;

  uint8_t *payload = nullptr;
  unsigned int payload_position = 0;
};

#endif /* __MVTX_POOL_H__ */
