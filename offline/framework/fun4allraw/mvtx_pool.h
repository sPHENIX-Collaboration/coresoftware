#ifndef __MVTX_POOL_H__
#define __MVTX_POOL_H__

#include "mvtx_decoder/PayLoadCont.h"
#include "mvtx_decoder/GBTLink.h"
#include "mvtx_decoder/StrobeData.h"

#include <set>
#include <unordered_map>
#include <vector>

#include <cstdint>

class Packet;

class mvtx_pool {

 public:
  mvtx_pool() = default;
  virtual ~mvtx_pool();

  virtual int addPacket(Packet* p);

  size_t get_feeidSet_size();

  int get_feeid(const uint16_t iLink);
  int get_hbfSet_size(const uint16_t iLink);
  int get_trgSet_size(const uint16_t iLink);
  int get_strbSet_size(const uint16_t iLink);

  int get_L1_IR_BC(const uint16_t iLnk, const uint32_t);
  int get_TRG_IR_BC(const uint16_t iLnk, const uint32_t);
  int get_TRG_DET_FIELD(const uint16_t iLnk, const uint32_t);
  int get_TRG_NR_HITS(const uint16_t iLnk, const uint32_t);

  long long int get_L1_IR_BCO(const uint16_t iLnk, const uint32_t);
  long long int get_TRG_IR_BCO(const uint16_t iLnk, const uint32_t);

  int iValue(const int, const char* what);
  int iValue(const int, const int, const char* what);
  int iValue(const int, const int, const int, const char* what);

  long long int lValue(const int, const int, const char* what);

  std::vector<mvtx::mvtx_hit*>& get_hits(const int feeId,
                                         const int i_strb);

  void set_verbosity(const int val) { verbosity = val; }
  int  get_verbosity() { return verbosity; }

 protected:
   uint32_t get_linkId(const uint16_t i);

 private:
  int mvtx_decode();
  bool m_is_decoded = false;

  int verbosity = 0;

  struct dumpEntry
  {
    int entry = -1;
  };

  void loadInput(Packet* p);
  void setupLinks();

  mvtx::PayLoadCont mBuffer = {};
  std::unordered_map<uint16_t, dumpEntry> mFeeId2LinkID; // link fee_id to GBTLinks
  std::vector<mvtx::GBTLink> mGBTLinks;
  std::set<uint16_t> feeid_set;

  uint8_t* payload = nullptr;
  unsigned int payload_position = 0;
};

#endif /* __MVTX_POOL_H__ */
