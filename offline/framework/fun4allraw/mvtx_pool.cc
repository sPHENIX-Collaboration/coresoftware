#include "mvtx_pool.h"

#include <string>
#include <cstdint>

#include <Event/packet.h>
#include "mvtx_decoder/RDH.h"

using namespace std;

//_________________________________________________
mvtx_pool::~mvtx_pool()
{
  if (get_verbosity() > 1)
  {
    std::cout << "LOG mvtx_pool::~mvtx_pool() called." << std::endl;
  }

  for (auto& link : mGBTLinks)
  {
    // clear data and the statistics
    link.clear(true, true);
  }

  feeid_set.clear();
}


//_________________________________________________
int mvtx_pool::addPacket(Packet* p)
{
  if (get_verbosity() > 1)
  {
    std::cout << "LOG mvtx_pool::addPacket(Packet *p) called." << std::endl;
  }

  if (! p)
  {
    return 0;
  }
  // force decode data after new packet
  m_is_decoded = false;

  for (auto& link : mGBTLinks)
  {
    // clear data but not the statistics
    link.clear(false, true);
  }
  feeid_set.clear();

  loadInput(p);
  setupLinks();

  mBuffer.movePtr(payload_position);

  if (mBuffer.isEmpty())
  {
    mBuffer.clear();
  }
  else
  {
    mBuffer.moveUnusedToHead();
  }

  return 0;
}


//_________________________________________________
void mvtx_pool::loadInput(Packet* p)
{
  // here begins the payload
  uint8_t* payload_start = reinterpret_cast<uint8_t*>(p->pValue(p->getIdentifier()));
  //padding is supposed to be in units of dwords, this assumes dwords
  unsigned int dlength = p->getDataLength() - p->getPadding();
  dlength *= 4;

  if ((dlength < mvtx_utils::FLXWordLength) || (dlength % mvtx_utils::FLXWordLength))
  {
    COUT << ENDL
         << "!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!! \n"
         << "DMA packet has incomplete FLX words, only "
         << dlength << " bytes(" << (dlength / mvtx_utils::FLXWordLength)
         << " FLX words), will be decoded. \n"
         << "!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!! \n"
         << ENDL;
    dlength -= (dlength % mvtx_utils::FLXWordLength);
  }

  // Add raw data from packet to buffer
  mBuffer.add(payload_start, dlength);

  // Repositioning pointer to unread data in the buffer
  payload = mBuffer.getPtr();
  payload_position = 0;

  return;
}

//_________________________________________________
void mvtx_pool::setupLinks()
{
  size_t dlength = mBuffer.getUnusedSize();
  do
  {
    // Skip FLX padding
    if (*(reinterpret_cast<uint16_t*>(&payload[payload_position] + 30)) == 0xFFFF)
    {
      payload_position += mvtx_utils::FLXWordLength;
    }
    // at least one combine FLX header and RDH words
    else if ((dlength - payload_position) >= static_cast<uint64_t>(2 * mvtx_utils::FLXWordLength))
    {
      if (*(reinterpret_cast<uint16_t*>(&payload[payload_position] + 30)) == 0xAB01)
      {
        const auto* rdhP = reinterpret_cast<const mvtx::RDH*>(&payload[payload_position]);
        if (get_verbosity() > 3)
        {
          mvtx::RDHUtils::printRDH(mvtx::RDHAny::voidify(*rdhP));
        }
        if (! mvtx::RDHUtils::checkRDH(mvtx::RDHAny::voidify(*rdhP), true, true))
        {
          // In case of corrupt RDH, skip felix word and continue to next
          payload_position += mvtx_utils::FLXWordLength;
          continue;
        }
        const size_t pageSizeInBytes = ((*rdhP).pageSize + 1ULL /*add Flx Hdr word*/) * mvtx_utils::FLXWordLength;
        if (pageSizeInBytes > (dlength - payload_position))
        {
          if (get_verbosity() > 1)
          {
            std::cout << "WARNING: "
                      << "Skipping Incomplete FELIX packet" << std::endl;
          }
          // skip incomplete felix packet, return to fetch more data
          break;
        }
        else
        {
          feeid_set.insert((*rdhP).feeId);
          auto& lnkref = mFeeId2LinkID[(*rdhP).feeId];
          if (lnkref.entry == -1)
          {
            lnkref.entry = mGBTLinks.size();
            mGBTLinks.emplace_back((*rdhP).flxId, (*rdhP).feeId);
          }
          auto& gbtLink = mGBTLinks[lnkref.entry];

          if (! (*rdhP).packetCounter) // start HB
          {
            // close previous HBF without stop Bit
            if (gbtLink.hbf_length)
            {
              log_error << "FLX: " << gbtLink.flxId << ", FeeId: " << gbtLink.feeId \
                << ". Found new HBF before stop previous HBF. Previous HBF will be ignored." << std::endl;
                gbtLink.cacheData(gbtLink.hbf_length, (gbtLink.hbf_error |= mvtx::PayLoadSG::HBF_ERRORS::Incomplete));
            }
            gbtLink.hbf_length = pageSizeInBytes;
            gbtLink.hbf_error = mvtx::PayLoadSG::HBF_ERRORS::NoError;
          }
          else
          {
            if ((*rdhP).packetCounter != gbtLink.prev_pck_cnt + 1)
            {
              log_error << "Incorrect pages count " << (*rdhP).packetCounter <<", previous page count was " \
                        << gbtLink.prev_pck_cnt << std::endl;
              gbtLink.hbf_length += pageSizeInBytes;
              gbtLink.hbf_error |= mvtx::PayLoadSG::HBF_ERRORS::Incomplete;
            }
            else
            {
              if (! gbtLink.hbf_length)
              {
                log_error << "FLX: " << gbtLink.flxId << ", FeeId: " << gbtLink.feeId
                          << ". Found continuous HBF before start new HBF. data will be ignored." << std::endl;
                gbtLink.hbf_error |= mvtx::PayLoadSG::HBF_ERRORS::Incomplete;
              }
              gbtLink.hbf_length += pageSizeInBytes;

              if ((*rdhP).stopBit) // found HB end
              {
                gbtLink.cacheData(gbtLink.hbf_length, gbtLink.hbf_error);
                gbtLink.hbf_length = 0;
              }
            }
          }
          gbtLink.prev_pck_cnt = (*rdhP).packetCounter;

          // move packet to buffer
          gbtLink.data.add((payload + payload_position), pageSizeInBytes);
          payload_position += pageSizeInBytes;
        }
      }
      else
      {
        // skip raw data without a initial FLX header
        // (YCM)TODO: OK for OM but error otherwise
        if (get_verbosity() > 0)
        {
          std::cout << "Felix header: " << std::hex << "0x" << std::setfill('0') << std::setw(4);
          std::cout << *(reinterpret_cast<uint16_t*>(&payload[payload_position] + 30)) << std::dec <<std::endl;
        }
        // move to next flx word and continue flx word loop
        payload_position += mvtx_utils::FLXWordLength;
        continue;
      }
    }
    else
    {
      break; // skip incomplete flx_header
    }
  } while (payload_position < dlength);

  return;
}


//_________________________________________________
int mvtx_pool::mvtx_decode()
{
  if (m_is_decoded)
  {
    return 0;
  }
  m_is_decoded = true;

  for (auto& link : mGBTLinks)
  {
    link.collectROFCableData();
  }

  return 0;
}

//_________________________________________________
uint32_t mvtx_pool::get_linkId(const uint16_t iLnk)
{
  if (mFeeId2LinkID.find(iLnk) == mFeeId2LinkID.cend())
  {
    log_error << "FeeId " << iLnk
              << " was not found in the feeId mapping for this packet"
              << std::endl;
    return std::numeric_limits<uint32_t>::quiet_NaN();
  }
  return mFeeId2LinkID[iLnk].entry;
}

//_________________________________________________
size_t mvtx_pool::get_feeidSet_size()
{
  mvtx_decode();
  return feeid_set.size();
}

//_________________________________________________
int mvtx_pool::get_feeid(const uint16_t iLnk)
{
  return (iLnk < feeid_set.size()) ? *(next(feeid_set.begin(), iLnk)) : -1;
}

//_________________________________________________
int mvtx_pool::get_hbfSet_size(const uint16_t iLnk)
{
  mvtx_decode();
  auto lnkId = get_linkId(iLnk);
  if (mGBTLinks[lnkId].rawData.getNPieces() != mGBTLinks[lnkId].hbf_count)
  {
    log_error << "Mismatch size for HBF from hbfData: "
              << mGBTLinks[lnkId].hbf_count << " and link rawData Pieces: "
              << mGBTLinks[lnkId].rawData.getNPieces() << std::endl;
    return -1;
  }
  return mGBTLinks[lnkId].hbf_count;
}

//_________________________________________________
int mvtx_pool::get_trgSet_size(const uint16_t iLnk)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mL1TrgTime.size();
}

//_________________________________________________
int mvtx_pool::get_strbSet_size(const uint16_t iLnk)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mTrgData.size();
}

//_________________________________________________
int mvtx_pool::get_L1_IR_BC(const uint16_t iLnk, const uint32_t index)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  // return (index < mGBTLinks[lnkId].mL1TrgTime.size()) ?
  //         mGBTLinks[lnkId].mL1TrgTime[index].bc : -1;
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mL1TrgTime[index].bc;
}

//_________________________________________________
int mvtx_pool::get_TRG_IR_BC(const uint16_t iLnk, const uint32_t index)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  // return (index < mGBTLinks[lnkId].mTrgData.size()) ?
  //         mGBTLinks[lnkId].mTrgData[index].ir.bc : -1;
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mTrgData[index].ir.bc;
}

//_________________________________________________
int mvtx_pool::get_TRG_DET_FIELD(const uint16_t iLnk, const uint32_t index)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  // return (index < mGBTLinks[lnkId].mTrgData.size()) ?
  //         mGBTLinks[lnkId].mTrgData[index].detectorField : -1;
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mTrgData[index].detectorField;
}

//_________________________________________________
int mvtx_pool::get_TRG_NR_HITS(const uint16_t iLnk, const uint32_t index)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  // return (index < mGBTLinks[lnkId].mTrgData.size()) ?
  //         mGBTLinks[lnkId].mTrgData[index].hit_vector.size() : -1;
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mTrgData[index].hit_vector.size();
}

//_________________________________________________
long long int mvtx_pool::get_L1_IR_BCO(const uint16_t iLnk, const uint32_t index)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  // return (index < mGBTLinks[lnkId].mL1TrgTime.size()) ?
  //         mGBTLinks[lnkId].mL1TrgTime[index].orbit : -1;
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mL1TrgTime[index].orbit;
}

//_________________________________________________
long long int mvtx_pool::get_TRG_IR_BCO(const uint16_t iLnk, const uint32_t index)
{
  // mvtx_decode();
  // auto lnkId = get_linkId(iLnk);
  // return (index < mGBTLinks[lnkId].mTrgData.size()) ?
  //         mGBTLinks[lnkId].mTrgData[index].ir.orbit : -1;
  return mGBTLinks[mFeeId2LinkID[iLnk].entry].mTrgData[index].ir.orbit;
}

//_________________________________________________
std::vector<mvtx::mvtx_hit *>& mvtx_pool::get_hits(const int feeId, const int i_strb)
{
  return mGBTLinks[mFeeId2LinkID[feeId].entry].mTrgData[i_strb].hit_vector;
}

//_________________________________________________
int mvtx_pool::iValue(const int n, const char *what)
{
  mvtx_decode();
  if (n == -1) // Global Information.
  {
    if (strcmp(what, "NR_LINKS") == 0)
    {
      return get_feeidSet_size();
    }
    else
    {
      std::cout << "Unknow option " << what << std::endl;
      return -1;
    }
  }

  unsigned int i = n;
  if (strcmp(what, "FEEID") == 0)
  {
    return get_feeid(i);
  }
  else
  {
    if ( strcmp(what, "NR_HBF") == 0 )
    {
      return get_hbfSet_size(i);
    }
    else if ( strcmp(what, "NR_PHYS_TRG") == 0 )
    {
      return get_trgSet_size(i);
    }
    else if ( strcmp(what, "NR_STROBES") == 0 )
    {
    return get_strbSet_size(i);
    }
    else if ( strcmp(what, "NR_HITS") == 0 )  // the number of datasets
    {
      return -1;
    }
    else
    {
      std::cout << "Unknow option " << what << std::endl;
      return -1;
    }
  }
  return 0;
}

//_________________________________________________
int mvtx_pool::iValue(const int i_feeid, const int idx, const char *what)
{
  mvtx_decode();
  uint32_t feeId = i_feeid;
  uint32_t index = idx;

  if ( strcmp(what, "L1_IR_BC") == 0 )
  {
    return get_L1_IR_BC(feeId, index);
  }
  else if ( strcmp(what, "TRG_IR_BC") == 0 )
  {
    return get_TRG_IR_BC(feeId, index);
  }
  else if ( strcmp(what, "TRG_DET_FIELD") == 0 )
  {
    return get_TRG_DET_FIELD(feeId, index);
  }
  else if ( strcmp(what, "TRG_NR_HITS") == 0)
  {
    return get_TRG_NR_HITS(feeId, index);
  }
  else
  {
    std::cout << "Unknow option " << what << std::endl;
    return -1;
  }
  return 0;
}


//_________________________________________________
int mvtx_pool::iValue(const int i_feeid, const int i_trg, const int i_hit, const char *what)
{
  mvtx_decode();

  uint32_t feeId = i_feeid;
  uint32_t trg = i_trg;
  uint32_t hit = i_hit;

  if (mFeeId2LinkID.find(feeId) == mFeeId2LinkID.cend())
  {
    log_error << "FeeId " << feeId << "was not found in the feeId mapping for this packet" << std::endl;
    assert(false);
  }
  uint32_t lnkId =  mFeeId2LinkID[feeId].entry;

  if ( strcmp(what, "HIT_CHIP_ID") == 0 )
  {
    return ( (i_hit >= 0) && (hit < mGBTLinks[lnkId].mTrgData[trg].hit_vector.size()) ) ? \
                     mGBTLinks[lnkId].mTrgData[trg].hit_vector[hit]->chip_id : -1;
  }
  else if ( strcmp(what, "HIT_BC") == 0 )
  {
    return ( (i_hit >= 0) && (hit < mGBTLinks[lnkId].mTrgData[trg].hit_vector.size()) ) ? \
                     mGBTLinks[lnkId].mTrgData[trg].hit_vector[hit]->bunchcounter : -1;
  }
  else if ( strcmp(what, "HIT_ROW") == 0 )
  {
    return ( (i_hit >= 0) && (hit < mGBTLinks[lnkId].mTrgData[trg].hit_vector.size()) ) ? \
                     mGBTLinks[lnkId].mTrgData[trg].hit_vector[hit]->row_pos : -1;
  }
  else if ( strcmp(what, "HIT_COL") == 0 )
  {
    return ( (i_hit >= 0) && (hit < mGBTLinks[lnkId].mTrgData[trg].hit_vector.size()) ) ? \
                     mGBTLinks[lnkId].mTrgData[trg].hit_vector[hit]->col_pos : -1;
  }
  else
  {
    std::cout << "Unknow option " << what << std::endl;
    return -1;
  }
  return 0;
}

//_________________________________________________
long long int mvtx_pool::lValue(const int i_feeid, const int idx, const char *what)
{
  mvtx_decode();

  uint32_t feeId = i_feeid;
  uint32_t index = idx;

  if ( strcmp(what, "L1_IR_BCO") == 0 )
  {
    return get_L1_IR_BCO(feeId, index);
  }
  else if ( strcmp(what, "TRG_IR_BCO") == 0 )
  {
    return get_TRG_IR_BCO(feeId, index);
  }
  else
  {
    std::cout << "Unknow option " << what << std::endl;
    return -1;
  }

  return 0;
}
