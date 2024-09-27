#include "mvtx_pool.h"

#include <string>
#include <cstdint>

#include <Event/packet.h>

using namespace std;

//_________________________________________________
mvtx_pool::~mvtx_pool()
{
  for ( auto&& link : mGBTLinks )
  {
    link.clear(true, true);
  }

  feeid_set.clear();
}


//_________________________________________________
int mvtx_pool::addPacket(Packet *p)
{
  if (!p)
  {
    return 0;
  }
  m_is_decoded = false;

  feeid_set.clear();

  for (auto& link : mGBTLinks)
  {
    link.clear(false, true); // clear data but not the statistics
  }

  loadInput(p);
  setupLinks();

  mBuffer.movePtr(payload_position);

  if ( mBuffer.isEmpty() )
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
void mvtx_pool::loadInput(Packet *p)
{
  uint8_t* payload_start = (uint8_t *) p->pValue(p->getIdentifier());  // here begins the payload
  unsigned int dlength = p->getDataLength() - p->getPadding(); //padding is supposed to be in units of dwords, this assumes dwords
  dlength *= 4;

  if ( dlength % mvtx_utils::FLXWordLength )
  {
    dlength -= dlength % mvtx_utils::FLXWordLength;
    COUT
      << ENDL
      << "!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!! \n"
      << "DMA packet has incomplete FLX words, only "
      << dlength << " bytes(" << (dlength / mvtx_utils::FLXWordLength)
      << " FLX words), will be decoded. \n"
      << "!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!! \n"
      << ENDL;
  }

  mBuffer.add(payload_start, dlength);
  payload = mBuffer.getPtr();
  payload_position = 0;

  return;
}


//_________________________________________________
void mvtx_pool::setupLinks()
{
  mvtx_utils::RdhExt_t rdh = {};
  size_t dlength = mBuffer.getUnusedSize();
  do
  {
    // Skip FLX padding
    if ( *(reinterpret_cast<uint16_t*>(&payload[payload_position] + 30)) == 0xFFFF )
    {
      payload_position += mvtx_utils::FLXWordLength;
    }
    else if ( (dlength - payload_position) >= static_cast<uint64_t>(2 * mvtx_utils::FLXWordLength) ) // at least FLX header and RDH
    {
      if ( *(reinterpret_cast<uint16_t*>(&payload[payload_position] + 30)) == 0xAB01 )
      {
        rdh.decode(&payload[payload_position]);
        if ( ! rdh.checkRDH(true) )
        {
          // In case of corrupt RDH, skip felix word and continue to next
          payload_position += mvtx_utils::FLXWordLength;
          continue;
        }
        const size_t pageSizeInBytes = static_cast<size_t>((rdh.pageSize + 1) * mvtx_utils::FLXWordLength);
        if ( pageSizeInBytes > (dlength - payload_position) )
        {
          break; // skip incomplete felix packet, return to fetch more data
        }
        else
        {
          feeid_set.insert(rdh.feeId);
          auto& lnkref = mFeeId2LinkID[rdh.feeId];
          if ( lnkref.entry == -1 )
          {
            lnkref.entry = mGBTLinks.size();
            mGBTLinks.emplace_back(rdh.flxId, rdh.feeId);
          }
          auto& gbtLink = mGBTLinks[lnkref.entry];

          if ( (rdh.packetCounter) && (rdh.packetCounter != gbtLink.prev_pck_cnt + 1) )
          {
            log_error << "Incorrect pages count " << rdh.packetCounter <<", previous page count was " \
              << gbtLink.prev_pck_cnt << std::endl;
            payload_position += pageSizeInBytes;
            continue;
          }
          gbtLink.prev_pck_cnt = rdh.packetCounter;

          gbtLink.data.add((payload + payload_position), pageSizeInBytes);

          if ( ! rdh.packetCounter ) // start HB
          {
            if ( gbtLink.hbf_length )
            {
              log_error << "FLX: " << gbtLink.flxId << ", FeeId: " << gbtLink.feeId \
                << ". Found new HBF before stop previous HBF. Previous HBF will be ignored." << std::endl;
                gbtLink.cacheData(gbtLink.hbf_length, true);
            }
            gbtLink.hbf_length = pageSizeInBytes;
            gbtLink.hbf_error = false;
          }
          else
          {
            if (! gbtLink.hbf_length )
            {
              log_error << "FLX: " << gbtLink.flxId << ", FeeId: " << gbtLink.feeId
              << ". Found continuous HBF before start new HBF. data will be ignored." << std::endl;
              gbtLink.hbf_error = true;
            }
            gbtLink.hbf_length += pageSizeInBytes;
            if ( rdh.stopBit ) // found HB end
            {
              gbtLink.cacheData(gbtLink.hbf_length, gbtLink.hbf_error);
              gbtLink.hbf_length = 0;
            }
          }
          payload_position += pageSizeInBytes;
        }
      }
      else
      {
        // skip raw data without a initial FLX header
        // (YCM)TODO: OK for OM but error otherwise
        if (false)
        {
          std::cout << "Felix header: " << std::hex << "0x" << std::setfill('0') << std::setw(4);
          std::cout << *(reinterpret_cast<uint16_t*>(&payload[payload_position] + 30)) << std::dec <<std::endl;
        }
        payload_position += mvtx_utils::FLXWordLength;
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

  for ( auto& link : mGBTLinks )
  {
    link.collectROFCableData();
  }

  return 0;
}


//_________________________________________________
int mvtx_pool::iValue(const int n, const char *what)
{
  mvtx_decode();
  if ( n == -1 ) // Global Information.
  {
    if ( strcmp(what, "NR_LINKS") == 0)
    {
      return feeid_set.size();
    }
    else
    {
      std::cout << "Unknow option " << what << std::endl;
      return -1;
    }
  }

  unsigned int i = n;
  if ( strcmp(what, "FEEID") == 0 )
  {
    return (i < feeid_set.size()) ? *(next(feeid_set.begin(), i)) : -1;
  }
  else
  {
    if (mFeeId2LinkID.find(i) == mFeeId2LinkID.cend())
    {
      log_error << "FeeId " << i << " was not found in the feeId mapping for this packet" << std::endl;
      assert(false);
    }
    uint32_t lnkId =  mFeeId2LinkID[i].entry;
    if ( strcmp(what, "NR_HBF") == 0 )
    {
      if ( mGBTLinks[lnkId].rawData.getNPieces() != mGBTLinks[lnkId].hbf_count)
      {
        log_error << "Mismatch size for HBF from hbfData: " << mGBTLinks[lnkId].hbf_count << " and link rawData Pieces: " \
          << mGBTLinks[lnkId].rawData.getNPieces() << std::endl;
        assert(false);
      }
      return mGBTLinks[lnkId].hbf_count;
    }
    else if ( strcmp(what, "NR_PHYS_TRG") == 0 )
    {
      return mGBTLinks[lnkId].mL1TrgTime.size();
    }
    else if ( strcmp(what, "NR_STROBES") == 0 )
    {
    return mGBTLinks[lnkId].mTrgData.size();
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

  if (mFeeId2LinkID.find(feeId) == mFeeId2LinkID.cend())
  {
    log_error << "FeeId " << feeId << " was not found in the feeId mapping for this packet" << std::endl;
    assert(false);
  }
  uint32_t lnkId =  mFeeId2LinkID[feeId].entry;

  if ( strcmp(what, "L1_IR_BC") == 0 )
  {
    return (index < mGBTLinks[lnkId].mL1TrgTime.size()) ? mGBTLinks[lnkId].mL1TrgTime[index].bc : -1;
  }
  else if ( strcmp(what, "TRG_IR_BC") == 0 )
  {
    return (index < mGBTLinks[lnkId].mTrgData.size()) ? mGBTLinks[lnkId].mTrgData[index].ir.bc : -1;
  }
  else if ( strcmp(what, "TRG_DET_FIELD") == 0 )
  {
    return (index < mGBTLinks[lnkId].mTrgData.size()) ? mGBTLinks[lnkId].mTrgData[index].detectorField : -1;
  }
  else if ( strcmp(what, "TRG_NR_HITS") == 0)
  {
    return (index < mGBTLinks[lnkId].mTrgData.size()) ? mGBTLinks[lnkId].mTrgData[index].hit_vector.size() : -1;
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
std::vector<mvtx::mvtx_hit *> &mvtx_pool::get_hits(const int feeId, const int i_strb)
{
  return mGBTLinks[mFeeId2LinkID[feeId].entry].mTrgData[i_strb].hit_vector;
}

//_________________________________________________
long long int mvtx_pool::lValue(const int i_feeid, const int idx, const char *what)
{
  mvtx_decode();

  uint32_t feeId = i_feeid;
  uint32_t index = idx;

  if (mFeeId2LinkID.find(feeId) == mFeeId2LinkID.cend())
  {
    log_error << "FeeId " << feeId << "was not found in the feeId mapping for this packet" << std::endl;
    assert(false);
  }
  uint32_t lnkId =  mFeeId2LinkID[feeId].entry;

  if ( strcmp(what, "L1_IR_BCO") == 0 )
  {
    return (index < mGBTLinks[lnkId].mL1TrgTime.size()) ? mGBTLinks[lnkId].mL1TrgTime[index].orbit : -1;
  }
  else if ( strcmp(what, "TRG_IR_BCO") == 0 )
  {
    return (index < mGBTLinks[lnkId].mTrgData.size()) ? mGBTLinks[lnkId].mTrgData[index].ir.orbit : -1;
  }
  else
  {
    std::cout << "Unknow option " << what << std::endl;
    return -1;
  }

  return 0;
}

