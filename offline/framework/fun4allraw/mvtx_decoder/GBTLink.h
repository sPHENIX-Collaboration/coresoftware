// @file GBTLink.h
// @brief Declarations of helper classes for the ITS/MFT raw data decoding
// @sa <O2/Detectors/ITSMFT/common/reconstruction/include/ITSMFTReconstruction/GBTLink.>
//     <760019308>

#ifndef MVTXDECODER_GBTLINK_H
#define MVTXDECODER_GBTLINK_H

#define _RAW_READER_ERROR_CHECKS_ // comment this to disable error checking

#include "mvtx_decoder/mvtx_utils.h"
#include "mvtx_decoder/PayLoadCont.h"
#include "mvtx_decoder/PayLoadSG.h"
#include "mvtx_decoder/DecodingStat.h"
#include "mvtx_decoder/GBTWord.h"
#include "mvtx_decoder/InteractionRecord.h"

//#include "MVTXDecoder/RUDecodeData.h"
//#include "MVTXDecoder/RUInfo.h"
//#include "MVTXDecoder/RAWDataHeader.h"
//#include "MVTXDecoder/RDHUtils.h"
//#include "MVTXDecoder/PhysTrigger.h"

#include <iostream>
#include <memory>
#include <iomanip>

#define GBTLINK_DECODE_ERRORCHECK(errRes, errEval)                            \
  errRes = errEval;                                                           \
  if ((errRes)&uint8_t(ErrorPrinted)) {                                       \
    ruPtr->linkHBFToDump[(uint64_t(subSpec) << 32) + hbfEntry] = irHBF.orbit; \
    errRes &= ~uint8_t(ErrorPrinted);                                         \
  }                                                                           \
  if ((errRes)&uint8_t(Abort)) {                                              \
    discardData();                                                            \
    return AbortedOnError;                                                    \
  }

namespace mvtx
{
  using namespace mvtx_utils;

  struct TRGData
  {
    InteractionRecord ir = {};
    bool hasCDW = false;
    GBTCalibDataWord calWord = {};
    size_t first_hit_pos = 0;
    size_t n_hits = 0;

    TRGData(uint64_t orb, uint16_t b) : ir(orb, b) {};

    void clear()
    {
      ir.clear();
      hasCDW = false;
      calWord = {};
      first_hit_pos = 0;
      n_hits = 0;
    }
  };

  typedef struct mvtx_hit
  {
    uint8_t chip_id = 0xf;
    uint16_t bunchcounter = 0xFFFF;
    uint16_t row_pos = 0xFFFF;
    uint16_t col_pos = 0xFFFF;
  } mvtx_hit;

/// support for the GBT single link data
struct GBTLink
{
//  enum RawDataDumps : int { DUMP_NONE, // no raw data dumps on error
//                            DUMP_HBF,  // dump HBF for FEEID with error
//                            DUMP_TF,   // dump whole TF at error
//                            DUMP_NTYPES };

  enum CollectedDataStatus : int8_t { None,
                                      AbortedOnError,
                                      StoppedOnEndOfData,
                                      DataSeen }; // None is set before starting collectROFCableData

//  enum ErrorType : uint8_t { NoError = 0x0,
//                             Warning = 0x1,
//                             Skip = 0x2,
//                             Abort = 0x4,
//                             ErrorPrinted = 0x1 << 7 };

  static constexpr int RawBufferMargin = 5000000;                      // keep uploaded at least this amount
  static constexpr int RawBufferSize = 10000000 + 2 * RawBufferMargin; // size in MB
  static constexpr uint8_t MaxCablesPerLink = 3;

  CollectedDataStatus status = None;

  uint16_t flxId = 0;     // FLX ID
  uint16_t feeId = 0;     // FEE ID

  PayLoadCont data; // data buffer for single feeeid
  std::array<PayLoadCont, MaxCablesPerLink> cableData;


  uint32_t hbfEntry = 0;      // entry of the current HBF page in the rawData SG list
  InteractionRecord ir = {};

  GBTLinkDecodingStat statistics; // link decoding statistics
  bool     hbf_found = false;

  uint32_t hbf_length = 0;
  uint32_t prev_pck_cnt = 0;
  uint32_t hbf_count = 0;

  PayLoadSG rawData;         // scatter-gatter buffer for cached CRU pages, each starting with RDH
  size_t dataOffset = 0;     //
  std::vector<InteractionRecord> mL1TrgTime;
  std::vector<TRGData> mTrgData;

  std::vector<mvtx_hit *> hit_vector = {};

  //------------------------------------------------------------------------
  GBTLink() = default;
  GBTLink(uint16_t _flx, uint16_t _fee);
  void clear(bool resetStat = true, bool resetTFRaw = false);

  CollectedDataStatus collectROFCableData();

  void cacheData(size_t start, size_t sz)
  {
    rawData.add(start, sz);
  }

  void clearCableData()
  {
    for ( auto&& _data : cableData )
    {
      _data.clear();
    }
  }

  int readFlxWord( GBTWord* gbtwords, uint16_t &w16 );
  int decode_lane( const uint8_t chipId, PayLoadCont& buffer );

  void getRowCol(const uint8_t reg, const uint16_t addr, uint16_t& row, uint16_t& col)
  {
    row = ( addr >> 0x1 ) & 0x1FF;
    col = ( reg << 5 | ( (addr >> 9) & 0x1E ) ) | ( (addr ^ addr >> 1) & 0x1 );
  }

  void addHit(const uint8_t laneId, const uint8_t bc, uint8_t reg, const uint16_t addr)
  {
    auto* hit = new mvtx_hit();
//	  memset(hit, 0, sizeof(*hit));

    hit->chip_id = laneId;
    hit->bunchcounter = bc;
    getRowCol(reg, addr, hit->row_pos, hit->col_pos);

	  hit_vector.push_back(hit);
  }

  void check_APE(const uint8_t& chipId, const uint8_t& dataC)
  {
    std::cerr << "Link: " << feeId << ", Chip: " << (int)chipId;
    switch (dataC)
    {
      case 0xF2:
        std::cerr << " APE_STRIP_START" << std::endl;
        break;
      case 0xF4:
        std::cerr << " APE_DET_TIMEOUT" << std::endl;
        break;
      case 0xF5:
        std::cerr << " APE_OOT" << std::endl;
        break;
      case 0xF6:
        std::cerr << " APE_PROTOCOL_ERROR" << std::endl;
        break;
      case 0xF7:
        std::cerr << " APE_LANE_FIFO_OVERFLOW_ERROR" << std::endl;
        break;
      case 0xF8:
        std::cerr << " APE_FSM_ERROR" << std::endl;
        break;
      case 0xF9:
        std::cerr << " APE_PENDING_DETECTOR_EVENT_LIMIT" << std::endl;
        break;
      case 0xFA:
        std::cerr << " APE_PENDING_LANE_EVENT_LIMIT" << std::endl;
        break;
      case 0xFB:
        std::cerr << " APE_O2N_ERROR" << std::endl;
        break;
      case 0xFC:
        std::cerr << " APE_RATE_MISSING_TRG_ERROR" << std::endl;
        break;
      case 0xFD:
        std::cerr << " APE_PE_DATA_MISSING" << std::endl;
        break;
      case 0xFE:
        std::cerr << " APE_OOT_DATA_MISSING" << std::endl;
        break;
      default:
        std::cerr << " Unknown APE code" << std::endl;
    }
    return;
  }

  void AlpideByteError(const uint8_t& chipId, PayLoadCont& buffer)
  {
    uint8_t dataC = 0;

    std::cerr << "Link: " << feeId << ", Chip: " << (int)chipId;
    std::cerr << " invalid byte 0x" << std::hex << (int)(dataC) << std::endl;
    while ( buffer.next(dataC) )
    {
      std::cerr << " " << std::hex << (int)(dataC) << " ";
    }
    std::cerr << std::endl;
    buffer.clear();
    return;
  }

  void PrintFlxWord(std::ostream& os, uint8_t* pos)
  {
    os  << std::setfill('0');
    for ( int i = 0; i < 32 ; i++)
    {
      os << std::hex << std::setw(2) << (int)pos[i] << " " << std::dec;
    }
    os  << std::setfill(' ') << std::endl;
  }

  void PrintBlock(std::ostream& os, uint8_t* pos, size_t n)
  {
    for (uint32_t i = 0; i < n; ++i)
    {
      PrintFlxWord(os, pos + 32 * i);
    }
  }

//  ClassDefNV(GBTLink, 1);
};

///_________________________________________________________________
/// collect cables data for single ROF, return number of real payload words seen,
/// -1 in case of critical error
inline GBTLink::CollectedDataStatus GBTLink::collectROFCableData(/*const Mapping& chmap*/)
{
  bool prev_evt_complete = false;
  bool header_found = false;
//  bool trailer_found = false;
  uint8_t* hbf_start = nullptr;

  status = None;

  auto currRawPiece = rawData.currentPiece();
  dataOffset = 0;
  while (currRawPiece)
  { // we may loop over multiple FLX page
    uint32_t n_no_continuation = 0;
    uint32_t n_packet_done = 0;

    if (dataOffset >= currRawPiece->size)
    {
      data.movePtr(currRawPiece->size);
      dataOffset = 0;
      // start of the RDH
      if ( ! (currRawPiece = rawData.nextPiece()) )
      { // fetch next CRU page
        break;                                     // Data chunk (TF?) is done
      }
    }

    if ( currRawPiece->hasError ) // Skip
    {
      dataOffset = currRawPiece->size;
      ++hbf_count;
      continue;
    }

    if ( !dataOffset )
    {
      hbf_start = data.getPtr();
    }

    // here we always start with the RDH
    RdhExt_t rdh = {};
    uint8_t* rdh_start = data.getPtr() + dataOffset;
    rdh.decode(rdh_start);

    size_t pagesize = (rdh.pageSize + 1) * FLXWordLength;
    const size_t nFlxWords = (pagesize - (2 * FLXWordLength)) / FLXWordLength;
    //Fill statistics
    if ( !rdh.packetCounter )
    {
      if ( dataOffset )
      {
        log_error << "Wrong dataOffset value " << dataOffset << " at the start of a HBF" << std::endl;
        assert(false);
      }
      statistics.clear();
      //TODO: initialize/clear alpide data buffer
      for ( uint32_t trg = GBTLinkDecodingStat::BitMaps::ORBIT; trg < GBTLinkDecodingStat::nBitMap; ++trg )
      {
        if  ( (rdh.trgType >> trg) & 1 )
        {
          statistics.trgBitCounts[trg]++;
        }
      }
      hbfEntry = rawData.currentPieceId(); // in case of problems with RDH, dump full TF
      ++hbf_count;
    }
    else if ( !rdh.stopBit )
    {
      if (prev_evt_complete)
      {
        log_error << "Previous event was already completed" << std::endl;
        assert(false);
      }
    }

    dataOffset += 2 * FLXWordLength;
    int prev_gbt_cnt = 3;
    GBTWord gbtWords[3];
    uint16_t w16 = 0;
    for ( size_t iflx = 0; iflx < nFlxWords; ++iflx )
    {
      readFlxWord(gbtWords, w16);
      int16_t n_gbt_cnt = (w16 & 0x3FF) - prev_gbt_cnt;
      prev_gbt_cnt = (w16 & 0x3FF);
      if (n_gbt_cnt < 1 || n_gbt_cnt > 3)
      {
        log_error << "Bad gbt counter in the flx packet. FLX: " << flxId << ", Feeid: " << feeId << ", n_gbt_cnt: " << n_gbt_cnt \
          << ", prev_gbt_cnt: " << prev_gbt_cnt << ", size: " << currRawPiece->size << ", dataOffset: " << dataOffset << std::endl;
        PrintBlock(std::cerr, rdh_start, nFlxWords + 2);
        std::cerr << "Full HBF" << std::endl;
        PrintBlock(std::cerr, hbf_start, (currRawPiece->size/32) );
        break;
      }
      for ( int i = 0; i < n_gbt_cnt; ++i )
      {
        auto &gbtWord = gbtWords[i];
        if ( gbtWord.isIHW() ) // ITS HEADER WORD
        {
          //TODO assert first word after RDH and active lanes
          if (! ( !gbtWord.activeLanes ||
                    ((gbtWord.activeLanes >> 0) & 0x7) == 0x7 ||
                    ((gbtWord.activeLanes >> 3) & 0x7) == 0x7 ||
                    ((gbtWord.activeLanes >> 6) & 0x7) == 0x7) )
          {
            log_error << "Expected all active lanes for links, but " << gbtWord.activeLanes << "found in HBF " << hbfEntry << ", " \
              << gbtWord.asString().data() << std::endl;
            assert(false);
          }
        }
        else if ( gbtWord.isTDH() ) // TRIGGER DATA HEADER (TDH)
        {
          header_found = true;
          ir.orbit = gbtWord.bco;
          ir.bc = gbtWord.bc;
          if ( gbtWord.bc ) //statistic trigger for first bc already filled on RDH
          {
            for ( uint32_t trg = GBTLinkDecodingStat::BitMaps::ORBIT; trg < GBTLinkDecodingStat::nBitMap; ++trg )
            {
              if ( trg == GBTLinkDecodingStat::BitMaps::FE_RST ) //  TDH save first 12 bits only
                break;
              if ( ((gbtWord.triggerType >> trg) & 1) )
              {
                statistics.trgBitCounts[trg]++;
              }
            }
          }

          if ( (gbtWord.triggerType >> GBTLinkDecodingStat::BitMaps::PHYSICS) & 0x1 )
          {
            mL1TrgTime.push_back(ir);
          }

          if ( !gbtWord.continuation && !gbtWord.noData)
          {
            n_no_continuation++;
            mTrgData.emplace_back(ir.orbit, ir.bc);
          } // end if not cont
        } // end TDH
        else if ( gbtWord.isCDW() ) // CALIBRATION DATA WORD
        {
          mTrgData.back().hasCDW = true;
          mTrgData.back().calWord = *(reinterpret_cast<GBTCalibDataWord*>(&gbtWord));
        }
        else if ( gbtWord.isTDT() )
        {
//          trailer_found = true;
          if ( gbtWord.packet_done )
          {
            n_packet_done++;
            if (n_packet_done < n_no_continuation)
            {
              log_error << "TDT packet done before TDH no continuation " << n_no_continuation \
                << " != " << n_packet_done << std::endl;
              assert(false);
            }
          }
          prev_evt_complete = gbtWord.packet_done;
          //TODO: YCM Add warning and counter for timeout and violation
        }
        else if ( gbtWord.isDDW() ) // DIAGNOSTIC DATA WORD (DDW)
        {
          if (! rdh.stopBit)
          {
            log_error << "" << std::endl;
            assert(false);
          }
        }
        else if ( gbtWord.isDiagnosticIB() ) // IB DIAGNOSTIC DATA
        {
            std::cout << "WARNING: IB Diagnostic word found." << std::endl;
            std::cout << "diagnostic_lane_id: " << (gbtWord.id >> 5);
            std::cout << " lane_error_id: " << gbtWord.lane_error_id;
            std::cout << " diasnotic_data: 0x" << std::hex << gbtWord.diagnostic_data << std::endl;
        }
        else if ( gbtWord.isData() ) //IS IB DATA
        {
          if (! header_found )
          {
            log_error << "Trigger header not found before chip data" << std::endl;
            assert(false);
          }
          auto lane = ( gbtWord.data8[9] & 0x1F ) % 3;
          cableData[lane].add(gbtWord.getW8(), 9);
        }

        if ( prev_evt_complete )
        {
          auto&& trgData = mTrgData.back();
          trgData.first_hit_pos = hit_vector.size();
          for( auto&& itr = cableData.begin(); itr != cableData.end(); ++itr)
          {
            if (!itr->isEmpty())
            {
              decode_lane(std::distance(cableData.begin(), itr), *itr);
            }
          }
          trgData.n_hits = hit_vector.size() - trgData.first_hit_pos;
          prev_evt_complete = false;
          header_found = false;
//          trailer_found = false;
          clearCableData();
        }
      }
    }
  }
  return  (status = StoppedOnEndOfData);
}

//_________________________________________________
inline int GBTLink::decode_lane( const uint8_t chipId, PayLoadCont& buffer)
{
  int ret = 0; // currently we just print stuff, but we will add stuff to our
               // structures and return a status later (that's why it's not a const function)

  if ( buffer.getSize() < 3 )
  {
    log_error << "chip data is too short: " << buffer.getSize() << std::endl;
    assert(false);
  }

  uint8_t dataC = 0;
  uint16_t dataS = 0;

  bool busy_on = false, busy_off = false;
  bool chip_header_found = false;
  bool chip_trailer_found = true;

  uint8_t laneId = 0xFF;
  uint8_t bc = 0xFF;
  uint8_t reg = 0xFF;

  if ( !( (buffer[0] & 0xF0) == 0xE0 || (buffer[0] & 0xF0) == 0xA0 ||\
          (buffer[0] == 0xF0) || (buffer[0] == 0xF1) || (buffer[0] & 0xF0) == 0xF0 ) )
  {
    AlpideByteError(chipId, buffer);
    return 0;
  }

  while ( buffer.next(dataC) )
  {
    if ( dataC == 0xF1 ) // BUSY ON
    {
      busy_on = true ;
    }
    else if ( dataC == 0xF0 ) // BUSY OFF
    {
      busy_off = true;
    }
    else if ( (dataC & 0xF0) == 0xF0) // APE
    {
      check_APE(chipId, dataC);
      chip_trailer_found = true;
      busy_on = busy_off = chip_header_found = 0;
    }
    else if ( (dataC & 0xF0) == 0xE0 ) // EMPTY
    {
      chip_header_found = false;
      chip_trailer_found = true;
      laneId = (dataC & 0x0F) % 3;
      if ( laneId != chipId )
      {
        log_error << "Error laneId " << laneId << " (" << (dataC & 0xF) << ") and chipId " << chipId << std::endl;
        assert(false);
      }
      buffer.next(bc);
      busy_on = busy_off = false;
    }
    else
    {
      if ( chip_header_found )
      {
        if ( (dataC & 0xE0) == 0xC0 ) // REGION HEADER
        {
          if ( buffer.getUnusedSize() < 2 )
          {
            log_error << "No data short would fit (at least a data short after region header!)" << std::endl;
            assert(false);
          }
          // TODO: move first region header out of loop, asserting its existence
          reg = dataC & 0x1F;
        }
        else if ( (dataC & 0xC0) == 0x40 ) // DATA SHORT
        {
          if( buffer.isEmpty() )
          {
            log_error << "data short do not fit" << std::endl;
            assert(false);
          }
          if ( reg == 0xFF )
          {
            log_error << "data short at " << buffer.getOffset() << " before region header" << std::endl;
            assert(false);
          }
          dataS = (dataC << 8);
          buffer.next(dataC);
          dataS |= dataC;
          addHit(laneId, bc, reg, (dataS & 0x3FFF));
        }
        else if ( (dataC & 0xC0) == 0x00) // DATA LONG
        {
          if ( buffer.getUnusedSize() < 3 )
          {
            log_error << "No data long would fit (at least a data short after region header!)" << std::endl;
            assert(false);
          }
          if ( reg == 0xFF )
          {
            log_error << "data short at " << buffer.getOffset() << " before region header" << std::endl;
            assert(false);
          }
          buffer.next(dataS);
          uint16_t addr = ((dataC & 0x3F) << 8) | ( (dataS >> 8) & 0xFF );
          addHit(laneId, bc, reg, addr);
          uint8_t hit_map = (dataS & 0xFF);
          if ( hit_map & 0x80 )
          {
            log_error << "Wrong bit before DATA LONG bit map" << std::endl;
            assert(false);
          }
          while( hit_map != 0x00 )
          {
            ++addr;
            if ( hit_map & 1 )
            {
              addHit(laneId, bc, reg, addr);
            }
            hit_map >>= 1;
          }
        }
        else if ( (dataC & 0xF0) == 0xB0 ) // CHIP TRAILER
        {
//          uint8_t flag = (dataC & 0x0F);
          //TODO: YCM add chipdata statistic
          chip_trailer_found = true;
          busy_on = busy_off = chip_header_found = false;
        }
        else // ERROR
        {
          AlpideByteError(chipId, buffer);
        }
      }
      else
      {
        if ( (dataC & 0xF0) == 0xA0 ) // CHIP HEADER
        {
          if (! chip_trailer_found)
          {
            log_error << "New chip header found before a previous chip trailer" << std::endl;
          }
          chip_header_found = true;
          chip_trailer_found = false;
          laneId = (dataC & 0x0F) % 3;
          if (laneId != chipId )
          {
            log_error << "Error laneId " << laneId << " (" << (dataC & 0xF) << ") and chipId " << chipId << std::endl;
            assert(false);
          }
          buffer.next(bc);
          reg = 0xFF;
        }
        else if ( dataC == 0x00 ) // PADDING
        {
          continue;
        }
        else
        { // ERROR
          AlpideByteError(chipId, buffer);
        } // else !chip_header_found
      }  // if chip_header_found
    } // busy_on, busy_off, chip_empty, other
    if (busy_on)
    {
    }
    if (busy_off)
    {
    }
  }  // while

  return ret;
}


} // namespace mvtx

#endif // _MVTX_DECODER_ITSMFT_GBTLINK_H_
