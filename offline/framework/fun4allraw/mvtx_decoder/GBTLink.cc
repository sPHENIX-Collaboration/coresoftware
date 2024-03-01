// @file GBTLink.cxx
// @brief Definitions of GBTLink class used for the ITS/MFT raw data decoding
// @sa <O2/Detectors/ITSMFT/common/reconstruction/src/GBTLink.cxx>
//     <e5b583efa>

#include "mvtx_decoder/GBTLink.h"

#include <bitset>

using namespace mvtx;

// using RDHUtils = mvtx::RDHUtils;
// using RDH = mvtx::RAWDataHeader;

///======================================================================
///                 GBT Link data decoding class
///======================================================================

///_________________________________________________________________
/// create link with given ids
GBTLink::GBTLink(uint16_t _flx, uint16_t _fee)
  : flxId(_flx)
  , feeId(_fee)
{
  //  chipStat.feeID = _fee;
  //  statistics.feeID = _fee;
  data.expand(RawBufferSize);
}

/////_________________________________________________________________
///// create string describing the link
// std::string GBTLink::describe() const
//{
//   std::string ss = fmt::format("link cruID:{:#06x}/lID{} feeID:{:#06x}", cruID, int(idInCRU), feeID);
//   if (lanes) {
//     ss += fmt::format(" lanes {}", std::bitset<28>(lanes).to_string());
//   }
//   return ss;
// }

///_________________________________________________________________
/// reset link
void GBTLink::clear(bool resetStat, bool resetTFRaw)
{
  if (data.isEmpty())
  {
    data.clear();
  }
  else
  {
    data.moveUnusedToHead();
    //    std::cout << "Link: " << feeID << " buffer Size: " << data.getSize() << std::endl;
  }

  if (resetTFRaw)
  {
    rawData.clear();
    mL1TrgTime.clear();
    mTrgData.clear();
    for (auto&& hit : hit_vector)
    {
      delete hit;
    }
    hit_vector.clear();
    dataOffset = 0;
    hbf_count = 0;
  }

  if (resetStat)
  {
    statistics.clear();
  }
  hbfEntry = 0;
  status = None;
}

///_________________________________________________________________
/// this function reads in 32 bytes  =  3 GBT words and 2 bytes
int GBTLink::readFlxWord(GBTWord* gbtwords, uint16_t& w16)
{
  for (uint8_t k = 0; k < 3; k++)
  {
    gbtwords[k] = *(reinterpret_cast<GBTWord*>(data.getPtr() + dataOffset));
    dataOffset += 10;
  }

  w16 = *(reinterpret_cast<uint16_t*>(data.getPtr() + dataOffset));
  dataOffset += 2;
  return 0;
}

//
/////_________________________________________________________________
// void GBTLink::printTrigger(const GBTTrigger* gbtTrg, int offs)
//{
//   std::bitset<12> trb(gbtTrg->triggerType);
//   LOG(info) << "Offs: " << offs << " Trigger : Orbit " << gbtTrg->orbit << " BC: " << gbtTrg->bc << " Trigger: " << trb << " noData:"
//             << gbtTrg->noData << " internal:" << gbtTrg->internal << " continuation:" << gbtTrg->continuation << " on " << describe();
//   gbtTrg->printX();
// }
//
/////_________________________________________________________________
// void GBTLink::printCalibrationWord(const GBTCalibration* gbtCal, int offs)
//{
//   LOGP(info, "Offs: {} Calibration word {:5} | user_data {:#08x} on {}", offs, gbtCal->calibCounter, gbtCal->calibUserField, describe());
//   gbtCal->printX();
// }
//
/////_________________________________________________________________
// void GBTLink::printHeader(const GBTDataHeader* gbtH, int offs)
//{
//   std::bitset<28> LA(gbtH->activeLanes);
//   LOG(info) << "Offs: " << offs << " Header : Active Lanes " << LA << " on " << describe();
//   gbtH->printX();
// }
//
/////_________________________________________________________________
////void GBTLink::printHeader(const GBTDataHeaderL* gbtH, int offs)
////{
////  std::bitset<28> LA(gbtH->activeLanesL);
////  LOG(info) << "Offs: " << offs << " HeaderL : Active Lanes " << LA << " on " << describe();
////  gbtH->printX(expectPadding);
////}
//
/////_________________________________________________________________
// void GBTLink::printTrailer(const GBTDataTrailer* gbtT, int offs)
//{
//   std::bitset<28> LT(gbtT->lanesTimeout), LS(gbtT->lanesStops); // RSTODO
//   LOG(info) << "Offs: " << offs << " Trailer: Done=" << gbtT->packetDone << " Lanes TO: " << LT << " | Lanes ST: " << LS << " on " << describe();
//   gbtT->printX();
// }
//
/////_________________________________________________________________
// void GBTLink::printDiagnostic(const GBTDiagnostic* gbtD, int offs)
//{
//   LOG(info) << "Offs: " << offs << " Diagnostic word on " << describe();
//   gbtD->printX();
// }
//
/////_________________________________________________________________
// void GBTLink::printCableDiagnostic(const GBTCableDiagnostic* gbtD)
//{
//   LOGP(info, "Diagnostic for {} Lane {} | errorID: {} data {:#018x} on {}", "IB", gbtD->getCableID(), gbtD->laneErrorID, gbtD->diagnosticData, describe());
//   gbtD->printX();
// }
//
/////_________________________________________________________________
// void GBTLink::printCableStatus(const GBTCableStatus* gbtS)
//{
//   LOGP(info, "Status data, not processed at the moment, on {}", describe());
//   gbtS->printX();
// }
//
/////====================================================================

#ifdef _RAW_READER_ERROR_CHECKS_

///_________________________________________________________________
/// Check RDH correctness
// uint8_t GBTLink::checkErrorsRDH(const RDH& rdh)
//{
//   uint8_t err = uint8_t(NoError);
//   if (!RDHUtils::checkRDH(rdh, true)) {
//     statistics.errorCounts[GBTLinkDecodingStat::ErrNoRDHAtStart]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrNoRDHAtStart])) {
//       err |= uint8_t(ErrorPrinted);
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrNoRDHAtStart];
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrNoRDHAtStart);
//     err |= uint8_t(Abort);
//     return err; // fatal error
//   }
//   /*
//   if ((RDHUtils::getPacketCounter(rdh) > packetCounter + 1) && packetCounter >= 0) {
//     if (irHBF.isDummy()) {
//       irHBF = RDHUtils::getHeartBeatIR(rdh);
//     }
//     statistics.errorCounts[GBTLinkDecodingStat::ErrPacketCounterJump]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrPacketCounterJump])) {
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrPacketCounterJump]
//                 << " : jump from " << int(packetCounter) << " to " << int(RDHUtils::getPacketCounter(rdh));
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrPacketCounterJump);
//     err |= uint8_t(Warning);
//   }
//   packetCounter = RDHUtils::getPacketCounter(rdh);
//   */
//   return err;
// }
//
///// Check encountered alignment padding word correctness //YCMTODO
// uint8_t GBTLink::checkErrorsAlignmentPadding()
//{
////  uint8_t err = uint8_t(NoError);
////  if (lastPageSize - dataOffset >= CRUPageAlignment) {
////    statistics.errorCounts[GBTLinkDecodingStat::ErrWrongAlignmentWord]++;
////    gbtErrStatUpadated = true;
////    if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrWrongAlignmentWord])) {
////      err |= uint8_t(ErrorPrinted);
////      LOG(info) << describe() << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrWrongAlignmentWord] << " at offset " << dataOffset << " for page size " << lastPageSize;
////    }
////    errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrWrongAlignmentWord);
////    err |= uint8_t(Warning);
////    return err; // fatal error
////  }
//  return err;
//}
//
/////_________________________________________________________________
///// Check RDH Stop correctness
// uint8_t GBTLink::checkErrorsRDHStop(const RDH& rdh)
//{
//   uint8_t err = uint8_t(NoError);
//   if (lastRDH && RDHUtils::getLHCBC(*lastRDH) != RDHUtils::getLHCBC(rdh)
//       && !RDHUtils::getStopBit(*lastRDH)) {
//     statistics.errorCounts[GBTLinkDecodingStat::ErrPageNotStopped]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrPageNotStopped])) {
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrPageNotStopped];
//       RDHUtils::printRDH(*lastRDH);
//       RDHUtils::printRDH(rdh);
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrPageNotStopped);
//     err |= uint8_t(Warning);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check if the RDH Stop page is empty
// uint8_t GBTLink::checkErrorsRDHStopPageEmpty(const RDH& rdh)
//{
//   uint8_t err = uint8_t(NoError);
//   if (RDHUtils::getStopBit(rdh) && RDHUtils::getPacketSize(rdh) != sizeof(RDH) + sizeof(FLXWord)) { // there could be only 1 diagnostic GBTWord after stop
//     statistics.errorCounts[GBTLinkDecodingStat::ErrStopPageNotEmpty]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrStopPageNotEmpty])) {
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrStopPageNotEmpty];
//       RDHUtils::printRDH(rdh);
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrStopPageNotEmpty);
//     err |= uint8_t(Warning);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check the GBT Trigger word correctness
// uint8_t GBTLink::checkErrorsTriggerWord(const GBTTrigger* gbtTrg)
//{
//   uint8_t err = uint8_t(NoError);
//   if (!gbtTrg->isTriggerWord()) { // check trigger word
//     statistics.errorCounts[GBTLinkDecodingStat::ErrMissingGBTTrigger]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrMissingGBTTrigger])) {
//       gbtTrg->printX();
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrMissingGBTTrigger];
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrMissingGBTTrigger);
//     err |= uint8_t(Abort);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check the GBT Calibration word correctness
// uint8_t GBTLink::checkErrorsCalibrationWord(const GBTCalibration* gbtCal)
//{
//   // at the moment do nothing
//   return uint8_t(NoError);
// }
//
/////_________________________________________________________________
///// Check the GBT Header word correctness
// uint8_t GBTLink::checkErrorsHeaderWord(const GBTDataHeader* gbtH)
//{
//   uint8_t err = uint8_t(NoError);
//   if (!gbtH->isDataHeader()) { // check header word
//     statistics.errorCounts[GBTLinkDecodingStat::ErrMissingGBTHeader]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrMissingGBTHeader])) {
//       gbtH->printX();
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrMissingGBTHeader];
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrMissingGBTHeader);
//     err |= uint8_t(Abort);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check active lanes status
// uint8_t GBTLink::checkErrorsActiveLanes(int cbl)
//{
//   uint8_t err = uint8_t(NoError);
//   if (~cbl & lanesActive) { // are there wrong lanes?
//     statistics.errorCounts[GBTLinkDecodingStat::ErrInvalidActiveLanes]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrInvalidActiveLanes])) {
//       std::bitset<32> expectL(cbl), gotL(lanesActive);
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrInvalidActiveLanes] << ' '
//                 << gotL << " vs " << expectL << " skip page";
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrInvalidActiveLanes);
//     err |= uint8_t(Warning);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check GBT Data word
// uint8_t GBTLink::checkErrorsGBTData(int cablePos)
//{
//   uint8_t err = uint8_t(NoError);
//   lanesWithData |= 0x1 << cablePos;    // flag that the data was seen on this lane
//   if (lanesStop & (0x1 << cablePos)) { // make sure stopped lanes do not transmit the data
//     statistics.errorCounts[GBTLinkDecodingStat::ErrDataForStoppedLane]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrDataForStoppedLane])) {
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrDataForStoppedLane] << cablePos;
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrDataForStoppedLane);
//     err |= uint8_t(Warning);
//   }
//
//   return err;
// }
//
/////_________________________________________________________________
///// Check GBT Data word ID: it might be diagnostic or status data
// uint8_t GBTLink::checkErrorsGBTDataID(const GBTData* gbtD)
//{
//   if (gbtD->isData()) {
//     return uint8_t(NoError);
//   }
//   uint8_t err = uint8_t(NoError);
//   statistics.errorCounts[GBTLinkDecodingStat::ErrGBTWordNotRecognized]++;
//   gbtErrStatUpadated = true;
//   if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrGBTWordNotRecognized])) {
//     if (gbtD->isCableDiagnostic()) {
//       printCableDiagnostic((GBTCableDiagnostic*)gbtD);
//     } else if (gbtD->isStatus()) {
//       printCableStatus((GBTCableStatus*)gbtD);
//     }
//     gbtD->printX();
//     LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrGBTWordNotRecognized];
//     err |= uint8_t(ErrorPrinted);
//   }
//   err |= uint8_t(Skip);
//   return err;
// }
//
/////_________________________________________________________________
///// Check the GBT Trailer word correctness
// uint8_t GBTLink::checkErrorsTrailerWord(const GBTDataTrailer* gbtT)
//{
//   uint8_t err = uint8_t(NoError);
//   if (!gbtT->isDataTrailer()) {
//     gbtT->printX();
//     statistics.errorCounts[GBTLinkDecodingStat::ErrMissingGBTTrailer]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrMissingGBTTrailer])) {
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrMissingGBTTrailer];
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrMissingGBTTrailer);
//     err |= uint8_t(Abort);
//     return err;
//   }
//   lanesStatus |= gbtT->lanesStatus; // register timeouts
//   return err;
// }
//
/////_________________________________________________________________
///// Check the Done status in GBT Trailer word
// uint8_t GBTLink::checkErrorsPacketDoneMissing(const GBTDataTrailer* gbtT, bool notEnd)
//{
//   uint8_t err = uint8_t(NoError);
//   if (!gbtT || (!gbtT->packetDone && notEnd)) { // Done may be missing only in case of carry-over to new CRU page
//     statistics.errorCounts[GBTLinkDecodingStat::ErrPacketDoneMissing]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrPacketDoneMissing])) {
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrPacketDoneMissing];
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrPacketDoneMissing);
//     err |= uint8_t(Warning);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check that all active lanes received their stop
// uint8_t GBTLink::checkErrorsLanesStops()
//{
//   // make sure all lane stops for finished page are received
//   uint8_t err = uint8_t(NoError);
//   if ((lanesActive & ~lanesStop)) {
//     if (RDHUtils::getTrgType(*lastRDH) != mvtx::TrgBitMap::SOT) { // only SOT trigger allows unstopped lanes?
//       statistics.errorCounts[GBTLinkDecodingStat::ErrUnstoppedLanes]++;
//       gbtErrStatUpadated = true;
//       if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrUnstoppedLanes])) {
//         std::bitset<32> active(lanesActive), stopped(lanesStop);
//         LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrUnstoppedLanes]
//                   << " | active: " << active << " stopped: " << stopped;
//         err |= uint8_t(ErrorPrinted);
//       }
//       errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrUnstoppedLanes);
//     }
//     err |= uint8_t(Warning);
//   }
//   // make sure all active lanes (except those in time-out) have sent some data
//   if ((~lanesWithData & lanesActive) != lanesTimeOut) {
//     statistics.errorCounts[GBTLinkDecodingStat::ErrNoDataForActiveLane]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrNoDataForActiveLane])) {
//       std::bitset<32> withData(lanesWithData), active(lanesActive), timeOut(lanesTimeOut);
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrNoDataForActiveLane]
//                 << " | with data: " << withData << " active: " << active << " timeOut: " << timeOut;
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrNoDataForActiveLane);
//     err |= uint8_t(Warning);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check diagnostic word
// uint8_t GBTLink::checkErrorsDiagnosticWord(const GBTDiagnostic* gbtD)
//{
//   uint8_t err = uint8_t(NoError);
//   if (RDHUtils::getPacketSize(lastRDH) != sizeof(RDH) + sizeof(FLXWord) || !gbtD->isDiagnosticWord()) { //
//     statistics.errorCounts[GBTLinkDecodingStat::ErrMissingDiagnosticWord]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrMissingDiagnosticWord])) {
//       gbtD->printX();
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrMissingDiagnosticWord];
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrMissingDiagnosticWord);
//     err |= uint8_t(Abort);
//   }
//   return err;
// }
//
/////_________________________________________________________________
///// Check cable ID validity
// uint8_t GBTLink::checkErrorsCableID(const GBTData* gbtD, uint8_t cableSW)
//{
//   uint8_t err = uint8_t(NoError);
//   if (cableSW == 0xff) {
//     statistics.errorCounts[GBTLinkDecodingStat::ErrWrongeCableID]++;
//     gbtErrStatUpadated = true;
//     if (needToPrintError(statistics.errorCounts[GBTLinkDecodingStat::ErrWrongeCableID])) {
//       gbtD->printX();
//       LOG(info) << describe() << ' ' << irHBF << ' ' << statistics.ErrNames[GBTLinkDecodingStat::ErrWrongeCableID] << ' ' << gbtD->getCableID();
//       err |= uint8_t(ErrorPrinted);
//     }
//     errorBits |= 0x1 << int(GBTLinkDecodingStat::ErrWrongeCableID);
//     err |= uint8_t(Skip);
//   }
//   return err;
// }

#endif
