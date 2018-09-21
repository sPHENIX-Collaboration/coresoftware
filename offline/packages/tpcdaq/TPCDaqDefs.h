#ifndef TPC_TPCDAQDEFS_H
#define TPC_TPCDAQDEFS_H

#include <utility>      // std::pair, std::make_pair
#include <stdint.h>

namespace TPCDaqDefs
{
//! TPC v1 FEE test stand decoder
namespace FEEv1
{

static const unsigned int kPACKET_ID = 1024;
static const unsigned int kPACKET_LENGTH = 137;
static const unsigned int kN_CHANNELS = 256;
static const unsigned int kSAMPLE_LENGTH = 128;

static const unsigned int kMaxPadX = 50;
static const unsigned int kMaxPadY = 12;
std::pair<int,int> SAMPAChan2PadXY(uint32_t fee_channel);

}  // namespace FEEv1

}  // namespace TPCDaqDefs

#endif  //TPC_TPCDAQDEFS_H
