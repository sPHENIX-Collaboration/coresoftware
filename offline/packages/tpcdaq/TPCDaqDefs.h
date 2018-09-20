#ifndef TPC_TPCDAQDEFS_H
#define TPC_TPCDAQDEFS_H

namespace TPCDaqDefs
{
//! TPC v1 FEE test stand decoder
namespace FEEv1
{

static const unsigned int kPACKET_ID = 1024;
static const unsigned int kPACKET_LENGTH = 137;
static const unsigned int kN_CHANNELS = 256;
static const unsigned int kSAMPLE_LENGTH = 128;

}  // namespace FEEv1

}  // namespace TPCDaqDefs

#endif  //TPC_TPCDAQDEFS_H
