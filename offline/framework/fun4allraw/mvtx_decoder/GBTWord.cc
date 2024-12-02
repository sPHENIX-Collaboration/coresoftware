// @file GBTWord.cxx
// @brief Classes for creation/interpretation of MVTX GBT data
// @sa <O2/Detectors/ITSMFT/common/reconstruction/src/GBTWord.cxx>
//     <1ecb2b0a2>

#include "mvtx_decoder/GBTWord.h"
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace mvtx;

//________________________________________________________________________________
std::string GBTWord::asString() const
{
  std::stringstream ss;
  for (int i = GBTWordLength; i--;)
  {
    ss << std::hex << std::setfill('0') << std::setw(2) << +data8[i] << " ";
  }
  return ss.str();
}

//________________________________________________________________________________
void GBTWord::printX() const
{
  /// print in right aligned hex format
  std::cout << asString() << std::endl;
}
