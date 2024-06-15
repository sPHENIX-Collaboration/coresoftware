// @file PixelData.cxx
// @brief Implementation for transient data of single pixel and set of pixels from current chip
// @sa <O2/Detectors/ITSMFT/common/reconstruction/src/PixelData.cxx>
//     <d44292025>

#include "mvtx_decoder/PixelData.h"
#include <iostream>
#include <sstream>

using namespace mvtx;

//________________________________________________________________________________
void ChipPixelData::print() const
{
  // print chip data
  std::cout << "Chip " << mChipID << " | " << mPixels.size() << " hits" << std::endl;
  for (std::size_t i = 0; i < mPixels.size(); i++)
  {
    std::cout << "#" << i << " C:" << mPixels[i].getCol() << " R: " << mPixels[i].getRow() << std::endl;
  }
}

//________________________________________________________________________________
std::string ChipPixelData::getErrorDetails(int pos) const
{
  // if possible, extract more detailed info about the error
  if (pos == int(ChipStat::RepeatingPixel))
  {
    std::stringstream ss;
    ss << ": " << (mErrorInfo & 0xffffU) << "/" << ((mErrorInfo >> 16U) & 0xffffU);
    return ss.str();
  }
  if (pos == int(ChipStat::UnknownWord))
  {
    std::string rbuf = ": 0x<";
    int nc = getNBytesInRawBuff();
    for (int i = 0; i < nc; i++)
    {
      std::stringstream ss;
      ss << (i ? " " : "") << std::hex << (int) getRawErrBuff()[i];
    }
    rbuf += '>';
    return rbuf;
  }
  return {};
}
