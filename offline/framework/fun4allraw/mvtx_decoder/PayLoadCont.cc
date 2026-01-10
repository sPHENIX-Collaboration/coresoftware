// @file PayLoadCont.cxx
// @brief Implementation of class for continuos buffer of ALPIDE data
// @sa <O2/Detectors/ITSMFT/common/reconstruction/src/PayLoadCont.cxx>
//     <03608ff89>

#include "PayLoadCont.h"

#include <algorithm>

//________________________________________________________________________________
mvtx::PayLoadCont::PayLoadCont(const PayLoadCont& src)
  :mBuffer(src.mBuffer)
{
  if (src.mPtr)
  {
    mPtr = mBuffer.data() + (src.mPtr - src.mBuffer.data());
  }
  if (src.mEnd)
  {
    mEnd = mBuffer.data() + (src.mEnd - src.mBuffer.data());
  }
}

//________________________________________________________________________________
mvtx::PayLoadCont& mvtx::PayLoadCont::operator=(const PayLoadCont& src)
{
  if (&src != this)
  {
    mBuffer = src.mBuffer;
    if (src.mPtr)
    {
      mPtr = mBuffer.data() + (src.mPtr - src.mBuffer.data());
    }
    if (src.mEnd)
    {
      mEnd = mBuffer.data() + (src.mEnd - src.mBuffer.data());
    }
  }
  return *this;
}

//________________________________________________________________________________
void mvtx::PayLoadCont::expand(size_t sz)
{
  ///< increase the buffer size
  auto* oldHead = mBuffer.data();
  sz = std::max(sz, MinCapacity);
  if (sz < mBuffer.size())
  {  // never decrease the size
    return;
  }
  mBuffer.resize(sz);
  if (oldHead)
  {  // fix the pointers to account for the reallocation
    int64_t diff = mBuffer.data() - oldHead;
    mPtr += diff;
    mEnd += diff;
  }
  else
  {  // new buffer
    clear();
  }
}
