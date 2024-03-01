// @file PayLoadSG.h
// @brief Declaration of class for scatter-gather buffer
// @author ruben.shahoyan@cern.ch
// @sa <O2/Detectors/ITSMFT/common/reconstruction/include/ITSMFTReconstruction/PayLoadSG.h>
//     <1c31a8d52>

#ifndef MVTXDECODER_PAYLOADSG_H
#define MVTXDECODER_PAYLOADSG_H

#include <cstdint>
#include <vector>

namespace mvtx
{

class PayLoadSG
{
  // scatter-gather buffer for the payload: base pointer + vector of references for pieces to collect
 public:

  PayLoadSG() = default;
  ~PayLoadSG() = default;

  ///< add n bytes to the buffer
  void add( size_t n, bool err )
  {
    if (n)
    {
      mBuffer.emplace_back(n, err);
    }
  }

  ///< move current pointer to the head
  void rewind()
  {
    mCurrentPieceId = 0;
  }

  ///< make buffer empty
  void clear()
  {
    mBuffer.clear();
    mCurrentPieceId = 0;
  }

  struct SGPiece
  {
    uint32_t size = 0;   // size of the piece
    bool hasError = false;
    SGPiece() = default;
    SGPiece(int n, bool err) : size(n), hasError(err) {}
  };

  void setDone() { mCurrentPieceId = mBuffer.size(); }

  size_t& currentPieceId() { return mCurrentPieceId; }
  size_t currentPieceId() const { return mCurrentPieceId; }

  const SGPiece* currentPiece() const { return mCurrentPieceId < mBuffer.size() ? &mBuffer[mCurrentPieceId] : nullptr; }

  const SGPiece* nextPiece()
  {
    // move to the next piece
    mCurrentPieceId++;
    return currentPiece();
  }

  const SGPiece* getPiece(int i) const { return &mBuffer[i]; }

  size_t getNPieces() const { return mBuffer.size(); }

 private:
  std::vector<SGPiece> mBuffer;   // list of pieces to fetch
  size_t mCurrentPieceId = 0;     // current piece

  //ClassDefNV(PayLoadSG, 1);
};

} // namespace mvtx

#endif
