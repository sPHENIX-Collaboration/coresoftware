// @file PayLoadSG.h
// @brief Declaration of class for scatter-gather buffer
// @author ruben.shahoyan@cern.ch
// @sa <O2/Detectors/ITSMFT/common/reconstruction/include/ITSMFTReconstruction/PayLoadSG.h>
//     <1c31a8d52>

#ifndef MVTXDECODER_PAYLOADSG_H
#define MVTXDECODER_PAYLOADSG_H

#include <Rtypes.h>

#include <cstdint>
#include <vector>

namespace mvtx {

// scatter-gather buffer for the payload: base pointer + vector of references for pieces to collect
class PayLoadSG {
 public:
  PayLoadSG() = default;
  ~PayLoadSG() = default;

  enum HBF_ERRORS : uint8_t {
    NoError = 0x0,
    Incomplete = 0x1,
    PacketWithError = 0x2
  };

  ///< add n bytes to the buffer
  void add(size_t n, uint8_t err)
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
    //mBuffer.clear();
    std::vector<SGPiece>().swap(mBuffer);
    mBuffer.shrink_to_fit();
    mCurrentPieceId = 0;
  }

  struct SGPiece
  {
    uint32_t size = 0;   //size of the piece
    uint8_t hasError = 0;
    SGPiece() = default;
    SGPiece(int n, uint8_t err) : size(n), hasError(err) {}
  };

  void setDone() { mCurrentPieceId = mBuffer.size(); }

  size_t& currentPieceId() { return mCurrentPieceId; }
  size_t currentPieceId() const { return mCurrentPieceId; }

  const SGPiece* currentPiece() const
  {
    return mCurrentPieceId < mBuffer.size() ? &mBuffer[mCurrentPieceId] : nullptr;
  }

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

  ClassDefNV(PayLoadSG, 1);
};

} // namespace mvtx

#endif