#ifndef FUN4ALLRAW_MVTXRAWEVTHEADER_H
#define FUN4ALLRAW_MVTXRAWEVTHEADER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <set>
#include <limits>

class  MvtxRawEvtHeader: public PHObject
{

public:
  MvtxRawEvtHeader() = default;
  virtual ~MvtxRawEvtHeader() = default;

  virtual void AddFeeId(const int&) { return; };
  virtual void AddL1Trg(const uint64_t&) { return; };

  virtual void AddFeeId(const std::set<uint16_t>&) { return; }
  virtual void AddL1Trg(const std::set<uint64_t>&) { return; }

private:
  ClassDefOverride(MvtxRawEvtHeader, 1)
};

#endif
