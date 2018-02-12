#ifndef __MvtxDefUtil_H__
#define __MvtxDefUtil_H__

#include <tracker/TrkrDefUtil.h>


class MvtxDefUtil : public TrkrDefUtil
{
public:

  /// ctor
  MvtxDefUtil() {};

  /// dtor
  ~MvtxDefUtil() {};

  /// get the stave id
  uint8_t get_staveid(TrkrDefs::hitsetkey key);
  uint8_t get_staveid(TrkrDefs::cluskey key);

  /// get the chip id
  uint8_t get_chipid(TrkrDefs::hitsetkey key);
  uint8_t get_chipid(TrkrDefs::cluskey key);

  /// generate hitsetkey (This is mvtx so tracker id is already known)
  TrkrDefs::hitsetkey gen_hitsetkey(const char lyr, const uint8_t stave, const uint8_t chip);

  /// generate cluskey
  TrkrDefs::cluskey gen_cluskey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid);


private:

  // hitsetkey layout:
  //  Mvtx specific lower 16 bits
  //   24 - 32  tracker id
  //   16 - 24  layer
  //   8  - 16  stave id
  //   0  -  8  chip id
  static const unsigned int bitshift_staveid = 8;  
  static const unsigned int bitshift_chipid = 0; 



};

#endif //__MvtxDefUtil_H__
