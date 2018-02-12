#include "MvtxDefUtil.h"


uint8_t
MvtxDefUtil::get_staveid(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_staveid);
  return tmp;
}

uint8_t
MvtxDefUtil::get_staveid(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_clusid);
  return get_staveid(tmp);
}

uint8_t 
MvtxDefUtil::get_chipid(TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_chipid);
  return tmp;
}

uint8_t 
MvtxDefUtil::get_chipid(TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_clusid);
  return get_chipid(tmp);
}

TrkrDefs::hitsetkey 
MvtxDefUtil::gen_hitsetkey(const char lyr, const uint8_t stave, const uint8_t chip)
{
  TrkrDefs::hitsetkey key = TrkrDefUtil::gen_hitsetkey(TrkrDefs::TRKRID::mvtx_id, lyr);
  TrkrDefs::hitsetkey tmp = stave;
  key |= (tmp << bitshift_staveid);
  tmp = chip;
  key |= (tmp << bitshift_chipid);
  return key;
}

TrkrDefs::cluskey 
MvtxDefUtil::gen_cluskey(const char lyr, const uint8_t stave, const uint8_t chip, const uint32_t clusid)
{
  TrkrDefs::cluskey tmp = gen_hitsetkey(lyr, stave, chip);
  TrkrDefs::cluskey key = (tmp << bitshift_clusid);
  key |= clusid;
  return key;
}



