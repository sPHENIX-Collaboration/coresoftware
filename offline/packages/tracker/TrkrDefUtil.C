#include "TrkrDefUtil.h"
#include <bitset>


void 
TrkrDefUtil::print_bits(const TrkrDefs::hitsetkey key, std::ostream& os)
{
  os << "key: " << std::bitset<32>(key) << std::endl;
}

void 
TrkrDefUtil::print_bits(const TrkrDefs::cluskey key, std::ostream& os)
{
  os << "key: " << std::bitset<64>(key) << std::endl;
}

// void 
// TrkrDefUtil::print_bits(const TrkrDefs::hitkey key, std::ostream& os = std::cout)
// {
//   os << "key: " << std::bitset<32>(key) << std::endl;
// }

uint8_t 
TrkrDefUtil::get_trackerid(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_trackerid);
  return tmp;
}

uint8_t 
TrkrDefUtil::get_trackerid(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_clusid);
  return get_trackerid(tmp);
}

uint8_t 
TrkrDefUtil::get_layer(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_layer);
  return tmp;
}

uint8_t 
TrkrDefUtil::get_layer(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> bitshift_clusid);
  return get_layer(tmp);
}

uint32_t 
TrkrDefUtil::get_index(const TrkrDefs::cluskey key)
{
  return key;
}

TrkrDefs::hitsetkey 
TrkrDefUtil::get_hitsetkeylo(const TrkrDefs::TRKRID trkr_id)
{

}

TrkrDefs::hitsetkey 
TrkrDefUtil::get_hitsetkeyhi(const TrkrDefs::TRKRID trkr_id)
{

}

TrkrDefs::hitsetkey 
TrkrDefUtil::get_hitsetkeylo(const TrkrDefs::TRKRID trkr_id, const char lyr)
{

}

TrkrDefs::hitsetkey 
TrkrDefUtil::get_hitsetkeyhi(const TrkrDefs::TRKRID trkr_id, const char lyr)
{

}

TrkrDefs::cluskey 
TrkrDefUtil::get_cluskeylo(const TrkrDefs::TRKRID trkr_id)
{

}

TrkrDefs::cluskey 
TrkrDefUtil::get_cluskeyhi(const TrkrDefs::TRKRID trkr_id)
{

}

TrkrDefs::cluskey 
TrkrDefUtil::get_cluskeylo(const TrkrDefs::TRKRID trkr_id, const char lyr)
{

}

TrkrDefs::cluskey 
TrkrDefUtil::get_cluskeyhi(const TrkrDefs::TRKRID trkr_id, const char lyr)
{

}

TrkrDefs::hitsetkey 
TrkrDefUtil::gen_hitsetkey(const TrkrDefs::TRKRID trkr_id, const char lyr)
{
  TrkrDefs::hitsetkey tmp = trkr_id;
  TrkrDefs::hitsetkey key = tmp << bitshift_trackerid;  // detector id
  tmp = lyr;
  key |= (tmp << bitshift_layer);  // layer
  return key;
}
