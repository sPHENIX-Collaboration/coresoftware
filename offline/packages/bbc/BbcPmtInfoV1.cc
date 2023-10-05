#include "BbcPmtInfoV1.h"

/*
BbcPmtInfoV1::BbcPmtInfoV1(TowerInfo& tower)
{
  _time = (tower.get_time());
  _energy = (tower.get_energy());
}
*/

void BbcPmtInfoV1::Reset()
{
  Clear();
}

void BbcPmtInfoV1::Clear(Option_t* )
{
  bpmt = -1;
  bq = NAN;
  btt = NAN;
  btq = NAN;
}

void BbcPmtInfoV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a BbcPmtInfoV1 object" << std::endl;
  out << "Pmt: " << bpmt << ", q: " << bq << ", btt: "
      << btt << ", btq: " << btq << std::endl;
}
