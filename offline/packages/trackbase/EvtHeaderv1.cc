#include "EvtHeaderv1.h"

void 
EvtHeaderv1::identify(std::ostream& os) const
{
  os << "EvtHeader version 1" << std::endl;
}

void
EvtHeaderv1::Reset()
{
  bco_map.clear();
}

void
EvtHeaderv1::AddBCO(TrkrDefs::TrkrId const& trkr_id, FEE_t const& fee, BCO_t const& bco)
{
  bco_map[std::make_pair(trkr_id, fee)] = bco;
}

EvtHeaderv1::ConstIterator
EvtHeaderv1::GetBCO(TrkrDefs::TrkrId const& trkr_id, FEE_t const& fee) const
{
  return bco_map.find(std::make_pair(trkr_id, fee));
}
