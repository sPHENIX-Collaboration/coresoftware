#include "EvtHeader.h"

namespace
{
  EvtHeader::Map dummy_map;
}

void 
EvtHeader::identify(std::ostream& os)
{
  os << "EvtHeader base class" << std::endl;
}

void
EvtHeader::Reset()
{
  std::cout << "EvtHeader::Reset()" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  gSystem->Exit(1);
  exit(1);
}

void
EvtHeader::AddBCO(TrkrDefs::TrkrId const& trkr_id, FEE_t const& fee, BCO_t const& bco)
{
  std::cout << "EvtHeader::AddBCO(TrkrDefs::TrkrId const& trkr_id, FEE_t const& fee, BCO_t const& bco)" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  gSystem->Exit(1);
  exit(1);
}

EvtHeader::ConstIterator
EvtHeader::GetBCO(TrkrDefs::TrkrId const& trkr_id, FEE_t const& fee, BCO_t const& bco) const
{
  std::cout << "EvtHeader::GetBCO(TrkrDefs::TrkrId const& trkr_id, FEE_t const& fee, BCO_t const& bco) const" << std::endl;
  std::cout << "\tUnimplemented (call to instance of base class)" << std::endl;
  std::cout << "\tExiting" << std::endl;

  gSystem->Exit(1);
  exit(1);

  return dummy_map.end();
}
