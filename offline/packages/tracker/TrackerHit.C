#include "TrackerHit.h"

TrackerHit::TrackerHit()
  : hitid(~0x0)
{
}

TrackerHit::TrackerHit(TrackerDefs::keytype id)
  : hitid(id)
{
}

void TrackerHit::identify(ostream& os) const
{
  std::cout << "Class " << this->ClassName() << std::endl;
  return;
}

void TrackerHit::Reset()
{
  std::cout << "Reset not implemented by daughter class" << std::endl;
  return;
}

void TrackerHit::Copy(TrackerHit const& hit)
{
  set_hitid(hit.get_hitid());
}
