#include "TrackerHit.h"

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
