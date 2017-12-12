#include "TrackerHit.h"

void
TrackerHit::identify(ostream& os) const
{
  cout << "Class " << this->ClassName() << endl;
  return;
}

void
TrackerHit::Reset()
{
  cout << "Reset not implemented by daughter class" << endl;
  return;
}

void 
TrackerHit::Copy(TrackerHit const &hit)
{
  set_hitid( hit.get_hitid() );
}
