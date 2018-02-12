#include "TrkrHitSet.h"

void TrkrHitSet::identify(ostream& os) const
{
  std::cout << "Class " << this->ClassName() << std::endl;
  return;
}

void TrkrHitSet::Reset()
{
  std::cout << "Reset not implemented by daughter class" << std::endl;
  return;
}

void TrkrHitSet::Copy(TrkrHitSet const& hit)
{
  set_hitid(hit.get_hitid());
  set_truthid(hit.get_truthid());
}
