#include "Gl1RawHit.h"

void Gl1RawHit::CopyFrom(Gl1RawHit *gl1hit)
{
  set_bco(gl1hit->get_bco());
  setEvtSequence(gl1hit->getEvtSequence());
  return;
}
