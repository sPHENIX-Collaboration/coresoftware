#include "BcoInfov1.h"

void BcoInfov1::Reset()
{
  bco.fill(0);
  evtno.fill(0);
  return;
}

void BcoInfov1::identify(std::ostream& out) const
{
  out << "identify yourself: I am an BcoInfov1 Object\n";
  out << "previous event: " << get_previous_evtno() << std::hex
      << " bco: 0x" << get_previous_bco() << "\n"
      << std::dec
      << "current  event:  " << get_current_evtno() << std::hex
      << " bco: 0x" << get_current_bco() << "\n"
      << std::dec
      << "future   event:  " << get_future_evtno() << std::hex
      << " bco: 0x" << get_future_bco() << std::dec
      << std::endl;
  return;
}

int BcoInfov1::isValid() const
{
  return (bco[2] ? 1 : 0);  // return 1 if future bco is not zero
}
