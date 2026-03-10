#include "BcoInfov1.h"

void BcoInfov1::Reset()
{
  bco.fill(0);
  return;
}

void BcoInfov1::identify(std::ostream& out) const
{
  out << "identify yourself: I am an BcoInfov1 Object\n";
  out << std::hex;
  out << "bco previous event: 0x" << get_previous_bco() << "\n"
      << "bco current  event: 0x" << get_current_bco() << "\n"
      << "bco future   event: 0x" << get_future_bco()
      << std::dec
  << std::endl;

  return;
}

int BcoInfov1::isValid() const
{
  return (bco[2] ? 1 : 0);  // return 1 if future bco is not zero
}


