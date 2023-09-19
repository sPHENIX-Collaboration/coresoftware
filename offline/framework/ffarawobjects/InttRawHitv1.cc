#include "InttRawHitv1.h"

InttRawHitv1::InttRawHitv1(InttRawHit *intthit)
{
  set_bco(intthit->get_bco());
  set_packetid(intthit->get_packetid());
  set_word(intthit->get_word());
  set_fee(intthit->get_fee());
  set_channel_id(intthit->get_channel_id());
  set_chip_id(intthit->get_chip_id());
  set_adc(intthit->get_adc());
  set_FPHX_BCO(intthit->get_FPHX_BCO());
  set_full_FPHX(intthit->get_full_FPHX());
  set_full_ROC(intthit->get_full_ROC());
  set_amplitude(intthit->get_amplitude());
}

void InttRawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}
