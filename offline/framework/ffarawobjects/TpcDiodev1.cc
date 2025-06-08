#include "TpcDiodev1.h"

TpcDiodev1::TpcDiodev1(TpcDiode *tpcdiode)
{
  // TpcDiodev1::set_bco(tpcdiode->get_bco());
  // TpcDiodev1::set_gtm_bco(tpcdiode->get_gtm_bco());
  TpcDiodev1::set_packetid(tpcdiode->get_packetid());
  // TpcDiodev1::set_fee(tpcdiode->get_fee());
  TpcDiodev1::set_channel(tpcdiode->get_channel());
  // TpcDiodev1::set_sampaaddress(tpcdiode->get_sampaaddress());
  // TpcDiodev1::set_sampachannel(tpcdiode->get_sampachannel());
  TpcDiodev1::set_samples(tpcdiode->get_samples());
  TpcDiodev1::set_maxadc(tpcdiode->get_maxadc());
  TpcDiodev1::set_maxbin(tpcdiode->get_maxbin());
  TpcDiodev1::set_integral(tpcdiode->get_integral());
  TpcDiodev1::set_pulsewidth(tpcdiode->get_pulsewidth());
  TpcDiodev1::set_nabovethreshold(tpcdiode->get_nabovethreshold());

  for (size_t i = 0; i < tpcdiode->get_samples(); ++i)
  {
    TpcDiodev1::set_adc(i, tpcdiode->get_adc(i));
  }
}

void TpcDiodev1::identify(std::ostream &os) const
{
  // os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}

void TpcDiodev1::Clear(Option_t * /*unused*/)
{
  adc = std::vector<uint16_t>();
}
