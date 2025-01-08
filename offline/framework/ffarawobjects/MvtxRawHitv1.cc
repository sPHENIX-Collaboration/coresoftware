#include "MvtxRawHitv1.h"

MvtxRawHitv1::MvtxRawHitv1(MvtxRawHit *mvtxhit)
{
  set_bco(mvtxhit->get_bco());
  set_strobe_bc(mvtxhit->get_strobe_bc());
  set_chip_bc(mvtxhit->get_chip_bc());
  set_layer_id(mvtxhit->get_layer_id());
  set_stave_id(mvtxhit->get_stave_id());
  set_chip_id(mvtxhit->get_chip_id());
  set_row(mvtxhit->get_row());
  set_col(mvtxhit->get_col());
}

void MvtxRawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "Strobe BC: " << strobe_bc << ", Chip BC: " << chip_bc << std::endl;
  os << "Layer: " << (unsigned) layer_id << ", stave: " << (unsigned) stave_id << ", chip: " << (unsigned) chip_id << std::endl;
  os << "Row: " << row << ", column: " << col << std::endl;
}
