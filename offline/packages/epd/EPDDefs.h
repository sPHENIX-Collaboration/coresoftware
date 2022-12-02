#ifndef EPD_EPDDEFS_H
#define EPD_EPDDEFS_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"

namespace EPDDefs
{
  //  key layout:
  //  10-31 empty
  //  9   arm
  //  5-8 sector
  //  0-4 tile
  static unsigned int get_arm(uint32_t key)
  {
    return (key >> 9) & 0x1;
  }
  static unsigned int get_sector(uint32_t key)
  {
    return (key >> 5) & 0xF;
  }
  static unsigned int get_tileid(uint32_t key)
  {
    return (key) &0x1F;
  }
  static uint32_t make_epd_key(uint32_t arm, uint32_t sector, uint32_t tile_id)
  {
    return (arm << 9U | sector << 5U | tile_id) & 0x3FF;
  }
}  // namespace EPDDefs

#pragma GCC diagnostic pop
#endif
