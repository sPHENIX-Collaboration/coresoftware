#ifndef FUN4ALLRAW_MVTXRAWHITV1_H
#define FUN4ALLRAW_MVTXRAWHITV1_H

#include "MvtxRawHit.h"

#include <limits>


class  MvtxRawHitv1: public MvtxRawHit
{


public:
  MvtxRawHitv1() {}
  MvtxRawHitv1(MvtxRawHit *mvtxhit);
  ~MvtxRawHitv1() override {};

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  uint64_t get_bco() const override {return bco;}
  void set_bco(const uint64_t val) override {bco = val;}

  uint32_t get_strobe_bc() const override {return strobe_bc;}
  void set_strobe_bc(const uint32_t val) override {strobe_bc = val;}

  uint32_t get_chip_bc() const override {return chip_bc;}
  void set_chip_bc(const uint32_t val) override {chip_bc = val;}

  uint8_t get_layer_id() const override {return layer_id;}
  void set_layer_id(uint8_t val) override {layer_id = val;}

  uint8_t get_stave_id() const override {return stave_id;}
  void set_stave_id(uint8_t val) override {stave_id = val;}

  uint8_t get_chip_id() const override {return chip_id;}
  void set_chip_id(uint8_t val) override {chip_id = val;}

  uint16_t get_row() const override {return row;}
  void set_row(uint16_t val) override {row = val;}

  uint16_t get_col() const override {return col;}
  void set_col(uint16_t val) override {col = val;}

protected:
    uint64_t bco = std::numeric_limits<uint64_t>::max();
    uint32_t strobe_bc = std::numeric_limits<uint32_t>::max();
    uint32_t chip_bc = std::numeric_limits<uint32_t>::max();
    uint8_t layer_id = std::numeric_limits<uint8_t>::max();
    uint8_t stave_id = std::numeric_limits<uint8_t>::max();
    uint8_t chip_id = std::numeric_limits<uint8_t>::max();
    uint16_t row = std::numeric_limits<uint16_t>::max();
    uint16_t col = std::numeric_limits<uint16_t>::max();

  ClassDefOverride(MvtxRawHitv1,1)
};

#endif 
