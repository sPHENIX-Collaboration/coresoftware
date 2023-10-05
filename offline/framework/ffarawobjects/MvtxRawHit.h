#ifndef FUN4ALLRAW_MVTXRAWTHIT_H
#define FUN4ALLRAW_MVTXRAWTHIT_H

#include <phool/PHObject.h>

#include <limits>


class  MvtxRawHit: public PHObject
  {


public:
    MvtxRawHit() = default;
  virtual ~MvtxRawHit() = default;

  virtual uint64_t get_bco() const {return std::numeric_limits<uint64_t>::max();}
  virtual void set_bco(const uint64_t) {return;}

  virtual uint32_t get_strobe_bc() const {return std::numeric_limits<uint32_t>::max();}
  virtual void set_strobe_bc(const uint32_t) {return;}

  virtual uint32_t get_chip_bc() const {return std::numeric_limits<uint32_t>::max();}
  virtual void set_chip_bc(const uint32_t) {return;}

  virtual uint8_t get_layer_id() const {return std::numeric_limits<uint8_t>::max();}
  virtual void set_layer_id(uint8_t) {return;}

  virtual uint8_t get_stave_id() const {return std::numeric_limits<uint8_t>::max();}
  virtual void set_stave_id(uint8_t) {return;}

  virtual uint8_t get_chip_id() const {return std::numeric_limits<uint8_t>::max();}
  virtual void set_chip_id(uint8_t) {return;}

  virtual uint16_t get_row() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_row(uint16_t) {return;}

  virtual uint16_t get_col() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_col(uint16_t) {return;}

private:
  ClassDefOverride(MvtxRawHit,1)
};

#endif 
