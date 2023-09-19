#ifndef FUN4ALLRAW_TPCRAWTHIT_H
#define FUN4ALLRAW_TPCRAWTHIT_H

#include <phool/PHObject.h>

#include <limits>


class  TpcRawHit: public PHObject
  {


public:
    TpcRawHit() = default;
  virtual ~TpcRawHit() = default;

  virtual uint64_t get_bco() const {return std::numeric_limits<uint64_t>::max();}
  virtual void set_bco(const uint64_t) {return;}

  virtual int32_t get_packetid() const {return std::numeric_limits<int32_t>::max();}
  virtual void set_packetid(const int32_t) {return;}

  virtual uint16_t get_fee() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_fee(uint16_t) {return;}

  virtual uint16_t get_channel() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_channel(uint16_t) {return;}

private:
  ClassDefOverride(TpcRawHit,1)
};

#endif 
