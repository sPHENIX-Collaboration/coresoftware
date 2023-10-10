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

  virtual uint64_t get_gtm_bco() const {return std::numeric_limits<uint64_t>::max();}
  virtual void set_gtm_bco(const uint64_t) {return;}

  virtual int32_t get_packetid() const {return std::numeric_limits<int32_t>::max();}
  virtual void set_packetid(const int32_t) {return;}

  virtual uint16_t get_fee() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_fee(const uint16_t) {return;}

  virtual uint16_t get_channel() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_channel(const uint16_t) {return;}

  virtual uint16_t get_sampaaddress() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_sampaaddress(const uint16_t) {return;}

  virtual uint16_t get_sampachannel() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_sampachannel(const uint16_t) {return;}

  virtual uint16_t get_samples() const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_samples(const uint16_t) {return;}
  
  virtual uint16_t get_adc(size_t /*sample*/ ) const {return std::numeric_limits<uint16_t>::max();}
  virtual void set_adc(size_t /*sample*/, const uint16_t) { return; }

private:
  ClassDefOverride(TpcRawHit,1)
};

#endif 
