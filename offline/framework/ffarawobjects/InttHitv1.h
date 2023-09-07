#ifndef FUN4ALLRAW_INTTHITV1_H
#define FUN4ALLRAW_INTTHITV1_H

#include "InttHit.h"

#include <limits>


class  InttHitv1: public InttHit  
{


public:
  InttHitv1() {}
  ~InttHitv1() override {};

  uint64_t get_bco() const override {return bco;}
  void set_bco(const uint64_t val) override {bco = val;}

  uint32_t get_word() const override {return word;}
  void set_word(uint32_t val) override {word = val;}

  uint16_t get_fee() const override {return fee;}
  void set_fee(uint16_t val) override {fee = val;}

  uint16_t get_channel_id() const override {return channel_id;}
  void set_channel_id(uint16_t val) override {channel_id = val;}

  uint16_t get_chip_id() const override {return chip_id;}
  void set_chip_id(uint16_t val) override {chip_id = val;}

  uint16_t get_adc() const override {return adc;}
  void set_adc(uint16_t val) override {adc = val;}

  uint16_t get_FPHX_BCO() const override {return FPHX_BCO;}
  void set_FPHX_BCO(uint16_t val) override {FPHX_BCO = val;}

  uint16_t get_full_FPHX() const override {return full_FPHX;}
  void set_full_FPHX(uint16_t val) override {full_FPHX = val;}

  uint16_t get_full_ROC() const override {return full_ROC;}
  void set_full_ROC(uint16_t val) override {full_ROC = val;}

  uint16_t get_amplitude() const override {return amplitude;}
  void set_amplitude(uint16_t val) override {amplitude = val;}

  uint16_t get_full_fphx() const override {return full_fphx;}
  void set_full_fphx(uint16_t val) override {full_fphx = val;}


protected:
    uint64_t bco = std::numeric_limits<uint64_t>::max();
    uint32_t word = std::numeric_limits<uint32_t>::max();
    uint16_t fee = std::numeric_limits<uint16_t>::max();
    uint16_t channel_id = std::numeric_limits<uint16_t>::max();
    uint16_t chip_id = std::numeric_limits<uint16_t>::max();
    uint16_t adc = std::numeric_limits<uint16_t>::max();
    uint16_t FPHX_BCO = std::numeric_limits<uint16_t>::max();
    uint16_t full_FPHX = std::numeric_limits<uint16_t>::max();
    uint16_t full_ROC = std::numeric_limits<uint16_t>::max();
    uint16_t amplitude = std::numeric_limits<uint16_t>::max();
    uint16_t full_fphx = std::numeric_limits<uint16_t>::max();

  ClassDefOverride(InttHitv1,1)    
};

 
#endif 
