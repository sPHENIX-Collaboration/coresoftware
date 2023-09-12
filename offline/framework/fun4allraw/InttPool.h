#ifndef FUN4ALLRAW_INTTPOOL_H__
#define FUN4ALLRAW_INTTPOOL_H__

#include <Event/packet.h>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <stdint.h>


class  InttPool  {


public:
  InttPool( const unsigned int depth=100);
  virtual ~InttPool() {};

  virtual int addPacket( Packet *p);

  virtual unsigned int  rawValue(const int fee, const int index);

  virtual int iValue(const int hit, const int field);


  virtual int iValue(const int hit,const char * what);

  virtual long long  lValue(const int hit, const int field);
  virtual long long  lValue(const int hit,const char * what);

  //void  dump ( std::ostream& os = std::cout);
  virtual unsigned int min_depth() const; // the lowest vector length
  virtual bool depth_ok() const;
  virtual
  int next();

  //int    iValue(const int , const int, const char * what);

protected:
  int intt_decode ();

  
  static const int MAX_FEECOUNT =16;
  
//  int _broken = 0;
  
  int _is_decoded = 0;
  
  unsigned int _depth;
  
  struct intt_hit
  {
    uint64_t bco;
    uint16_t fee;
    uint16_t channel_id;
    uint16_t chip_id;
    uint16_t adc;
    uint16_t FPHX_BCO;
    uint16_t full_FPHX;
    uint16_t full_ROC;
    uint16_t amplitude;
    uint16_t full_fphx;
    uint32_t word;
  };

    
  std::vector<unsigned int> fee_data[MAX_FEECOUNT];
  std::vector<intt_hit *> intt_hits;

};

 
#endif 
