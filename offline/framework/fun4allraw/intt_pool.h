#ifndef __INTT_POOL_H__
#define __INTT_POOL_H__

#include <Event/packet.h>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <stdint.h>


class  intt_pool  {


public:
  intt_pool( const unsigned int required_depth=1000, const unsigned int low_mark =100);
  virtual ~intt_pool() {};

  virtual int addPacket( Packet *p);

  virtual void drain() { _low_mark = 0;};

  virtual unsigned int  rawValue(const int fee, const int index);

  virtual int iValue(const int hit, const int field);


  virtual int iValue(const int hit,const char * what);

  virtual long long  lValue(const int hit, const int field);
  virtual long long  lValue(const int hit,const char * what);

  //void  dump ( std::ostream& os = std::cout);
  virtual unsigned int min_depth() const; // the lowest vector length
  virtual bool depth_ok() const;
  virtual int next();

  virtual void  dump ( OSTREAM& os = std::cout);
  virtual int getIdentifier() const {return _myPacketid;};

  
  //int    iValue(const int , const int, const char * what);

protected:
  int intt_decode ();

  int intt_decode_hitlist (std::vector<unsigned int> & /*hitlist*/ , const int /*fee*/);

  static const int MAX_FEECOUNT =16;
  
  int _broken;
  
  int _is_decoded;
  int _last_accessed_hit;
  
  unsigned int _required_depth;
  unsigned int _low_mark;
  int _myPacketid;

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
    uint32_t event_counter;
    uint32_t word;
  };

    
  std::vector<unsigned int> fee_data[MAX_FEECOUNT];
  std::vector<intt_hit *> intt_hits;

  unsigned int last_index[MAX_FEECOUNT];

  
};

 
#endif 
