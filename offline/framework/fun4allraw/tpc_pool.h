#ifndef __TPC_POOL_H__
#define __TPC_POOL_H__

#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <iostream>

class Packet;

class  tpc_pool  {

public:
  tpc_pool( const unsigned int required_depth, const unsigned int low_mark =1000000);
  virtual ~tpc_pool();

  virtual int addPacket( Packet *);

  int next();

  
  //! SAMPA waveform interfaces
  int    iValue(const int ch, const int sample);
  int    iValue(const int ,const char * what);

  //! Expose the Level 1 trigger and endat taggers
  long long  lValue(const int channel, const char *what) ;

  void  dump ( std::ostream& os = std::cout) ;

  virtual bool depth_ok() const;
  virtual void drain() { _low_mark = 0;};
  virtual unsigned int min_depth() const { return _thebuffer.size(); }; 


protected:
  int tpc_decode ();

  static const unsigned short  MAGIC_KEY_0 = 0xfe;
  static const unsigned short  MAGIC_KEY_1 = 0x00;

  static const unsigned short FEE_MAGIC_KEY = 0xba00;
  static const unsigned short GTM_MAGIC_KEY = 0xbb00;
  static const unsigned short GTM_LVL1_ACCEPT_MAGIC_KEY = 0xbbf0;
  static const unsigned short GTM_ENDAT_MAGIC_KEY = 0xbbf1;

  static const unsigned short  MAX_FEECOUNT = 26;   // that many FEEs
  static const unsigned short  MAX_CHANNELS   = 8*32; // that many channels per FEE
  static const unsigned short  HEADER_LENGTH  = 7;
  
  unsigned short reverseBits(const unsigned short x) const;
  unsigned short crc16(const unsigned int fee, const unsigned int index, const int  l) const;

  //  int find_header ( std::vector<unsigned short>::const_iterator &itr,  const std::vector<unsigned short> &orig);
  int find_header ( const unsigned int xx,  const std::vector<unsigned short> &orig);
  int decode_gtm_data(unsigned short gtm[16]);
  
//  int _broken {0};
  
  int _is_decoded;

  unsigned int _required_depth;
  unsigned int _low_mark;

  unsigned int _progress_index;



  
  std::vector<unsigned short> _thebuffer;
  
  struct sampa_waveform {
    unsigned short fee;
    unsigned short pkt_length;
    unsigned short channel;
    unsigned short sampa_channel;
    unsigned short sampa_address;
    unsigned int bx_timestamp;
    std::vector<unsigned short> waveform;
    unsigned short adc_length;
    unsigned short checksum;
    bool     valid;
  };

  struct gtm_payload {
      unsigned short pkt_type;
      bool is_endat;
      bool is_lvl1;
      unsigned long long bco;
      unsigned int lvl1_count;
      unsigned int endat_count;
      unsigned long long last_bco;
      unsigned char modebits;
  };
  
  // once vector per possible channel 16 cards * 256 channels
  //std::vector<sampa_waveform *> waveform_vector[MAX_FEECOUNT * MAX_CHANNELS];

  // our sort functional

struct bco_compare {
    bool operator() (const sampa_waveform *lhs, const sampa_waveform *rhs) const
    {
      return  ( lhs->bx_timestamp <= rhs->bx_timestamp );
    }
};


  
  typedef std::multiset< sampa_waveform* , bco_compare> waveform_set;
  //typedef waveform_set::iterator wf_iter;
  
  waveform_set waveforms;
  //  waveform_set waveforms[MAX_FEECOUNT * MAX_CHANNELS];

  int cacheIterator(const int n);
  
  waveform_set::iterator _cached_iter;
  int _last_requested_element;
  
  std::vector<unsigned short> fee_data[MAX_FEECOUNT];

  std::vector<gtm_payload *> gtm_data;

};

 
#endif 
