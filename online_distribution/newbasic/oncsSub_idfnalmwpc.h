#ifndef __ONCSSUB_IDFNALMWPC_H__
#define __ONCSSUB_IDFNALMWPC_H__

#include <vector>

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idfnalmwpc : public  oncsSubevent_w4 {
#else
class oncsSub_idfnalmwpc : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_idfnalmwpc( subevtdata_ptr);
  ~oncsSub_idfnalmwpc();


  int    iValue(const int n, const char *);


  // for a given trigger, information about the TDC
  int    iValue(const int trigger, const int tdc, const char * what);

  // this returns the wire number
  int    iValue(const int trigger, const int tdc, const int index);

  // this is the workhorse interface
  int    iValue(const int trigger, const int tdc, const int index, const char *what);


  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);

  int _decoded;

#define MAXNROFTDCS 16


  typedef struct {
    int wire;
    int timestamp;
  } TDC_hit;


  typedef struct {
    int hits;
    short words;
    short TDC;
    short EventStatus;
    short trigger_nr;
    short triggertype;
    unsigned short evt_timestamp;
    unsigned short local_ts_upper;
    unsigned short local_ts_lower;
    unsigned long long absolute_time;
    std::vector<TDC_hit *> TDCHitlist;
  } TDCData;


  typedef struct {
    int trigger_nr;
    unsigned short evt_timestamp;
    unsigned short local_ts_upper;
    unsigned short local_ts_lower;
    unsigned long long absolute_time;
    TDCData td[MAXNROFTDCS];
  } TDCEvent;

  typedef struct {
    int spillwords;
    int spillcounter;
    short year;
    short month;
    short day;
    short hour;
    short minute;
    short second;
    int triggercount;
    short TDC_Status_Bits;
    short Spill_Link_Status_Bits;
  } SpillInfo;


  typedef struct {
    short spillwords;
    short TDC;
    short spilltriggercount;
    short spillstatus;
  } TDCspillheader;


  int tdcmask;
  int n_tdcs;
  int length;

  SpillInfo spillinfo;

  TDCspillheader tsh[MAXNROFTDCS];


  std::vector<TDCEvent *> TDCEventVector;

};



#endif 
