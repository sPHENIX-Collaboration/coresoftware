#ifndef __PACKET_CDEVMADCH_H__
#define __PACKET_CDEVMADCH_H__

#include <packet_w124.h>

/**
   This is the packet decoding the CDEV MADCH data.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevmadch : public Packet_w4 {
#else
class  Packet_cdevmadch : public Packet_w4 {
#endif

public:
  Packet_cdevmadch(PACKET_ptr);

  /**

  
   long   m_avgOrbTimeStamp;
  float  m_avgOrbPosition;
  float  m_avgOrbVariance;
  float  m_avgOrbStat;
  */

  int    iValue(const int channel,const char *what);
  double    dValue(const int channel,const char *what);
  //  float    rValue(const int channel,const int y);


  void  dump ( OSTREAM& ) ;
  
protected:
  virtual int *decode (int *);
  struct cdevMadchData  *ps;
  int no_structures;

};

#endif /* __PACKET_CDEVMADCH_H__ */
