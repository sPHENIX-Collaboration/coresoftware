#ifndef __PACKET_CDEVWCM_H__
#define __PACKET_CDEVWCM_H__

#include <packet_w124.h>

/**
   This is the packet decoding the CDEV WCM data.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevwcm : public Packet_w4 {
#else
class  Packet_cdevwcm : public Packet_w4 {
#endif

public:
  Packet_cdevwcm(PACKET_ptr);

  /**

  */
 


  int    iValue(const int channel,const char *what);
  float    rValue(const int channel,const char *what);
  float    rValue(const int channel,const int y);


  void  dump ( OSTREAM& ) ;
  
protected:
  virtual int *decode (int *);
  int numberofreadings;
  struct cdevWCMHistory *ps;


};

#endif /* __PACKET_CDEVWCM_H__ */
