#ifndef __PACKET_ID4SCALER_H__
#define __PACKET_ID4SCALER_H__


#include "packet_w124.h"

/**
   This is the packet which deals with data in ID4SCALER format.
   It inherits from Packet\_w4 because the data are 32bit entities.


Since there is a variety of data contained in that packet, the main interaction with this packet 
is through the iValue (channel,"STRING") interface.

The packet contains a num,ber of trigger masks which were valid for this partition.


  iValue(0,"NUMBERMASKS") returns the number of masks
  iValue(i,"TRIGGERMASK") returns the triggermask i, where 0 <= i < iValue(0,"NUMBERMASKS")

Similarly, 

  iValue(0,"NUMBERSCALERS") returns the number of scalers valid for this partition.

Then you have 3 sets of scalers.

  iValue(i,"RAWSCALERS")
  iValue(i,"LIFESCALERS")
  iValue(i,"SCALEDSCALERS")

where 0 <= i <  iValue(0,"NUMBERSCALERS")

which return the respective scaler values. 

Finally, a more complicated one is the iValue(i,"TIMESTRING") interface. Similar to the IDCSTR
packet. It returns a character string losely modeled after the standard UNIX getc(FILE* fp) interface.

When the string is ended, it returns EOF, so if you are familiar with the construct

  int c; 
  while ( ( c = getc(fp) ) !=EOF ) *str++ = c;

you'll see the similarity with 

  int i = 0;
  while (  ( c = p->iValue(i++,"TIMESTRING") ) != EOF ) *str++ = c;

(p is the pointer to this packet) 

In the end, you will have a regular 0-terminated C character string.

*/



#ifndef __CINT__
class WINDOWSEXPORT Packet_id4scaler : public Packet_w4 {
#else
class  Packet_id4scaler : public Packet_w4 {
#endif

public:
  Packet_id4scaler(PACKET_ptr);
  int    iValue(const int channel,const char *what);
  long long  lValue(const int channel, const char *what);

  virtual void  dump ( OSTREAM& ) ;

protected:
  virtual int *decode (int *);
};

#endif /* __PACKET_ID4SCALER_H__ */





