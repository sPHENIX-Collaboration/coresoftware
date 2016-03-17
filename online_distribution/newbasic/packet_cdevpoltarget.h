#ifndef __PACKET_CDEVPOLTARGET_H__
#define __PACKET_CDEVPOLTARGET_H__

#include <packet_w124.h>

/**
   This is the packet decoding the CDEV poltarget data.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevpoltarget : public Packet_w4 {
#else
class  Packet_cdevpoltarget : public Packet_w4 {
#endif

public:
  Packet_cdevpoltarget(PACKET_ptr);

  /**

The only supported interface of the PolarimeterTarget packet is the 
\begin{verbatim}
 packet->iValue(i,"WHAT") 
\end{verbatim}
call, which returns a int value.
The supported "WHAT" values are:
\begin{verbatim}
packet->iValue(i,"positionEncLinear")  ; returns 
packet->iValue(i,"positionEncRot")  ; returns

\end{verbatim}


  */
 

  int      iValue(const int channel,const char *what);


  void  dump ( OSTREAM& ) ;
  
protected:
  virtual int *decode (int *);
  struct cdevPolTargetData *ps;


};

#endif /* __PACKET_CDEVPOLTARGET_H__ */
