#ifndef __PACKET_MUTRG_DCM0_H__
#define __PACKET_MUTRG_DCM0_H__


#include "packet_w124.h"

/**
   This is the packet which deals with data in MUTRG\_DCM0 format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_mutrg_dcm0 : public Packet_w4 {
#else
class  Packet_mutrg_dcm0 : public Packet_w4 {
#endif

public:
  Packet_mutrg_dcm0();
  Packet_mutrg_dcm0(PACKET_ptr);
  ~Packet_mutrg_dcm0();

/** with the "what" parameter you can decide which aspect of
the data is made available. This class is one of those which have
several different "kinds" of data; we use this to bring up 
the misc. items in the FEM headers and trailers.


In addition, there is 
\begin{verbatim}
  packet->iValue(0,"EVTNR")
  packet->iValue(0,"FLAG")
  packet->iValue(0,"DETID")
  packet->iValue(0,"MODADDR")
  packet->iValue(0,"BCLCK")
  packet->iValue(0,"USERWORD")
  packet->iValue(0,"PARITY")

\end{verbatim}
*/


  int    iValue(const int channel);
  int    iValue(const int channel,const char *what);
  void  dump ( OSTREAM& ) ;


protected:
  int *decode (int *);
  int *decode_misc (int *);
};

#endif /* __PACKET_MUTRG_DCM0_H__ */



