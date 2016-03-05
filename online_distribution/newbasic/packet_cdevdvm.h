#ifndef __PACKET_CDEVDVM_H__
#define __PACKET_CDEVDVM_H__


#include <packet_w124.h>

/**
   This is the packet which deals with data in CDEVDVM format.
   It inherits from Packet\_w4 because the data are 32bit entities.


*/
/** with the "what" parameter you can decide which aspect of
the data is made available. This class is one of those which have
several different "kinds" of data; we use this to bring up 
the misc. .

\begin{verbatim}
packet->dValue(i, "beamCurrent") gives you the  beam Current    - double
packet->dValue(i, "beamLifeTime") you get the   beam lifetime   - double
packet->fValue(i, "beamLifeTime") you get the   beam lifetime   - float

\end{verbatim}

*/



#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevdvm : public Packet_w4 {
#else
class  Packet_cdevdvm : public Packet_w4 {
#endif

public:
  Packet_cdevdvm(PACKET_ptr);
  virtual void  dump ( OSTREAM& ) ;
  virtual double dValue(const int channel,const char *what);
  virtual float  fValue(const int channel,const char *what); // ?? TKH ??
  virtual float  rValue(const int channel,const char *what);

protected:
  virtual int *decode (int *);
  struct cdevDvmData *ps;
  int decoded;
};

#endif /* __PACKET_CDEVDVM_H__ */
