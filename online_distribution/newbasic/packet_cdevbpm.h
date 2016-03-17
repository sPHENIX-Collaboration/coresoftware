#ifndef __PACKET_CDEVBPM_H__
#define __PACKET_CDEVBPM_H__

#include <packet_w124.h>

/**
   This is the packet decoding the CDEV BPM data.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevbpm : public Packet_w4 {
#else
class  Packet_cdevbpm : public Packet_w4 {
#endif

public:
  Packet_cdevbpm(PACKET_ptr);
/** with the "what" parameter you can decide which aspect of
the data is made available. This class is one of those which have
several different "kinds" of data; we use this to bring up 
the misc. .

\begin{verbatim}
packet->iValue(i,"NOREADINGS")       returns the number of bpm devices
packet->iValue(i, "avgOrbTimeStamp") you get the avgOrbTimeStamp - int
packet->rValue(i, "avgOrbPosition") gives you the avgOrbPosition - float
packet->rValue(i, "avgOrbVariance") you get the avgOrbVariance   - float
packet->rValue(i, "avgOrbStat")     you get the avgOrbStat       - float
\end{verbatim}

*/


  int    iValue(const int channel,const char *what);
  float    rValue(const int channel,const char *what);
  //  float    rValue(const int channel,const int y);


  void  dump ( OSTREAM& ) ;
  
protected:
  virtual int *decode (int *);
  struct cdevBPMData *ps;
  int no_structures;

};

#endif /* __PACKET_CDEVBPM_H__ */
