#ifndef __PACKET_CDEVRINGNOPOL_H__
#define __PACKET_CDEVRINGNOPOL_H__


#include <packet_w124.h>

/**
   This is the packet which deals with data in CDEVRING format.
   It inherits from Packet\_w4 because the data are 32bit entities.


*/



#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevringnopol : public Packet_w4 {
#else
class  Packet_cdevringnopol : public Packet_w4 {
#endif

public:
  Packet_cdevringnopol(PACKET_ptr);
  virtual void  dump ( OSTREAM& ) ;
  virtual double  dValue(const int channel,const char *what);
  virtual int  iValue(const int channel,const char *what);
/** with the "what" parameter you can decide which aspect of
the data is made available. This class is one of those which have
several different "kinds" of data; we use this to bring up 
the misc. .



iValue(i, "fillNumber")  returns the fillNumber

\begin{verbatim}
dValue(i,"beamEnergy")  ; gets the beamEnergy
dValue(i,"gamma")       ; get  gamma
dValue(i,"betaPMHorizontal");
dValue(i,"betaPMVertical");
dValule(i,"chromaticityVertical");
dValue(i,"chromaticityHorizontal");
\end{verbatim}

In addition, there is 
\begin{verbatim}
 packet->iValue(0,"stoneType")    The stoneType
 packet->iValue(0,"timeOfFillStart")   The fill Start Time
 packet->iValue(0,"timeOfLuminosityStart")     The luminosity start time
 

\end{verbatim}
*/
protected:
  virtual int *decode (int *);
  struct cdevRingNoPolData *ps;
  int hasfilldata;
  int decoded;
};

#endif /* __PACKET_CDEVIR_H__ */
