#ifndef __PACKET_CDEVBUCKETS_H__
#define __PACKET_CDEVBUCKETS_H__


#include <packet_w124.h>

/**
   This is the packet which deals with data in CDEVBUCKETS format.
   It inherits from Packet\_w4 because the data are 32bit entities.


*/



#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevbuckets : public Packet_w4 {
#else
class  Packet_cdevbuckets : public Packet_w4 {
#endif

public:
  Packet_cdevbuckets(PACKET_ptr);
  virtual void  dump ( OSTREAM& ) ;
  virtual double  dValue(const int channel,const char *what);
  virtual int  iValue(const int channel,const char *what);
/** with the "what" parameter you can decide which aspect of
the data is made available. This class is one of those which have
several different "kinds" of data; we use this to bring up 
the misc. .


With iValue(i, "measuredFillPattern")     returns the measuredFillPattern, and
iValue(i, "intendedFillPattern")          returns the intendedFillPattern.
With iValue(i, "polarizationFillPattern") returns the polarizationFillPattern,


\begin{verbatim}
dValue(i,"bunchLength")                ; gets the bunchLength
dValue(i,"fillPatternThreshold")       ; get  fillPatternThreshold
dValue(i,"bunchOneRelativePhase")      ; gets the bunch one relative phase
\end{verbatim}


*/
protected:
  virtual int *decode (int *);
  struct cdevBucketsData *ps;
  int decoded;
};

#endif /* __PACKET_CDEVBUCKETS_H__ */
