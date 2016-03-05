#ifndef __PACKET_CDEVPOLARIMETER_H__
#define __PACKET_CDEVPOLARIMETER_H__

#include <packet_w124.h>

/**
   This is the packet decoding the CDEV polarimeter data.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevpolarimeter : public Packet_w4 {
#else
class  Packet_cdevpolarimeter : public Packet_w4 {
#endif

public:
  Packet_cdevpolarimeter(PACKET_ptr);

  /**
Because the polarimeter data are a mixture of all kind of data,
the only supported interface is the 
\begin{verbatim}
 packet->rValue(i,"WHAT") 
\end{verbatim}
call, which returns a float value (most of the parameters are delivered as floats
from RHIC).

As a rule of thumb, if a field is just a single value, such as "avgAsymXS", then
the "i" parameter is ignored, so just p->rValue(0,"avgAsymXS") will do.

For the arrays, such as "encoderPositionS" (array of 2)  or "bunchAsymXS" (array of 360),
the index has to be in the right range, 0 or 1, or 0<= i < 360, respectively.

If an array is out of bounds, or the keyword is not recognized, the call returns 0.

Here are the keywords, which are just the (somewhat cryptic at times) names of the 
fields assigned by CA:

\begin{verbatim}

 "m_cdevCaptureTimeStamp"
 "runIdS"
 "startTimeS"
 "stopTimeS"
 "encoderPositionS"       (array of 2)
 "statusS"
 "totalCountsS"
 "upCountsS"
 "downCountsS"
 "unpolCountsS"
 "avgAsymXS"
 "avgAsymX45S"
 "avgAsymX90S"
 "avgAsymYS"
 "avgAsymErrorXS"
 "avgAsymErrorX45S"
 "avgAsymErrorX90S"
 "avgAsymErrorYS"
 "beamEnergyS"
 "analyzingPowerS"
 "analyzingPowerErrorS"
 "numberEventsS"
 "maxTimeS"
 "polarizationM"
 "bunchAsymXS"               (array of 360)
 "bunchAsymYS"               (array of 360)
 "bunchAsymErrorXS"          (array of 360)
 "bunchAsymErrorYS"          (array of 360)
 "countsUpLeftS"             (array of 360)
 "countsLeftS"               (array of 360)
 "countsDownLeftS"           (array of 360)
 "countsDownRightS"          (array of 360)
 "countsRightS"              (array of 360)
 "countsUpRightS"            (array of 360)
\end{verbatim}


  */
 


  double   dValue(const int channel,const char *what);
  float    rValue(const int channel,const char *what);
  int      iValue(const int channel,const char *what);


  void  dump ( OSTREAM& ) ;
  
protected:
  virtual int *decode (int *);
  struct cdevPolarimeterData *ps;
  int haspoldata;

};

#endif /* __PACKET_CDEVPOLARIMETER_H__ */
