#ifndef __QPILEUPPOINT_H__
#define __QPILEUPPOINT_H__
#include "QPileUp.h"

class QPileUpPoint:public QPileUp {
 public:
  QPileUpPoint(float,float,float,float);
  virtual void Make();
 protected:
  float fRper;
  float fPper;
  float fZper;
  float fQabs;
};

#endif /* __QPILEUPPOINT_H__ */
