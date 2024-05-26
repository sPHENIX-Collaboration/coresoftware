#ifndef INTTVERTEXUTIL_H
#define INTTVERTEXUTIL_H

#include <TMath.h>
#include <TF1.h>

namespace InttVertexUtil {
  double gaus_func(double *x, double *par);
  std::vector<float> sigmaEff_avg (std::vector<float> v, float threshold);
};

#endif

