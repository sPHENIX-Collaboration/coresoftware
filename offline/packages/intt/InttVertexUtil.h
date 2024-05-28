#ifndef INTTVERTEXUTIL_H
#define INTTVERTEXUTIL_H

#include <TF1.h>
#include <TMath.h>

namespace InttVertexUtil {
  double gaus_func(double *x, double *par);
  std::vector<float> sigmaEff_avg (const std::vector<float>& v, float threshold);
};

#endif

