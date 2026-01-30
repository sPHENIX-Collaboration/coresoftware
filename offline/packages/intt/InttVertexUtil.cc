#include "InttVertexUtil.h"

#include <TF1.h>
#include <TMath.h>

#include <algorithm>
#include <numeric>
#include <vector>

namespace InttVertexUtil
{
  double gaus_func(double* x, double* par)
  {
    // note : par[0] : size
    // note : par[1] : mean
    // note : par[2] : width
    // note : par[3] : offset
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3];
  }

  std::vector<float> sigmaEff(std::vector<float> v, float threshold)
  {
    sort(v.begin(), v.end());

    int total = v.size();
    int max = (int) (threshold * total);

    std::vector<float> start;
    std::vector<float> stop;
    std::vector<float> width;

    unsigned i = 0;
    while (i != v.size() - 1)
    {
      int count = 0;
      unsigned j = i;
      while (j != v.size() - 1 && count < max)
      {
        ++count;
        ++j;
      }

      if (j != v.size() - 1)
      {
        start.push_back(v[i]);
        stop.push_back(v[j]);
        width.push_back(v[j] - v[i]);
      }
      ++i;
    }

    float minwidth = *min_element(width.begin(), width.end());

    unsigned pos = min_element(width.begin(), width.end()) - width.begin();

    float xmin = start[pos];
    float xmax = stop[pos];

    // cout<<"sigEffi test return width : "<<minwidth<<" edgel - edger : "<<xmin-xmax<<endl;

    return {minwidth, xmin, xmax};
  }

  float vector_average(std::vector<float> input_vector)
  {
    return accumulate(input_vector.begin(), input_vector.end(), 0.0) / double(input_vector.size());
  }

  std::vector<float> sigmaEff_avg(const std::vector<float>& v, float threshold)
  {
    std::vector<float> sigmaEff_vec = sigmaEff(v, threshold);

    std::vector<float> v_range;
    v_range.clear();

    // for (unsigned int i = 0; i < v.size(); i++){
    for (const auto& vi : v)
    {
      if (vi >= sigmaEff_vec[1] && vi <= sigmaEff_vec[2])
      {
        v_range.push_back(vi);
      }
    }

    // note : return the
    // note : 0. the avg among the sigeff range
    // note : 1. xmin
    // note : 2. xmax
    return {vector_average(v_range), sigmaEff_vec[1], sigmaEff_vec[2]};
  }

};  // namespace InttVertexUtil
