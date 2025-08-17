#include "AfterburnerAlgo.h"

#include <iostream>
#include <cmath>
#include <cstdlib>  // for exit


AfterburnerAlgo::AfterburnerAlgo( flowAfterburnerAlgorithm algorithm)
    : m_algorithm(algorithm)
{
  for (unsigned int i = 0; i < 6; ++i)
  {
    m_vn[i] = 0.0f;
    m_vn_scalefactors[i] = 1.0f;
  }
}

void AfterburnerAlgo::set_single_scale_N( const unsigned int n, const float scale )
{
  if (n < 1 || n > 6)
  {
    std::cerr << "AfterburnerAlgo::set_single_scale_N: n must be between 1 and 6, got " << n << std::endl;
    return;
  }
  m_vn_scalefactors[n - 1] = scale;
}

void AfterburnerAlgo::set_scale_all( const float scale )
{
  for (unsigned int i = 0; i < 6; ++i)
  {
    m_vn_scalefactors[i] = scale;
  }
}

float AfterburnerAlgo::get_vn(unsigned int n) const
{
  if (n < 1 || n > 6)
  {
    std::cerr << "AfterburnerAlgo::get_vn: n must be between 1 and 6, got " << n << std::endl;
    return 0.0f;
  }
  return m_vn[n - 1];
}


float AfterburnerAlgo::calc_v2(double b, double eta, double pt)
{
    
    float a1 = 0.4397 * exp(-(b - 4.526) * (b - 4.526) / 72.0) + 0.636;
    float a2 = 1.916 / (b + 2) + 0.1;
    float a3 = 4.79 * 0.0001 * (b - 0.621) * (b - 10.172) * (b - 23) + 1.2;  // this is >0 for b>0
    float a4 = 0.135 * exp(-0.5 * (b - 10.855) * (b - 10.855) / 4.607 / 4.607) + 0.0120;

    float temp1 = pow(pt, a1) / (1 + exp((pt - 3.0) / a3));
    float temp2 = pow(pt + 0.1, -a2) / (1 + exp(-(pt - 4.5) / a3));
    float temp3 = 0.01 / (1 + exp(-(pt - 4.5) / a3));

    // Adjust flow rapidity dependence to better match PHOBOS 200 GeV Au+Au data
    // JGL 9/9/2019
    // See JS ToG talk at https://indico.bnl.gov/event/6764/
    float val = (a4 * (temp1 + temp2) + temp3) * exp(-0.5 * eta * eta / 3.43 / 3.43);
    return val;
}

void AfterburnerAlgo::calc_flow(double eta, double pt)
{
    
    float v1 = 0, v2 = 0, v3 = 0, v4 = 0, v5 = 0, v6 = 0;
    if ( m_algorithm == custom_algorithm )
    { // Custom flow parameters
        v1 = 0.0000f;
        v2 = 0.0500f;
        v3 = 0.0280f;
        v4 = 0.0130f;
        v5 = 0.0045f;
        v6 = 0.0015f;
    }
    else if (m_algorithm == minbias_algorithm)
    { // all other algorithms need to calculate the flow
        v1 = 0;
        v2 = AfterburnerAlgo::calc_v2(m_impact_parameter, eta, pt);
        
        float fb = (
            0.97 
            + (1.06 * exp(-0.5 * m_impact_parameter * m_impact_parameter / 3.2 / 3.2))
         ) * sqrt(v2);
        float gb = (
            1.096 
            + (1.36 * exp(-0.5 * m_impact_parameter * m_impact_parameter / 3.0 / 3.0))
         ) * sqrt(v2);

        v3 = pow( fb , 3);
        v4 = pow( gb, 4);
        v5 = pow( gb, 5);
        v6 = pow( gb, 6);
    }
    else if (m_algorithm == minbias_v2_algorithm)
    { // only v2 is calculated
        v1 = 0;
        v2 = AfterburnerAlgo::calc_v2(m_impact_parameter, eta, pt);
        v3 = 0;
        v4 = 0;
        v5 = 0;
        v6 = 0;
    }

    m_vn[0] = v1 * m_vn_scalefactors[0];
    m_vn[1] = v2 * m_vn_scalefactors[1];
    m_vn[2] = v3 * m_vn_scalefactors[2];
    m_vn[3] = v4 * m_vn_scalefactors[3];
    m_vn[4] = v5 * m_vn_scalefactors[4];
    m_vn[5] = v6 * m_vn_scalefactors[5];

    if ( _do_fluctuations )
    {
        // not implemented yet
    }

    return;
}

std::string AfterburnerAlgo::getAlgoName(flowAfterburnerAlgorithm algo)
{
  switch (algo)
  {
    case minbias_algorithm:
      return "MINBIAS";
    case minbias_v2_algorithm:
      return "MINBIAS_V2_ONLY";
    case custom_algorithm:
      return "CUSTOM";
    default:
      std::cerr << "AfterburnerAlgo::getAlgoName: Unknown algorithm type, returning MINBIAS" << std::endl;
      return "MINBIAS";
  }
}

AfterburnerAlgo::flowAfterburnerAlgorithm AfterburnerAlgo::getAlgoFromName(const std::string &name)
{
  if (name == "MINBIAS")
  {
    return minbias_algorithm;
  }
  else if (name == "MINBIAS_V2_ONLY")
  {
    return minbias_v2_algorithm;
  }
  else if (name == "CUSTOM")
  {
    return custom_algorithm;
  }
  else
  {
    std::cerr << "AfterburnerAlgo::getAlgoFromName: Unknown algorithm name: " << name << std::endl;
    return minbias_algorithm; // default
  }
}