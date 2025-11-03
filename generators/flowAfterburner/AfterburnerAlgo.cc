#include "AfterburnerAlgo.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>  // for exit
#include <array>

#include <CLHEP/Random/RandomEngine.h>

AfterburnerAlgo::AfterburnerAlgo( flowAfterburnerAlgorithm algorithm)
    : m_algorithm(algorithm)
{
  for (unsigned int i = 0; i < 6; ++i)
  {
    m_vn[i] = 0.0;
    m_vn_scalefactors[i] = 1.0;
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
  for (float & m_vn_scalefactor : m_vn_scalefactors)
  {
    m_vn_scalefactor = scale;
  }
}

float AfterburnerAlgo::get_vn(unsigned int n) const
{
  if (n < 1 || n > 6)
  {
    std::cerr << "AfterburnerAlgo::get_vn: n must be between 1 and 6, got " << n << std::endl;
    return 0.0;
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

void AfterburnerAlgo::calc_flow(double eta, double pt, CLHEP::HepRandomEngine* engine)
{
    
    float v1 = 0;
    float v2 = 0;
    float v3 = 0;
    float v4 = 0;
    float v5 = 0;
    float v6 = 0;
    if ( m_algorithm == custom_algorithm )
    { // Custom flow parameters
        v1 = 0.0000;
        v2 = 0.0500;
        v3 = 0.0280;
        v4 = 0.0130;
        v5 = 0.0045;
        v6 = 0.0015;
    }
    else if (m_algorithm == minbias_algorithm)
    { // all other algorithms need to calculate the flow
        v1 = 0;
        v2 = AfterburnerAlgo::calc_v2(m_impact_parameter, eta, pt); 
        float v2_sqrt = std::sqrt(v2);

        float fb = (
            0.97 
            + (1.06 * exp(-0.5 * m_impact_parameter * m_impact_parameter / 3.2 / 3.2))
         ) * v2_sqrt;
        float gb = (
            1.096 
            + (1.36 * exp(-0.5 * m_impact_parameter * m_impact_parameter / 3.0 / 3.0))
         ) * v2_sqrt;

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

    if ( _do_fluctuations )
    {
      flucatate(engine, v1, v2, v3, v4, v5, v6);
    }

    m_vn[0] = v1 * m_vn_scalefactors[0];
    m_vn[1] = v2 * m_vn_scalefactors[1];
    m_vn[2] = v3 * m_vn_scalefactors[2];
    m_vn[3] = v4 * m_vn_scalefactors[3];
    m_vn[4] = v5 * m_vn_scalefactors[4];
    m_vn[5] = v6 * m_vn_scalefactors[5];

    
    

    return;
}

void AfterburnerAlgo::flucatate(CLHEP::HepRandomEngine* engine, float &v1, float &v2, float &v3, float &v4, float &v5, float &v6) const
{
  // calc polynomial
  // parameters are from fit of sigma to centrality, sampled over HIJING b
  const std::array<float, 8> coeffs = {
      -7.00411e-09, 4.24567e-07, 
      -9.87748e-06, 0.000112689,
      -0.000694686, 0.002413930, 
      -0.00324709, 0.0107906};
  float sigma = 0.0;
  for (float coeff : coeffs) 
  {
      sigma = sigma * m_impact_parameter + coeff;
  }
  sigma = std::max<double>(sigma, 0.0);
  const float fb = (
          0.97 
          + (1.06 * exp(-0.5 * m_impact_parameter * m_impact_parameter / 3.2 / 3.2))
        ) * std::sqrt(v2);
  const float gb = (
      1.096 
      + (1.36 * exp(-0.5 * m_impact_parameter * m_impact_parameter / 3.0 / 3.0))
    ) * std::sqrt(v2);

  const float s1 = 0; // Not implemented, v1 is always 0
  const float s2 = sigma;
  const float s3 = 1.5*std::pow(fb,3) * std::sqrt(v2) * sigma;
  const float s4 = 2.0*std::pow(gb,4) * v2 * sigma;
  const float s5 = 2.5*std::pow(gb,5) * std::sqrt(v2*v2*v2) * sigma;
  const float s6 = 3.0*std::pow(gb,6) * v2*v2 * sigma;

  // Marsaglia polar: two N(0,1) 
  float u;
  float v;
  float s;
  do 
  {
      u = 2.0*engine->flat() - 1.0;
      v = 2.0*engine->flat() - 1.0;
      s = u*u + v*v;
  } while (s >= 1.0 || s == 0.0);

  const float mul = std::sqrt(-2.0*std::log(s)/s);
  const float z0  = u*mul;
  const float z1  = v*mul;


  if ( v1 != 0.0 )
  {
    float _v1 = std::hypot(v1 + s1*z0,  s1*z1);
    v1 = std::clamp(_v1, 0.0F, 1.0F);
  }
  if ( v2 != 0.0 )
  {
    float _v2 = std::hypot(v2 + s2*z0,  s2*z1);
    v2 = std::clamp(_v2, 0.0F, 1.0F);
  }
  if ( v3 != 0.0 )
  {
    float _v3 = std::hypot(v3 + s3*z0,  s3*z1);
    v3 = std::clamp(_v3, 0.0F, 1.0F);
  }
  if ( v4 != 0.0 )
  {
    float _v4 = std::hypot(v4 + s4*z0,  s4*z1);
    v4 = std::clamp(_v4, 0.0F, 1.0F);
  }
  if ( v5 != 0.0 )
  {
    float _v5 = std::hypot(v5 + s5*z0,  s5*z1);
    v5 = std::clamp(_v5, 0.0F, 1.0F);
  }
  if ( v6 != 0.0 )
  {
    float _v6 = std::hypot(v6 + s6*z0,  s6*z1);
    v6 = std::clamp(_v6, 0.0F, 1.0F);
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
  if (name == "MINBIAS_V2_ONLY")
  {
    return minbias_v2_algorithm;
  }
  if (name == "CUSTOM")
  {
    return custom_algorithm;
  }
 
  std::cerr << "AfterburnerAlgo::getAlgoFromName: Unknown algorithm name: " << name << std::endl;
  return minbias_algorithm; // default
}

void AfterburnerAlgo::print(std::ostream &os) const
{
  os << "AfterburnerAlgo parameters:" << std::endl;
  os << "Algorithm: " << getAlgoName(m_algorithm) << std::endl;
  os << "Scale factors: ";
  for (const auto &scale : m_vn_scalefactors)
  {
    os << scale << " ";
  }
  os << std::endl;
  os << "Flow fluctuations enabled: " << (_do_fluctuations ? "Yes" : "No") << std::endl;

  return;
}