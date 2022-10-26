#include "FourHitSeedFinder.h"
#include <cmath>
#include <iostream>
#include <algorithm>


using namespace std;


bool FourHitSeedFinder::breakRecursion(const vector<SimpleHit3D>& hits, const HelixRange& range)
{
//   if(0.003*Bfield/range.max_k > 1.0){return false;}
//   else
//   {
//     if( ( (range.max_phi - range.min_phi)/0.025 + (range.max_d - range.min_d)/0.02 ) < 4. ){return true;}
//   }
  return false;
}


float FourHitSeedFinder::phiError(SimpleHit3D& hit, float min_k, float max_k, float min_dzdl, float max_dzdl)
{
  if(3.33333333333333314e+02*max_k*Bfield_inv < 1.0){return 0.;}
  
  float p_inv = 3.33333333333333314e+02*max_k*Bfield_inv*sqrt(1. - max_dzdl*max_dzdl);
  float scatter = 0.;
  for(int l=0;l<hit.layer;++l)
  {
    scatter += detector_scatter[l]*detector_scatter[l];
  }
  scatter = p_inv*sqrt(scatter)*0.5;
  return scatter;
}


float FourHitSeedFinder::dzdlError(SimpleHit3D& hit, float min_k, float max_k, float min_dzdl, float max_dzdl)
{
  if(3.33333333333333314e+02*max_k*Bfield_inv < 1.0){return 0.;}
  
  float p_inv = 3.33333333333333314e+02*max_k*Bfield_inv*sqrt(1. - max_dzdl*max_dzdl);
  float scatter = 0.;
  for(int l=0;l<hit.layer;++l)
  {
    scatter += detector_scatter[l]*detector_scatter[l];
  }
  scatter = p_inv*sqrt(scatter)*0.5;
  return scatter;
}
