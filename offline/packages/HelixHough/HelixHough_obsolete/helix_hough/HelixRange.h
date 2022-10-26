#ifndef __HELIXRANGE__
#define __HELIXRANGE__

class HelixRange
{
public:
  HelixRange(){}
  HelixRange(float min_cAngle, float max_cAngle, float min_dca, float max_dca, float min_invR, float max_invR, float min_theta, float max_theta, float min_z, float max_z) : min_k(min_invR), max_k(max_invR), min_phi(min_cAngle), max_phi(max_cAngle), min_d(min_dca), max_d(max_dca), min_dzdl(min_theta), max_dzdl(max_theta), min_z0(min_z), max_z0(max_z) {}
  ~HelixRange() {}
  
  float min_k, max_k;
  float min_phi, max_phi;
  float min_d, max_d;
  float min_dzdl, max_dzdl;
  float min_z0, max_z0;
};

#endif
