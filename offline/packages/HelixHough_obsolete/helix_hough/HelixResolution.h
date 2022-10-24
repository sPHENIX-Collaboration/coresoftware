#ifndef __HELIXRESOLUTION__
#define __HELIXRESOLUTION__

class HelixResolution
{
  public:
    HelixResolution(float cAngle, float dca, float invR, float theta, float z) : phi_res(cAngle), d_res(dca), k_res(invR), dzdl_res(theta), z0_res(z) {}
    ~HelixResolution(){}
    
    float phi_res, d_res, k_res, dzdl_res, z0_res;
};



#endif
