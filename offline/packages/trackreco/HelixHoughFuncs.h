#ifndef G4HOUGH_HELIXHOUGHFUNCS_H
#define G4HOUGH_HELIXHOUGHFUNCS_H

#include "HelixHoughSpace_v1.h"

#include <phool/PHObject.h>

#include <climits>
#include <cmath>

class HelixHoughSpace;

class HelixHoughFuncs : public PHObject {

public :
  virtual ~HelixHoughFuncs() {}


  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "HelixHoughFuncs base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int  isValid() const 			{return 0;}
  virtual PHObject* CloneMe() const 	{return nullptr;}

  // Define Hough space for helical tracks 
  virtual void set_current_zoom(unsigned int cur_zoom)		{}
  virtual void set_hough_space(HelixHoughSpace* hough_space)    {}
  virtual void calculate_dzdl_range(float* hitpos3d, std::vector<float>& z0_range, std::vector<float>& kappa_phi_d_ranges, float* dzdl_range) {};
  virtual void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, float* phi_r_range, float* phi_l_range) {};
  virtual void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_next_range){};
  virtual void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_prev_range, float* phi_next_range){};

/*
  virtual void set_z0_max(float z0_max)         {}
  virtual unsigned int get_z0_max() const       {return UINT_MAX;}

  virtual unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins){return -9999;}
*/
protected:
  HelixHoughFuncs(){};
  ClassDef(HelixHoughFuncs,1);
};

#endif
