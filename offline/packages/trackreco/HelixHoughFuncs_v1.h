#ifndef G4HOUGH_HELIXHOUGHFUNCSV1_H
#define G4HOUGH_HELIXHOUGHFUNCSV1_H

#include "HelixHoughFuncs.h"

#include <iostream>  // for cout, ostream
#include <vector>    // for vector

class HelixHoughSpace;
class PHObject;

class HelixHoughFuncs_v1 : public HelixHoughFuncs
{
 public:
  HelixHoughFuncs_v1();
  HelixHoughFuncs_v1(const HelixHoughFuncs_v1& hough_funcs);
  ~HelixHoughFuncs_v1() override{};

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override{};
  void Reset() override {}
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new HelixHoughFuncs_v1(*this); }

  void set_current_zoom(unsigned int cur_zoom) override { _cur_zoom = cur_zoom; }
  void set_hough_space(HelixHoughSpace* hough_space) override;

  void calculate_dzdl_range(float* hitpos3d, std::vector<float>& z0_range, std::vector<float>& kappa_phi_d_ranges, float* dzdl_range) override;
  void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, float* phi_r_range, float* phi_l_range) override;
  void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_next_range) override;
  void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_prev_range, float* phi_next_range) override;

  /*
  unsigned int get_kappa_min() const	{return _para_min[0];}
  unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins);
*/
 private:
  /*
  float _para_min[5];
  float _para_max[5];

  unsigned int _zoom_profile[ZOOMLEVEL_MAX][5];
  unsigned int _max_zoom;
*/
  HelixHoughSpace* _hough_space;
  unsigned int _cur_zoom;

  ClassDefOverride(HelixHoughFuncs_v1, 1)
};

#endif
