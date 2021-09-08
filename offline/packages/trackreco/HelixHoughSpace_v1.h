#ifndef G4HOUGH_HELIXHOUGHSPACEV1_H
#define G4HOUGH_HELIXHOUGHSPACEV1_H

#include "HelixHoughSpace.h"

#include <iostream>           // for cout, ostream
#include <vector>             // for vector

class PHObject;

class HelixHoughSpace_v1 : public HelixHoughSpace {

public:

  HelixHoughSpace_v1();
  HelixHoughSpace_v1(const HelixHoughSpace_v1& hough_space);
  ~HelixHoughSpace_v1() override {};


  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const override {};
  void Reset() override {}
  int  isValid() const override {return 1;}
  PHObject* CloneMe() const override {return new HelixHoughSpace_v1(*this);}

  void add_one_zoom(std::vector<unsigned int>& one_zoom) override;
  unsigned int get_max_zoom() override;
  void print_zoom_profile() override;
  void print_para_range() override;

  void set_kappa_min(float kappa_min)  override{_para_min[0] = kappa_min;}
  float get_kappa_min() const override		{return _para_min[0];}
  void set_kappa_max(float kappa_max) override 	{_para_max[0] = kappa_max;}
  float get_kappa_max() const  override		{return _para_max[0];}
  void set_phi_min(float phi_min)  override	{_para_min[1] = phi_min;}
  float get_phi_min() const override		{return _para_min[1];}
  void set_phi_max(float phi_max)  override	{_para_max[1] = phi_max;}
  float get_phi_max() const  override		{return _para_max[1];}
  void set_d_min(float d_min)  override		{_para_min[2] = d_min;}
  float get_d_min() const  override		{return _para_min[2];}
  void set_d_max(float d_max) override		{_para_max[2] = d_max;}
  float get_d_max() const  override		{return _para_max[2];}
  void set_dzdl_min(float dzdl_min) override	{_para_min[3] = dzdl_min;}
  float get_dzdl_min() const  override		{return _para_min[3];}
  void set_dzdl_max(float dzdl_max) override 	{_para_max[3] = dzdl_max;}
  float get_dzdl_max() const  override		{return _para_max[3];}
  void set_z0_min(float z0_min)   override       {_para_min[4] = z0_min;}
  float get_z0_min() const  override		{return _para_min[4];}
  void set_z0_max(float z0_max)  override        {_para_max[4] = z0_max;}
  float get_z0_max() const  override		{return _para_max[4];}

  unsigned int get_n_kappa_bins(unsigned int zoomlevel) const override {return _zoom_profile[zoomlevel][0];}
  unsigned int get_n_phi_bins(unsigned int zoomlevel) const override {return _zoom_profile[zoomlevel][1];}
  unsigned int get_n_d_bins(unsigned int zoomlevel) const override {return _zoom_profile[zoomlevel][2];}
  unsigned int get_n_dzdl_bins(unsigned int zoomlevel) const override {return _zoom_profile[zoomlevel][3];}
  unsigned int get_n_z0_bins(unsigned int zoomlevel) const override {return _zoom_profile[zoomlevel][4];}

  float get_kappa_bin_size(unsigned int zoomlevel) const override;
  float get_phi_bin_size(unsigned int zoomlevel) const override;
  float get_d_bin_size(unsigned int zoomlevel) const override;
  float get_dzdl_bin_size(unsigned int zoomlevel) const override;
  float get_z0_bin_size(unsigned int zoomlevel) const override;
/*
  float get_kappa_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ik) const;
  float get_phi_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ip) const;
  float get_d_center(unsigned int zoomlevel, std::vector<unsigned int>& v_id) const;
  float get_dzdl_center(unsigned int zoomlevel, std::vector<unsigned int>& v_il) const;
  float get_z0_center(unsigned int zoomlevel, std::vector<unsigned int>& v_iz) const;
*/
  unsigned int get_kappa_bin(unsigned int zoomlevel, float kappa) const override;
  unsigned int get_phi_bin(unsigned int zoomlevel, float phi) const override;
  unsigned int get_d_bin(unsigned int zoomlevel, float d) const override;
  unsigned int get_dzdl_bin(unsigned int zoomlevel, float dzdl) const override;
  unsigned int get_z0_bin(unsigned int zoomlevel, float z0) const override;

  unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins) const override;

private:

  float _para_min[5];
  float _para_max[5];

  unsigned int _zoom_profile[ZOOMLEVEL_MAX][5];
  unsigned int _max_zoom;

  ClassDefOverride(HelixHoughSpace_v1,1)
};

#endif
