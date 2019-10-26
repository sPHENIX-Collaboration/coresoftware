#ifndef G4HOUGH_HELIXHOUGHSPACE_H
#define G4HOUGH_HELIXHOUGHSPACE_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>

#ifndef ZOOMLEVEL_MAX
#define ZOOMLEVEL_MAX           5
#endif

class HelixHoughSpace : public PHObject {

public :
  virtual ~HelixHoughSpace() {}


  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "HelixHough base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int  isValid() const 			{return 0;}
  virtual PHObject* CloneMe() const 	{return nullptr;}

  // Define Hough space for helical tracks 
  virtual void add_one_zoom(std::vector<unsigned int>& one_zoom) {};
  virtual unsigned int get_max_zoom() 		{return UINT_MAX;}
  virtual void print_zoom_profile()		{}
  virtual void print_para_range()		{}

  virtual void set_kappa_min(float kappa_min)   {}
  virtual float get_kappa_min() const    	{return NAN;}
  virtual void set_kappa_max(float kappa_max)   {}
  virtual float get_kappa_max() const    	{return NAN;}
  virtual void set_phi_min(float phi_min)       {}
  virtual float get_phi_min() const      	{return NAN;} 
  virtual void set_phi_max(float phi_max)       {}
  virtual float get_phi_max() const      	{return NAN;}
  virtual void set_d_min(float d_min)           {}
  virtual float get_d_min() const        	{return NAN;}
  virtual void set_d_max(float d_max)           {}
  virtual float get_d_max() const        	{return NAN;}
  virtual void set_dzdl_min(float dzdl_min)     {}
  virtual float get_dzdl_min() const     	{return NAN;}
  virtual void set_dzdl_max(float dzdl_max)     {}
  virtual float get_dzdl_max() const     	{return NAN;}
  virtual void set_z0_min(float z0_min)         {}
  virtual float get_z0_min() const       	{return NAN;}
  virtual void set_z0_max(float z0_max)         {}
  virtual float get_z0_max() const       	{return NAN;}


  virtual unsigned int get_n_kappa_bins(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual unsigned int get_n_phi_bins(unsigned int zoomlevel) const   {return UINT_MAX;}
  virtual unsigned int get_n_d_bins(unsigned int zoomlevel) const     {return UINT_MAX;}
  virtual unsigned int get_n_dzdl_bins(unsigned int zoomlevel) const  {return UINT_MAX;}
  virtual unsigned int get_n_z0_bins(unsigned int zoomlevel) const    {return UINT_MAX;}

  virtual float get_kappa_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_phi_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_d_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_dzdl_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_z0_bin_size(unsigned int zoomlevel) const      {return NAN;}

/*
  virtual float get_kappa_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ik) const {return NAN;}
  virtual float get_phi_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ip) const   {return NAN;}
  virtual float get_d_center(unsigned int zoomlevel, std::vector<unsigned int>& v_id) const     {return NAN;}
  virtual float get_dzdl_center(unsigned int zoomlevel, std::vector<unsigned int>& v_il) const  {return NAN;}
  virtual float get_z0_center(unsigned int zoomlevel, std::vector<unsigned int>& v_iz) const    {return NAN;}
*/
  virtual unsigned int get_kappa_bin(unsigned int zoomlevel, float kappa) const {return UINT_MAX;} 
  virtual unsigned int get_phi_bin(unsigned int zoomlevel, float phi) const   {return UINT_MAX;}
  virtual unsigned int get_d_bin(unsigned int zoomlevel, float d) const     {return UINT_MAX;} 
  virtual unsigned int get_dzdl_bin(unsigned int zoomlevel, float dzdl) const  {return UINT_MAX;}
  virtual unsigned int get_z0_bin(unsigned int zoomlevel, float z0) const    {return UINT_MAX;}


  virtual unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins) const {return UINT_MAX;}

protected:
  HelixHoughSpace(){};
  ClassDef(HelixHoughSpace,1);
};

#endif
