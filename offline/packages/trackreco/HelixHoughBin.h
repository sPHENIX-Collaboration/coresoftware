#ifndef G4HOUGH_HELIXHOUGHBIN_H
#define G4HOUGH_HELIXHOUGHBIN_H


#include "HelixHoughSpace.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <set>

class HelixHoughSpace;

// bin has hierarchy, 5 bins for helix parameters

class HelixHoughBin : public PHObject {

public:

  typedef std::set<unsigned int> ClusterSet;
  typedef std::set<unsigned int>::const_iterator ConstClusterIter;
  typedef std::set<unsigned int>::iterator       ClusterIter;


  virtual ~HelixHoughBin() {}

  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "HelixHoughBin base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const {return 0;}
  virtual PHObject* CloneMe() const {return nullptr;}

  virtual void init() {}

//     get_cluster_IDs() {}  
  virtual void add_cluster_ID(unsigned int cluster_ID) {}
  virtual unsigned int get_count() const {return UINT_MAX;}
  virtual void clear_clusters() {}
  virtual bool empty_clusters() {return true;}
  virtual size_t           erase_cluster(unsigned int cluster_id)       {return 0;}
  virtual ConstClusterIter begin_clusters() const                       {return ClusterSet().end();}
  virtual ConstClusterIter find_cluster(unsigned int cluster_id) const  {return ClusterSet().end();}
  virtual ConstClusterIter end_clusters() const                         {return ClusterSet().end();}
  virtual ClusterIter      begin_clusters()                             {return ClusterSet().end();}
  virtual ClusterIter      find_cluster(unsigned int cluster_id)        {return ClusterSet().end();}
  virtual ClusterIter      end_clusters()                               {return ClusterSet().end();}

  virtual unsigned int get_bin(unsigned int zoomlevel) const    {return UINT_MAX;}
  virtual void set_bin(unsigned int zoomlevel, unsigned int bin){}

  virtual unsigned int get_zoomlevel() const {return UINT_MAX;}
  virtual void set_zoomlevel(unsigned int zoomlevel) {}

  virtual unsigned int get_kappa_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_kappa_bin(unsigned int zoomlevel, unsigned int kappa_bin) {}
  virtual unsigned int get_phi_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_phi_bin(unsigned int zoomlevel, unsigned int phi_bin) {}
  virtual unsigned int get_phi_high_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_phi_high_bin(unsigned int zoomlevel) {}
  virtual void set_phi_high_bin(unsigned int zoomlevel, unsigned int phi_high_bin) {}
  virtual unsigned int get_phi_low_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_phi_low_bin(unsigned int zoomlevel) {}
  virtual void set_phi_low_bin(unsigned int zoomlevel, unsigned int phi_low_bin) {}
  virtual unsigned int get_d_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_d_bin(unsigned int zoomlevel, unsigned int d_bin) {}
  virtual unsigned int get_dzdl_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_dzdl_bin(unsigned int zoomlevel, unsigned int dzdl_bin) {}
  virtual unsigned int get_dzdl_high_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_dzdl_high_bin(unsigned int zoomlevel) {}
  virtual void set_dzdl_high_bin(unsigned int zoomlevel, unsigned int dzdl_high_bin) {}
  virtual unsigned int get_dzdl_low_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_dzdl_low_bin(unsigned int zoomlevel) {}
  virtual void set_dzdl_low_bin(unsigned int zoomlevel, unsigned int dzdl_low_bin) {}
  virtual unsigned int  get_z0_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_z0_bin(unsigned int zoomlevel, unsigned int z0_bin) {}

  virtual  void set_hough_space(HelixHoughSpace* hough_space) {};
  virtual  void set_bins(unsigned int zoomlevel, unsigned int bin) {};

  virtual unsigned int get_global_bin(unsigned int zoomlevel) {return UINT_MAX;}
  virtual void set_global_bin(unsigned int zoomlevel) {}

  virtual unsigned int get_neighbors_global_bin(unsigned int zoomlevel, unsigned int var, unsigned int bit_sign) {return UINT_MAX;}

  virtual float get_kappa_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_phi_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_d_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_dzdl_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_z0_center(unsigned int zoomlevel) {return 999.;}

protected:
  HelixHoughBin() {}

  ClassDef(HelixHoughBin,1);
};

#endif

