#ifndef G4HOUGH_HELIXHOUGHBINV1_H
#define G4HOUGH_HELIXHOUGHBINV1_H

#include "HelixHoughBin.h"
#include "HelixHoughSpace.h"

#include <stddef.h>           // for size_t
#include <iostream>           // for cout, ostream

class PHObject;

class HelixHoughBin_v1 : public HelixHoughBin {

public:

  HelixHoughBin_v1(unsigned int bin);
  virtual ~HelixHoughBin_v1() {}

  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const;
  void Reset() {};
  int  isValid() const {return 1;}
  PHObject* CloneMe() const {return new HelixHoughBin_v1(*this);}

  void init();

//     get_cluster_IDs() {}  
  void 		   add_cluster_ID(unsigned int cluster_ID) 	{_cluster_IDs.insert(cluster_ID);} 
  unsigned int 	   get_count() const				{return _cluster_IDs.size();}
  void		   clear_clusters()				{_cluster_IDs.clear();}
  bool 		   empty_clusters()				{return _cluster_IDs.empty();}
  size_t       	   erase_cluster(unsigned int cluster_id)       {return _cluster_IDs.erase(cluster_id);}

  ConstClusterIter begin_clusters() const                    	{return _cluster_IDs.begin();}
  ConstClusterIter find_cluster(unsigned int cluster_id) const  {return _cluster_IDs.find(cluster_id);}
  ConstClusterIter end_clusters() const                      	{return _cluster_IDs.end();}
  ClusterIter      begin_clusters()                          	{return _cluster_IDs.begin();}
  ClusterIter      find_cluster(unsigned int cluster_id)        {return _cluster_IDs.find(cluster_id);}
  ClusterIter      end_clusters()                            	{return _cluster_IDs.end();}

  unsigned int get_bin(unsigned int zoomlevel) const 	{return _bin[zoomlevel];}
  void set_bin(unsigned int zoomlevel, unsigned int bin){_bin[zoomlevel]=bin;}

  unsigned int get_zoomlevel() const 	{return _zoomlevel;}
  void set_zoomlevel(unsigned int zoomlevel) {_zoomlevel = zoomlevel;}

  unsigned int get_kappa_bin(unsigned int zoomlevel) const {return _kappa_bins[zoomlevel];}
  void set_kappa_bin(unsigned int zoomlevel, unsigned int kappa_bin) {_kappa_bins[zoomlevel] = kappa_bin;}
  unsigned int get_phi_bin(unsigned int zoomlevel) const {return _phi_bins[zoomlevel];}
  void set_phi_bin(unsigned int zoomlevel, unsigned int phi_bin) {_phi_bins[zoomlevel] = phi_bin;}
  unsigned int get_phi_high_bin(unsigned int zoomlevel) const {return _phi_high_bins[zoomlevel];}
  void set_phi_high_bin(unsigned int zoomlevel, unsigned int phi_high_bin) {_phi_high_bins[zoomlevel] = phi_high_bin;}
  unsigned int get_phi_low_bin(unsigned int zoomlevel) const {return _phi_low_bins[zoomlevel];}
  void set_phi_low_bin(unsigned int zoomlevel, unsigned int phi_low_bin) {_phi_low_bins[zoomlevel] = phi_low_bin;}
  unsigned int get_d_bin(unsigned int zoomlevel) const {return _d_bins[zoomlevel];}
  void set_d_bin(unsigned int zoomlevel, unsigned int d_bin) {_d_bins[zoomlevel] = d_bin;}
  unsigned int get_dzdl_bin(unsigned int zoomlevel) const {return _dzdl_bins[zoomlevel];}
  void set_dzdl_bin(unsigned int zoomlevel, unsigned int dzdl_bin) {_dzdl_bins[zoomlevel] = dzdl_bin;}
  unsigned int get_dzdl_high_bin(unsigned int zoomlevel) const {return _dzdl_high_bins[zoomlevel];}
  void set_dzdl_high_bin(unsigned int zoomlevel, unsigned int dzdl_high_bin) {_dzdl_high_bins[zoomlevel] = dzdl_high_bin;}
  unsigned int get_dzdl_low_bin(unsigned int zoomlevel) const {return _dzdl_low_bins[zoomlevel];}
  void set_dzdl_low_bin(unsigned int zoomlevel, unsigned int dzdl_low_bin) {_dzdl_low_bins[zoomlevel] = dzdl_low_bin;}
  unsigned int get_z0_bin(unsigned int zoomlevel) const {return _z0_bins[zoomlevel];}
  void set_z0_bin(unsigned int zoomlevel, unsigned int z0_bin) {_z0_bins[zoomlevel] = z0_bin;}


  void set_hough_space(HelixHoughSpace* hough_space);   
  void set_bins(unsigned int zoomlevel, unsigned int bin);

  unsigned int get_global_bin(unsigned int zoomlevel);
  void set_global_bin(unsigned int zoomlevel);
  unsigned int get_neighbors_global_bin(unsigned int zoomlevel, unsigned int var, unsigned int bit_sign);

  float get_kappa_center(unsigned int zoomlevel);
  float get_phi_center(unsigned int zoomlevel);
  float get_d_center(unsigned int zoomlevel);
  float get_dzdl_center(unsigned int zoomlevel);
  float get_z0_center(unsigned int zoomlevel);

private:

  ClusterSet _cluster_IDs;// hits who voted for this bin

  unsigned int _global_bin;
  unsigned int _bin[ZOOMLEVEL_MAX];
  unsigned int _kappa_bins[ZOOMLEVEL_MAX];//={3,2,4} 3rd in most coarse bins, 4th in the narrowest bins 
  unsigned int _phi_bins[ZOOMLEVEL_MAX];
  unsigned int _phi_high_bins[ZOOMLEVEL_MAX];
  unsigned int _phi_low_bins[ZOOMLEVEL_MAX];
  unsigned int _d_bins[ZOOMLEVEL_MAX];
  unsigned int _dzdl_bins[ZOOMLEVEL_MAX];
  unsigned int _dzdl_high_bins[ZOOMLEVEL_MAX];
  unsigned int _dzdl_low_bins[ZOOMLEVEL_MAX];
  unsigned int _z0_bins[ZOOMLEVEL_MAX];
  unsigned int _zoomlevel;


  HelixHoughSpace* _hough_space;

  ClassDef(HelixHoughBin_v1,1);
};

#endif

