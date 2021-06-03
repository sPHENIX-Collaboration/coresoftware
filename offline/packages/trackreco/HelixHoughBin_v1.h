#ifndef G4HOUGH_HELIXHOUGHBINV1_H
#define G4HOUGH_HELIXHOUGHBINV1_H

#include "HelixHoughBin.h"
#include "HelixHoughSpace.h"

#include <stddef.h>  // for size_t
#include <iostream>  // for cout, ostream

class PHObject;

class HelixHoughBin_v1 : public HelixHoughBin
{
 public:
  HelixHoughBin_v1(unsigned int bin);
  ~HelixHoughBin_v1() override {}

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override{};
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new HelixHoughBin_v1(*this); }

  void init() override;

  //     get_cluster_IDs() {}
  void add_cluster_ID(unsigned int cluster_ID) override { _cluster_IDs.insert(cluster_ID); }
  unsigned int get_count() const override { return _cluster_IDs.size(); }
  void clear_clusters() override { _cluster_IDs.clear(); }
  bool empty_clusters() override { return _cluster_IDs.empty(); }
  size_t erase_cluster(unsigned int cluster_id) override { return _cluster_IDs.erase(cluster_id); }

  ConstClusterIter begin_clusters() const override { return _cluster_IDs.begin(); }
  ConstClusterIter find_cluster(unsigned int cluster_id) const override { return _cluster_IDs.find(cluster_id); }
  ConstClusterIter end_clusters() const override { return _cluster_IDs.end(); }
  ClusterIter begin_clusters() override { return _cluster_IDs.begin(); }
  ClusterIter find_cluster(unsigned int cluster_id) override { return _cluster_IDs.find(cluster_id); }
  ClusterIter end_clusters() override { return _cluster_IDs.end(); }

  unsigned int get_bin(unsigned int zoomlevel) const override { return _bin[zoomlevel]; }
  void set_bin(unsigned int zoomlevel, unsigned int bin) override { _bin[zoomlevel] = bin; }

  unsigned int get_zoomlevel() const override { return _zoomlevel; }
  void set_zoomlevel(unsigned int zoomlevel) override { _zoomlevel = zoomlevel; }

  unsigned int get_kappa_bin(unsigned int zoomlevel) const override { return _kappa_bins[zoomlevel]; }
  void set_kappa_bin(unsigned int zoomlevel, unsigned int kappa_bin) override { _kappa_bins[zoomlevel] = kappa_bin; }
  unsigned int get_phi_bin(unsigned int zoomlevel) const override { return _phi_bins[zoomlevel]; }
  void set_phi_bin(unsigned int zoomlevel, unsigned int phi_bin) override { _phi_bins[zoomlevel] = phi_bin; }
  unsigned int get_phi_high_bin(unsigned int zoomlevel) const override { return _phi_high_bins[zoomlevel]; }
  void set_phi_high_bin(unsigned int zoomlevel, unsigned int phi_high_bin) override { _phi_high_bins[zoomlevel] = phi_high_bin; }
  unsigned int get_phi_low_bin(unsigned int zoomlevel) const override { return _phi_low_bins[zoomlevel]; }
  void set_phi_low_bin(unsigned int zoomlevel, unsigned int phi_low_bin) override { _phi_low_bins[zoomlevel] = phi_low_bin; }
  unsigned int get_d_bin(unsigned int zoomlevel) const override { return _d_bins[zoomlevel]; }
  void set_d_bin(unsigned int zoomlevel, unsigned int d_bin) override { _d_bins[zoomlevel] = d_bin; }
  unsigned int get_dzdl_bin(unsigned int zoomlevel) const override { return _dzdl_bins[zoomlevel]; }
  void set_dzdl_bin(unsigned int zoomlevel, unsigned int dzdl_bin) override { _dzdl_bins[zoomlevel] = dzdl_bin; }
  unsigned int get_dzdl_high_bin(unsigned int zoomlevel) const override { return _dzdl_high_bins[zoomlevel]; }
  void set_dzdl_high_bin(unsigned int zoomlevel, unsigned int dzdl_high_bin) override { _dzdl_high_bins[zoomlevel] = dzdl_high_bin; }
  unsigned int get_dzdl_low_bin(unsigned int zoomlevel) const override { return _dzdl_low_bins[zoomlevel]; }
  void set_dzdl_low_bin(unsigned int zoomlevel, unsigned int dzdl_low_bin) override { _dzdl_low_bins[zoomlevel] = dzdl_low_bin; }
  unsigned int get_z0_bin(unsigned int zoomlevel) const override { return _z0_bins[zoomlevel]; }
  void set_z0_bin(unsigned int zoomlevel, unsigned int z0_bin) override { _z0_bins[zoomlevel] = z0_bin; }

  void set_hough_space(HelixHoughSpace* hough_space) override;
  void set_bins(unsigned int zoomlevel, unsigned int bin) override;

  unsigned int get_global_bin(unsigned int zoomlevel) override;
  void set_global_bin(unsigned int zoomlevel) override;
  unsigned int get_neighbors_global_bin(unsigned int zoomlevel, unsigned int var, unsigned int bit_sign) override;

  float get_kappa_center(unsigned int zoomlevel) override;
  float get_phi_center(unsigned int zoomlevel) override;
  float get_d_center(unsigned int zoomlevel) override;
  float get_dzdl_center(unsigned int zoomlevel) override;
  float get_z0_center(unsigned int zoomlevel) override;

 private:
  ClusterSet _cluster_IDs;  // hits who voted for this bin

  unsigned int _global_bin;
  unsigned int _bin[ZOOMLEVEL_MAX];
  unsigned int _kappa_bins[ZOOMLEVEL_MAX];  //={3,2,4} 3rd in most coarse bins, 4th in the narrowest bins
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

  ClassDefOverride(HelixHoughBin_v1, 1);
};

#endif
