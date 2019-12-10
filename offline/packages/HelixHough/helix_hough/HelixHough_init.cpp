#include "HelixHough.h"
#include "HelixKalmanState.h"
#include "HelixRange.h"
#include "HelixResolution.h"
#include "SimpleHit3D.h"
#include "SimpleTrack3D.h"

#include <iostream>
#include <utility>
#include <vector>

using namespace std;


HelixHough::HelixHough(unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, HelixResolution& min_resolution, HelixResolution& max_resolution, HelixRange& range) : vote_time(0.), xy_vote_time(0.), z_vote_time(0.), print_timings(false), separate_by_helicity(true), helicity(false), check_layers(false), req_layers(0), bin_scale(1.), z_bin_scale(1.), remove_hits(false), only_one_helicity(false), start_zoom(0), max_hits_pairs(0), cluster_start_bin(2), layers_at_a_time(4), n_layers(6), smooth_back(false), cull_input_hits(false), iterate_clustering(false)
{
  initHelixHough(n_phi, n_d, n_k, n_dzdl, n_z0, min_resolution, max_resolution, range);
  hit_used = new vector<unsigned int>;
}


HelixHough::HelixHough(vector<vector<unsigned int> >& zoom_profile, unsigned int minzoom, HelixRange& range) : vote_time(0.), xy_vote_time(0.), z_vote_time(0.), print_timings(false), separate_by_helicity(true), helicity(false), check_layers(false), req_layers(0), bin_scale(1.), z_bin_scale(1.), remove_hits(false), only_one_helicity(false), start_zoom(0), max_hits_pairs(0), cluster_start_bin(2), layers_at_a_time(4), n_layers(6), layer_start(-1), layer_end(-1), smooth_back(false), cull_input_hits(false), iterate_clustering(false)
{
  for(unsigned int i=0;i<hits_vec.size();i++){delete hits_vec[i];}
  hits_vec.clear();
  for(unsigned int i=0;i<pairs_vec.size();i++)
    {
      (*pairs_vec[i]).clear();
      delete pairs_vec[i];
    }
  pairs_vec.clear();
  for(unsigned int i=0;i<bins_vec.size();i++){delete bins_vec[i];}
  bins_vec.clear();
  for(unsigned int i=0;i<seeds_vec.size();i++){delete seeds_vec[i];}
  seeds_vec.clear();
  for(unsigned int i=0;i<clusters_vec.size();i++){delete clusters_vec[i];}
  clusters_vec.clear();
  num_clusters.clear();
  
  for(unsigned int i=0;i<seeds_vec.size();i++){delete seeds_vec[i];}
  seeds_vec.clear();
  
  
  top_range = range;
  max_zoom = zoom_profile.size()-1;
  min_zoom = minzoom;
  for(unsigned int i=0;i<=(max_zoom);++i)
  {
    unsigned int phibins=zoom_profile[i][0];
    unsigned int dbins=zoom_profile[i][1];
    unsigned int kbins=zoom_profile[i][2];
    unsigned int dzdlbins=zoom_profile[i][3];
    unsigned int z0bins=zoom_profile[i][4];
    
    n_phi_bins.push_back(phibins);
    n_d_bins.push_back(dbins);
    n_k_bins.push_back(kbins);
    n_dzdl_bins.push_back(dzdlbins);
    n_z0_bins.push_back(z0bins);
    
    hits_vec.push_back(new vector<SimpleHit3D>);
    pairs_vec.push_back(new vector<pair<unsigned int,unsigned int> >);
    bins_vec.push_back(new vector<BinEntryPair5D>);
    seeds_vec.push_back(new vector<SimpleTrack3D>);
    clusters_vec.push_back(new vector<ParameterCluster>);
    num_clusters.push_back(0);
  }
  
  hit_used = new vector<unsigned int>;
}


void HelixHough::initHelixHough(unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, HelixResolution& min_resolution, HelixResolution& max_resolution, HelixRange& range)
{
  for(unsigned int i=0;i<hits_vec.size();i++){delete hits_vec[i];}
  hits_vec.clear();
  for(unsigned int i=0;i<pairs_vec.size();i++){delete pairs_vec[i];}
  pairs_vec.clear();
  for(unsigned int i=0;i<bins_vec.size();i++){delete bins_vec[i];}
  bins_vec.clear();
  for(unsigned int i=0;i<clusters_vec.size();i++){delete clusters_vec[i];}
  clusters_vec.clear();
  num_clusters.clear();
  
  for(unsigned int i=0;i<seeds_vec.size();i++){delete seeds_vec[i];}
  seeds_vec.clear();
  
  top_range = range;
  
  unsigned int zoom_level = 0;
  float k_res = (range.max_k - range.min_k);
  float phi_res = (range.max_phi - range.min_phi);
  float d_res = (range.max_d - range.min_d);
  float z0_res = (range.max_z0 - range.min_z0);
  float dzdl_res = (range.max_dzdl - range.min_dzdl);
  min_zoom=0;
  bool min_zoom_set=false;
  
  unsigned int nzoom_k=0;
  unsigned int nzoom_d=0;
  unsigned int nzoom_phi=0;
  unsigned int nzoom_dzdl=0;
  unsigned int nzoom_z0=0;
  // find the number of iterations it takes to get to the maximum zooming for each dimension
  while(true)
  {
    if(zoom_level != 0)
    {
      if(k_res<=max_resolution.k_res && phi_res<=max_resolution.phi_res && d_res<=max_resolution.d_res && z0_res<=max_resolution.z0_res && dzdl_res<=max_resolution.dzdl_res)
      {
        max_zoom = zoom_level;
        break;
      }
    }
    
    unsigned int kbins=1;
    unsigned int phibins=1;
    unsigned int dbins=1;
    unsigned int z0bins=1;
    unsigned int dzdlbins=1;
    if(k_res>max_resolution.k_res)
    {
      nzoom_k+=1;
      kbins = (unsigned int)(ceil(k_res/max_resolution.k_res));
      if(kbins > n_k){kbins = n_k;}
      k_res/=((float)kbins);
    }
    if(phi_res>max_resolution.phi_res)
    {
      nzoom_phi+=1;
      phibins = (unsigned int)(ceil(phi_res/max_resolution.phi_res));
      if(phibins > n_phi){phibins = n_phi;}
      phi_res/=((float)phibins);
    }
    if(d_res>max_resolution.d_res)
    {
      nzoom_d+=1;
      dbins = (unsigned int)(ceil(d_res/max_resolution.d_res));
      if(dbins > n_d){dbins = n_d;}
      d_res/=((float)dbins);
    }
    if(z0_res>max_resolution.z0_res)
    {
      nzoom_z0+=1;
      z0bins = (unsigned int)(ceil(z0_res/max_resolution.z0_res));
      if(z0bins > n_z0){z0bins = n_z0;}
      z0_res/=((float)z0bins);
    }
    if(dzdl_res>max_resolution.dzdl_res)
    {
      nzoom_dzdl+=1;
      dzdlbins = (unsigned int)(ceil(dzdl_res/max_resolution.dzdl_res));
      if(dzdlbins > n_dzdl){dzdlbins = n_dzdl;}
      dzdl_res/=((float)dzdlbins);
    }
    zoom_level+=1;
  }
  
  k_res = (range.max_k - range.min_k);
  phi_res = (range.max_phi - range.min_phi);
  d_res = (range.max_d - range.min_d);
  z0_res = (range.max_z0 - range.min_z0);
  dzdl_res = (range.max_dzdl - range.min_dzdl);
  for(unsigned int zoom=0;zoom<max_zoom;++zoom)
  {
    if(k_res<=min_resolution.k_res && phi_res<=min_resolution.phi_res && d_res<=min_resolution.d_res && z0_res<=min_resolution.z0_res && dzdl_res<=min_resolution.dzdl_res && min_zoom_set==false)
    {
      min_zoom = zoom;
      min_zoom_set=true;
    }
    //if we need to zoom more, we determine the binning of the next zoom level.
    unsigned int kbins=1;
    unsigned int phibins=1;
    unsigned int dbins=1;
    unsigned int z0bins=1;
    unsigned int dzdlbins=1;
    if((k_res>max_resolution.k_res) && ( (max_zoom - zoom) <= nzoom_k ))
    {
      kbins = (unsigned int)(ceil(k_res/max_resolution.k_res));
      if( kbins > n_k ){kbins = n_k;}
      k_res/=((float)kbins);
    }
    if((phi_res>max_resolution.phi_res) && ( (max_zoom - zoom) <= nzoom_phi ))
    {
      phibins = (unsigned int)(ceil(phi_res/max_resolution.phi_res));
      if( phibins > n_phi ){phibins = n_phi;}
      phi_res/=((float)phibins);
    }
    if((d_res>max_resolution.d_res) && ( (max_zoom - zoom) <= nzoom_d ))
    {
      dbins = (unsigned int)(ceil(d_res/max_resolution.d_res));
      if( dbins > n_d ){dbins = n_d;}
      d_res/=((float)dbins);
    }
    if((z0_res>max_resolution.z0_res) && ( (max_zoom - zoom) <= nzoom_z0 ))
    {
      z0bins = (unsigned int)(ceil(z0_res/max_resolution.z0_res));
      if( z0bins > n_z0 ){z0bins = n_z0;}
      z0_res/=((float)z0bins);
    }
    if((dzdl_res>max_resolution.dzdl_res) && ( (max_zoom - zoom) <= nzoom_dzdl ))
    {
      dzdlbins = (unsigned int)(ceil(dzdl_res/max_resolution.dzdl_res));
      if( dzdlbins > n_dzdl ){dzdlbins = n_dzdl;}
      dzdl_res/=((float)dzdlbins);
    }
    n_phi_bins.push_back(phibins);
    n_d_bins.push_back(dbins);
    n_k_bins.push_back(kbins);
    n_dzdl_bins.push_back(dzdlbins);
    n_z0_bins.push_back(z0bins);
    
    hits_vec.push_back(new vector<SimpleHit3D>);
    pairs_vec.push_back(new vector<pair<unsigned int,unsigned int> >);
    bins_vec.push_back(new vector<BinEntryPair5D>);
    seeds_vec.push_back(new vector<SimpleTrack3D>);
    clusters_vec.push_back(new vector<ParameterCluster>);
    num_clusters.push_back(0);
  }
  max_zoom -= 1;
  if(min_zoom_set==false){min_zoom=max_zoom;}
  cout<<"min_zoom = "<<min_zoom<<endl;
}


HelixHough::~HelixHough()
{
  for(unsigned int i=0;i<hits_vec.size();i++){delete hits_vec[i];}
  for(unsigned int i=0;i<pairs_vec.size();i++){delete pairs_vec[i];}
  for(unsigned int i=0;i<bins_vec.size();i++){delete bins_vec[i];}
  for(unsigned int i=0;i<seeds_vec.size();i++){delete seeds_vec[i];}
  for(unsigned int i=0;i<clusters_vec.size();i++){delete clusters_vec[i];}
  delete hit_used;
}
