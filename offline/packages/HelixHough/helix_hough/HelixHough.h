#ifndef __HELIX_HOUGH_H__
#define __HELIX_HOUGH_H__

#include "fastvec.h"
#include "HelixRange.h"
#include "HelixResolution.h"
#include "SimpleHit3D.h"
#include "SimpleTrack3D.h"
#include "HelixKalmanState.h"
#include <xmmintrin.h>
#include <emmintrin.h>
#include <cmath>
#include <vector>
#include <set>
#include <map>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif





class ParRange
{
public:
  ParRange(){}
  ParRange(unsigned int min_cAngle, unsigned int max_cAngle, unsigned int min_dca, unsigned int max_dca, unsigned int min_invR, unsigned int max_invR, unsigned int min_theta, unsigned int max_theta, unsigned int min_z, unsigned int max_z) : min_k(min_invR), max_k(max_invR), min_phi(min_cAngle), max_phi(max_cAngle), min_d(min_dca), max_d(max_dca), min_dzdl(min_theta), max_dzdl(max_theta), min_z0(min_z), max_z0(max_z) {}
  ~ParRange() {}
  
  void mergeRange(unsigned int phi, unsigned int d, unsigned int k, unsigned int dzdl, unsigned int z0)
  {
    min_phi = (((phi < min_phi)-1)&min_phi) ^ (((phi >= min_phi)-1)&phi);
    max_phi = (((phi > max_phi)-1)&max_phi) ^ (((phi <= max_phi)-1)&phi);
    min_d = (((d < min_d)-1)&min_d) ^ (((d >= min_d)-1)&d);
    max_d = (((d > max_d)-1)&max_d) ^ (((d <= max_d)-1)&d);
    min_k = (((k < min_k)-1)&min_k) ^ (((k >= min_k)-1)&k);
    max_k = (((k > max_k)-1)&max_k) ^ (((k <= max_k)-1)&k);
    min_dzdl = (((dzdl < min_dzdl)-1)&min_dzdl) ^ (((dzdl >= min_dzdl)-1)&dzdl);
    max_dzdl = (((dzdl > max_dzdl)-1)&max_dzdl) ^ (((dzdl <= max_dzdl)-1)&dzdl);
    min_z0 = (((z0 < min_z0)-1)&min_z0) ^ (((z0 >= min_z0)-1)&z0);
    max_z0 = (((z0 > max_z0)-1)&max_z0) ^ (((z0 <= max_z0)-1)&z0);
  }
  
  void mergeRange(ParRange const& other)
  {
    min_phi = (((other.min_phi < min_phi)-1)&min_phi) ^ (((other.min_phi >= min_phi)-1)&other.min_phi);
    max_phi = (((other.max_phi > max_phi)-1)&max_phi) ^ (((other.max_phi <= max_phi)-1)&other.max_phi);
    min_d = (((other.min_d < min_d)-1)&min_d) ^ (((other.min_d >= min_d)-1)&other.min_d);
    max_d = (((other.max_d > max_d)-1)&max_d) ^ (((other.max_d <= max_d)-1)&other.max_d);
    min_k = (((other.min_k < min_k)-1)&min_k) ^ (((other.min_k >= min_k)-1)&other.min_k);
    max_k = (((other.max_k > max_k)-1)&max_k) ^ (((other.max_k <= max_k)-1)&other.max_k);
    min_dzdl = (((other.min_dzdl < min_dzdl)-1)&min_dzdl) ^ (((other.min_dzdl >= min_dzdl)-1)&other.min_dzdl);
    max_dzdl = (((other.max_dzdl > max_dzdl)-1)&max_dzdl) ^ (((other.max_dzdl <= max_dzdl)-1)&other.max_dzdl);
    min_z0 = (((other.min_z0 < min_z0)-1)&min_z0) ^ (((other.min_z0 >= min_z0)-1)&other.min_z0);
    max_z0 = (((other.max_z0 > max_z0)-1)&max_z0) ^ (((other.max_z0 <= max_z0)-1)&other.max_z0);
  }
  
  unsigned int min_k, max_k;
  unsigned int min_phi, max_phi;
  unsigned int min_d, max_d;
  unsigned int min_dzdl, max_dzdl;
  unsigned int min_z0, max_z0;
};


class ParameterCluster
{
  public:
    ParameterCluster(){}
    ~ParameterCluster(){}
    
    ParRange range;
    std::vector<unsigned int> hit_indexes;
};

class BinEntryPair5D
{
public:
  BinEntryPair5D(unsigned int b=0, unsigned int e=0, bool is=false) : bin(b), entry(e), is_seed(is) {}
  ~BinEntryPair5D(){}
  
  bool operator<(const BinEntryPair5D& other) const
  {
    return ( bin < other.bin );
  }
  
  static unsigned int linearBin(unsigned int nbins2, unsigned int nbins3, unsigned int nbins4, unsigned int nbins5, unsigned int bin1, unsigned int bin2, unsigned int bin3, unsigned int bin4, unsigned int bin5)
  {
    return bin5 + nbins5*(bin4 + nbins4*(bin3 + nbins3*( bin2 + nbins2*bin1 )));
  }
  
  void bin5D(unsigned int nbins2, unsigned int nbins3, unsigned int nbins4, unsigned int nbins5, unsigned int& bin1, unsigned int& bin2, unsigned int& bin3, unsigned int& bin4, unsigned int& bin5) const
  {
    unsigned int temp1 = nbins2*nbins3*nbins4*nbins5;
    bin1 = bin/temp1;
    unsigned int temp2 = bin - bin1*temp1;
    temp1 = nbins3*nbins4*nbins5;
    bin2 = temp2/temp1;
    temp2 -= bin2*temp1;
    temp1 = nbins4*nbins5;
    bin3 = temp2/temp1;
    temp2 -= bin3*temp1;
    temp1 = nbins5;
    bin4 = temp2/temp1;
    temp2 -= bin4*temp1;
    bin5 = temp2;
  }
  
  unsigned int bin;
  unsigned int entry;
  bool is_seed;
};



class HelixHough
{
  public:
    HelixHough(unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, HelixResolution& min_resolution, HelixResolution& max_resolution, HelixRange& range);
    HelixHough(std::vector<std::vector<unsigned int> >& zoom_profile, unsigned int minzoom, HelixRange& range);
    virtual ~HelixHough();
    
    void initHelixHough(unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, HelixResolution& min_resolution, HelixResolution& max_resolution, HelixRange& range);
    
    void findHelices(std::vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks=0);
    
    void findHelices(unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks, unsigned int zoomlevel);
    
    void findSeededHelices(std::vector<SimpleTrack3D>& seeds, std::vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks=0);
    void findSeededHelices_run(std::vector<SimpleTrack3D>& seeds, std::vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks=0);
    void findSeededHelices(unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks, unsigned int zoomlevel);
    
    void vote(unsigned int zoomlevel);
    static void phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi_1, float* max_phi_1, float* min_phi_2, float* max_phi_2);
    static void phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi, float* max_phi, float hel, __m128& phi_3_out, __m128& phi_4_out);
    static void phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* max_k, float* min_phi, float* max_phi, float hel, __m128& phi_3, __m128& phi_4, __m128& phi_3_out, __m128& phi_4_out);
    static void dzdlRange_sse(float* x_a, float* y_a, float* z_a, float cosphi1, float sinphi1, float cosphi2, float sinphi2, float min_k, float max_k, float min_d, float max_d, float* min_z0, float* max_z0, float* min_dzdl_a, float* max_dzdl_a);
    static void phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi, float* max_phi, float* min_phi_2, float* max_phi_2, float hel, __m128& phi_3_out, __m128& phi_4_out, float* hit_x_2, float* hit_y_2, __m128& phi_3_out_2, __m128& phi_4_out_2);
    static void phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi, float* max_phi, float* min_phi_2, float* max_phi_2, float hel, __m128& phi_3, __m128& phi_4, __m128& phi_3_out, __m128& phi_4_out, float* hit_x_2, float* hit_y_2, __m128& phi_3_2, __m128& phi_4_2, __m128& phi_3_out_2, __m128& phi_4_out_2);
    
    void setPrintTimings(bool pt){print_timings=pt;}
    
    virtual void finalize(std::vector<SimpleTrack3D>& input, std::vector<SimpleTrack3D>& output){}
    virtual void findTracks(std::vector<SimpleHit3D>& hits, std::vector<SimpleTrack3D>& tracks, const HelixRange& range) = 0;
    virtual void initEvent(std::vector<SimpleHit3D>& hits, unsigned int min_hits){}
    virtual void findSeededTracks(std::vector<SimpleTrack3D>& seeds, std::vector<SimpleHit3D>& hits, std::vector<SimpleTrack3D>& tracks, const HelixRange& range){}
    
    // return true if we want to go straight into the user-written findTracks function instead of possibly continuing the Hough Transform
    virtual bool breakRecursion(const std::vector<SimpleHit3D>& hits, const HelixRange& range){return false;}
    
    // additional phi error to add to a hit in the voting stage, given an bin with the k and dzdl parameters passed.  This is useful 
    // to take into account multiple scattering
    virtual float phiError(SimpleHit3D& hit, float min_k, float max_k, float min_d, float max_d, float min_z0, float max_z0, float min_dzdl, float max_dzdl, bool pairvoting=false){return 0.;}
    // same as above, but for additional error on dzdl
    virtual float dzdlError(SimpleHit3D& hit, float min_k, float max_k, float min_d, float max_d, float min_z0, float max_z0, float min_dzdl, float max_dzdl, bool pairvoting=false){return 0.;}
    
    virtual void initSeeding(std::vector<SimpleTrack3D>& seeds){}
    
    virtual void setSeparateByHelicity(bool sbh){separate_by_helicity=sbh;}
    virtual void setOnlyOneHelicity(bool ooh){only_one_helicity=ooh;}
    void setHelicity(bool hel){helicity=hel;}
    
    virtual void requireLayers(unsigned int nl)
    {
      check_layers = true;
      req_layers = nl;
    }
    
    virtual void setBinScale(float b_scl){bin_scale = b_scl;}
    virtual void setZBinScale(float b_scl){z_bin_scale = b_scl;}
    virtual void setRemoveHits(bool rh){remove_hits=rh;}
    
    virtual void setRangeFromSeed(HelixRange& range, SimpleTrack3D& seed);
    
    virtual void setTopRange(HelixRange& tr){top_range = tr;}
    void splitIntoBins(unsigned int min_hits, unsigned int max_hits, std::vector<HelixRange>& ranges, std::vector<std::vector<SimpleHit3D> >& split_hits, unsigned int zoomlevel);
    
    virtual void clear(){}
    
    virtual void setStartZoom(unsigned int sz){start_zoom=sz;}
    
    virtual void setMaxHitsPairs(unsigned int mhp){max_hits_pairs=mhp;}
    
    void setClusterStartBin(unsigned int csb){cluster_start_bin=csb;}
    
    void setSeedStates(std::vector<HelixKalmanState>& states){seed_states = states;}
    
    std::vector<HelixKalmanState>& getKalmanStates(){return track_states;}
    
    void setLayersAtATime(unsigned int l){layers_at_a_time = l;}
    
    void setSmoothBack(bool sb){smooth_back=sb;}

    void setCullInputHits( bool cih ){ cull_input_hits = cih; }
    void setIterateClustering( bool icl ){ iterate_clustering = icl; }
    
  protected:
    bool remove_hits;
    std::vector<unsigned int>* hit_used;
    
    void setRange(const BinEntryPair5D& bp, HelixRange& range1, HelixRange& range2, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0);
    
    //number of bins for each zoom level
    std::vector<unsigned int> n_phi_bins;
    std::vector<unsigned int> n_d_bins;
    std::vector<unsigned int> n_k_bins;
    std::vector<unsigned int> n_dzdl_bins;
    std::vector<unsigned int> n_z0_bins;
    //top level hits
    std::vector<SimpleHit3D>* base_hits;
    // vector of hits used for isolating hits in a single zoom level
    std::vector<std::vector<SimpleHit3D>* > hits_vec;
    // for each zoomlevel, a vector of pair indexes, with the index mapping to the entry position in hits_vec
    std::vector<std::vector<std::pair<unsigned int,unsigned int> >* > pairs_vec;
    // stores external index values
    std::vector<unsigned int> index_mapping;
    /// vector of BinEntryPairs which list the parameter bins each hit gets voted into
    std::vector<std::vector<BinEntryPair5D>* > bins_vec;
    
    std::vector<std::vector<SimpleTrack3D>* > seeds_vec;
    
    HelixRange top_range;
    HelixRange current_range;
    std::vector<HelixRange> zoomranges;
    std::vector<std::vector<ParameterCluster>*> clusters_vec;
    std::vector<unsigned int> num_clusters;
    // vector on which to perform unique sorting
    std::vector<unsigned int> C_clus;
    // vector to contain merged hits from a cluster and candidate bin
    std::vector<unsigned int> temp_merged_clus;
    
    unsigned int max_zoom;//maximum number of times to zoom in on the parameter array
    unsigned int min_zoom;//minimum number of times to zoom in on the parameter array
    
    bool using_vertex;
    
    unsigned int max_tracks;
    
    double vote_time;
    double xy_vote_time;
    double z_vote_time;
    double cluster_time;
    
    bool print_timings;
    
    bool separate_by_helicity,helicity,only_one_helicity;
    
    bool check_layers;
    unsigned int req_layers;
    
    float bin_scale; // proportion of the bin size in each direction which is actually used for voting
    float z_bin_scale; // proportion of the bin size in the z-direction in each direction which is actually used for voting
    
    unsigned int start_zoom; // top-level zoomlevel : defaults to zero
    
    
    static void allButKappaRange_sse(float* x1_a,float* x2_a,float* y1_a,float* y2_a,float* z1_a,float* z2_a, float* min_k_a,float* max_k_a, float* min_phi_1_a,float* max_phi_1_a,float* min_phi_2_a,float* max_phi_2_a, float* min_d_1_a,float* max_d_1_a,float* min_d_2_a,float* max_d_2_a, float* min_dzdl_a,float* max_dzdl_a, float* min_z0_1_a,float* max_z0_1_a,float* min_z0_2_a,float* max_z0_2_a);
    
    
    void fillBins(unsigned int total_bins, unsigned int hit_counter, float* min_phi_a, float* max_phi_a, std::vector<SimpleHit3D>& four_hits, fastvec2d& z_bins, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, unsigned int d_bin, unsigned int k_bin, unsigned int n_phi, unsigned int zoomlevel, float low_phi, float high_phi, float inv_phi_range, fastvec& vote_array);
    
    void makeClusters(unsigned int zoomlevel, unsigned int MAX, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, unsigned int min_hits, std::vector<ParameterCluster>& clusters, bool& use_clusters, bool& is_super_bin);
    
    bool attemptClusterMerge(unsigned int zoomlevel, unsigned int MAX, unsigned int ca, unsigned int d, unsigned int r, unsigned int th, unsigned int zz0, unsigned int bin, unsigned int newbin, std::vector<unsigned char>& good_bins, unsigned int volume, float cluster_size_cut, float overlap_cut, std::vector<ParameterCluster>& clusters, unsigned int* bins_start, unsigned int* bins_end, std::vector<unsigned int>& map_clus, std::vector<unsigned char>& too_big, std::vector<unsigned int>& temp_merged, std::vector<unsigned int>& C);
    
    void vote_z(unsigned int zoomlevel, unsigned int n_phi, unsigned int n_d, unsigned int
    n_k, unsigned int n_dzdl, unsigned int n_z0, fastvec2d& z_bins);
    
    
    void vote_pairs(unsigned int zoomlevel);
    void fillBins(unsigned int total_bins, unsigned int pair_counter, unsigned int* pair_index, float* min_phi, float* max_phi, float* min_d, float* max_d, float* min_dzdl, float* max_dzdl, float* min_z0, float* max_z0, std::vector<std::vector<SimpleHit3D> > & four_pairs, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, unsigned int k_bin, unsigned int n_phi, unsigned int zoomlevel, float low_phi, float high_phi, float low_d, float high_d, float low_z0, float high_z0, float low_dzdl, float high_dzdl, float inv_phi_range, float inv_d_range, float inv_z0_range, float inv_dzdl_range, fastvec& vote_array);
    void findHelicesByPairsBegin(unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks, unsigned int zoomlevel);
    void findHelicesByPairs(unsigned int min_hits, unsigned int max_hits, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks, unsigned int zoomlevel);
    
    unsigned int max_hits_pairs;
    
    std::vector<std::pair<unsigned int, unsigned int> > temp_pairs;
    std::set<unsigned int> new_hits;
    std::map<unsigned int, unsigned int> old_to_new;
    
    unsigned int cluster_start_bin;
    
    unsigned int bins_start[1<<12];
    unsigned int bins_end[1<<12];
    
    unsigned int layers_at_a_time;
    
    std::vector<HelixKalmanState> track_states;
    std::vector<HelixKalmanState> seed_states;
    
    unsigned int n_layers;
    int layer_start;
    int layer_end;
    bool smooth_back;
    bool cull_input_hits;
    bool iterate_clustering;
};

#endif

