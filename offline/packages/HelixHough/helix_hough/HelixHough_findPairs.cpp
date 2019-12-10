#include "HelixHough.h"
#include "HelixRange.h"

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <sys/time.h>
#include <utility>
#include <vector>

class SimpleTrack3D;

using namespace std;

static inline void setClusterRange(HelixRange& range1, HelixRange& range2,
                                   ParRange& prange, unsigned int n_phi,
                                   unsigned int n_d, unsigned int n_k,
                                   unsigned int n_dzdl, unsigned int n_z0) {
  float dzdl_size = (range1.max_dzdl - range1.min_dzdl) / ((float)n_dzdl);
  float z0_size = (range1.max_z0 - range1.min_z0) / ((float)n_z0);
  float phi_size = (range1.max_phi - range1.min_phi) / ((float)n_phi);
  float d_size = (range1.max_d - range1.min_d) / ((float)n_d);
  float k_size = (range1.max_k - range1.min_k) / ((float)n_k);

  range2.min_phi = range1.min_phi + phi_size * ((float)(prange.min_phi));
  range2.max_phi = range1.min_phi + phi_size * ((float)(prange.max_phi + 1));
  range2.min_d = range1.min_d + d_size * ((float)(prange.min_d));
  range2.max_d = range1.min_d + d_size * ((float)(prange.max_d + 1));
  range2.min_k = range1.min_k + k_size * ((float)(prange.min_k));
  range2.max_k = range1.min_k + k_size * ((float)(prange.max_k + 1));
  range2.min_dzdl = range1.min_dzdl + dzdl_size * ((float)(prange.min_dzdl));
  range2.max_dzdl =
      range1.min_dzdl + dzdl_size * ((float)(prange.max_dzdl + 1));
  range2.min_z0 = range1.min_z0 + z0_size * ((float)(prange.min_z0));
  range2.max_z0 = range1.min_z0 + z0_size * ((float)(prange.max_z0 + 1));
}

static inline void setRange(const BinEntryPair5D& bp, HelixRange& range1,
                            HelixRange& range2, unsigned int n_phi,
                            unsigned int n_d, unsigned int n_k,
                            unsigned int n_dzdl, unsigned int n_z0) {
  float dzdl_size = (range1.max_dzdl - range1.min_dzdl) / ((float)n_dzdl);
  float z0_size = (range1.max_z0 - range1.min_z0) / ((float)n_z0);
  float phi_size = (range1.max_phi - range1.min_phi) / ((float)n_phi);
  float d_size = (range1.max_d - range1.min_d) / ((float)n_d);
  float k_size = (range1.max_k - range1.min_k) / ((float)n_k);

  unsigned int z0_bin = 0;
  unsigned int dzdl_bin = 0;
  unsigned int k_bin = 0;
  unsigned int d_bin = 0;
  unsigned int phi_bin = 0;

  bp.bin5D(n_d, n_k, n_dzdl, n_z0, phi_bin, d_bin, k_bin, dzdl_bin, z0_bin);
  range2.min_phi = range1.min_phi + phi_size * ((float)(phi_bin));
  range2.max_phi = range2.min_phi + phi_size;
  range2.min_d = range1.min_d + d_size * ((float)(d_bin));
  range2.max_d = range2.min_d + d_size;
  range2.min_k = range1.min_k + k_size * ((float)(k_bin));
  range2.max_k = range2.min_k + k_size;
  range2.min_dzdl = range1.min_dzdl + dzdl_size * ((float)(dzdl_bin));
  range2.max_dzdl = range2.min_dzdl + dzdl_size;
  range2.min_z0 = range1.min_z0 + z0_size * ((float)(z0_bin));
  range2.max_z0 = range2.min_z0 + z0_size;
}

void HelixHough::findHelicesByPairsBegin(unsigned int min_hits,
                                         unsigned int max_hits,
                                         vector<SimpleTrack3D>& tracks,
                                         unsigned int maxtracks,
                                         unsigned int zoomlevel) {
  pairs_vec[zoomlevel]->clear();
  pair<unsigned int, unsigned int> onepair;
  // make list of all pairs from the hits in this zoomlevel
  for (unsigned int i = 0; i < (hits_vec[zoomlevel]->size()); ++i) {
    //     if( (*(hits_vec[zoomlevel]))[i].get_layer() != 0 ){continue;}
    onepair.first = i;
    for (unsigned int j = 0; j < (hits_vec[zoomlevel]->size()); ++j) {
      if (((*(hits_vec[zoomlevel]))[j].get_layer() <=
           (*(hits_vec[zoomlevel]))[i].get_layer())) {
        continue;
      }

      //       if( (*(hits_vec[zoomlevel]))[i].get_layer() != 4 ){continue;}
      onepair.second = j;
      pairs_vec[zoomlevel]->push_back(onepair);
    }
  }
  findHelicesByPairs(min_hits, max_hits, tracks, maxtracks, zoomlevel);
}

void HelixHough::findHelicesByPairs(unsigned int min_hits,
                                    unsigned int max_hits,
                                    vector<SimpleTrack3D>& tracks,
                                    unsigned int maxtracks,
                                    unsigned int zoomlevel) {
  if ((maxtracks != 0) && (tracks.size() >= max_tracks)) {
    return;
  }
  unsigned int tracks_at_start = tracks.size();

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;
  if (print_timings == true) {
    gettimeofday(&t1, NULL);
  }
  vote_pairs(zoomlevel);
  if (print_timings == true) {
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    vote_time += (time2 - time1);
  }

  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];

  unsigned int n_entries = bins_vec[zoomlevel]->size();
  if (n_entries == 0) {
    return;
  }

  if (print_timings == true) {
    gettimeofday(&t1, NULL);
  }
  num_clusters[zoomlevel] = 0;
  bool use_clusters = true;
  bool is_super_bin = false;
  makeClusters(zoomlevel, pairs_vec[zoomlevel]->size(), n_phi, n_d, n_k, n_dzdl,
               n_z0, min_hits, *(clusters_vec[zoomlevel]), use_clusters,
               is_super_bin);
  if (print_timings) {
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    cluster_time += (time2 - time1);
  }

  pair<unsigned int, unsigned int> onepair;
  for (unsigned int i = 0, size = num_clusters[zoomlevel]; i < size; ++i) {
    temp_pairs.clear();
    new_hits.clear();
    old_to_new.clear();
    hits_vec[zoomlevel + 1]->clear();
    pairs_vec[zoomlevel + 1]->clear();
    vector<unsigned int>::iterator index_iter;
    for (index_iter = (*(clusters_vec[zoomlevel]))[i].hit_indexes.begin();
         index_iter != (*(clusters_vec[zoomlevel]))[i].hit_indexes.end();
         ++index_iter) {
      if (remove_hits == true) {
        if (((*hit_used)[(*(hits_vec[zoomlevel]))
                             [(*(pairs_vec[zoomlevel]))[*index_iter].first]
                                 .get_id()] == false) &&
            ((*hit_used)[(*(hits_vec[zoomlevel]))
                             [(*(pairs_vec[zoomlevel]))[*index_iter].second]
                                 .get_id()] == false)) {
          temp_pairs.push_back((*(pairs_vec[zoomlevel]))[*index_iter]);
          new_hits.insert((temp_pairs.back()).first);
          new_hits.insert((temp_pairs.back()).second);
        }
      } else {
        temp_pairs.push_back((*(pairs_vec[zoomlevel]))[*index_iter]);
        new_hits.insert((temp_pairs.back()).first);
        new_hits.insert((temp_pairs.back()).second);
      }
    }

    unsigned int temp_index = 0;
    set<unsigned int>::iterator it;
    for (it = new_hits.begin(); it != new_hits.end(); ++it) {
      hits_vec[zoomlevel + 1]->push_back((*(hits_vec[zoomlevel]))[*it]);
      old_to_new[*it] = temp_index;
      temp_index += 1;
    }
    for (unsigned int j = 0, size = temp_pairs.size(); j != size; ++j) {
      onepair.first = old_to_new[temp_pairs[j].first];
      onepair.second = old_to_new[temp_pairs[j].second];
      pairs_vec[zoomlevel + 1]->push_back(onepair);
    }
    setClusterRange(zoomranges[zoomlevel], zoomranges[zoomlevel + 1],
                    (*(clusters_vec[zoomlevel]))[i].range, n_phi, n_d, n_k,
                    n_dzdl, n_z0);
    if ((breakRecursion(*(hits_vec[zoomlevel + 1]),
                        zoomranges[zoomlevel + 1]) == true)) {
    } else if (((zoomlevel + 1) == max_zoom) ||
               (hits_vec[zoomlevel + 1]->size() <= max_hits)) {
      //           findTracksByPairs(*(hits_vec[zoomlevel+1]),
      //           *(pairs_vec[zoomlevel+1]), tracks, zoomranges[zoomlevel+1]);
      findTracks(*(hits_vec[zoomlevel + 1]), tracks, zoomranges[zoomlevel + 1]);
    } else {
      findHelicesByPairs(min_hits, max_hits, tracks, maxtracks, zoomlevel + 1);
    }
  }
}
