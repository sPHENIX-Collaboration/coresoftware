#include "HelixHough.h"
#include "HelixRange.h"
#include "SimpleHit3D.h"
#include "SimpleTrack3D.h"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <memory>
#include <sys/time.h>
#include <utility>
#include <vector>

using namespace std;

static inline void in_place_counting_sort_unique(vector<unsigned int>& A,
                                                 vector<unsigned int>& C,
                                                 unsigned int MAX) {
  unsigned int SIZE = A.size();
  for (unsigned int i = 0; i < SIZE; ++i) {
    ++C[A[i]];
  }
  unsigned int current = 0;
  for (unsigned int i = 0; i < MAX; ++i) {
    A[current] = ((((C[i] != 0) - 1) & (A[current])) ^ (((C[i] == 0) - 1) & i));
    current += (C[i] != 0);
  }
  A.resize(current);
}

static inline void in_place_counting_unique(vector<unsigned int>& A,
                                            vector<unsigned int>& C) {
  unsigned int SIZE = A.size();
  unsigned int current = 0;
  for (unsigned int i = 0; i < SIZE; ++i) {
    C[A[i]] += 1;
    if (C[A[i]] == 1) {
      A[current] = A[i];
      current += 1;
    }
  }
  A.resize(current);
}

void HelixHough::setRange(const BinEntryPair5D& bp, HelixRange& range1,
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

void HelixHough::findHelices(vector<SimpleHit3D>& hits_init, unsigned int min_hits,
                             unsigned int max_hits,
                             vector<SimpleTrack3D>& tracks,
                             unsigned int maxtracks) {
  vote_time = 0.;
  xy_vote_time = 0.;
  z_vote_time = 0.;
  cluster_time = 0.;

  vector<SimpleHit3D> hits;

  if (layer_end > 0 && layer_start > 0){
    for (unsigned int in = 0; in < hits_init.size(); in++) {
    if(hits_init[in].get_layer() > layer_start-1 && hits_init[in].get_layer() < layer_end+1 )
    hits.push_back(hits_init[in]);
    }
  } else if ( layer_end < 0 && layer_start > 0) {
    for (unsigned int in = 0; in < hits_init.size(); in++) {
    if(hits_init[in].get_layer() > layer_start-1)
    hits.push_back(hits_init[in]);
    }
  } else if ( layer_end > 0 && layer_start < 0) {
    for (unsigned int in = 0; in < hits_init.size(); in++) {
    if( hits_init[in].get_layer() < layer_end+1 )
    hits.push_back(hits_init[in]);
    }
  } else {
    for (unsigned int in = 0; in < hits_init.size(); in++) {
    hits.push_back(hits_init[in]);
    }
  }

  index_mapping.clear();
  index_mapping.resize(hits.size(), 0);
  hit_used->clear();
  for (unsigned int i = 0; i < hits.size(); i++) {
    index_mapping[i] = hits[i].get_id();
    hits[i].set_id(i);
  }
  (*hit_used).assign(hits.size(), false);

  initEvent(hits, min_hits);

  max_tracks = maxtracks;
  base_hits = &hits;
  (*(hits_vec[start_zoom])) = hits;
  current_range = top_range;
  zoomranges.clear();
  for (unsigned int z = 0; z <= max_zoom; z++) {
    zoomranges.push_back(top_range);
  }
  vector<SimpleTrack3D> temp_tracks;
  if ((separate_by_helicity == true) && (only_one_helicity == false)) {
    helicity = true;
    findHelices(min_hits, max_hits, temp_tracks, maxtracks, start_zoom);
    helicity = false;
    findHelices(min_hits, max_hits, temp_tracks, maxtracks, start_zoom);
  } else {
    findHelices(min_hits, max_hits, temp_tracks, maxtracks, start_zoom);
  }

  vector<SimpleHit3D> tr_hits;
  if (cull_input_hits == true) {
    for (unsigned int i = 0; i < hits.size(); i++) {
      if ((*hit_used)[hits[i].get_id()] == false) {
        tr_hits.push_back(hits[i]);
        tr_hits.back().set_id(index_mapping[i]);
      }
    }
  }

  for (unsigned int i = 0; i < hits.size(); i++) {
    hits[i].set_id(index_mapping[i]);
  }
  for (unsigned int t = 0; t < temp_tracks.size(); t++) {
    for (unsigned int h = 0; h < temp_tracks[t].hits.size(); h++) {
      if (temp_tracks[t].hits[h].get_id() != (unsigned)-1) {
        temp_tracks[t].hits[h].set_id(index_mapping[temp_tracks[t].hits[h].get_id()]);
      }
    }
  }

  finalize(temp_tracks, tracks);

  if (cull_input_hits == true) {
    hits = tr_hits;
  }

  if (print_timings == true) {
    cout << "vote time = " << vote_time << endl;
    cout << "xy vote time = " << xy_vote_time << endl;
    cout << "z vote time = " << z_vote_time << endl;
    cout << "cluster time = " << cluster_time << endl;
  }
}

static bool mergeClusters(ParameterCluster& clus1,
                          ParameterCluster const& clus2, unsigned int MAX,
                          float overlap_cut) {
  vector<unsigned int> old_indexes = clus1.hit_indexes;
  for (unsigned int i = 0; i < clus2.hit_indexes.size(); ++i) {
    clus1.hit_indexes.push_back(clus2.hit_indexes[i]);
  }
  vector<unsigned int> C;
  C.assign(MAX + 1, 0);
  in_place_counting_unique(clus1.hit_indexes, C);

  float size_diff_1 = ((float)(clus1.hit_indexes.size() - old_indexes.size()) /
                       ((float)(old_indexes.size())));
  float size_diff_2 =
      ((float)(clus1.hit_indexes.size() - clus2.hit_indexes.size()) /
       ((float)(clus2.hit_indexes.size())));

  if ((size_diff_1 < overlap_cut) || (size_diff_2 < overlap_cut)) {
    clus1.range.mergeRange(clus2.range);
    return true;
  } else {
    clus1.hit_indexes = old_indexes;
    return false;
  }
}

bool HelixHough::attemptClusterMerge(
    unsigned int zoomlevel, unsigned int MAX, unsigned int ca, unsigned int d,
    unsigned int r, unsigned int th, unsigned int zz0, unsigned int bin,
    unsigned int newbin, vector<unsigned char>& good_bins, unsigned int volume,
    float cluster_size_cut, float overlap_cut,
    vector<ParameterCluster>& clusters, unsigned int* bins_start,
    unsigned int* bins_end, vector<unsigned int>& map_clus,
    vector<unsigned char>& too_big, vector<unsigned int>& temp_merged,
    vector<unsigned int>& C) {
  if ((good_bins[newbin] == 1) && (map_clus[newbin] != 4294967295)) {
    if (too_big[map_clus[newbin]] == 0) {
      unsigned int tempsize = clusters[map_clus[newbin]].hit_indexes.size();
      for (unsigned int ind = bins_start[bin]; ind <= bins_end[bin]; ++ind) {
        clusters[map_clus[newbin]].hit_indexes.push_back(
            (*(bins_vec[zoomlevel]))[ind].entry);
      }
      C.clear();
      C.assign(MAX + 1, 0);
      in_place_counting_unique(clusters[map_clus[newbin]].hit_indexes, C);
      unsigned int size_diff =
          clusters[map_clus[newbin]].hit_indexes.size() - tempsize;
      unsigned int size_diff_2 = clusters[map_clus[newbin]].hit_indexes.size() -
                                 (1 + bins_end[bin] - bins_start[bin]);
      if ((((float)size_diff) / ((float)(1 + bins_end[bin] - bins_start[bin])) <
           overlap_cut) ||
          (((float)size_diff_2) / ((float)(tempsize)) < overlap_cut)) {
        clusters[map_clus[newbin]].range.mergeRange(ca, d, r, th, zz0);
        ParameterCluster* cluster = &(clusters[map_clus[newbin]]);
        unsigned int cluster_volume =
            ((cluster->range.max_phi - cluster->range.min_phi) + 1) *
            ((cluster->range.max_d - cluster->range.min_d) + 1) *
            ((cluster->range.max_k - cluster->range.min_k) + 1) *
            ((cluster->range.max_dzdl - cluster->range.min_dzdl) + 1) *
            ((cluster->range.max_z0 - cluster->range.min_z0) + 1);
        if (((float)cluster_volume) / ((float)volume) > cluster_size_cut) {
          too_big[map_clus[newbin]] = 1;
        }
        map_clus[bin] = map_clus[newbin];
        return true;
      } else {
        clusters[map_clus[newbin]].hit_indexes.resize(tempsize);
      }
    }
  }
  return false;
}

void HelixHough::makeClusters(unsigned int zoomlevel, unsigned int MAX,
                              unsigned int n_phi, unsigned int n_d,
                              unsigned int n_k, unsigned int n_dzdl,
                              unsigned int n_z0, unsigned int min_hits,
                              vector<ParameterCluster>& clusters,
                              bool& use_clusters, bool& is_super_bin) {
  unsigned int volume = n_phi * n_d * n_k * n_dzdl * n_z0;
  float cluster_size_cut = 1.0;
  float bin_size_cut = 0.75;
  float overlap_cut = 0.1;
  is_super_bin = false;

  vector<unsigned int> map_clus(volume, 4294967295);
  vector<unsigned char> good_bins(volume, 0);
  vector<unsigned char> too_big(volume, 0);

  for (unsigned int ca = 0; ca < n_phi; ++ca) {
    for (unsigned int d = 0; d < n_d; ++d) {
      for (unsigned int r = 0; r < n_k; ++r) {
        for (unsigned int th = 0; th < n_dzdl; ++th) {
          for (unsigned int zz0 = 0; zz0 < n_z0; ++zz0) {
            unsigned int bin = BinEntryPair5D::linearBin(n_d, n_k, n_dzdl, n_z0,
                                                         ca, d, r, th, zz0);

            if (bins_end[bin] == 4294967295) {
              continue;
            }
            if ((1 + bins_end[bin] - bins_start[bin]) >= min_hits) {
              if (check_layers == true) {
                unsigned int layer_mask[4] = {0, 0, 0, 0};
                for (unsigned int i = bins_start[bin]; i <= bins_end[bin];
                     ++i) {
                  if ((*(hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[i]
                                                   .entry]
                          .get_layer() < 32) {
                    layer_mask[0] =
                        layer_mask[0] |
                        (1 << (*(
                             hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[i]
                                                       .entry]
                                  .get_layer());
                  } else if ((*(hits_vec[zoomlevel]))
                                 [(*(bins_vec[zoomlevel]))[i].entry]
                                     .get_layer() < 64) {
                    layer_mask[1] =
                        layer_mask[1] |
                        (1 << ((*(hits_vec[zoomlevel]))
                                   [(*(bins_vec[zoomlevel]))[i].entry]
                                       .get_layer() -
                               32));
                  } else if ((*(hits_vec[zoomlevel]))
                                 [(*(bins_vec[zoomlevel]))[i].entry]
                                     .get_layer() < 96) {
                    layer_mask[2] =
                        layer_mask[2] |
                        (1 << ((*(hits_vec[zoomlevel]))
                                   [(*(bins_vec[zoomlevel]))[i].entry]
                                       .get_layer() -
                               64));
                  } else if ((*(hits_vec[zoomlevel]))
                                 [(*(bins_vec[zoomlevel]))[i].entry]
                                     .get_layer() < 128) {
                    layer_mask[3] =
                        layer_mask[3] |
                        (1 << ((*(hits_vec[zoomlevel]))
                                   [(*(bins_vec[zoomlevel]))[i].entry]
                                       .get_layer() -
                               96));
                  }
                }
                unsigned int nlayers = __builtin_popcount(layer_mask[0]) +
                                       __builtin_popcount(layer_mask[1]) +
                                       __builtin_popcount(layer_mask[2]) +
                                       __builtin_popcount(layer_mask[3]);
                if (nlayers >= req_layers) {
                  good_bins[bin] = 1;
                }
              } else {
                good_bins[bin] = 1;
              }
            } else {
              continue;
            }

            if (good_bins[bin] == 0) {
              continue;
            }

            if (ca > 0) {
              unsigned int newbin = bin - n_d * n_k * n_dzdl * n_z0;
              if (attemptClusterMerge(
                      zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins,
                      volume, cluster_size_cut, overlap_cut, clusters,
                      (unsigned int*)(bins_start), (unsigned int*)(bins_end),
                      map_clus, too_big, temp_merged_clus, C_clus) == true) {
                continue;
              }
            }

            if (d > 0) {
              unsigned int newbin = bin - n_k * n_dzdl * n_z0;
              if (attemptClusterMerge(
                      zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins,
                      volume, cluster_size_cut, overlap_cut, clusters,
                      (unsigned int*)(bins_start), (unsigned int*)(bins_end),
                      map_clus, too_big, temp_merged_clus, C_clus) == true) {
                continue;
              }
            }

            if (r > 0) {
              unsigned int newbin = bin - n_dzdl * n_z0;
              if (attemptClusterMerge(
                      zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins,
                      volume, cluster_size_cut, overlap_cut, clusters,
                      (unsigned int*)(bins_start), (unsigned int*)(bins_end),
                      map_clus, too_big, temp_merged_clus, C_clus) == true) {
                continue;
              }
            }

            if (th > 0) {
              unsigned int newbin = bin - n_z0;
              if (attemptClusterMerge(
                      zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins,
                      volume, cluster_size_cut, overlap_cut, clusters,
                      (unsigned int*)(bins_start), (unsigned int*)(bins_end),
                      map_clus, too_big, temp_merged_clus, C_clus) == true) {
                continue;
              }
            }

            if (zz0 > 0) {
              unsigned int newbin = bin - 1;
              if (attemptClusterMerge(
                      zoomlevel, MAX, ca, d, r, th, zz0, bin, newbin, good_bins,
                      volume, cluster_size_cut, overlap_cut, clusters,
                      (unsigned int*)(bins_start), (unsigned int*)(bins_end),
                      map_clus, too_big, temp_merged_clus, C_clus) == true) {
                continue;
              }
            }
            if (num_clusters[zoomlevel] >= clusters.size()) {
              clusters.push_back(ParameterCluster());
            }
            clusters[num_clusters[zoomlevel]].range =
                ParRange(ca, ca, d, d, r, r, th, th, zz0, zz0);
            map_clus[bin] = num_clusters[zoomlevel];
            clusters[num_clusters[zoomlevel]].hit_indexes.clear();
            for (unsigned int ind = bins_start[bin]; ind <= bins_end[bin];
                 ++ind) {
              clusters[num_clusters[zoomlevel]].hit_indexes.push_back(
                  (*(bins_vec[zoomlevel]))[ind].entry);
            }

            num_clusters[zoomlevel] += 1;
          }
        }
      }
    }
  }
  if (iterate_clustering == false) {
    return;
  }
  if (num_clusters[zoomlevel] == 0) {
    return;
  }
  vector<ParameterCluster> in_clusters;
  for (unsigned int i = 0; i < num_clusters[zoomlevel]; ++i) {
    in_clusters.push_back(clusters[i]);
  }
  vector<ParameterCluster> out_clusters;
  out_clusters.push_back(in_clusters[0]);
  for (unsigned int i = 1; i < in_clusters.size(); ++i) {
    bool merged = false;
    for (unsigned int j = 0; j < out_clusters.size(); ++j) {
      merged = mergeClusters(out_clusters[j], in_clusters[i], MAX, overlap_cut);
      if (merged == true) {
        merged = true;
        break;
      }
    }
    if (merged == false) {
      out_clusters.push_back(in_clusters[i]);
    }
  }
  clusters = out_clusters;
  num_clusters[zoomlevel] = out_clusters.size();
}

void HelixHough::findHelices(unsigned int min_hits, unsigned int max_hits,
                             vector<SimpleTrack3D>& tracks,
                             unsigned int maxtracks, unsigned int zoomlevel) {
  unsigned int tracks_at_start = tracks.size();

  if ((maxtracks != 0) && (tracks.size() >= max_tracks)) {
    return;
  }

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;
  if (print_timings == true) {
    gettimeofday(&t1, NULL);
  }
  vote(zoomlevel);
  if (print_timings == true) {
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    vote_time += (time2 - time1);
  }

  unsigned int n_entries = bins_vec[zoomlevel]->size();
  if (n_entries == 0) {
    return;
  }

  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];

  num_clusters[zoomlevel] = 0;
  bool use_clusters = true;
  bool is_super_bin = false;
  if (zoomlevel >= cluster_start_bin) {
    if (print_timings == true) {
      gettimeofday(&t1, NULL);
    }
    makeClusters(zoomlevel, hits_vec[zoomlevel]->size(), n_phi, n_d, n_k,
                 n_dzdl, n_z0, min_hits, *(clusters_vec[zoomlevel]),
                 use_clusters, is_super_bin);

    if (print_timings) {
      gettimeofday(&t2, NULL);
      time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
      time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
      cluster_time += (time2 - time1);
    }
  } else {
    use_clusters = false;
  }

  if (use_clusters == false) {
    unsigned int count = 0;
    hits_vec[zoomlevel + 1]->clear();
    setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel],
             zoomranges[zoomlevel + 1], n_phi, n_d, n_k, n_dzdl, n_z0);
    // scan over the bins in 5-D hough space
    while (count < n_entries) {
      if (remove_hits == true) {
        if ((*hit_used)[(*(
                hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]
                            .get_id()] == false) {
          hits_vec[zoomlevel + 1]->push_back(
              (*(hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]);
        }
      } else {
        hits_vec[zoomlevel + 1]->push_back(
            (*(hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]);
      }

      count += 1;
      // we have collected all hits from this bin. now zoom again or find tracks
      // with user routine
      if ((count == n_entries) || ((*(bins_vec[zoomlevel]))[count].bin !=
                                   (*(bins_vec[zoomlevel]))[count - 1].bin)) {
        if (hits_vec[zoomlevel + 1]->size() >= min_hits) {
          if (breakRecursion(*(hits_vec[zoomlevel + 1]),
                             zoomranges[zoomlevel + 1]) == true) {
          } else if (((zoomlevel + 1) == max_zoom)) {
            findTracks(*(hits_vec[zoomlevel + 1]), tracks,
                       zoomranges[zoomlevel + 1]);
          } else if (((zoomlevel + 1) >= min_zoom) &&
                     ((hits_vec[zoomlevel + 1]->size() <= max_hits))) {
            findTracks(*(hits_vec[zoomlevel + 1]), tracks,
                       zoomranges[zoomlevel + 1]);
          } else if ((hits_vec[zoomlevel + 1]->size() <= max_hits_pairs) &&
                     ((zoomlevel + 1) >= min_zoom)) {
            findHelicesByPairsBegin(min_hits, max_hits, tracks, maxtracks,
                                    zoomlevel + 1);
          } else {
            findHelices(min_hits, max_hits, tracks, maxtracks, zoomlevel + 1);
          }
          if (maxtracks != 0) {
            double phi_proportion =
                (zoomranges[zoomlevel].max_phi -
                 zoomranges[zoomlevel].min_phi) /
                ((zoomranges[0].max_phi - zoomranges[0].min_phi));
            double d_proportion =
                (zoomranges[zoomlevel].max_d - zoomranges[zoomlevel].min_d) /
                ((zoomranges[0].max_d - zoomranges[0].min_d));
            double k_proportion =
                (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k) /
                ((zoomranges[0].max_k - zoomranges[0].min_k));
            unsigned int expected_tracks =
                (unsigned int)(fabs(((double)maxtracks) * phi_proportion *
                                    d_proportion * k_proportion)) +
                1;

            if ((tracks.size() - tracks_at_start) > expected_tracks) {
              return;
            }
          }
        }
        if (count == n_entries) {
          break;
        }
        hits_vec[zoomlevel + 1]->clear();
        setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel],
                 zoomranges[zoomlevel + 1], n_phi, n_d, n_k, n_dzdl, n_z0);
      }
    }
  } else {
    if (clusters_vec[zoomlevel]->size() == 0) {
      return;
    }
    // for each cluster, eiter perform the Hough again or break into
    // user-defined routine
    for (unsigned int i = 0, size = num_clusters[zoomlevel]; i < size; ++i) {
      hits_vec[zoomlevel + 1]->clear();
      vector<unsigned int>::iterator index_iter;
      for (index_iter = (*(clusters_vec[zoomlevel]))[i].hit_indexes.begin();
           index_iter != (*(clusters_vec[zoomlevel]))[i].hit_indexes.end();
           ++index_iter) {
        if (remove_hits == true) {
          if ((*hit_used)[(*(hits_vec[zoomlevel]))[*index_iter].get_id()] ==
              false) {
            hits_vec[zoomlevel + 1]->push_back(
                (*(hits_vec[zoomlevel]))[*index_iter]);
          }
        } else {
          hits_vec[zoomlevel + 1]->push_back(
              (*(hits_vec[zoomlevel]))[*index_iter]);
        }
      }

      setClusterRange(zoomranges[zoomlevel], zoomranges[zoomlevel + 1],
                      (*(clusters_vec[zoomlevel]))[i].range, n_phi, n_d, n_k,
                      n_dzdl, n_z0);
      if ((breakRecursion(*(hits_vec[zoomlevel + 1]),
                          zoomranges[zoomlevel + 1]) == true)) {
      } else if ((zoomlevel + 1) == max_zoom) {
        findTracks(*(hits_vec[zoomlevel + 1]), tracks,
                   zoomranges[zoomlevel + 1]);
      } else if (((zoomlevel + 1) >= min_zoom) &&
                 (hits_vec[zoomlevel + 1]->size() <= max_hits)) {
        findTracks(*(hits_vec[zoomlevel + 1]), tracks,
                   zoomranges[zoomlevel + 1]);
      } else if ((hits_vec[zoomlevel + 1]->size() <= max_hits_pairs) &&
                 ((zoomlevel + 1) >= min_zoom)) {
        findHelicesByPairsBegin(min_hits, max_hits, tracks, maxtracks,
                                zoomlevel + 1);
      } else {
        findHelices(min_hits, max_hits, tracks, maxtracks, zoomlevel + 1);
      }
    }
  }
}

void HelixHough::findSeededHelices(vector<SimpleTrack3D>& seeds,
                                   vector<SimpleHit3D>& hits,
                                   unsigned int min_hits, unsigned int max_hits,
                                   vector<SimpleTrack3D>& tracks,
                                   unsigned int maxtracks) {
  unsigned int n_layers_orig = n_layers;
  unsigned int req_layers_orig = req_layers;

  bool smooth_back_orig = smooth_back;

  int max_layer = 0;
  int min_layer = 999999;
  for (unsigned int i = 0; i < hits.size(); ++i) {
    if (hits[i].get_layer() > max_layer) {
      max_layer = hits[i].get_layer();
    }
    if (hits[i].get_layer() < min_layer) {
      min_layer = hits[i].get_layer();
    }
    if (hits[i].get_layer() < 0) {
      findSeededHelices_run(seeds, hits, min_hits, max_hits, tracks, maxtracks);
      return;
    }
  }

  if (min_layer == 999999) {
    findSeededHelices_run(seeds, hits, min_hits, max_hits, tracks, maxtracks);
    return;
  }

  n_layers = min_layer;

  float req_prop = ((float)req_layers) / ((float)((1 + max_layer - min_layer)));

  int layers_left = (1 + max_layer - min_layer);
  unsigned int cur_layer = min_layer;
  vector<SimpleTrack3D> cur_seeds = seeds;
  vector<SimpleTrack3D> cur_tracks;
  vector<SimpleHit3D> cur_hits;
  smooth_back = false;
  while (layers_left > 0) {
    unsigned int layers_iter = layers_at_a_time;
    if ((cur_layer + layers_at_a_time - 1) > max_layer) {
      layers_iter = (max_layer - (cur_layer - 1));
    }

    cur_tracks.clear();
    cur_hits.clear();
    for (unsigned int i = 0; i < hits.size(); ++i) {
      if ((hits[i].get_layer() >= cur_layer) &&
          (hits[i].get_layer() <= (cur_layer + layers_at_a_time - 1))) {
        cur_hits.push_back(hits[i]);
      }
    }

    n_layers += layers_iter;

    unsigned int nreq = floor(((float)(layers_iter)) * req_prop);
    if (nreq == 0) {
      nreq = 1;
    }
    requireLayers(nreq);

    clear();
    req_layers = nreq;

    if (layers_left <= layers_at_a_time) {
      smooth_back = smooth_back_orig;
    }

    findSeededHelices_run(cur_seeds, cur_hits, nreq, max_hits, cur_tracks,
                          maxtracks);
    setSeedStates(track_states);
    cur_seeds = cur_tracks;
    layers_left -= layers_at_a_time;
    cur_layer += layers_at_a_time;
  }
  tracks = cur_tracks;

  n_layers = n_layers_orig;
  req_layers = req_layers_orig;
}

void HelixHough::findSeededHelices_run(vector<SimpleTrack3D>& seeds,
                                       vector<SimpleHit3D>& hits,
                                       unsigned int min_hits,
                                       unsigned int max_hits,
                                       vector<SimpleTrack3D>& tracks,
                                       unsigned int maxtracks) {
  vote_time = 0.;
  xy_vote_time = 0.;
  z_vote_time = 0.;

  initEvent(hits, min_hits);
  initSeeding(seeds);

  hit_used->clear();
  index_mapping.clear();
  index_mapping.resize(hits.size(), 0);
  for (unsigned int i = 0; i < hits.size(); i++) {
    index_mapping[i] = hits[i].get_id();
    hits[i].set_id(i);
  }
  (*hit_used).assign(hits.size(), false);

  max_tracks = maxtracks;
  base_hits = &hits;
  (*(hits_vec[0])) = hits;
  (*(seeds_vec[0])) = seeds;
  current_range = top_range;
  zoomranges.clear();
  for (unsigned int z = 0; z <= max_zoom; z++) {
    zoomranges.push_back(top_range);
  }
  vector<SimpleTrack3D> temp_tracks;

  if (separate_by_helicity == true) {
    helicity = true;
    findSeededHelices(min_hits, max_hits, temp_tracks, maxtracks, 0);
    helicity = false;
    findSeededHelices(min_hits, max_hits, temp_tracks, maxtracks, 0);
  } else {
    findSeededHelices(min_hits, max_hits, temp_tracks, maxtracks, 0);
  }

  for (unsigned int i = 0; i < hits.size(); i++) {
    hits[i].set_id(index_mapping[i]);
  }

  finalize(temp_tracks, tracks);

  if (print_timings == true) {
    cout << "vote time = " << vote_time << endl;
    cout << "xy vote time = " << xy_vote_time << endl;
    cout << "z vote time = " << z_vote_time << endl;
    cout << endl;
  }
}

class floatBin {
 public:
  floatBin(float l, float h) : low(l), high(h) {}
  ~floatBin() {}
  float low;
  float high;
  bool operator<(const floatBin& other) const { return (high < other.low); }
};

// return which bin in bins that val belongs to, or -1 if it doesn't belong to
// any
static int seed_bin(vector<floatBin>& bins, float val) {
  floatBin bin(val, val);
  if ((bin < bins[0]) || (bins.back() < bin)) {
    return -1;
  }
  pair<vector<floatBin>::iterator, vector<floatBin>::iterator> bounds;
  bounds = equal_range(bins.begin(), bins.end(), bin);
  return ((int)(bounds.first - bins.begin()));
}

void HelixHough::setRangeFromSeed(HelixRange& range, SimpleTrack3D& seed) {
  range.min_phi = seed.phi - 0.001;
  range.max_phi = seed.phi + 0.001;
  if (range.min_phi < 0.) {
    range.min_phi = 0.;
  }
  if (range.max_phi > 2. * M_PI) {
    range.max_phi = 2. * M_PI;
  }
  range.min_d = seed.d - 0.001;
  range.max_d = seed.d + 0.001;
  range.min_k = seed.kappa - 0.00001;
  range.max_k = seed.kappa + 0.00001;
  if (range.min_k < 0.) {
    range.min_k = 0.;
  }

  range.min_k = range.min_k * range.min_k;
  range.max_k = range.max_k * range.max_k;

  range.min_dzdl = seed.dzdl - 0.05;
  range.max_dzdl = seed.dzdl + 0.05;
  range.min_z0 = seed.z0 - 0.01;
  range.max_z0 = seed.z0 + 0.01;
}

void HelixHough::findSeededHelices(unsigned int min_hits, unsigned int max_hits,
                                   vector<SimpleTrack3D>& tracks,
                                   unsigned int maxtracks,
                                   unsigned int zoomlevel) {
  unsigned int tracks_at_start = tracks.size();

  if ((maxtracks != 0) && (tracks.size() >= max_tracks)) {
    return;
  }

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;
  gettimeofday(&t1, NULL);
  vote(zoomlevel);
  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
  vote_time += (time2 - time1);

  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];

  float low = 0.;
  float high = 0.;
  float delta = 0.;

  vector<floatBin> phibins;
  low = zoomranges[zoomlevel].min_phi;
  delta = (zoomranges[zoomlevel].max_phi - zoomranges[zoomlevel].min_phi) /
          ((float)(n_phi));
  high = low + delta;
  for (unsigned int i = 0; i < n_phi; ++i) {
    phibins.push_back(floatBin(low, high));
    low += delta;
    high += delta;
  }

  vector<floatBin> dbins;
  low = zoomranges[zoomlevel].min_d;
  delta = (zoomranges[zoomlevel].max_d - zoomranges[zoomlevel].min_d) /
          ((float)(n_d));
  high = low + delta;
  for (unsigned int i = 0; i < n_d; ++i) {
    dbins.push_back(floatBin(low, high));
    low += delta;
    high += delta;
  }

  vector<floatBin> kbins;
  low = zoomranges[zoomlevel].min_k;
  delta = (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k) /
          ((float)(n_k));
  high = low + delta;
  for (unsigned int i = 0; i < n_k; ++i) {
    kbins.push_back(floatBin(low, high));
    low += delta;
    high += delta;
  }

  vector<floatBin> z0bins;
  low = zoomranges[zoomlevel].min_z0;
  delta = (zoomranges[zoomlevel].max_z0 - zoomranges[zoomlevel].min_z0) /
          ((float)(n_z0));
  high = low + delta;
  for (unsigned int i = 0; i < n_z0; ++i) {
    z0bins.push_back(floatBin(low, high));
    low += delta;
    high += delta;
  }

  vector<floatBin> dzdlbins;
  low = zoomranges[zoomlevel].min_dzdl;
  delta = (zoomranges[zoomlevel].max_dzdl - zoomranges[zoomlevel].min_dzdl) /
          ((float)(n_dzdl));
  high = low + delta;
  for (unsigned int i = 0; i < n_dzdl; ++i) {
    dzdlbins.push_back(floatBin(low, high));
    low += delta;
    high += delta;
  }

  // voting for the seeds
  for (unsigned int i = 0; i < seeds_vec[zoomlevel]->size(); ++i) {
    float d = (*(seeds_vec[zoomlevel]))[i].d;
    float phi = (*(seeds_vec[zoomlevel]))[i].phi;
    float kappa = pow((*(seeds_vec[zoomlevel]))[i].kappa, 1.);
    float z0 = (*(seeds_vec[zoomlevel]))[i].z0;
    float dzdl = (*(seeds_vec[zoomlevel]))[i].dzdl;

    int tempbin = 0;

    unsigned int phi_bin = 0;
    tempbin = seed_bin(phibins, phi);
    if (tempbin < 0) {
      continue;
    } else {
      phi_bin = (unsigned int)(tempbin);
    }

    unsigned int d_bin = 0;
    tempbin = seed_bin(dbins, d);
    if (tempbin < 0) {
      continue;
    } else {
      d_bin = (unsigned int)(tempbin);
    }

    unsigned int kappa_bin = 0;
    tempbin = seed_bin(kbins, kappa);
    if (tempbin < 0) {
      continue;
    } else {
      kappa_bin = (unsigned int)(tempbin);
    }

    unsigned int z0_bin = 0;
    tempbin = seed_bin(z0bins, z0);
    if (tempbin < 0) {
      continue;
    } else {
      z0_bin = (unsigned int)(tempbin);
    }

    unsigned int dzdl_bin = 0;
    tempbin = seed_bin(dzdlbins, dzdl);
    if (tempbin < 0) {
      continue;
    } else {
      dzdl_bin = (unsigned int)(tempbin);
    }

    unsigned int bin = BinEntryPair5D::linearBin(
        n_d, n_k, n_dzdl, n_z0, phi_bin, d_bin, kappa_bin, dzdl_bin, z0_bin);

    bins_vec[zoomlevel]->push_back(BinEntryPair5D(bin, i, true));
  }
  sort(bins_vec[zoomlevel]->begin(), bins_vec[zoomlevel]->end());

  unsigned int n_entries = bins_vec[zoomlevel]->size();
  if (n_entries == 0) {
    return;
  }

  unsigned int count = 0;
  hits_vec[zoomlevel + 1]->clear();
  seeds_vec[zoomlevel + 1]->clear();
  setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel],
           zoomranges[zoomlevel + 1], n_phi, n_d, n_k, n_dzdl, n_z0);
  // scan over the bins in 5-D hough space
  while (count < n_entries) {
    if ((*(bins_vec[zoomlevel]))[count].is_seed == false) {
      if (remove_hits == true) {
        if ((*hit_used)[(*(
                hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]
                            .get_id()] == false) {
          hits_vec[zoomlevel + 1]->push_back(
              (*(hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]);
        }
      } else {
        hits_vec[zoomlevel + 1]->push_back(
            (*(hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]);
      }
    } else {
      seeds_vec[zoomlevel + 1]->push_back(
          (*(seeds_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]);
    }

    count++;

    // we have collected all hits from this bin. now zoom again or find tracks
    // with user routine
    if ((count == n_entries) || ((*(bins_vec[zoomlevel]))[count].bin !=
                                 (*(bins_vec[zoomlevel]))[count - 1].bin)) {
      if ((hits_vec[zoomlevel + 1]->size() >= min_hits) &&
          (seeds_vec[zoomlevel + 1]->size() != 0)) {
        if (breakRecursion(*(hits_vec[zoomlevel + 1]),
                           zoomranges[zoomlevel + 1]) == true) {
          findSeededTracks(*(seeds_vec[zoomlevel + 1]),
                           *(hits_vec[zoomlevel + 1]), tracks,
                           zoomranges[zoomlevel + 1]);
        } else if (((zoomlevel + 1) >= min_zoom) &&
                   ((hits_vec[zoomlevel + 1]->size() <= max_hits) ||
                    (zoomlevel + 1) == max_zoom)) {
          findSeededTracks(*(seeds_vec[zoomlevel + 1]),
                           *(hits_vec[zoomlevel + 1]), tracks,
                           zoomranges[zoomlevel + 1]);
        } else {
          findSeededHelices(min_hits, max_hits, tracks, maxtracks,
                            zoomlevel + 1);
        }
        if (maxtracks != 0) {
          double phi_proportion =
              (zoomranges[zoomlevel].max_phi - zoomranges[zoomlevel].min_phi) /
              ((zoomranges[0].max_phi - zoomranges[0].min_phi));
          double d_proportion =
              (zoomranges[zoomlevel].max_d - zoomranges[zoomlevel].min_d) /
              ((zoomranges[0].max_d - zoomranges[0].min_d));
          double k_proportion =
              (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k) /
              ((zoomranges[0].max_k - zoomranges[0].min_k));
          unsigned int expected_tracks =
              (unsigned int)(fabs(((double)maxtracks) * phi_proportion *
                                  d_proportion * k_proportion)) +
              1;

          if ((tracks.size() - tracks_at_start) > expected_tracks) {
            return;
          }
        }
      }
      if (count == n_entries) {
        break;
      }
      hits_vec[zoomlevel + 1]->clear();
      seeds_vec[zoomlevel + 1]->clear();
      setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel],
               zoomranges[zoomlevel + 1], n_phi, n_d, n_k, n_dzdl, n_z0);
    }
  }
}
