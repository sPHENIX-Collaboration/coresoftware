#include "HelixHough.h"

#include "fastvec.h"
#include "SimpleHit3D.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace std;

static void counting_sort(unsigned int MAX, std::vector<BinEntryPair5D>& A) {
  unsigned int SIZE = A.size();
  std::vector<BinEntryPair5D> B(SIZE);
  std::vector<unsigned int> C(MAX + 1, 0);

  for (unsigned int i = 0; i < SIZE; ++i) {
    C[A[i].bin] += 1;
  }

  for (unsigned int i = 1; i <= MAX; ++i) {
    C[i] += C[i - 1];
  }

  for (long i = SIZE - 1; i >= 0; --i) {
    B[C[A[i].bin] - 1] = A[i];
    --C[A[i].bin];
  }

  for (unsigned int i = 0; i < SIZE; ++i) {
    A[i] = B[i];
  }
}

static void counting_sort(unsigned int SIZE, unsigned int MAX, unsigned int* A,
                          unsigned int* B, unsigned int* C) {
  for (unsigned int i = 0; i < SIZE; ++i) {
    C[(A[i]) >> 20] += 1;
  }
  for (unsigned int i = 1; i < MAX; ++i) {
    C[i] += C[i - 1];
  }
  for (int i = (SIZE - 1); i >= 0; --i) {
    B[C[A[i] >> 20] - 1] = A[i];
    --C[A[i] >> 20];
  }
  for (unsigned int i = 0; i < SIZE; ++i) {
    A[i] = B[i];
  }
}

void HelixHough::vote_pairs(unsigned int zoomlevel) {
  bins_vec[zoomlevel]->clear();
  fastvec vote_array;

  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];

  unsigned int total_bins = n_phi * n_d * n_k * n_dzdl * n_z0;

  float low_phi = zoomranges[zoomlevel].min_phi;
  float high_phi = zoomranges[zoomlevel].max_phi;
  float low_d = zoomranges[zoomlevel].min_d;
  float high_d = zoomranges[zoomlevel].max_d;
  float low_z0 = zoomranges[zoomlevel].min_d;
  float high_z0 = zoomranges[zoomlevel].max_d;
  float low_dzdl = zoomranges[zoomlevel].min_dzdl;
  float high_dzdl = zoomranges[zoomlevel].max_dzdl;
  float inv_phi_range = 1. / (high_phi - low_phi);
  float inv_d_range = 1. / (high_d - low_d);
  float inv_z0_range = 1. / (high_z0 - low_z0);
  float inv_dzdl_range = 1. / (high_dzdl - low_dzdl);
  float k_size = (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k) /
                 ((float)n_k);

  unsigned int pair_counter = 0;
  vector<vector<SimpleHit3D> > four_pairs;
  vector<SimpleHit3D> onepair;
  unsigned int pair_index[4];
  onepair.assign(2, SimpleHit3D());
  four_pairs.assign(4, onepair);
  float x1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float z1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float x2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float z2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_phi_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_phi_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_phi_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_phi_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_d_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_d_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_d_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_d_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_dzdl_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_dzdl_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_z0_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_z0_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_z0_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_z0_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_k = zoomranges[zoomlevel].min_k;
  float max_k = min_k + k_size;
  float error_scale = 2.;
  for (unsigned int k_bin = 0; k_bin < n_k; ++k_bin) {
    float pwr = 1.0;
    float min_kappa = pow(min_k, pwr);
    float max_kappa = pow(max_k, pwr);
    float avg = 0.5 * (max_kappa + min_kappa);
    float width = 0.5 * (max_kappa - min_kappa);
    max_kappa = avg + width * bin_scale;
    min_kappa = avg - width * bin_scale;

    float min_k_a[4] __attribute__((aligned(16))) = {min_kappa, min_kappa,
                                                     min_kappa, min_kappa};
    float max_k_a[4] __attribute__((aligned(16))) = {max_kappa, max_kappa,
                                                     max_kappa, max_kappa};

    pair_counter = 0;
    for (unsigned int i = 0; i < pairs_vec[zoomlevel]->size(); i++) {
      pair_index[pair_counter] = i;
      four_pairs[pair_counter][0] =
          ((*(hits_vec[zoomlevel]))[(*(pairs_vec[zoomlevel]))[i].first]);
      four_pairs[pair_counter][1] =
          ((*(hits_vec[zoomlevel]))[(*(pairs_vec[zoomlevel]))[i].second]);
      x1_a[pair_counter] = four_pairs[pair_counter][0].get_x();
      y1_a[pair_counter] = four_pairs[pair_counter][0].get_y();
      z1_a[pair_counter] = four_pairs[pair_counter][0].get_z();
      x2_a[pair_counter] = four_pairs[pair_counter][1].get_x();
      y2_a[pair_counter] = four_pairs[pair_counter][1].get_y();
      z2_a[pair_counter] = four_pairs[pair_counter][1].get_z();
      pair_counter += 1;
      if (pair_counter == 4) {
        HelixHough::allButKappaRange_sse(
            x1_a, x2_a, y1_a, y2_a, z1_a, z2_a, min_k_a, max_k_a, min_phi_1_a,
            max_phi_1_a, min_phi_2_a, max_phi_2_a, min_d_1_a, max_d_1_a,
            min_d_2_a, max_d_2_a, min_dzdl_a, max_dzdl_a, min_z0_1_a,
            max_z0_1_a, min_z0_2_a, max_z0_2_a);
        for (unsigned int h = 0; h < pair_counter; ++h) {

          float dr0 = sqrt((2.0*sqrt(four_pairs[h][0].get_size(0,0))) *
			   (2.0*sqrt(four_pairs[h][0].get_size(0,0))) +
                           (2.0*sqrt(four_pairs[h][0].get_size(1,1))) *
			   (2.0*sqrt(four_pairs[h][0].get_size(1,1))));
	  
          float r0 = sqrt(four_pairs[h][0].get_x() * four_pairs[h][0].get_x() +
                          four_pairs[h][0].get_y() * four_pairs[h][0].get_y());
	  
          float dr1 = sqrt((2.0*sqrt(four_pairs[h][1].get_size(0,0))) *
			   (2.0*sqrt(four_pairs[h][1].get_size(0,0))) +
                           (2.0*sqrt(four_pairs[h][1].get_size(1,1))) *
			   (2.0*sqrt(four_pairs[h][1].get_size(1,1))));

          float r1 = sqrt(four_pairs[h][1].get_x() * four_pairs[h][1].get_x() +
                          four_pairs[h][1].get_y() * four_pairs[h][1].get_y());
          float r1r0_inv = 1. / (r1 - r0);

          // phi error from hit error
          float dphi = (dr0 + dr1) * r1r0_inv;
          // phi error from m.s. or whatever else
          float phi_scatt = phiError(
              four_pairs[h][1], min_kappa, max_kappa,
              zoomranges[zoomlevel].min_d, zoomranges[zoomlevel].max_d,
              zoomranges[zoomlevel].min_z0, zoomranges[zoomlevel].max_z0,
              zoomranges[zoomlevel].min_dzdl, zoomranges[zoomlevel].max_dzdl,
              true);
          dphi += phi_scatt;
          dphi *= error_scale;
          min_phi_1_a[h] -= dphi;
          min_phi_2_a[h] -= dphi;
          max_phi_1_a[h] += dphi;
          max_phi_2_a[h] += dphi;

          // d error from hit error
          float dd = r0 * (dr0 + dr1) * r1r0_inv + dr0;
          dd += phi_scatt * r1;
          dd *= error_scale;
          min_d_1_a[h] -= dd;
          min_d_2_a[h] -= dd;
          max_d_1_a[h] += dd;
          max_d_2_a[h] += dd;

          // dzdl error from hit error
          float ddzdl = ((2.0*sqrt(four_pairs[h][0].get_size(2,2))) +
			 (2.0*sqrt(four_pairs[h][1].get_size(2,2)))) * r1r0_inv;
          // dzdl error from m.s. or whatever else
          float dzdl_scatt = dzdlError(
              four_pairs[h][1], min_kappa, max_kappa,
              zoomranges[zoomlevel].min_d, zoomranges[zoomlevel].max_d,
              zoomranges[zoomlevel].min_z0, zoomranges[zoomlevel].max_z0,
              zoomranges[zoomlevel].min_dzdl, zoomranges[zoomlevel].max_dzdl,
              true);
          ddzdl += dzdl_scatt;
          ddzdl *= error_scale;
          min_dzdl_a[h] -= ddzdl;
          max_dzdl_a[h] += ddzdl;

          // z0 error from hit error
          float dz0 =
              r0 * ((2.0*sqrt(four_pairs[h][0].get_size(2,2))) +
		    (2.0*sqrt(four_pairs[h][1].get_size(2,2)))) * r1r0_inv;
          dz0 += dzdl_scatt * r1;
          dz0 *= error_scale;
          min_z0_1_a[h] -= dz0;
          min_z0_2_a[h] -= dz0;
          max_z0_1_a[h] += dz0;
          max_z0_2_a[h] += dz0;
        }
        if (separate_by_helicity == false) {
          fillBins(total_bins, pair_counter, pair_index, min_phi_1_a,
                   max_phi_1_a, min_d_1_a, max_d_1_a, min_dzdl_a, max_dzdl_a,
                   min_z0_1_a, max_z0_1_a, four_pairs, n_d, n_k, n_dzdl, n_z0,
                   k_bin, n_phi, zoomlevel, low_phi, high_phi, low_d, high_d,
                   low_z0, high_z0, low_dzdl, high_dzdl, inv_phi_range,
                   inv_d_range, inv_z0_range, inv_dzdl_range, vote_array);
          fillBins(total_bins, pair_counter, pair_index, min_phi_2_a,
                   max_phi_2_a, min_d_2_a, max_d_2_a, min_dzdl_a, max_dzdl_a,
                   min_z0_2_a, max_z0_2_a, four_pairs, n_d, n_k, n_dzdl, n_z0,
                   k_bin, n_phi, zoomlevel, low_phi, high_phi, low_d, high_d,
                   low_z0, high_z0, low_dzdl, high_dzdl, inv_phi_range,
                   inv_d_range, inv_z0_range, inv_dzdl_range, vote_array);
        } else {
          if (helicity == false) {
            fillBins(total_bins, pair_counter, pair_index, min_phi_1_a,
                     max_phi_1_a, min_d_1_a, max_d_1_a, min_dzdl_a, max_dzdl_a,
                     min_z0_1_a, max_z0_1_a, four_pairs, n_d, n_k, n_dzdl, n_z0,
                     k_bin, n_phi, zoomlevel, low_phi, high_phi, low_d, high_d,
                     low_z0, high_z0, low_dzdl, high_dzdl, inv_phi_range,
                     inv_d_range, inv_z0_range, inv_dzdl_range, vote_array);
          } else {
            fillBins(total_bins, pair_counter, pair_index, min_phi_2_a,
                     max_phi_2_a, min_d_2_a, max_d_2_a, min_dzdl_a, max_dzdl_a,
                     min_z0_2_a, max_z0_2_a, four_pairs, n_d, n_k, n_dzdl, n_z0,
                     k_bin, n_phi, zoomlevel, low_phi, high_phi, low_d, high_d,
                     low_z0, high_z0, low_dzdl, high_dzdl, inv_phi_range,
                     inv_d_range, inv_z0_range, inv_dzdl_range, vote_array);
          }
        }
        pair_counter = 0;
      }
    }
    if (pair_counter != 0) {
      HelixHough::allButKappaRange_sse(
          x1_a, x2_a, y1_a, y2_a, z1_a, z2_a, min_k_a, max_k_a, min_phi_1_a,
          max_phi_1_a, min_phi_2_a, max_phi_2_a, min_d_1_a, max_d_1_a,
          min_d_2_a, max_d_2_a, min_dzdl_a, max_dzdl_a, min_z0_1_a, max_z0_1_a,
          min_z0_2_a, max_z0_2_a);
      for (unsigned int h = 0; h < pair_counter; ++h) {
	
	float dr0 = sqrt((2.0*sqrt(four_pairs[h][0].get_size(0,0))) *
			 (2.0*sqrt(four_pairs[h][0].get_size(0,0))) +
			 (2.0*sqrt(four_pairs[h][0].get_size(1,1))) *
			 (2.0*sqrt(four_pairs[h][0].get_size(1,1))));

        float r0 = sqrt(four_pairs[h][0].get_x() * four_pairs[h][0].get_x() +
                        four_pairs[h][0].get_y() * four_pairs[h][0].get_y());

	float dr1 = sqrt((2.0*sqrt(four_pairs[h][1].get_size(0,0))) *
			 (2.0*sqrt(four_pairs[h][1].get_size(0,0))) +
			 (2.0*sqrt(four_pairs[h][1].get_size(1,1))) *
			 (2.0*sqrt(four_pairs[h][1].get_size(1,1))));

        float r1 = sqrt(four_pairs[h][1].get_x() * four_pairs[h][1].get_x() +
                        four_pairs[h][1].get_y() * four_pairs[h][1].get_y());
        float r1r0_inv = 1. / (r1 - r0);

        // phi error from hit error
        float dphi = (dr0 + dr1) * r1r0_inv;
        // phi error from m.s. or whatever else
        float phi_scatt = phiError(
            four_pairs[h][1], min_kappa, max_kappa, zoomranges[zoomlevel].min_d,
            zoomranges[zoomlevel].max_d, zoomranges[zoomlevel].min_z0,
            zoomranges[zoomlevel].max_z0, zoomranges[zoomlevel].min_dzdl,
            zoomranges[zoomlevel].max_dzdl, true);
        dphi += phi_scatt;
        dphi *= error_scale;
        min_phi_1_a[h] -= dphi;
        min_phi_2_a[h] -= dphi;
        max_phi_1_a[h] += dphi;
        max_phi_2_a[h] += dphi;

        // d error from hit error
        float dd = r0 * (dr0 + dr1) * r1r0_inv + dr0;
        dd += phi_scatt * r1;
        dd *= error_scale;
        min_d_1_a[h] -= dd;
        min_d_2_a[h] -= dd;
        max_d_1_a[h] += dd;
        max_d_2_a[h] += dd;

        // dzdl error from hit error
        float ddzdl = ((2.0*sqrt(four_pairs[h][0].get_size(2,2))) +
		       (2.0*sqrt(four_pairs[h][1].get_size(2,2)))) * r1r0_inv;
        // dzdl error from m.s. or whatever else
        float dzdl_scatt = dzdlError(
            four_pairs[h][1], min_kappa, max_kappa, zoomranges[zoomlevel].min_d,
            zoomranges[zoomlevel].max_d, zoomranges[zoomlevel].min_z0,
            zoomranges[zoomlevel].max_z0, zoomranges[zoomlevel].min_dzdl,
            zoomranges[zoomlevel].max_dzdl, true);
        ddzdl += dzdl_scatt;
        ddzdl *= error_scale;
        min_dzdl_a[h] -= ddzdl;
        max_dzdl_a[h] += ddzdl;

        // z0 error from hit error
        float dz0 = r0 * ((2.0*sqrt(four_pairs[h][0].get_size(2,2))) +
			  (2.0*sqrt(four_pairs[h][1].get_size(2,2)))) * r1r0_inv;
        dz0 += dzdl_scatt * r1;
        dz0 *= error_scale;
        min_z0_1_a[h] -= dz0;
        min_z0_2_a[h] -= dz0;
        max_z0_1_a[h] += dz0;
        max_z0_2_a[h] += dz0;
      }
      if (separate_by_helicity == false) {
        fillBins(total_bins, pair_counter, pair_index, min_phi_1_a, max_phi_1_a,
                 min_d_1_a, max_d_1_a, min_dzdl_a, max_dzdl_a, min_z0_1_a,
                 max_z0_1_a, four_pairs, n_d, n_k, n_dzdl, n_z0, k_bin, n_phi,
                 zoomlevel, low_phi, high_phi, low_d, high_d, low_z0, high_z0,
                 low_dzdl, high_dzdl, inv_phi_range, inv_d_range, inv_z0_range,
                 inv_dzdl_range, vote_array);
        fillBins(total_bins, pair_counter, pair_index, min_phi_2_a, max_phi_2_a,
                 min_d_2_a, max_d_2_a, min_dzdl_a, max_dzdl_a, min_z0_2_a,
                 max_z0_2_a, four_pairs, n_d, n_k, n_dzdl, n_z0, k_bin, n_phi,
                 zoomlevel, low_phi, high_phi, low_d, high_d, low_z0, high_z0,
                 low_dzdl, high_dzdl, inv_phi_range, inv_d_range, inv_z0_range,
                 inv_dzdl_range, vote_array);
      } else {
        if (helicity == false) {
          fillBins(total_bins, pair_counter, pair_index, min_phi_1_a,
                   max_phi_1_a, min_d_1_a, max_d_1_a, min_dzdl_a, max_dzdl_a,
                   min_z0_1_a, max_z0_1_a, four_pairs, n_d, n_k, n_dzdl, n_z0,
                   k_bin, n_phi, zoomlevel, low_phi, high_phi, low_d, high_d,
                   low_z0, high_z0, low_dzdl, high_dzdl, inv_phi_range,
                   inv_d_range, inv_z0_range, inv_dzdl_range, vote_array);
        } else {
          fillBins(total_bins, pair_counter, pair_index, min_phi_2_a,
                   max_phi_2_a, min_d_2_a, max_d_2_a, min_dzdl_a, max_dzdl_a,
                   min_z0_2_a, max_z0_2_a, four_pairs, n_d, n_k, n_dzdl, n_z0,
                   k_bin, n_phi, zoomlevel, low_phi, high_phi, low_d, high_d,
                   low_z0, high_z0, low_dzdl, high_dzdl, inv_phi_range,
                   inv_d_range, inv_z0_range, inv_dzdl_range, vote_array);
        }
      }
      pair_counter = 0;
    }
    min_k += k_size;
    max_k += k_size;
  }
  if (vote_array.size == 0) {
    return;
  }
  if (vote_array.size < (1 << 14)) {
    unsigned int B[1 << 14];
    unsigned int C[(1 << 12) + 1] = {0};
    counting_sort(vote_array.size, total_bins, vote_array.arr, B, C);
    bins_vec[zoomlevel]->resize(vote_array.size, BinEntryPair5D(0, 0));
    for (unsigned int i = 0; i < vote_array.size; ++i) {
      (*(bins_vec[zoomlevel]))[i].bin = (vote_array.arr[i]) >> 20;
      (*(bins_vec[zoomlevel]))[i].entry =
          (vote_array.arr[i]) & ((unsigned int)1048575);
      (*(bins_vec[zoomlevel]))[i].is_seed = false;
    }
  } else {
    bins_vec[zoomlevel]->resize(vote_array.size, BinEntryPair5D(0, 0));
    for (unsigned int i = 0; i < vote_array.size; ++i) {
      (*(bins_vec[zoomlevel]))[i].bin = (vote_array[i]) >> 20;
      (*(bins_vec[zoomlevel]))[i].entry =
          (vote_array[i]) & ((unsigned int)1048575);
      (*(bins_vec[zoomlevel]))[i].is_seed = false;
    }
    if (bins_vec[zoomlevel]->size() > total_bins) {
      counting_sort(total_bins, *(bins_vec[zoomlevel]));
    } else {
      sort(bins_vec[zoomlevel]->begin(), bins_vec[zoomlevel]->end());
    }
  }
}

void HelixHough::fillBins(
    unsigned int total_bins, unsigned int pair_counter,
    unsigned int* pair_index, float* min_phi, float* max_phi, float* min_d,
    float* max_d, float* min_dzdl, float* max_dzdl, float* min_z0,
    float* max_z0, vector<vector<SimpleHit3D> >& four_pairs, unsigned int n_d,
    unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0,
    unsigned int k_bin, unsigned int n_phi, unsigned int zoomlevel,
    float low_phi, float high_phi, float low_d, float high_d, float low_z0,
    float high_z0, float low_dzdl, float high_dzdl, float inv_phi_range,
    float inv_d_range, float inv_z0_range, float inv_dzdl_range,
    fastvec& vote_array) {
  unsigned int buffer[(1 << 10)];
  unsigned int bufnum = 0;

  for (unsigned int i = 0; i < pair_counter; ++i) {
    if ((max_d[i] < low_d) || (min_d[i] > high_d)) {
      continue;
    }
    if ((max_z0[i] < low_z0) || (min_z0[i] > high_z0)) {
      continue;
    }
    if ((max_dzdl[i] < low_dzdl) || (min_dzdl[i] > high_dzdl)) {
      continue;
    }

    unsigned int index = pair_index[i];

    unsigned int low_d_bin = 0;
    unsigned int high_d_bin = (n_d - 1);
    unsigned int low_z0_bin = 0;
    unsigned int high_z0_bin = (n_z0 - 1);
    unsigned int low_dzdl_bin = 0;
    unsigned int high_dzdl_bin = (n_dzdl - 1);
    if (min_d[i] > low_d) {
      low_d_bin =
          (unsigned int)(((min_d[i] - low_d) * inv_d_range) * ((float)(n_d)));
    }
    if (max_d[i] < high_d) {
      high_d_bin =
          (unsigned int)(((max_d[i] - low_d) * inv_d_range) * ((float)(n_d)));
    }
    if (min_z0[i] > low_z0) {
      low_z0_bin = (unsigned int)(((min_z0[i] - low_z0) * inv_z0_range) *
                                  ((float)(n_z0)));
    }
    if (max_z0[i] < high_z0) {
      high_z0_bin = (unsigned int)(((max_z0[i] - low_z0) * inv_z0_range) *
                                   ((float)(n_z0)));
    }
    if (min_dzdl[i] > low_dzdl) {
      low_dzdl_bin =
          (unsigned int)(((min_dzdl[i] - low_dzdl) * inv_dzdl_range) *
                         ((float)(n_dzdl)));
    }
    if (max_dzdl[i] < high_dzdl) {
      high_dzdl_bin =
          (unsigned int)(((max_dzdl[i] - low_dzdl) * inv_dzdl_range) *
                         ((float)(n_dzdl)));
    }

    if (min_phi[i] >= 0.) {
      if ((max_phi[i] < low_phi) || (min_phi[i] > high_phi)) {
        continue;
      }
      unsigned int low_phi_bin = 0;
      unsigned int high_phi_bin = (n_phi - 1);
      if (min_phi[i] > low_phi) {
        low_phi_bin = (unsigned int)(((min_phi[i] - low_phi) * inv_phi_range) *
                                     ((float)(n_phi)));
      }
      if (max_phi[i] < high_phi) {
        high_phi_bin = (unsigned int)(((max_phi[i] - low_phi) * inv_phi_range) *
                                      ((float)(n_phi)));
      }
      for (unsigned int phib = low_phi_bin; phib <= high_phi_bin; ++phib) {
        for (unsigned int db = low_d_bin; db <= high_d_bin; ++db) {
          for (unsigned int z0b = low_z0_bin; z0b <= high_z0_bin; ++z0b) {
            for (unsigned int dzdlb = low_dzdl_bin; dzdlb <= high_dzdl_bin;
                 ++dzdlb) {
              unsigned int bin = BinEntryPair5D::linearBin(
                  n_d, n_k, n_dzdl, n_z0, phib, db, k_bin, dzdlb, z0b);
              buffer[bufnum] = ((bin << 20) ^ index);
              bufnum += 1;
            }
          }
        }
      }
    }
  }
  vote_array.push_back(buffer, bufnum);
}
