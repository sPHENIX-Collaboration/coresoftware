#include "HelixHough.h"

#include "fastvec.h"
#include "SimpleHit3D.h"
#include "vector_math_inline.h"

#include <cmath>
#include <cstddef>
#include <emmintrin.h>
#include <memory>
#include <sys/time.h>
#include <vector>
#include <xmmintrin.h>

using namespace std;

static const unsigned int digits = 2;
static const unsigned int r = 7;
static const unsigned int radix = 1 << r;
static const unsigned int mask = radix - 1;

static void radix_sort(std::vector<BinEntryPair5D>& A) {
  unsigned int SIZE = A.size();
  std::vector<BinEntryPair5D> B(SIZE);
  std::vector<unsigned int> cnt(radix);

  for (unsigned int i = 0, shift = 0; i < digits; i++, shift += r) {
    for (unsigned int j = 0; j < radix; ++j) {
      cnt[j] = 0;
    }

    for (unsigned int j = 0; j < SIZE; ++j) {
      ++cnt[(A[j].bin >> shift) & mask];
    }

    for (unsigned int j = 1; j < radix; ++j) {
      cnt[j] += cnt[j - 1];
    }

    for (long j = SIZE - 1; j >= 0; --j) {
      B[--cnt[(A[j].bin >> shift) & mask]] = A[j];
    }

    for (unsigned int j = 0; j < SIZE; ++j) {
      A[j] = B[j];
    }
  }
}

static void radix_sort(unsigned int* A, unsigned int* B, unsigned int SIZE) {
  unsigned int cnt[1 << 7];

  for (unsigned int i = 0, shift = 0; i < digits; i++, shift += r) {
    for (unsigned int j = 0; j < radix; ++j) {
      cnt[j] = 0;
    }

    for (unsigned int j = 0; j < SIZE; ++j) {
      ++cnt[((A[j] >> 20) >> shift) & mask];
    }

    for (unsigned int j = 1; j < radix; ++j) {
      cnt[j] += cnt[j - 1];
    }

    for (long j = SIZE - 1; j >= 0; --j) {
      B[--cnt[((A[j] >> 20) >> shift) & mask]] = A[j];
    }

    for (unsigned int j = 0; j < SIZE; ++j) {
      A[j] = B[j];
    }
  }
}

static void counting_sort(unsigned int MAX, std::vector<BinEntryPair5D>& A,
                          unsigned int* C, unsigned int* bins_start,
                          unsigned int* bins_end) {
  unsigned int SIZE = A.size();
  std::vector<BinEntryPair5D> B(SIZE);

  for (unsigned int i = 0; i < SIZE; ++i) {
    C[A[i].bin] += 1;
  }

  bins_start[0] = 0;
  bins_end[0] = C[0] - 1;
  for (unsigned int i = 1; i < MAX; ++i) {
    bins_start[i] = C[i - 1];
    C[i] += C[i - 1];
    bins_end[i] = C[i] - 1;
  }

  for (int i = SIZE - 1; i >= 0; --i) {
    B[C[A[i].bin] - 1] = A[i];
    --C[A[i].bin];
  }

  for (unsigned int i = 0; i < SIZE; ++i) {
    A[i] = B[i];
  }
}

static void counting_sort(unsigned int SIZE, unsigned int MAX, unsigned int* A,
                          unsigned int* B, unsigned int* C,
                          unsigned int* bins_start, unsigned int* bins_end) {
  for (unsigned int i = 0; i < SIZE; ++i) {
    C[(A[i]) >> 20] += 1;
  }
  bins_start[0] = 0;
  bins_end[0] = C[0] - 1;
  for (unsigned int i = 1; i < MAX; ++i) {
    bins_start[i] = C[i - 1];
    C[i] += C[i - 1];
    bins_end[i] = C[i] - 1;
  }
  for (int i = (SIZE - 1); i >= 0; --i) {
    B[C[A[i] >> 20] - 1] = A[i];
    --C[A[i] >> 20];
  }
  for (unsigned int i = 0; i < SIZE; ++i) {
    A[i] = B[i];
  }
}

void fillBins4_sse(float* min_phi_a, float* max_phi_a, float n_phi_val,
                   float inv_phi_range_val, float low_phi_a, float high_phi_a,
                   unsigned int* philow1, unsigned int* philow2,
                   unsigned int* phihi1, unsigned int* phihi2) {
  const unsigned int zero_int_a[4] __attribute__((aligned(16))) = {
      0x00000000, 0x00000000, 0x00000000, 0x00000000};
  __m128i zero_int = _mm_load_si128((__m128i*)zero_int_a);
  const unsigned int neg1_int_a[4] __attribute__((aligned(16))) = {
      0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF};
  __m128i neg1_int = _mm_load_si128((__m128i*)neg1_int_a);
  const unsigned int one_int_a[4] __attribute__((aligned(16))) = {
      0x00000001, 0x00000001, 0x00000001, 0x00000001};
  __m128i one_int = _mm_load_si128((__m128i*)one_int_a);

  __m128 min_phi = _mm_load_ps(min_phi_a);
  __m128 max_phi = _mm_load_ps(max_phi_a);
  __m128 n_phi = _mm_load1_ps(&(n_phi_val));
  __m128i n_phi_min1 = _mm_cvtps_epi32(n_phi);
  n_phi_min1 = _mm_sub_epi32(n_phi_min1, one_int);
  __m128 inv_phi_range = _mm_load1_ps(&(inv_phi_range_val));
  __m128 low_phi = _mm_load1_ps(&(low_phi_a));
  __m128 high_phi = _mm_load1_ps(&(high_phi_a));

  __m128i low_phi_bin_1 =
      _mm_cvttps_epi32((((min_phi - low_phi) * inv_phi_range) * ((n_phi))));
  __m128i high_phi_bin_1 =
      _mm_cvttps_epi32((((max_phi - low_phi) * inv_phi_range) * ((n_phi))));
  __m128i low_phi_bin_2 = _mm_cvttps_epi32(
      (((min_phi + twopi - low_phi) * inv_phi_range) * ((n_phi))));

  __m128i cmp1 = (__m128i)(_mm_cmplt_ps(min_phi, low_phi));
  __m128i cmp2 = (__m128i)(_mm_cmplt_ps(high_phi, max_phi));
  __m128i cmp3 = (__m128i)(_mm_cmplt_ps(min_phi + twopi, low_phi));
  __m128i cmp4 = (__m128i)(_mm_cmplt_ps(zero, min_phi));
  //   ( high_phi < ( min_phi_a[i] + 2.*M_PI) ) && ( low_phi > max_phi_a[i] )
  __m128i cmp5 = (__m128i)(_mm_cmplt_ps(high_phi, min_phi + twopi));
  __m128i cmptmp1 = (__m128i)(_mm_cmplt_ps(max_phi, low_phi));
  cmp5 = _mm_and_si128(cmp5, cmptmp1);
  // ( high_phi >= ( min_phi_a[i] + 2.*M_PI) ) && ( low_phi > max_phi_a[i] )
  __m128i cmp5_1 = (__m128i)(_mm_cmplt_ps(min_phi + twopi, high_phi));
  cmptmp1 = (__m128i)(_mm_cmplt_ps(max_phi, low_phi));
  cmp5_1 = _mm_and_si128(cmp5_1, cmptmp1);
  // ( high_phi < ( min_phi_a[i] + 2.*M_PI) ) && ( low_phi <= max_phi_a[i] )
  __m128i cmp5_2 = (__m128i)(_mm_cmplt_ps(high_phi, min_phi + twopi));
  cmptmp1 = (__m128i)(_mm_cmplt_ps(low_phi, max_phi));
  cmp5_2 = _mm_and_si128(cmp5_2, cmptmp1);
  // ( high_phi >= ( min_phi_a[i] + 2.*M_PI) ) && ( low_phi <= max_phi_a[i] )
  __m128i cmp5_3 = (__m128i)(_mm_cmplt_ps(min_phi + twopi, high_phi));
  cmptmp1 = (__m128i)(_mm_cmplt_ps(low_phi, max_phi));
  cmp5_3 = _mm_and_si128(cmp5_3, cmptmp1);
  cmp5_1 = _mm_or_si128(cmp5_1, cmp5_3);
  cmp5_2 = _mm_or_si128(cmp5_2, cmp5_3);
  // ( max_phi_a[i] < low_phi ) || ( min_phi_a[i] > high_phi )
  __m128i cmp6 = (__m128i)(_mm_cmplt_ps(max_phi, low_phi));
  cmptmp1 = (__m128i)(_mm_cmplt_ps(high_phi, min_phi));
  cmp6 = _mm_or_si128(cmp6, cmptmp1);

  // split into two cases due to 0,2pi wraparound

  // low bin :

  // choose between zero and low_phi_bin_2
  __m128i tmp1 = _mm_and_si128(cmp3, zero_int);
  __m128i tmp2 = _mm_andnot_si128(cmp3, low_phi_bin_2);
  __m128i lowphi_sel_1 = _mm_xor_si128(tmp1, tmp2);
  // set result to negative based on cmp5
  tmp1 = _mm_andnot_si128(cmp5_1, neg1_int);
  tmp2 = _mm_and_si128(cmp5_1, lowphi_sel_1);
  lowphi_sel_1 = _mm_xor_si128(tmp1, tmp2);
  // choose between above result and ((low_phi_bin_1 or 0), or neg if cmp6)
  tmp1 = _mm_and_si128(cmp1, zero_int);
  tmp2 = _mm_andnot_si128(cmp1, low_phi_bin_1);
  __m128i tmp4 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp6, neg1_int);
  tmp2 = _mm_andnot_si128(cmp6, tmp4);
  __m128i tmp3 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp4, tmp3);
  tmp2 = _mm_andnot_si128(cmp4, lowphi_sel_1);
  lowphi_sel_1 = _mm_xor_si128(tmp1, tmp2);

  // high bin :
  tmp1 = _mm_and_si128(cmp2, n_phi_min1);
  tmp2 = _mm_andnot_si128(cmp2, high_phi_bin_1);
  tmp4 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp6, neg1_int);
  tmp2 = _mm_andnot_si128(cmp6, tmp4);
  tmp4 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp5, neg1_int);
  tmp2 = _mm_andnot_si128(cmp5, n_phi_min1);
  tmp3 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp4, tmp4);
  tmp2 = _mm_andnot_si128(cmp4, tmp3);
  __m128i highphi_sel_1 = _mm_xor_si128(tmp1, tmp2);

  tmp1 = _mm_and_si128(cmp5_2, zero_int);
  tmp2 = _mm_andnot_si128(cmp5_2, neg1_int);
  tmp4 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp4, neg1_int);
  tmp2 = _mm_andnot_si128(cmp4, tmp4);
  __m128i lowphi_sel_2 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp2, n_phi_min1);
  tmp2 = _mm_andnot_si128(cmp2, high_phi_bin_1);
  tmp4 = _mm_xor_si128(tmp1, tmp2);
  tmp1 = _mm_and_si128(cmp4, neg1_int);
  tmp2 = _mm_andnot_si128(cmp4, tmp4);
  __m128i highphi_sel_2 = _mm_xor_si128(tmp1, tmp2);

  _mm_store_si128((__m128i*)philow1, lowphi_sel_1);
  _mm_store_si128((__m128i*)philow2, lowphi_sel_2);
  _mm_store_si128((__m128i*)phihi1, highphi_sel_1);
  _mm_store_si128((__m128i*)phihi2, highphi_sel_2);
}

void HelixHough::fillBins(unsigned int total_bins, unsigned int hit_counter,
                          float* min_phi_a, float* max_phi_a,
                          vector<SimpleHit3D>& four_hits, fastvec2d& z_bins,
                          unsigned int n_d, unsigned int n_k,
                          unsigned int n_dzdl, unsigned int n_z0,
                          unsigned int d_bin, unsigned int k_bin,
                          unsigned int n_phi, unsigned int zoomlevel,
                          float low_phi, float high_phi, float inv_phi_range,
                          fastvec& vote_array) {
  unsigned int buffer[(1 << 10)];
  unsigned int bufnum = 0;
  unsigned int zbuffer[1 << 10];
  unsigned int zbufnum[8];
  unsigned int size2 = n_z0 * n_dzdl;

  unsigned int philow1[4] __attribute__((aligned(16))) = {
      0x00000000, 0x00000000, 0x00000000, 0x00000000};
  unsigned int philow2[4] __attribute__((aligned(16))) = {
      0x00000000, 0x00000000, 0x00000000, 0x00000000};
  unsigned int phihi1[4] __attribute__((aligned(16))) = {
      0x00000000, 0x00000000, 0x00000000, 0x00000000};
  unsigned int phihi2[4] __attribute__((aligned(16))) = {
      0x00000000, 0x00000000, 0x00000000, 0x00000000};

  z_bins.fetch(four_hits[0].get_id(), four_hits[hit_counter - 1].get_id(), zbuffer,
               zbufnum);

  unsigned int zoff = n_z0 * n_dzdl * (k_bin + n_k * d_bin);
  unsigned int binprod = n_z0 * n_dzdl * n_k * n_d;

  unsigned int count = hit_counter;
  unsigned int offset = 0;
  unsigned int cur = 4;
  if (count < 4) {
    cur = count;
  }
  while (true) {
    float minphi_a[4] __attribute__((aligned(16)));
    float maxphi_a[4] __attribute__((aligned(16)));
    for (unsigned int i = 0; i < cur; ++i) {
      minphi_a[i] = min_phi_a[i + offset];
      maxphi_a[i] = max_phi_a[i + offset];
    }

    fillBins4_sse(minphi_a, maxphi_a, (float)n_phi, inv_phi_range, low_phi,
                  high_phi, philow1, philow2, phihi1, phihi2);

    for (unsigned int i = 0; i < cur; ++i) {
      unsigned int index = four_hits[i + offset].get_id();

      unsigned int pos = (i + offset) * size2;
      for (unsigned int zbin = 0; zbin < zbufnum[i + offset]; ++zbin) {
        unsigned int zint = zbuffer[pos];
        pos += 1;
        unsigned int zstart = (zint & 65535) + n_z0 * (zint >> 16) + zoff;

        for (unsigned int b = philow1[i]; b <= phihi1[i]; ++b) {
          if (philow1[i] == 4294967295) {
            break;
          }
          if (b >= n_phi) {
            break;
          }
          unsigned int zadd = binprod * b;
          unsigned int bin = zstart + zadd;
          buffer[bufnum] = ((bin << 20) ^ index);
          bufnum += 1;
        }
        for (unsigned int b = philow2[i]; b <= phihi2[i]; ++b) {
          if (philow2[i] == 4294967295) {
            break;
          }
          if (b >= n_phi) {
            break;
          }
          unsigned int zadd = binprod * b;
          unsigned int bin = zstart + zadd;
          buffer[bufnum] = ((bin << 20) ^ index);
          bufnum += 1;
        }
      }
    }
    if (count == cur) {
      break;
    }
    count -= cur;
    offset += cur;
    cur = 4;
    if (count < 4) {
      cur = count;
    }
  }
  vote_array.push_back(buffer, bufnum);
}

void HelixHough::vote_z(unsigned int zoomlevel, unsigned int n_phi,
                        unsigned int n_d, unsigned int n_k, unsigned int n_dzdl,
                        unsigned int n_z0, fastvec2d& z_bins) {
  float z0_size =
      (zoomranges[zoomlevel].max_z0 - zoomranges[zoomlevel].min_z0) /
      ((float)n_z0);
  float dzdl_size =
      (zoomranges[zoomlevel].max_dzdl - zoomranges[zoomlevel].min_dzdl) /
      ((float)n_dzdl);
  float low_phi = zoomranges[zoomlevel].min_phi;
  float high_phi = zoomranges[zoomlevel].max_phi;
  float low_dzdl = zoomranges[zoomlevel].min_dzdl;
  float high_dzdl = zoomranges[zoomlevel].max_dzdl;
  float dzdl_bin_size_inv = 1. / dzdl_size;

  // cache cosine and sine calculations
  float min_cos = cos(zoomranges[zoomlevel].min_phi);
  float max_cos = cos(zoomranges[zoomlevel].max_phi);
  float min_sin = sin(zoomranges[zoomlevel].min_phi);
  float max_sin = sin(zoomranges[zoomlevel].max_phi);
  if ((high_phi - low_phi) > M_PI) {
    min_cos = 1.;
    max_cos = -1.;
    min_sin = 0.;
    max_sin = 0.;
  }

  unsigned int one_z_bin;
  float pwr = 1.0;
  float min_kappa = pow(zoomranges[zoomlevel].min_k, pwr);
  float max_kappa = pow(zoomranges[zoomlevel].max_k, pwr);

  unsigned int hit_counter = 0;
  float x_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float z_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_dzdl_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_dzdl_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_z0_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_z0_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dz_a[4] = {0., 0., 0., 0.};
  vector<SimpleHit3D> four_hits;
  SimpleHit3D temphit;
  four_hits.assign(4, temphit);
  unsigned int temp_zcount[4];
  unsigned buffer[4][1 << 8];
  for (unsigned int i = 0; i < hits_vec[zoomlevel]->size(); i++) {
    x_a[hit_counter] = (*(hits_vec[zoomlevel]))[i].get_x();
    y_a[hit_counter] = (*(hits_vec[zoomlevel]))[i].get_y();
    z_a[hit_counter] = (*(hits_vec[zoomlevel]))[i].get_z();

    dz_a[hit_counter] = (2.0*sqrt((*(hits_vec[zoomlevel]))[i].get_size(2,2)));
    four_hits[hit_counter] = (*(hits_vec[zoomlevel]))[i];

    hit_counter++;

    if (hit_counter == 4) {
      for (unsigned int h = 0; h < hit_counter; ++h) {
        temp_zcount[h] = 0.;
      }

      for (unsigned int zz = 0; zz < n_z0; ++zz) {
        float min_z0 = zoomranges[zoomlevel].min_z0 + ((float)(zz)) * z0_size;
        float max_z0 = min_z0 + z0_size;

        float avg = 0.5 * (min_z0 + max_z0);
        float width = 0.5 * (max_z0 - min_z0);
        max_z0 = avg + width * z_bin_scale;
        min_z0 = avg - width * z_bin_scale;

        for (unsigned int h = 0; h < hit_counter; ++h) {
          float dz = dz_a[h];
          min_z0_a[h] = min_z0 - dz;
          max_z0_a[h] = max_z0 + dz;
        }

        dzdlRange_sse(x_a, y_a, z_a, min_cos, min_sin, max_cos, max_sin,
                      min_kappa, max_kappa, zoomranges[zoomlevel].min_d,
                      zoomranges[zoomlevel].max_d, min_z0_a, max_z0_a,
                      min_dzdl_a, max_dzdl_a);

        unsigned int low_bin = 0;
        unsigned int high_bin = 0;

        for (unsigned int h = 0; h < hit_counter; ++h) {
          float d_dzdl = dzdlError(
              four_hits[h], min_kappa, max_kappa, zoomranges[zoomlevel].min_d,
              zoomranges[zoomlevel].max_d, min_z0, max_z0,
              zoomranges[zoomlevel].min_dzdl, zoomranges[zoomlevel].max_dzdl);

          float min_dzdl = min_dzdl_a[h] - d_dzdl;
          float max_dzdl = max_dzdl_a[h] + d_dzdl;

          if ((min_dzdl > high_dzdl) || (max_dzdl < low_dzdl)) {
            low_bin = 1;
            high_bin = 0;
          } else {
            if (min_dzdl > low_dzdl) {
              low_bin =
                  (unsigned int)((min_dzdl - low_dzdl) * dzdl_bin_size_inv);
            } else {
              low_bin = 0;
            }

            if (max_dzdl < high_dzdl) {
              high_bin =
                  (unsigned int)((max_dzdl - low_dzdl) * dzdl_bin_size_inv);
            } else {
              high_bin = n_dzdl - 1;
            }
          }
          for (unsigned int bb = low_bin; bb <= high_bin; bb++) {
            if (bb >= n_dzdl) {
              continue;
            }
            one_z_bin = (zz ^ (bb << 16));
            buffer[h][temp_zcount[h]] = one_z_bin;
            temp_zcount[h] += 1.;
          }
        }
      }
      for (unsigned int h = 0; h < hit_counter; ++h) {
        z_bins.fill(buffer[h], temp_zcount[h]);
      }
      hit_counter = 0;
    }
  }
  if (hit_counter != 0) {
    for (unsigned int h = 0; h < hit_counter; ++h) {
      temp_zcount[h] = 0.;
    }

    for (unsigned int zz = 0; zz < n_z0; ++zz) {
      float min_z0 = zoomranges[zoomlevel].min_z0 + ((float)(zz)) * z0_size;
      float max_z0 = min_z0 + z0_size;

      float avg = 0.5 * (min_z0 + max_z0);
      float width = 0.5 * (max_z0 - min_z0);
      max_z0 = avg + width * z_bin_scale;
      min_z0 = avg - width * z_bin_scale;

      for (unsigned int h = 0; h < hit_counter; ++h) {
        float dz = dz_a[h];
        min_z0_a[h] = min_z0 - dz;
        max_z0_a[h] = max_z0 + dz;
      }

      dzdlRange_sse(x_a, y_a, z_a, min_cos, min_sin, max_cos, max_sin,
                    min_kappa, max_kappa, zoomranges[zoomlevel].min_d,
                    zoomranges[zoomlevel].max_d, min_z0_a, max_z0_a, min_dzdl_a,
                    max_dzdl_a);

      unsigned int low_bin = 0;
      unsigned int high_bin = 0;

      for (unsigned int h = 0; h < hit_counter; ++h) {
        float d_dzdl = dzdlError(
            four_hits[h], min_kappa, max_kappa, zoomranges[zoomlevel].min_d,
            zoomranges[zoomlevel].max_d, min_z0, max_z0,
            zoomranges[zoomlevel].min_dzdl, zoomranges[zoomlevel].max_dzdl);

        float min_dzdl = min_dzdl_a[h] - d_dzdl;
        float max_dzdl = max_dzdl_a[h] + d_dzdl;

        if ((min_dzdl > high_dzdl) || (max_dzdl < low_dzdl)) {
          low_bin = 1;
          high_bin = 0;
        } else {
          if (min_dzdl > low_dzdl) {
            low_bin = (unsigned int)((min_dzdl - low_dzdl) * dzdl_bin_size_inv);
          } else {
            low_bin = 0;
          }

          if (max_dzdl < high_dzdl) {
            high_bin =
                (unsigned int)((max_dzdl - low_dzdl) * dzdl_bin_size_inv);
          } else {
            high_bin = n_dzdl - 1;
          }
        }
        for (unsigned int bb = low_bin; bb <= high_bin; bb++) {
          if (bb >= n_dzdl) {
            continue;
          }
          one_z_bin = (zz ^ (bb << 16));
          buffer[h][temp_zcount[h]] = one_z_bin;
          temp_zcount[h] += 1.;
        }
      }
    }
    for (unsigned int h = 0; h < hit_counter; ++h) {
      z_bins.fill(buffer[h], temp_zcount[h]);
    }
    hit_counter = 0;
  }
}

void HelixHough::vote(unsigned int zoomlevel) {
  bins_vec[zoomlevel]->clear();
  fastvec vote_array;

  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];

  fastvec2d z_bins(n_dzdl * n_z0);

  unsigned int total_bins = n_phi * n_d * n_k * n_dzdl * n_z0;

  float d_size = (zoomranges[zoomlevel].max_d - zoomranges[zoomlevel].min_d) /
                 ((float)n_d);
  float k_size = (zoomranges[zoomlevel].max_k - zoomranges[zoomlevel].min_k) /
                 ((float)n_k);
  float low_phi = zoomranges[zoomlevel].min_phi;
  float high_phi = zoomranges[zoomlevel].max_phi;
  float inv_phi_range = 1. / (high_phi - low_phi);

  float pwr = 1.0;
  float min_kappa = pow(zoomranges[zoomlevel].min_k, pwr);
  float max_kappa = pow(zoomranges[zoomlevel].max_k, pwr);

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;
  if (print_timings == true) {
    gettimeofday(&t1, NULL);
  }
  vote_z(zoomlevel, n_phi, n_d, n_k, n_dzdl, n_z0, z_bins);
  if (print_timings == true) {
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    z_vote_time += (time2 - time1);
  }

  // now vote in xy
  __m128 phi_3_in;
  __m128 phi_4_in;
  __m128 phi_3_out;
  __m128 phi_4_out;
  __m128 phi_3_in_2;
  __m128 phi_4_in_2;
  __m128 phi_3_out_2;
  __m128 phi_4_out_2;
  unsigned int hit_counter = 0;
  float x_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float x_a_2[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y_a_2[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  vector<SimpleHit3D> four_hits;
  SimpleHit3D temphit;
  four_hits.assign(4, temphit);
  vector<SimpleHit3D> four_hits_2;
  four_hits_2.assign(4, temphit);
  vector<SimpleHit3D> eight_hits;
  eight_hits.assign(8, temphit);
  float min_phi_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_phi_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_phi_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float max_phi_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float min_phi_8[8];
  float max_phi_8[8];
  if (print_timings == true) {
    gettimeofday(&t1, NULL);
  }
  float min_k_array[1 << 8];
  float max_k_array[1 << 8];
  min_k_array[0] = zoomranges[zoomlevel].min_k;
  max_k_array[0] = zoomranges[zoomlevel].min_k + k_size;
  for (unsigned int k_bin = 1; k_bin < n_k; ++k_bin) {
    min_k_array[k_bin] = min_k_array[k_bin - 1] + k_size;
    max_k_array[k_bin] = max_k_array[k_bin - 1] + k_size;
  }
  for (unsigned int k_bin = 0; k_bin < n_k; ++k_bin) {
    min_k_array[k_bin] = pow(min_k_array[k_bin], pwr);
    max_k_array[k_bin] = pow(max_k_array[k_bin], pwr);
  }
  float min_d_array[1 << 8];
  float max_d_array[1 << 8];
  min_d_array[0] = zoomranges[zoomlevel].min_d;
  max_d_array[0] = zoomranges[zoomlevel].min_d + d_size;
  for (unsigned int d_bin = 1; d_bin < n_d; ++d_bin) {
    min_d_array[d_bin] = min_d_array[d_bin - 1] + d_size;
    max_d_array[d_bin] = max_d_array[d_bin - 1] + d_size;
  }

  for (unsigned int k_bin = 0; k_bin < n_k; ++k_bin) {
    float avg = 0.5 * (max_k_array[k_bin] + min_k_array[k_bin]);
    float width = 0.5 * (max_k_array[k_bin] - min_k_array[k_bin]);
    max_k_array[k_bin] = avg + width * bin_scale;
    min_k_array[k_bin] = avg - width * bin_scale;
  }
  for (unsigned int d_bin = 0; d_bin < n_d; ++d_bin) {
    float avg = 0.5 * (max_d_array[d_bin] + min_d_array[d_bin]);
    float width = 0.5 * (max_d_array[d_bin] - min_d_array[d_bin]);
    max_d_array[d_bin] = avg + width * bin_scale;
    min_d_array[d_bin] = avg - width * bin_scale;
  }

  for (unsigned int i = 0; i < hits_vec[zoomlevel]->size(); i++) {
    if (hit_counter < 4) {
      four_hits[hit_counter] = ((*(hits_vec[zoomlevel]))[i]);
      x_a[hit_counter] = four_hits[hit_counter].get_x();
      y_a[hit_counter] = four_hits[hit_counter].get_y();
      four_hits[hit_counter].set_id(i);
      eight_hits[hit_counter] = four_hits[hit_counter];
    } else {
      four_hits_2[hit_counter - 4] = ((*(hits_vec[zoomlevel]))[i]);
      x_a_2[hit_counter - 4] = four_hits_2[hit_counter - 4].get_x();
      y_a_2[hit_counter - 4] = four_hits_2[hit_counter - 4].get_y();
      four_hits_2[hit_counter - 4].set_id(i);
      eight_hits[hit_counter] = four_hits_2[hit_counter - 4];
    }
    hit_counter++;
    if ((hit_counter == 8) && (separate_by_helicity == true)) {
      for (unsigned int d_bin = 0; d_bin < n_d; ++d_bin) {
        float min_d_a[4] __attribute__((aligned(16))) = {
            min_d_array[d_bin], min_d_array[d_bin], min_d_array[d_bin],
            min_d_array[d_bin]};
        float max_d_a[4] __attribute__((aligned(16))) = {
            max_d_array[d_bin], max_d_array[d_bin], max_d_array[d_bin],
            max_d_array[d_bin]};

        for (unsigned int k_bin = 0; k_bin < n_k; ++k_bin) {
          float min_k_a[4] __attribute__((aligned(16))) = {
              min_k_array[k_bin], min_k_array[k_bin], min_k_array[k_bin],
              min_k_array[k_bin]};
          float max_k_a[4] __attribute__((aligned(16))) = {
              max_k_array[k_bin], max_k_array[k_bin], max_k_array[k_bin],
              max_k_array[k_bin]};
          float hel = -1.;
          if (helicity == true) {
            hel = 1.;
          }
          if (k_bin == 0) {
            HelixHough::phiRange_sse(
                x_a, y_a, min_d_a, max_d_a, min_k_a, max_k_a, min_phi_1_a,
                max_phi_1_a, min_phi_2_a, max_phi_2_a, hel, phi_3_out,
                phi_4_out, x_a_2, y_a_2, phi_3_out_2, phi_4_out_2);

            phi_3_in = phi_3_out;
            phi_4_in = phi_4_out;
            phi_3_in_2 = phi_3_out_2;
            phi_4_in_2 = phi_4_out_2;
          } else {
            HelixHough::phiRange_sse(
                x_a, y_a, min_d_a, max_d_a, min_k_a, max_k_a, min_phi_1_a,
                max_phi_1_a, min_phi_2_a, max_phi_2_a, hel, phi_3_in, phi_4_in,
                phi_3_out, phi_4_out, x_a_2, y_a_2, phi_3_in_2, phi_4_in_2,
                phi_3_out_2, phi_4_out_2);

            phi_3_in = phi_3_out;
            phi_4_in = phi_4_out;
            phi_3_in_2 = phi_3_out_2;
            phi_4_in_2 = phi_4_out_2;
          }
          for (unsigned int h = 0; h < hit_counter; ++h) {
            if (h < 4) {
	      
              float dphi = sqrt(((2.0*sqrt(four_hits[h].get_size(0,0))) *
				 (2.0*sqrt(four_hits[h].get_size(0,0))) +
                                 (2.0*sqrt(four_hits[h].get_size(1,1))) *
				 (2.0*sqrt(four_hits[h].get_size(1,1)))) /
                                (four_hits[h].get_x() * four_hits[h].get_x() +
                                 four_hits[h].get_y() * four_hits[h].get_y()));
              dphi += phiError(
                  four_hits[h], min_kappa, max_kappa, min_d_array[d_bin],
                  max_d_array[d_bin], zoomranges[zoomlevel].min_z0,
                  zoomranges[zoomlevel].max_z0, zoomranges[zoomlevel].min_dzdl,
                  zoomranges[zoomlevel].max_dzdl);

              min_phi_1_a[h] -= dphi;
              max_phi_1_a[h] += dphi;
              min_phi_8[h] = min_phi_1_a[h];
              max_phi_8[h] = max_phi_1_a[h];
            } else {	      
              float dphi =
                  sqrt(((2.0*sqrt(four_hits_2[h - 4].get_size(0,0))) *
			(2.0*sqrt(four_hits_2[h - 4].get_size(0,0))) +
                        (2.0*sqrt(four_hits_2[h - 4].get_size(1,1))) *
			(2.0*sqrt(four_hits_2[h - 4].get_size(1,1)))) /
                       (four_hits_2[h - 4].get_x() * four_hits_2[h - 4].get_x() +
                        four_hits_2[h - 4].get_y() * four_hits_2[h - 4].get_y()));
              dphi += phiError(
                  four_hits[h - 4], min_kappa, max_kappa, min_d_array[d_bin],
                  max_d_array[d_bin], zoomranges[zoomlevel].min_z0,
                  zoomranges[zoomlevel].max_z0, zoomranges[zoomlevel].min_dzdl,
                  zoomranges[zoomlevel].max_dzdl);

              min_phi_2_a[h - 4] -= dphi;
              max_phi_2_a[h - 4] += dphi;
              min_phi_8[h] = min_phi_2_a[h - 4];
              max_phi_8[h] = max_phi_2_a[h - 4];
            }
          }
          fillBins(total_bins, 8, min_phi_8, max_phi_8, eight_hits, z_bins, n_d,
                   n_k, n_dzdl, n_z0, d_bin, k_bin, n_phi, zoomlevel, low_phi,
                   high_phi, inv_phi_range, vote_array);
        }
      }
      hit_counter = 0;
    }
    if ((hit_counter == 4) && (((hits_vec[zoomlevel]->size() - (i + 1)) < 4) ||
                               (separate_by_helicity == false))) {
      for (unsigned int d_bin = 0; d_bin < n_d; ++d_bin) {
        float min_d_a[4] __attribute__((aligned(16))) = {
            min_d_array[d_bin], min_d_array[d_bin], min_d_array[d_bin],
            min_d_array[d_bin]};
        float max_d_a[4] __attribute__((aligned(16))) = {
            max_d_array[d_bin], max_d_array[d_bin], max_d_array[d_bin],
            max_d_array[d_bin]};

        for (unsigned int k_bin = 0; k_bin < n_k; ++k_bin) {
          float min_k_a[4] __attribute__((aligned(16))) = {
              min_k_array[k_bin], min_k_array[k_bin], min_k_array[k_bin],
              min_k_array[k_bin]};
          float max_k_a[4] __attribute__((aligned(16))) = {
              max_k_array[k_bin], max_k_array[k_bin], max_k_array[k_bin],
              max_k_array[k_bin]};
          if (separate_by_helicity == true) {
            float hel = -1.;
            if (helicity == true) {
              hel = 1.;
            }
            if (k_bin == 0) {
              HelixHough::phiRange_sse(x_a, y_a, min_d_a, max_d_a, min_k_a,
                                       max_k_a, min_phi_1_a, max_phi_1_a, hel,
                                       phi_3_out, phi_4_out);
              phi_3_in = phi_3_out;
              phi_4_in = phi_4_out;
            } else {
              HelixHough::phiRange_sse(x_a, y_a, min_d_a, max_d_a, max_k_a,
                                       min_phi_1_a, max_phi_1_a, hel, phi_3_in,
                                       phi_4_in, phi_3_out, phi_4_out);
              phi_3_in = phi_3_out;
              phi_4_in = phi_4_out;
            }
          } else {
            HelixHough::phiRange_sse(x_a, y_a, min_d_a, max_d_a, min_k_a,
                                     max_k_a, min_phi_1_a, max_phi_1_a,
                                     min_phi_2_a, max_phi_2_a);
          }
          for (unsigned int h = 0; h < hit_counter; ++h) {
	    
            float dphi = sqrt((2.0*sqrt(four_hits[h].get_size(0,0))) *
			      (2.0*sqrt(four_hits[h].get_size(0,0))) +
                              (2.0*sqrt(four_hits[h].get_size(1,1))) *
			      (2.0*sqrt(four_hits[h].get_size(1,1)))) /
                              (four_hits[h].get_x() * four_hits[h].get_x() +
                               four_hits[h].get_y() * four_hits[h].get_y());
            dphi += phiError(
                four_hits[h], min_kappa, max_kappa, min_d_array[d_bin],
                max_d_array[d_bin], zoomranges[zoomlevel].min_z0,
                zoomranges[zoomlevel].max_z0, zoomranges[zoomlevel].min_dzdl,
                zoomranges[zoomlevel].max_dzdl);

            min_phi_1_a[h] -= dphi;
            min_phi_2_a[h] -= dphi;
            max_phi_1_a[h] += dphi;
            max_phi_2_a[h] += dphi;
          }
          if (separate_by_helicity == true) {
            fillBins(total_bins, hit_counter, min_phi_1_a, max_phi_1_a,
                     four_hits, z_bins, n_d, n_k, n_dzdl, n_z0, d_bin, k_bin,
                     n_phi, zoomlevel, low_phi, high_phi, inv_phi_range,
                     vote_array);
          } else {
            fillBins(total_bins, hit_counter, min_phi_1_a, max_phi_1_a,
                     four_hits, z_bins, n_d, n_k, n_dzdl, n_z0, d_bin, k_bin,
                     n_phi, zoomlevel, low_phi, high_phi, inv_phi_range,
                     vote_array);
            fillBins(total_bins, hit_counter, min_phi_2_a, max_phi_2_a,
                     four_hits, z_bins, n_d, n_k, n_dzdl, n_z0, d_bin, k_bin,
                     n_phi, zoomlevel, low_phi, high_phi, inv_phi_range,
                     vote_array);
          }
        }
      }
      hit_counter = 0;
    }
  }
  if (hit_counter != 0) {
    for (unsigned int d_bin = 0; d_bin < n_d; ++d_bin) {
      float min_d_a[4] __attribute__((aligned(16))) = {
          min_d_array[d_bin], min_d_array[d_bin], min_d_array[d_bin],
          min_d_array[d_bin]};
      float max_d_a[4] __attribute__((aligned(16))) = {
          max_d_array[d_bin], max_d_array[d_bin], max_d_array[d_bin],
          max_d_array[d_bin]};

      for (unsigned int k_bin = 0; k_bin < n_k; ++k_bin) {
        float min_k_a[4] __attribute__((aligned(16))) = {
            min_k_array[k_bin], min_k_array[k_bin], min_k_array[k_bin],
            min_k_array[k_bin]};
        float max_k_a[4] __attribute__((aligned(16))) = {
            max_k_array[k_bin], max_k_array[k_bin], max_k_array[k_bin],
            max_k_array[k_bin]};
        if (separate_by_helicity == true) {
          float hel = -1.;
          if (helicity == true) {
            hel = 1.;
          }
          if (k_bin == 0) {
            HelixHough::phiRange_sse(x_a, y_a, min_d_a, max_d_a, min_k_a,
                                     max_k_a, min_phi_1_a, max_phi_1_a, hel,
                                     phi_3_out, phi_4_out);

            phi_3_in = phi_3_out;
            phi_4_in = phi_4_out;
          } else {
            HelixHough::phiRange_sse(x_a, y_a, min_d_a, max_d_a, max_k_a,
                                     min_phi_1_a, max_phi_1_a, hel, phi_3_in,
                                     phi_4_in, phi_3_out, phi_4_out);
            phi_3_in = phi_3_out;
            phi_4_in = phi_4_out;
          }
        } else {
          HelixHough::phiRange_sse(x_a, y_a, min_d_a, max_d_a, min_k_a, max_k_a,
                                   min_phi_1_a, max_phi_1_a, min_phi_2_a,
                                   max_phi_2_a);
        }
        for (unsigned int h = 0; h < hit_counter; ++h) {
	  
          float dphi = sqrt(((2.0*sqrt(four_hits[h].get_size(0,0))) *
			     (2.0*sqrt(four_hits[h].get_size(0,0))) +
                             (2.0*sqrt(four_hits[h].get_size(1,1))) *
			     (2.0*sqrt(four_hits[h].get_size(1,1)))) /
                            (four_hits[h].get_x() * four_hits[h].get_x() +
                             four_hits[h].get_y() * four_hits[h].get_y()));
          dphi += phiError(
              four_hits[h], min_kappa, max_kappa, min_d_array[d_bin],
              max_d_array[d_bin], zoomranges[zoomlevel].min_z0,
              zoomranges[zoomlevel].max_z0, zoomranges[zoomlevel].min_dzdl,
              zoomranges[zoomlevel].max_dzdl);

          min_phi_1_a[h] -= dphi;
          min_phi_2_a[h] -= dphi;
          max_phi_1_a[h] += dphi;
          max_phi_2_a[h] += dphi;
        }
        if (separate_by_helicity == true) {
          fillBins(total_bins, hit_counter, min_phi_1_a, max_phi_1_a, four_hits,
                   z_bins, n_d, n_k, n_dzdl, n_z0, d_bin, k_bin, n_phi,
                   zoomlevel, low_phi, high_phi, inv_phi_range, vote_array);
        } else {
          fillBins(total_bins, hit_counter, min_phi_1_a, max_phi_1_a, four_hits,
                   z_bins, n_d, n_k, n_dzdl, n_z0, d_bin, k_bin, n_phi,
                   zoomlevel, low_phi, high_phi, inv_phi_range, vote_array);
          fillBins(total_bins, hit_counter, min_phi_2_a, max_phi_2_a, four_hits,
                   z_bins, n_d, n_k, n_dzdl, n_z0, d_bin, k_bin, n_phi,
                   zoomlevel, low_phi, high_phi, inv_phi_range, vote_array);
        }
      }
    }
    hit_counter = 0;
  }
  if (print_timings == true) {
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    xy_vote_time += (time2 - time1);
  }
  if (vote_array.size == 0) {
    return;
  }
  unsigned int C[(1 << 12)] = {0};
  if (vote_array.size < (1 << 14)) {
    unsigned int B[1 << 14];
    counting_sort(vote_array.size, total_bins, vote_array.arr, B, C, bins_start,
                  bins_end);
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
    counting_sort(total_bins, *(bins_vec[zoomlevel]), C, bins_start, bins_end);
    //     if(bins_vec[zoomlevel]->size() >
    //     total_bins){counting_sort(total_bins, *(bins_vec[zoomlevel]), C,
    //     bins_start, bins_end);}
    //     else{sort(bins_vec[zoomlevel]->begin(), bins_vec[zoomlevel]->end());}
  }
}
