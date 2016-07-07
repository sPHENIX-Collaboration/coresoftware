#include "sPHENIXTracker.h"
#include <sys/time.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "vector_math_inline.h"

using namespace std;
using namespace Eigen;

void sPHENIXTracker::calculateKappaTangents(
    float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a,
    float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a,
    float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a,
    float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a,
    float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a,
    float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a) {
  static const __m128 two = {2., 2., 2., 2.};

  __m128 x1 = _mm_load_ps(x1_a);
  __m128 x2 = _mm_load_ps(x2_a);
  __m128 x3 = _mm_load_ps(x3_a);
  __m128 y1 = _mm_load_ps(y1_a);
  __m128 y2 = _mm_load_ps(y2_a);
  __m128 y3 = _mm_load_ps(y3_a);
  __m128 z1 = _mm_load_ps(z1_a);
  __m128 z2 = _mm_load_ps(z2_a);
  __m128 z3 = _mm_load_ps(z3_a);

  __m128 dx1 = _mm_load_ps(dx1_a);
  __m128 dx2 = _mm_load_ps(dx2_a);
  __m128 dx3 = _mm_load_ps(dx3_a);
  __m128 dy1 = _mm_load_ps(dy1_a);
  __m128 dy2 = _mm_load_ps(dy2_a);
  __m128 dy3 = _mm_load_ps(dy3_a);
  __m128 dz1 = _mm_load_ps(dz1_a);
  __m128 dz2 = _mm_load_ps(dz2_a);
  __m128 dz3 = _mm_load_ps(dz3_a);

  __m128 D12 = _mm_sub_ps(x2, x1);
  D12 = _mm_mul_ps(D12, D12);
  __m128 tmp1 = _mm_sub_ps(y2, y1);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D12 = _mm_add_ps(D12, tmp1);
  D12 = _vec_sqrt_ps(D12);

  __m128 D23 = _mm_sub_ps(x3, x2);
  D23 = _mm_mul_ps(D23, D23);
  tmp1 = _mm_sub_ps(y3, y2);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D23 = _mm_add_ps(D23, tmp1);
  D23 = _vec_sqrt_ps(D23);

  __m128 D31 = _mm_sub_ps(x1, x3);
  D31 = _mm_mul_ps(D31, D31);
  tmp1 = _mm_sub_ps(y1, y3);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D31 = _mm_add_ps(D31, tmp1);
  D31 = _vec_sqrt_ps(D31);

  __m128 k = _mm_mul_ps(D12, D23);
  k = _mm_mul_ps(k, D31);
  k = _vec_rec_ps(k);
  tmp1 = (D12 + D23 + D31) * (D23 + D31 - D12) * (D12 + D31 - D23) *
         (D12 + D23 - D31);
  tmp1 = _vec_sqrt_ps(tmp1);
  k *= tmp1;

  __m128 tmp2 = _mm_cmpgt_ps(tmp1, zero);
  tmp1 = _mm_and_ps(tmp2, k);
  tmp2 = _mm_andnot_ps(tmp2, zero);
  k = _mm_xor_ps(tmp1, tmp2);

  _mm_store_ps(kappa_a, k);
  __m128 k_inv = _vec_rec_ps(k);

  __m128 D12_inv = _vec_rec_ps(D12);
  __m128 D23_inv = _vec_rec_ps(D23);
  __m128 D31_inv = _vec_rec_ps(D31);

  __m128 dr1 = dx1 * dx1 + dy1 * dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2 * dx2 + dy2 * dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3 * dx3 + dy3 * dy3;
  dr3 = _vec_sqrt_ps(dr3);

  __m128 dk1 = (dr1 + dr2) * D12_inv * D12_inv;
  __m128 dk2 = (dr2 + dr3) * D23_inv * D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);

  __m128 ux12 = (x2 - x1) * D12_inv;
  __m128 uy12 = (y2 - y1) * D12_inv;
  __m128 ux23 = (x3 - x2) * D23_inv;
  __m128 uy23 = (y3 - y2) * D23_inv;
  __m128 ux13 = (x3 - x1) * D31_inv;
  __m128 uy13 = (y3 - y1) * D31_inv;

  __m128 cosalpha = ux12 * ux13 + uy12 * uy13;
  __m128 sinalpha = ux13 * uy12 - ux12 * uy13;

  __m128 ux_mid = ux23 * cosalpha - uy23 * sinalpha;
  __m128 uy_mid = ux23 * sinalpha + uy23 * cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);

  __m128 ux_end = ux23 * cosalpha + uy23 * sinalpha;
  __m128 uy_end = uy23 * cosalpha - ux23 * sinalpha;

  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);

  // asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s2 = _vec_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s2);
  tmp1 = _mm_andnot_ps(tmp1, D23);
  s2 = _mm_xor_ps(tmp1, tmp2);

  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m128 del_z_2 = z3 - z2;
  __m128 dzdl_2 = s2 * s2 + del_z_2 * del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3) * D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);

  sinalpha = ux13 * uy23 - ux23 * uy13;
  v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s1 = _vec_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s1);
  tmp1 = _mm_andnot_ps(tmp1, D12);
  s1 = _mm_xor_ps(tmp1, tmp2);

  __m128 del_z_1 = z2 - z1;
  __m128 dzdl_1 = s1 * s1 + del_z_1 * del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2) * D12_inv;
  _mm_store_ps(dzdl_1_a, dzdl_1);
  _mm_store_ps(ddzdl_1_a, ddzdl_1);
}

void sPHENIXTracker::calculateKappaTangents(
    float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a,
    float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a,
    float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a,
    float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a,
    float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a,
    float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a,
    float sinang_cut, float cosang_diff_inv, float* cur_kappa_a,
    float* cur_dkappa_a, float* cur_ux_a, float* cur_uy_a, float* cur_chi2_a,
    float* chi2_a) {
  static const __m128 two = {2., 2., 2., 2.};

  __m128 x1 = _mm_load_ps(x1_a);
  __m128 x2 = _mm_load_ps(x2_a);
  __m128 x3 = _mm_load_ps(x3_a);
  __m128 y1 = _mm_load_ps(y1_a);
  __m128 y2 = _mm_load_ps(y2_a);
  __m128 y3 = _mm_load_ps(y3_a);
  __m128 z1 = _mm_load_ps(z1_a);
  __m128 z2 = _mm_load_ps(z2_a);
  __m128 z3 = _mm_load_ps(z3_a);

  __m128 dx1 = _mm_load_ps(dx1_a);
  __m128 dx2 = _mm_load_ps(dx2_a);
  __m128 dx3 = _mm_load_ps(dx3_a);
  __m128 dy1 = _mm_load_ps(dy1_a);
  __m128 dy2 = _mm_load_ps(dy2_a);
  __m128 dy3 = _mm_load_ps(dy3_a);
  __m128 dz1 = _mm_load_ps(dz1_a);
  __m128 dz2 = _mm_load_ps(dz2_a);
  __m128 dz3 = _mm_load_ps(dz3_a);

  __m128 D12 = _mm_sub_ps(x2, x1);
  D12 = _mm_mul_ps(D12, D12);
  __m128 tmp1 = _mm_sub_ps(y2, y1);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D12 = _mm_add_ps(D12, tmp1);
  D12 = _vec_sqrt_ps(D12);

  __m128 D23 = _mm_sub_ps(x3, x2);
  D23 = _mm_mul_ps(D23, D23);
  tmp1 = _mm_sub_ps(y3, y2);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D23 = _mm_add_ps(D23, tmp1);
  D23 = _vec_sqrt_ps(D23);

  __m128 D31 = _mm_sub_ps(x1, x3);
  D31 = _mm_mul_ps(D31, D31);
  tmp1 = _mm_sub_ps(y1, y3);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D31 = _mm_add_ps(D31, tmp1);
  D31 = _vec_sqrt_ps(D31);

  __m128 k = _mm_mul_ps(D12, D23);
  k = _mm_mul_ps(k, D31);
  k = _vec_rec_ps(k);
  tmp1 = (D12 + D23 + D31) * (D23 + D31 - D12) * (D12 + D31 - D23) *
         (D12 + D23 - D31);
  tmp1 = _vec_sqrt_ps(tmp1);
  k *= tmp1;

  __m128 tmp2 = _mm_cmpgt_ps(tmp1, zero);
  tmp1 = _mm_and_ps(tmp2, k);
  tmp2 = _mm_andnot_ps(tmp2, zero);
  k = _mm_xor_ps(tmp1, tmp2);

  _mm_store_ps(kappa_a, k);
  __m128 k_inv = _vec_rec_ps(k);

  __m128 D12_inv = _vec_rec_ps(D12);
  __m128 D23_inv = _vec_rec_ps(D23);
  __m128 D31_inv = _vec_rec_ps(D31);

  __m128 dr1 = dx1 * dx1 + dy1 * dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2 * dx2 + dy2 * dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3 * dx3 + dy3 * dy3;
  dr3 = _vec_sqrt_ps(dr3);

  __m128 dk1 = (dr1 + dr2) * D12_inv * D12_inv;
  __m128 dk2 = (dr2 + dr3) * D23_inv * D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);

  __m128 ux12 = (x2 - x1) * D12_inv;
  __m128 uy12 = (y2 - y1) * D12_inv;
  __m128 ux23 = (x3 - x2) * D23_inv;
  __m128 uy23 = (y3 - y2) * D23_inv;
  __m128 ux13 = (x3 - x1) * D31_inv;
  __m128 uy13 = (y3 - y1) * D31_inv;

  __m128 cosalpha = ux12 * ux13 + uy12 * uy13;
  __m128 sinalpha = ux13 * uy12 - ux12 * uy13;

  __m128 ux_mid = ux23 * cosalpha - uy23 * sinalpha;
  __m128 uy_mid = ux23 * sinalpha + uy23 * cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);

  __m128 ux_end = ux23 * cosalpha + uy23 * sinalpha;
  __m128 uy_end = uy23 * cosalpha - ux23 * sinalpha;

  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);

  // asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s2 = _vec_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s2);
  tmp1 = _mm_andnot_ps(tmp1, D23);
  s2 = _mm_xor_ps(tmp1, tmp2);

  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m128 del_z_2 = z3 - z2;
  __m128 dzdl_2 = s2 * s2 + del_z_2 * del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3) * D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);

  sinalpha = ux13 * uy23 - ux23 * uy13;
  v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s1 = _vec_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s1);
  tmp1 = _mm_andnot_ps(tmp1, D12);
  s1 = _mm_xor_ps(tmp1, tmp2);

  __m128 del_z_1 = z2 - z1;
  __m128 dzdl_1 = s1 * s1 + del_z_1 * del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2) * D12_inv;
  _mm_store_ps(dzdl_1_a, dzdl_1);
  _mm_store_ps(ddzdl_1_a, ddzdl_1);

  __m128 c_dk = _mm_load_ps(cur_dkappa_a);
  __m128 c_k = _mm_load_ps(cur_kappa_a);
  __m128 c_ux = _mm_load_ps(cur_ux_a);
  __m128 c_uy = _mm_load_ps(cur_uy_a);
  __m128 c_chi2 = _mm_load_ps(cur_chi2_a);
  __m128 sinang = _mm_load1_ps(&sinang_cut);
  __m128 cosdiff = _mm_load1_ps(&cosang_diff_inv);

  __m128 kdiff = c_k - k;
  __m128 n_dk = c_dk + dk + sinang * k;
  __m128 chi2_k = kdiff * kdiff / (n_dk * n_dk);
  __m128 cos_scatter = c_ux * ux_mid + c_uy * uy_mid;
  __m128 chi2_ang =
      (one - cos_scatter) * (one - cos_scatter) * cosdiff * cosdiff;
  tmp1 = dzdl_1 * sinang;
  _vec_fabs_ps(tmp1);
  __m128 chi2_dzdl = (dzdl_1 - dzdl_2) / (ddzdl_1 + ddzdl_2 + tmp1);
  chi2_dzdl *= chi2_dzdl;
  chi2_dzdl *= one_o_2;

  __m128 n_chi2 = c_chi2 + chi2_ang + chi2_k + chi2_dzdl;
  _mm_store_ps(chi2_a, n_chi2);
}

struct TempComb {
  TempComb() {}
  HelixKalmanState state;
  SimpleTrack3D track;

  inline bool operator<(const TempComb& other) const {
    if (track.hits.size() > other.track.hits.size()) {
      return true;
    } else if (track.hits.size() == other.track.hits.size()) {
      return state.chi2 < other.state.chi2;
    } else {
      return false;
    }
  }
};

void sPHENIXTracker::initDummyHits(vector<SimpleHit3D>& dummies,
                                   const HelixRange& range,
                                   HelixKalmanState& init_state) {
  SimpleTrack3D dummy_track;
  dummy_track.hits.push_back(SimpleHit3D());
  dummy_track.kappa = 0.5 * (range.min_k + range.max_k);
  dummy_track.phi = 0.5 * (range.min_phi + range.max_phi);
  dummy_track.d = 0.5 * (range.min_d + range.max_d);
  dummy_track.dzdl = 0.5 * (range.min_dzdl + range.max_dzdl);
  dummy_track.z0 = 0.5 * (range.min_z0 + range.max_z0);

  init_state.kappa = dummy_track.kappa;
  init_state.nu = sqrt(dummy_track.kappa);
  init_state.phi = dummy_track.phi;
  init_state.d = dummy_track.d;
  init_state.dzdl = dummy_track.dzdl;
  init_state.z0 = dummy_track.z0;

  init_state.C = Matrix<float, 5, 5>::Zero(5, 5);
  init_state.C(0, 0) = pow(range.max_phi - range.min_phi, 2.);
  init_state.C(1, 1) = pow(range.max_d - range.min_d, 2.);
  init_state.C(2, 2) = pow(10. * sqrt(range.max_k - range.min_k), 2.);
  init_state.C(3, 3) = pow(range.max_z0 - range.min_z0, 2.);
  init_state.C(4, 4) = pow(range.max_dzdl - range.min_dzdl, 2.);
  init_state.chi2 = 0.;
  init_state.position = 0;
  init_state.x_int = 0.;
  init_state.y_int = 0.;
  init_state.z_int = 0.;

  for (unsigned int i = 0; i < n_layers; ++i) {
    float x, y, z;
    projectToLayer(dummy_track, i, x, y, z);
    dummies[i].set_x(x);
    dummies[i].set_ex(5.);
    dummies[i].set_y(x);
    dummies[i].set_ey(5.);
    dummies[i].set_z(x);
    dummies[i].set_ez(5.);
    dummies[i].set_layer(i);
  }
}

static bool next_combo_n(vector<int> const& lsizes, vector<int>& comb_n) {
  unsigned int n = lsizes.size() / 2;
  for (int l = 0; l < n; ++l) {
    if (comb_n[l] == (lsizes[l] - 1)) {
      comb_n[l] = 0;
    } else {
      comb_n[l] += 1;
      return true;
    }
  }
  return false;
}

static SimpleHit3D& get_hit(vector<SimpleHit3D>& hits,
                            vector<SimpleHit3D>& dummies, int index) {
  if (index >= 0) {
    return hits[index];
  } else {
    return dummies[(-index) - 1];
  }
}

class hit_triplet {
 public:
  hit_triplet(unsigned int h1, unsigned int h2, unsigned int h3, unsigned int t,
              float c)
      : hit1(h1), hit2(h2), hit3(h3), track(t), chi2(c) {}
  ~hit_triplet() {}

  bool operator<(const hit_triplet& other) const {
    return (hit1 < other.hit1) ||
           ((hit2 < other.hit2) && (hit1 == other.hit1)) ||
           ((hit3 < other.hit3) && (hit1 == other.hit1) &&
            (hit2 == other.hit2));
  }

  bool operator==(const hit_triplet& other) const {
    return ((hit1 == other.hit1) && (hit2 == other.hit2) &&
            (hit3 == other.hit3));
  }

  unsigned int hit1, hit2, hit3, track;
  float chi2;
};

static void triplet_rejection(vector<SimpleTrack3D>& input,
                              vector<float>& chi2s, vector<bool>& usetrack) {
  vector<hit_triplet> trips;
  for (unsigned int i = 0; i < input.size(); ++i) {
    for (unsigned int h1 = 0; h1 < input[i].hits.size(); ++h1) {
      for (unsigned int h2 = (h1 + 1); h2 < input[i].hits.size(); ++h2) {
        for (unsigned int h3 = (h2 + 1); h3 < input[i].hits.size(); ++h3) {
          trips.push_back(hit_triplet(input[i].hits[h1].get_id(),
                                      input[i].hits[h2].get_id(),
                                      input[i].hits[h3].get_id(), i, chi2s[i]));
        }
      }
    }
  }
  if (trips.size() == 0) {
    return;
  }
  sort(trips.begin(), trips.end());
  unsigned int pos = 0;
  unsigned int cur_h1 = trips[pos].hit1;
  unsigned int cur_h2 = trips[pos].hit2;
  while (pos < trips.size()) {
    unsigned int next_pos = pos + 1;
    if (next_pos >= trips.size()) {
      break;
    }
    while (trips[pos] == trips[next_pos]) {
      next_pos += 1;
      if (next_pos >= trips.size()) {
        break;
      }
    }
    if ((next_pos - pos) > 1) {
      float best_chi2 = trips[pos].chi2;
      float next_chi2 = trips[pos + 1].chi2;
      unsigned int best_pos = pos;
      for (unsigned int i = (pos + 1); i < next_pos; ++i) {
        if (input[trips[i].track].hits.size() <
            input[trips[best_pos].track].hits.size()) {
          continue;
        } else if ((input[trips[i].track].hits.size() >
                    input[trips[best_pos].track].hits.size()) ||
                   (input[trips[i].track].hits.back().get_layer() >
                    input[trips[best_pos].track].hits.back().get_layer())) {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
          continue;
        }
        if ((trips[i].chi2 < best_chi2) ||
            (usetrack[trips[best_pos].track] == false)) {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
        } else if (trips[i].chi2 < next_chi2) {
          next_chi2 = trips[i].chi2;
        }
      }
      for (unsigned int i = pos; i < next_pos; ++i) {
        if (i != best_pos) {
          usetrack[trips[i].track] = false;
        }
      }
    }
    pos = next_pos;
    cur_h1 = trips[pos].hit1;
    cur_h2 = trips[pos].hit2;
  }
}

static bool remove_bad_hits(SimpleTrack3D& track, float cut) {
  SimpleTrack3D temp_track = track;
  float fit_chi2 = 0.;
  vector<float> chi2_hit;
  vector<float> temp_hits;
  while (true) {
    temp_track = track;
    fit_chi2 = sPHENIXTracker::fitTrack(temp_track, chi2_hit);
    bool all_good = true;
    track.hits.clear();
    for (int h = 0; h < temp_track.hits.size(); h += 1) {
      if (chi2_hit[h] < cut) {
        track.hits.push_back(temp_track.hits[h]);
      } else {
        all_good = false;
      }
    }
    if (track.hits.size() < 3) {
      return false;
    }
    if (all_good == true) {
      return true;
    }
  }
}

static void initial_combos(int nhits, vector<SimpleHit3D>& hits,
                           vector<vector<int> >& layer_indexes,
                           vector<TempComb>& cur_comb, float CHI2_CUT,
                           int n_layers, CylinderKalman* kalman) {
  vector<int> lsizes;
  for (int i = 0; i < nhits; ++i) {
    lsizes.push_back(layer_indexes[i].size());
  }
  vector<int> comb_n;
  comb_n.assign(nhits, 0);
  vector<TempComb> comb;
  while (true) {
    SimpleTrack3D temp_track;
    TempComb tc;
    float sqrt12_inv = 1. / sqrt(12.);
    for (unsigned int h = 0; h < nhits; ++h) {
      if (layer_indexes[h][comb_n[h]] < 0) {
        continue;
      }
      SimpleHit3D& hit = hits[layer_indexes[h][comb_n[h]]];
      temp_track.hits.push_back(hit);
      temp_track.hits.back().set_ex( temp_track.hits.back().get_ex() * sqrt12_inv);
      temp_track.hits.back().set_ey( temp_track.hits.back().get_ex() * sqrt12_inv);
      temp_track.hits.back().set_ez( temp_track.hits.back().get_ex() * sqrt12_inv);
    }
    if (temp_track.hits.size() >= 3) {
      float init_chi2 = sPHENIXTracker::fitTrack(temp_track);
      comb.push_back(TempComb());
      TempComb& curcomb = comb.back();
      if (init_chi2 / (2. * (temp_track.hits.size()) - 5.) > CHI2_CUT) {
        comb.pop_back();
      } else {
        curcomb.state.phi = temp_track.phi;
        if (curcomb.state.phi < 0.) {
          curcomb.state.phi += 2. * M_PI;
        }
        curcomb.state.d = temp_track.d;
        curcomb.state.kappa = temp_track.kappa;
        curcomb.state.nu = sqrt(curcomb.state.kappa);
        curcomb.state.z0 = temp_track.z0;
        curcomb.state.dzdl = temp_track.dzdl;
        curcomb.state.C = Matrix<float, 5, 5>::Zero(5, 5);
        curcomb.state.C(0, 0) = pow(0.01, 2.);
        curcomb.state.C(1, 1) = pow(0.5, 2.);
        curcomb.state.C(2, 2) = pow(0.3 * curcomb.state.nu, 2.);
        curcomb.state.C(3, 3) = pow(0.5, 2.);
        curcomb.state.C(4, 4) = pow(0.05, 2.);
        curcomb.state.chi2 = 0.;
        curcomb.state.position = n_layers;
        curcomb.state.x_int = 0.;
        curcomb.state.y_int = 0.;
        curcomb.state.z_int = 0.;

        curcomb.track = temp_track;

        for (int h = 0; h < curcomb.track.hits.size(); h += 1) {
          kalman->addHit(curcomb.track.hits[h], curcomb.state);
        }
      }
    }
    if (next_combo_n(lsizes, comb_n) == false) {
      break;
    }
  }

  vector<SimpleTrack3D> comb_tracks;
  vector<TempComb> comb2;
  vector<float> chi2s;
  for (int t = 0; t < comb.size(); t += 1) {
    if (remove_bad_hits(comb[t].track, CHI2_CUT + 2.) == false) {
      continue;
    }
    comb_tracks.push_back(comb[t].track);
    chi2s.push_back(comb[t].state.chi2);
    comb2.push_back(comb[t]);
  }
  vector<bool> usetrack(comb2.size(), true);

  triplet_rejection(comb_tracks, chi2s, usetrack);

  for (int t = 0; t < comb_tracks.size(); t += 1) {
    if (usetrack[t]) {
      if (comb2[t].state.chi2 / (2. * (comb2[t].track.hits.size()) - 5.) <
          CHI2_CUT) {
        cur_comb.push_back(comb2[t]);
      }
    }
  }

  sort(cur_comb.begin(), cur_comb.end());
  if (cur_comb.size() > 4) {
    cur_comb.resize(4, TempComb());
  }
}

static void extend_combos(vector<SimpleHit3D>& hits,
                          vector<vector<int> >& layer_indexes,
                          vector<TempComb>& cur_comb, float CHI2_CUT,
                          CylinderKalman* kalman, int layer_begin,
                          int layer_end) {
  vector<TempComb> comb;
  vector<TempComb> comb2;
  vector<TempComb> comb3;
  for (int c = 0; c < cur_comb.size(); c += 1) {
    comb2.clear();
    comb2.push_back(cur_comb[c]);
    for (int l = layer_begin; l <= layer_end; l += 1) {
      comb3.clear();
      for (int i = 0; i < comb2.size(); i += 1) {
        comb3.push_back(comb2[i]);
        for (int h = 0; h < (layer_indexes[l].size() - 1); h += 1) {
          comb3.push_back(comb2[i]);
          kalman->addHit(hits[layer_indexes[l][h]], comb3.back().state);
          if (comb3.back().state.chi2 >
              CHI2_CUT * (2. * (comb3.back().track.hits.size()) - 5.)) {
            comb3.pop_back();
          } else {
            comb3.back().track.hits.push_back(hits[layer_indexes[l][h]]);
          }
        }
      }
      comb2 = comb3;
    }
    for (int i = 0; i < comb2.size(); i += 1) {
      comb.push_back(comb2[i]);
    }
  }
  cur_comb.clear();

  vector<SimpleTrack3D> comb_tracks;
  comb2.clear();
  vector<float> chi2s;
  for (int t = 0; t < comb.size(); t += 1) {
    if (remove_bad_hits(comb[t].track, CHI2_CUT + 2.) == false) {
      continue;
    }
    comb_tracks.push_back(comb[t].track);
    chi2s.push_back(comb[t].state.chi2);
    comb2.push_back(comb[t]);
  }
  vector<bool> usetrack(comb2.size(), true);

  triplet_rejection(comb_tracks, chi2s, usetrack);

  for (int t = 0; t < comb_tracks.size(); t += 1) {
    if (usetrack[t]) {
      if (comb2[t].state.chi2 / (2. * (comb2[t].track.hits.size()) - 5.) <
          CHI2_CUT) {
        cur_comb.push_back(comb2[t]);
      }
    }
  }

  sort(cur_comb.begin(), cur_comb.end());
  if (cur_comb.size() > 4) {
    cur_comb.resize(4, TempComb());
  }
}

void sPHENIXTracker::findTracksByCombinatorialKalman(
    vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks,
    const HelixRange& range) {
  float CHI2_CUT = chi2_cut + 2.;

  vector<vector<int> > layer_indexes;
  layer_indexes.assign(n_layers, vector<int>());
  for (unsigned int i = 0; i < hits.size(); ++i) {
    layer_indexes[hits[i].get_layer()].push_back(i);
  }
  for (int i = 0; i < (int)n_layers; ++i) {
    layer_indexes[i].push_back(-(i + 1));
  }
  vector<TempComb> cur_comb;

  initial_combos(12, hits, layer_indexes, cur_comb, CHI2_CUT, n_layers, kalman);
  if (cur_comb.size() == 0) {
    return;
  }

  int cur_layer = 12;
  int layer_iter = 1;
  while (true) {
    bool finished = false;
    int layer_end = cur_layer + layer_iter - 1;
    if (layer_end >= (n_layers - 1)) {
      finished = true;
      layer_end = n_layers - 1;
    }
    extend_combos(hits, layer_indexes, cur_comb, CHI2_CUT, kalman, cur_layer,
                  layer_end);
    if (finished == true) {
      break;
    }
    cur_layer = layer_end + 1;

    vector<TempComb> tcomb;
    for (int i = 0; i < cur_comb.size(); i += 1) {
      if (cur_comb[i].track.hits.size() > ((cur_layer + 1) / 2)) {
        if (cur_comb[i].state.chi2 <
            chi2_cut * (2. * (cur_comb[i].track.hits.size()) - 5.)) {
          tcomb.push_back(cur_comb[i]);
        }
      }
    }
    cur_comb = tcomb;

    if (cur_comb.size() == 0) {
      return;
    }
  }

  for (int i = 0; i < cur_comb.size(); i += 1) {
    if (!(cur_comb[i].state.chi2 == cur_comb[i].state.chi2)) {
      continue;
    }
    tracks.push_back(cur_comb[i].track);
    track_states.push_back(cur_comb[i].state);
    if (remove_hits == true) {
      for (unsigned int i = 0; i < tracks.back().hits.size(); ++i) {
        (*hit_used)[tracks.back().hits[i].get_id()] = true;
      }
    }
  }
}

void sPHENIXTracker::findTracksBySegments(vector<SimpleHit3D>& hits,
                                          vector<SimpleTrack3D>& tracks,
                                          const HelixRange& range) {
  if (n_layers > 20) {
    findTracksByCombinatorialKalman(hits, tracks, range);
    return;
  }

  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;
  unsigned int curseg_size = 0;
  unsigned int nextseg_size = 0;

  vector<TrackSegment> complete_segments;

  unsigned int allowed_missing = n_layers - req_layers;

  for (unsigned int l = 0; l < n_layers; ++l) {
    layer_sorted[l].clear();
  }
  for (unsigned int i = 0; i < hits.size(); ++i) {
    unsigned int min = (hits[i].get_layer() - allowed_missing);
    if (allowed_missing > hits[i].get_layer()) {
      min = 0;
    }
    for (unsigned int l = min; l <= hits[i].get_layer(); l += 1) {
      layer_sorted[l].push_back(hits[i]);
      SimpleHit3D& hit = layer_sorted[l].back();
      float err_scale = 1.;
      int layer = hit.get_layer();
      if ((layer >= 0) && (layer < (int)(hit_error_scale.size()))) {
        err_scale *= 3. * hit_error_scale[layer];
      }
      hit.set_ex( hit.get_ex() * err_scale);
      hit.set_ey( hit.get_ey() * err_scale);
      hit.set_ez( hit.get_ez() * err_scale);
    }
  }
  for (unsigned int l = 0; l < n_layers; ++l) {
    if (layer_sorted[l].size() == 0) {
      return;
    }
  }

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;

  gettimeofday(&t1, NULL);

  float cosang_diff = 1. - cosang_cut;
  float cosang_diff_inv = 1. / cosang_diff;
  float sinang_cut = sqrt(1. - cosang_cut * cosang_cut);
  float easy_chi2_cut = ca_chi2_cut;

  vector<float> inv_layer;
  inv_layer.assign(n_layers, 1.);
  for (unsigned int l = 3; l < n_layers; ++l) {
    inv_layer[l] = 1. / (((float)l) - 2.);
  }

  unsigned int hit_counter = 0;
  float x1_a[4] __attribute__((aligned(16)));
  float x2_a[4] __attribute__((aligned(16)));
  float x3_a[4] __attribute__((aligned(16)));
  float y1_a[4] __attribute__((aligned(16)));
  float y2_a[4] __attribute__((aligned(16)));
  float y3_a[4] __attribute__((aligned(16)));
  float z1_a[4] __attribute__((aligned(16)));
  float z2_a[4] __attribute__((aligned(16)));
  float z3_a[4] __attribute__((aligned(16)));

  float dx1_a[4] __attribute__((aligned(16)));
  float dx2_a[4] __attribute__((aligned(16)));
  float dx3_a[4] __attribute__((aligned(16)));
  float dy1_a[4] __attribute__((aligned(16)));
  float dy2_a[4] __attribute__((aligned(16)));
  float dy3_a[4] __attribute__((aligned(16)));
  float dz1_a[4] __attribute__((aligned(16)));
  float dz2_a[4] __attribute__((aligned(16)));
  float dz3_a[4] __attribute__((aligned(16)));

  float kappa_a[4] __attribute__((aligned(16)));
  float dkappa_a[4] __attribute__((aligned(16)));

  float ux_mid_a[4] __attribute__((aligned(16)));
  float uy_mid_a[4] __attribute__((aligned(16)));
  float ux_end_a[4] __attribute__((aligned(16)));
  float uy_end_a[4] __attribute__((aligned(16)));

  float dzdl_1_a[4] __attribute__((aligned(16)));
  float dzdl_2_a[4] __attribute__((aligned(16)));
  float ddzdl_1_a[4] __attribute__((aligned(16)));
  float ddzdl_2_a[4] __attribute__((aligned(16)));

  float cur_kappa_a[4] __attribute__((aligned(16)));
  float cur_dkappa_a[4] __attribute__((aligned(16)));
  float cur_ux_a[4] __attribute__((aligned(16)));
  float cur_uy_a[4] __attribute__((aligned(16)));
  float cur_chi2_a[4] __attribute__((aligned(16)));
  float chi2_a[4] __attribute__((aligned(16)));

  unsigned int hit1[4];
  unsigned int hit2[4];
  unsigned int hit3[4];

  TrackSegment temp_segment;
  temp_segment.hits.assign(n_layers, 0);
  // make segments out of first 3 layers
  for (unsigned int i = 0, sizei = layer_sorted[0].size(); i < sizei; ++i) {
    for (unsigned int j = 0, sizej = layer_sorted[1].size(); j < sizej; ++j) {
      for (unsigned int k = 0, sizek = layer_sorted[2].size(); k < sizek; ++k) {
        if ((layer_sorted[0][i].get_layer() >= layer_sorted[1][j].get_layer()) ||
            (layer_sorted[1][j].get_layer() >= layer_sorted[2][k].get_layer())) {
          continue;
        }

        x1_a[hit_counter] = layer_sorted[0][i].get_x();
        y1_a[hit_counter] = layer_sorted[0][i].get_y();
        z1_a[hit_counter] = layer_sorted[0][i].get_z();
        dx1_a[hit_counter] = layer_sorted[0][i].get_ex();
        dy1_a[hit_counter] = layer_sorted[0][i].get_ey();
        dz1_a[hit_counter] = layer_sorted[0][i].get_ez();

        x2_a[hit_counter] = layer_sorted[1][j].get_x();
        y2_a[hit_counter] = layer_sorted[1][j].get_y();
        z2_a[hit_counter] = layer_sorted[1][j].get_z();
        dx2_a[hit_counter] = layer_sorted[1][j].get_ex();
        dy2_a[hit_counter] = layer_sorted[1][j].get_ey();
        dz2_a[hit_counter] = layer_sorted[1][j].get_ez();

        x3_a[hit_counter] = layer_sorted[2][k].get_x();
        y3_a[hit_counter] = layer_sorted[2][k].get_y();
        z3_a[hit_counter] = layer_sorted[2][k].get_z();
        dx3_a[hit_counter] = layer_sorted[2][k].get_ex();
        dy3_a[hit_counter] = layer_sorted[2][k].get_ey();
        dz3_a[hit_counter] = layer_sorted[2][k].get_ez();

        hit1[hit_counter] = i;
        hit2[hit_counter] = j;
        hit3[hit_counter] = k;

        hit_counter += 1;

        if (hit_counter == 4) {
          calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
                                 z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
                                 dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
                                 ux_mid_a, uy_mid_a, ux_end_a, uy_end_a,
                                 dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);

          for (unsigned int h = 0; h < hit_counter; ++h) {
            temp_segment.chi2 =
                (dzdl_1_a[h] - dzdl_2_a[h]) /
                (ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
            temp_segment.chi2 *= temp_segment.chi2;
            if (temp_segment.chi2 > 2.0) {
              continue;
            }
            temp_segment.ux = ux_end_a[h];
            temp_segment.uy = uy_end_a[h];
            temp_segment.kappa = kappa_a[h];
            if (temp_segment.kappa > top_range.max_k) {
              continue;
            }
            temp_segment.dkappa = dkappa_a[h];
            temp_segment.hits[0] = hit1[h];
            temp_segment.hits[1] = hit2[h];
            temp_segment.hits[2] = hit3[h];
            temp_segment.n_hits = 3;
            unsigned int outer_layer =
                layer_sorted[2][temp_segment.hits[2]].get_layer();
            if ((outer_layer - 2) > allowed_missing) {
              continue;
            }
            if ((n_layers - 3) <= allowed_missing) {
              complete_segments.push_back(temp_segment);
            }
            if (next_seg->size() == nextseg_size) {
              next_seg->push_back(temp_segment);
              nextseg_size += 1;
            } else {
              (*next_seg)[nextseg_size] = temp_segment;
              nextseg_size += 1;
            }
          }

          hit_counter = 0;
        }
      }
    }
  }
  if (hit_counter != 0) {
    calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a,
                           dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a,
                           dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a,
                           ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a,
                           ddzdl_2_a);

    for (unsigned int h = 0; h < hit_counter; ++h) {
      temp_segment.chi2 =
          (dzdl_1_a[h] - dzdl_2_a[h]) /
          (ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
      temp_segment.chi2 *= temp_segment.chi2;
      if (temp_segment.chi2 > 2.0) {
        continue;
      }
      temp_segment.ux = ux_end_a[h];
      temp_segment.uy = uy_end_a[h];
      temp_segment.kappa = kappa_a[h];
      if (temp_segment.kappa > top_range.max_k) {
        continue;
      }
      temp_segment.dkappa = dkappa_a[h];
      temp_segment.hits[0] = hit1[h];
      temp_segment.hits[1] = hit2[h];
      temp_segment.hits[2] = hit3[h];
      temp_segment.n_hits = 3;
      unsigned int outer_layer = layer_sorted[2][temp_segment.hits[2]].get_layer();
      if ((outer_layer - 2) > allowed_missing) {
        continue;
      }
      if ((n_layers - 3) <= allowed_missing) {
        complete_segments.push_back(temp_segment);
      }
      if (next_seg->size() == nextseg_size) {
        next_seg->push_back(temp_segment);
        nextseg_size += 1;
      } else {
        (*next_seg)[nextseg_size] = temp_segment;
        nextseg_size += 1;
      }
    }

    hit_counter = 0;
  }

  swap(cur_seg, next_seg);
  swap(curseg_size, nextseg_size);

  // add hits to segments layer-by-layer, cutting out bad segments
  unsigned int whichseg[4];
  for (unsigned int l = 3; l < n_layers; ++l) {
    if (l == (n_layers - 1)) {
      easy_chi2_cut *= 0.25;
    }
    nextseg_size = 0;
    for (unsigned int i = 0, sizei = curseg_size; i < sizei; ++i) {
      for (unsigned int j = 0, sizej = layer_sorted[l].size(); j < sizej; ++j) {
        if ((layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_layer() >=
             layer_sorted[l][j].get_layer())) {
          continue;
        }

        x1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_x();
        y1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_y();
        z1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_z();
        x2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_x();
        y2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_y();
        z2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_z();
        x3_a[hit_counter] = layer_sorted[l][j].get_x();
        y3_a[hit_counter] = layer_sorted[l][j].get_y();
        z3_a[hit_counter] = layer_sorted[l][j].get_z();

        dx1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_ex();
        dy1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_ey();
        dz1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_ez();
        dx2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_ex();
        dy2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_ey();
        dz2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_ez();
        dx3_a[hit_counter] = layer_sorted[l][j].get_ex();
        dy3_a[hit_counter] = layer_sorted[l][j].get_ey();
        dz3_a[hit_counter] = layer_sorted[l][j].get_ez();

        cur_kappa_a[hit_counter] = (*cur_seg)[i].kappa;
        cur_dkappa_a[hit_counter] = (*cur_seg)[i].dkappa;
        cur_ux_a[hit_counter] = (*cur_seg)[i].ux;
        cur_uy_a[hit_counter] = (*cur_seg)[i].uy;
        cur_chi2_a[hit_counter] = (*cur_seg)[i].chi2;

        whichseg[hit_counter] = i;
        hit1[hit_counter] = j;

        hit_counter += 1;
        if (hit_counter == 4) {
          calculateKappaTangents(
              x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a,
              dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a,
              dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a,
              dzdl_2_a, ddzdl_1_a, ddzdl_2_a, sinang_cut, cosang_diff_inv,
              cur_kappa_a, cur_dkappa_a, cur_ux_a, cur_uy_a, cur_chi2_a,
              chi2_a);

          for (unsigned int h = 0; h < hit_counter; ++h) {
            if ((chi2_a[h]) * inv_layer[l] < easy_chi2_cut) {
              temp_segment.chi2 = chi2_a[h];
              temp_segment.ux = ux_end_a[h];
              temp_segment.uy = uy_end_a[h];
              temp_segment.kappa = kappa_a[h];
              if (temp_segment.kappa > top_range.max_k) {
                continue;
              }
              temp_segment.dkappa = dkappa_a[h];
              for (unsigned int ll = 0; ll < l; ++ll) {
                temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
              }
              temp_segment.hits[l] = hit1[h];
              unsigned int outer_layer =
                  layer_sorted[l][temp_segment.hits[l]].get_layer();
              temp_segment.n_hits = l + 1;
              if ((n_layers - (l + 1)) <= allowed_missing) {
                complete_segments.push_back(temp_segment);
              }
              if ((outer_layer - l) > allowed_missing) {
                continue;
              }
              if (next_seg->size() == nextseg_size) {
                next_seg->push_back(temp_segment);
                nextseg_size += 1;
              } else {
                (*next_seg)[nextseg_size] = temp_segment;
                nextseg_size += 1;
              }
            }
          }
          hit_counter = 0;
        }
      }
    }
    if (hit_counter != 0) {
      calculateKappaTangents(
          x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a,
          dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
          ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a,
          ddzdl_2_a, sinang_cut, cosang_diff_inv, cur_kappa_a, cur_dkappa_a,
          cur_ux_a, cur_uy_a, cur_chi2_a, chi2_a);

      for (unsigned int h = 0; h < hit_counter; ++h) {
        if ((chi2_a[h]) * inv_layer[l] < easy_chi2_cut) {
          temp_segment.chi2 = chi2_a[h];
          temp_segment.ux = ux_end_a[h];
          temp_segment.uy = uy_end_a[h];
          temp_segment.kappa = kappa_a[h];
          if (temp_segment.kappa > top_range.max_k) {
            continue;
          }

          temp_segment.dkappa = dkappa_a[h];
          for (unsigned int ll = 0; ll < l; ++ll) {
            temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
          }
          temp_segment.hits[l] = hit1[h];
          unsigned int outer_layer =
              layer_sorted[l][temp_segment.hits[l]].get_layer();
          temp_segment.n_hits = l + 1;
          if ((n_layers - (l + 1)) <= allowed_missing) {
            complete_segments.push_back(temp_segment);
          }
          if ((outer_layer - l) > allowed_missing) {
            continue;
          }
          if (next_seg->size() == nextseg_size) {
            next_seg->push_back(temp_segment);
            nextseg_size += 1;
          } else {
            (*next_seg)[nextseg_size] = temp_segment;
            nextseg_size += 1;
          }
        }
      }
      hit_counter = 0;
    }
    swap(cur_seg, next_seg);
    swap(curseg_size, nextseg_size);
  }

  for (unsigned int i = 0; i < complete_segments.size(); ++i) {
    if (cur_seg->size() == curseg_size) {
      cur_seg->push_back(complete_segments[i]);
      curseg_size += 1;
    } else {
      (*cur_seg)[curseg_size] = complete_segments[i];
      curseg_size += 1;
    }
  }

  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
  CAtime += (time2 - time1);
  SimpleTrack3D temp_track;
  temp_track.hits.assign(n_layers, SimpleHit3D());
  vector<SimpleHit3D> temp_hits;
  for (unsigned int i = 0, sizei = curseg_size; i < sizei; ++i) {
    temp_track.hits.assign((*cur_seg)[i].n_hits, SimpleHit3D());

    temp_comb.assign((*cur_seg)[i].n_hits, 0);
    for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
      temp_comb[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].get_id();
    }
    sort(temp_comb.begin(), temp_comb.end());
    set<vector<unsigned int> >::iterator it = combos.find(temp_comb);
    if (it != combos.end()) {
      continue;
    }
    if (combos.size() > 10000) {
      combos.clear();
    }
    combos.insert(temp_comb);

    for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
      temp_track.hits[l] = layer_sorted[l][(*cur_seg)[i].hits[l]];
    }

    gettimeofday(&t1, NULL);

    // unsigned int layer_out = temp_track.hits.size()-1;
    // unsigned int layer_mid = layer_out/2;
    // temp_track_3hits.hits[0] = temp_track.hits[0];
    // temp_track_3hits.hits[1] = temp_track.hits[layer_mid];
    // temp_track_3hits.hits[2] = temp_track.hits[layer_out];

    float init_chi2 = fitTrack(temp_track);

    if (init_chi2 > fast_chi2_cut_max) {
      if (init_chi2 > fast_chi2_cut_par0 +
                          fast_chi2_cut_par1 / kappaToPt(temp_track.kappa)) {
        gettimeofday(&t2, NULL);
        time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
        time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
        KALtime += (time2 - time1);
        continue;
      }
    }
    HelixKalmanState state;
    state.phi = temp_track.phi;
    if (state.phi < 0.) {
      state.phi += 2. * M_PI;
    }
    state.d = temp_track.d;
    state.kappa = temp_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track.z0;
    state.dzdl = temp_track.dzdl;
    state.C = Matrix<float, 5, 5>::Zero(5, 5);
    state.C(0, 0) = pow(0.01, 2.);
    state.C(1, 1) = pow(0.01, 2.);
    state.C(2, 2) = pow(0.01 * state.nu, 2.);
    state.C(3, 3) = pow(0.05, 2.);
    state.C(4, 4) = pow(0.05, 2.);
    state.chi2 = 0.;
    state.position = 0;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;

    for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
      kalman->addHit(temp_track.hits[h], state);
      nfits += 1;
    }

    // fudge factor for non-gaussian hit sizes
    state.C *= 3.;
    state.chi2 *= 6.;

    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    KALtime += (time2 - time1);

    if (!(temp_track.kappa == temp_track.kappa)) {
      continue;
    }
    if (temp_track.kappa > top_range.max_k) {
      continue;
    }
    if (!(state.chi2 == state.chi2)) {
      continue;
    }
    if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > chi2_cut) {
      continue;
    }

    if (cut_on_dca == true) {
      if (fabs(temp_track.d) > dca_cut) {
        continue;
      }
      if (fabs(temp_track.z0) > dca_cut) {
        continue;
      }
    }

    tracks.push_back(temp_track);
    track_states.push_back(state);
    if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
        (temp_track.hits.size() >= n_removal_hits)) {
      for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
        (*hit_used)[temp_track.hits[i].get_id()] = true;
      }
    }
  }
}

// // void sPHENIXTracker::findSeededTracksbySegments(vector<SimpleTrack3D>& seeds,
// //                                                 vector<SimpleHit3D>& hits,
// //                                                 vector<SimpleTrack3D>& tracks,
// //                                                 const HelixRange& range) {
// //   unsigned int allowed_missing = n_layers - seed_layer - req_layers;

// //   for (unsigned int l = 0; l < n_layers; ++l) {
// //     layer_sorted[l].clear();
// //   }
// //   for (unsigned int i = 0; i < hits.size(); ++i) {
// //     unsigned int min = (hits[i].get_layer() - allowed_missing);
// //     if (allowed_missing > hits[i].get_layer()) {
// //       min = 0;
// //     }
// //     for (unsigned int l = min; l <= hits[i].get_layer(); l += 1) {
// //       layer_sorted[l].push_back(hits[i]);
// //     }
// //   }

// //   findSeededTracksbySegments_run(seeds, hits, tracks);
// // }

// // void sPHENIXTracker::findSeededTracksbySegments_run(
// //     vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits,
// //     vector<SimpleTrack3D>& tracks) {
// //   if (seeds.size() == 0) {
// //     return;
// //   }

// //   timeval t1, t2;
// //   double time1 = 0.;
// //   double time2 = 0.;

// //   gettimeofday(&t1, NULL);

// //   unsigned int first_new_layer = seed_layer;

// //   unsigned int allowed_missing = n_layers - seed_layer - req_layers;

// //   vector<TrackSegment>* cur_seg = &segments1;
// //   vector<TrackSegment>* next_seg = &segments2;
// //   unsigned int curseg_size = 0;
// //   unsigned int nextseg_size = 0;

// //   vector<TrackSegment> complete_segments;

// //   float cosang_diff = 1. - cosang_cut;
// //   float cosang_diff_inv = 1. / cosang_diff;
// //   float sinang_cut = sqrt(1. - cosang_cut * cosang_cut);
// //   float easy_chi2_cut = ca_chi2_cut;

// //   unsigned int hit_counter = 0;
// //   float x1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float x2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float x3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float y1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float y2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float y3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float z1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float z2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float z3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

// //   float dx1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dx2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dx3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dy1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dy2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dy3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dz1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dz2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dz3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

// //   float kappa_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dkappa_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

// //   float ux_mid_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float uy_mid_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float ux_end_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float uy_end_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

// //   float dzdl_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float dzdl_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float ddzdl_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
// //   float ddzdl_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

// //   TrackSegment temp_segment;
// //   temp_segment.hits.assign(n_layers, 0);
// //   unsigned int whichhit[4];
// //   unsigned int whichseed[4];
// //   for (unsigned int seed = 0, seedsize = seeds.size(); seed < seedsize;
// //        ++seed) {
// //     unsigned int firsthit = seeds[seed].hits.size();

// //     x1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_x();
// //     y1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_y();
// //     z1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_z();
// //     x2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_x();
// //     y2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_y();
// //     z2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_z();
// //     x3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_x();
// //     y3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_y();
// //     z3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_z();

// //     dx1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_ex();
// //     dy1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_ey();
// //     dz1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_ez();
// //     dx2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_ex();
// //     dy2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_ey();
// //     dz2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_ez();
// //     dx3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_ex();
// //     dy3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_ey();
// //     dz3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_ez();

// //     whichseed[hit_counter] = seed;
// //     hit_counter += 1;
// //     if (hit_counter == 4) {
// //       calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
// //                              z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
// //                              dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a,
// //                              uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a,
// //                              ddzdl_1_a, ddzdl_2_a);

// //       for (unsigned int h = 0; h < hit_counter; ++h) {
// //         temp_segment.chi2 = 0.;
// //         temp_segment.ux = ux_end_a[h];
// //         temp_segment.uy = uy_end_a[h];
// //         temp_segment.kappa = kappa_a[h];
// //         temp_segment.dkappa = dkappa_a[h];
// //         temp_segment.seed = whichseed[h];
// //         if (next_seg->size() == nextseg_size) {
// //           next_seg->push_back(temp_segment);
// //           nextseg_size += 1;
// //         } else {
// //           (*next_seg)[nextseg_size] = temp_segment;
// //           nextseg_size += 1;
// //         }
// //       }

// //       hit_counter = 0;
// //     }
// //   }
// //   if (hit_counter != 0) {
// //     calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a,
// //                            dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a,
// //                            dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a,
// //                            ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a,
// //                            ddzdl_2_a);

// //     for (unsigned int h = 0; h < hit_counter; ++h) {
// //       temp_segment.chi2 = 0.;
// //       temp_segment.ux = ux_end_a[h];
// //       temp_segment.uy = uy_end_a[h];
// //       temp_segment.kappa = kappa_a[h];
// //       temp_segment.dkappa = dkappa_a[h];
// //       temp_segment.seed = whichseed[h];
// //       if (next_seg->size() == nextseg_size) {
// //         next_seg->push_back(temp_segment);
// //         nextseg_size += 1;
// //       } else {
// //         (*next_seg)[nextseg_size] = temp_segment;
// //         nextseg_size += 1;
// //       }
// //     }

// //     hit_counter = 0;
// //   }
// //   swap(cur_seg, next_seg);
// //   swap(curseg_size, nextseg_size);
// //   unsigned int whichseg[4];
// //   for (unsigned int l = first_new_layer; l < n_layers; ++l) {
// //     nextseg_size = 0;
// //     for (unsigned int j = 0, sizej = curseg_size; j < sizej; ++j) {
// //       unsigned int firsthit = seeds[(*cur_seg)[j].seed].hits.size();
// //       for (unsigned int i = 0, sizei = layer_sorted[l].size(); i < sizei; ++i) {
// //         if ((l - 2) <= (first_new_layer - 1)) {
// //           x1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 2 - (first_new_layer - firsthit)]
// // 	    .get_x();
// //           y1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 2 - (first_new_layer - firsthit)]
// // 	    .get_y();
// //           z1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 2 - (first_new_layer - firsthit)]
// // 	    .get_z();
// //           dx1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 2 - (first_new_layer - firsthit)]
// // 	    .get_ex();
// //           dy1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 2 - (first_new_layer - firsthit)]
// // 	    .get_ey();
// //           dz1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 2 - (first_new_layer - firsthit)]
// // 	    .get_ez();
// //         } else {
// //           x1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_x();
// //           y1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_y();
// //           z1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_z();
// //           dx1_a[hit_counter] =
// //               layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_ex();
// //           dy1_a[hit_counter] =
// //               layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_ey();
// //           dz1_a[hit_counter] =
// //               layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_ez();
// //         }
// //         if ((l - 1) <= (first_new_layer - 1)) {
// //           x2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 1 - (first_new_layer - firsthit)]
// // 	    .get_x();
// //           y2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 1 - (first_new_layer - firsthit)]
// // 	    .get_y();
// //           z2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 1 - (first_new_layer - firsthit)]
// // 	    .get_z();
// //           dx2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 1 - (first_new_layer - firsthit)]
// // 	    .get_ex();
// //           dy2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 1 - (first_new_layer - firsthit)]
// // 	    .get_ey();
// //           dz2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
// // 	    .hits[l - 1 - (first_new_layer - firsthit)]
// // 	    .get_ez();
// //         } else {
// //           x2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_x();
// //           y2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_y();
// //           z2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_z();
// //           dx2_a[hit_counter] =
// //               layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_ex();
// //           dy2_a[hit_counter] =
// //               layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_ey();
// //           dz2_a[hit_counter] =
// //               layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_ez();
// //         }
// //         x3_a[hit_counter] = layer_sorted[l][i].get_x();
// //         y3_a[hit_counter] = layer_sorted[l][i].get_y();
// //         z3_a[hit_counter] = layer_sorted[l][i].get_z();
// //         dx3_a[hit_counter] = layer_sorted[l][i].get_ex();
// //         dy3_a[hit_counter] = layer_sorted[l][i].get_ey();
// //         dz3_a[hit_counter] = layer_sorted[l][i].get_ez();

// //         whichhit[hit_counter] = i;
// //         whichseg[hit_counter] = j;
// //         hit_counter += 1;
// //         if (hit_counter == 4) {
// //           calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
// //                                  z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
// //                                  dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
// //                                  ux_mid_a, uy_mid_a, ux_end_a, uy_end_a,
// //                                  dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);

// //           for (unsigned int h = 0; h < hit_counter; ++h) {
// //             float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
// //             float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
// //             dk += sinang_cut * kappa_a[h];
// //             float chi2_k = kdiff * kdiff / (dk * dk);
// //             float cos_scatter = (*cur_seg)[whichseg[h]].ux * ux_mid_a[h] +
// //                                 (*cur_seg)[whichseg[h]].uy * uy_mid_a[h];
// //             float chi2_ang = (1. - cos_scatter) * (1. - cos_scatter) *
// //                              cosang_diff_inv * cosang_diff_inv;
// //             float chi2_dzdl =
// //                 (dzdl_1_a[h] - dzdl_2_a[h]) /
// //                 (ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
// //             chi2_dzdl *= chi2_dzdl;
// //             if (((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl) /
// //                     ((float)l - 2.) <
// //                 easy_chi2_cut) {
// //               temp_segment.chi2 =
// //                   (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
// //               temp_segment.ux = ux_end_a[h];
// //               temp_segment.uy = uy_end_a[h];
// //               temp_segment.kappa = kappa_a[h];
// //               temp_segment.dkappa = dkappa_a[h];
// //               for (unsigned int ll = 0; ll < l; ++ll) {
// //                 temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
// //               }
// //               temp_segment.hits[l] = whichhit[h];
// //               temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
// //               unsigned int outer_layer =
// //                   layer_sorted[l][temp_segment.hits[l]].get_layer();
// //               temp_segment.n_hits = l + 1;
// //               if ((n_layers - (l + 1)) <= allowed_missing) {
// //                 complete_segments.push_back(temp_segment);
// //               }
// //               if ((outer_layer - l) > allowed_missing) {
// //                 continue;
// //               }
// //               if (next_seg->size() == nextseg_size) {
// //                 next_seg->push_back(temp_segment);
// //                 nextseg_size += 1;
// //               } else {
// //                 (*next_seg)[nextseg_size] = temp_segment;
// //                 nextseg_size += 1;
// //               }
// //             }
// //           }
// //           hit_counter = 0;
// //         }
// //       }
// //     }
// //     if (hit_counter != 0) {
// //       calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
// //                              z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
// //                              dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a,
// //                              uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a,
// //                              ddzdl_1_a, ddzdl_2_a);

// //       for (unsigned int h = 0; h < hit_counter; ++h) {
// //         float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
// //         float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
// //         dk += sinang_cut * kappa_a[h];
// //         float chi2_k = kdiff * kdiff / (dk * dk);
// //         float cos_scatter = (*cur_seg)[whichseg[h]].ux * ux_mid_a[h] +
// //                             (*cur_seg)[whichseg[h]].uy * uy_mid_a[h];
// //         float chi2_ang = (1. - cos_scatter) * (1. - cos_scatter) *
// //                          cosang_diff_inv * cosang_diff_inv;
// //         float chi2_dzdl =
// //             (dzdl_1_a[h] - dzdl_2_a[h]) /
// //             (ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
// //         chi2_dzdl *= chi2_dzdl;
// //         if (((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl) /
// //                 ((float)l - 2.) <
// //             easy_chi2_cut) {
// //           temp_segment.chi2 =
// //               (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
// //           temp_segment.ux = ux_end_a[h];
// //           temp_segment.uy = uy_end_a[h];
// //           temp_segment.kappa = kappa_a[h];
// //           temp_segment.dkappa = dkappa_a[h];
// //           for (unsigned int ll = 0; ll < l; ++ll) {
// //             temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
// //           }
// //           temp_segment.hits[l] = whichhit[h];
// //           temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
// //           unsigned int outer_layer =
// //               layer_sorted[l][temp_segment.hits[l]].get_layer();
// //           temp_segment.n_hits = l + 1;
// //           if ((n_layers - (l + 1)) <= allowed_missing) {
// //             complete_segments.push_back(temp_segment);
// //           }
// //           if ((outer_layer - l) > allowed_missing) {
// //             continue;
// //           }
// //           if (next_seg->size() == nextseg_size) {
// //             next_seg->push_back(temp_segment);
// //             nextseg_size += 1;
// //           } else {
// //             (*next_seg)[nextseg_size] = temp_segment;
// //             nextseg_size += 1;
// //           }
// //         }
// //       }
// //       hit_counter = 0;
// //     }

// //     swap(cur_seg, next_seg);
// //     swap(curseg_size, nextseg_size);
// //   }

// //   for (unsigned int i = 0; i < complete_segments.size(); ++i) {
// //     if (cur_seg->size() == curseg_size) {
// //       cur_seg->push_back(complete_segments[i]);
// //       curseg_size += 1;
// //     } else {
// //       (*cur_seg)[curseg_size] = complete_segments[i];
// //       curseg_size += 1;
// //     }
// //   }

// //   gettimeofday(&t2, NULL);
// //   time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
// //   time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
// //   CAtime += (time2 - time1);

// //   gettimeofday(&t1, NULL);

// //   SimpleTrack3D temp_track;
// //   temp_track.hits.assign(n_layers, SimpleHit3D());
// //   for (unsigned int i = 0, sizei = curseg_size; i < sizei; ++i) {
// //     temp_track.hits.clear();
// //     for (unsigned int l = 0; l < seeds[(*cur_seg)[i].seed].hits.size(); ++l) {
// //       temp_track.hits.push_back(seeds[(*cur_seg)[i].seed].hits[l]);
// //     }
// //     for (unsigned int l = seed_layer; l < (*cur_seg)[i].n_hits; ++l) {
// //       temp_track.hits.push_back(layer_sorted[l][(*cur_seg)[i].hits[l]]);
// //     }

// //     HelixKalmanState state;

// //     fitTrack(temp_track);
// //     state.C = Matrix<float, 5, 5>::Zero(5, 5);
// //     state.phi = temp_track.phi;
// //     if (state.phi < 0.) {
// //       state.phi += 2. * M_PI;
// //     }
// //     state.d = temp_track.d;
// //     state.kappa = temp_track.kappa;
// //     state.nu = sqrt(state.kappa);
// //     state.z0 = temp_track.z0;
// //     state.dzdl = temp_track.dzdl;
// //     state.C(0, 0) = pow(0.1 * state.phi, 2.);
// //     state.C(1, 1) = pow(0.1 * state.d, 2.);
// //     state.C(2, 2) = pow(0.1 * state.nu, 2.);
// //     state.C(3, 3) = pow(0.1 * state.z0, 2.);
// //     state.C(4, 4) = pow(0.1 * state.dzdl, 2.);

// //     state.chi2 = 0.;
// //     state.position = 0;
// //     state.x_int = 0.;
// //     state.y_int = 0.;
// //     state.z_int = 0.;
// //     for (unsigned int l = 0; l < temp_track.hits.size(); ++l) {
// //       kalman->addHit(temp_track.hits[l], state);
// //       nfits += 1;
// //     }

// //     if (state.chi2 / (2. * ((float)(temp_track.hits.size() + 1)) - 5.) >
// //         chi2_cut) {
// //       continue;
// //     }
// //     if (state.chi2 != state.chi2) {
// //       continue;
// //     }

// //     temp_track.phi = state.phi;
// //     if (temp_track.phi < 0.) {
// //       temp_track.phi += 2. * M_PI;
// //     }
// //     if (temp_track.phi > 2. * M_PI) {
// //       temp_track.phi -= 2. * M_PI;
// //     }
// //     temp_track.d = state.d;
// //     temp_track.kappa = state.kappa;
// //     temp_track.z0 = state.z0;
// //     temp_track.dzdl = state.dzdl;

// //     if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
// //         (temp_track.hits.size() >= n_removal_hits)) {
// //       for (unsigned int l = (seeds[(*cur_seg)[i].seed].hits.size());
// //            l < temp_track.hits.size(); ++l) {
// //         (*hit_used)[temp_track.hits[l].get_id()] = true;
// //         temp_track.hits[l].set_id(index_mapping[temp_track.hits[l].get_id()]);
// //       }
// //     } else {
// //       for (unsigned int l = (seeds[(*cur_seg)[i].seed].hits.size());
// //            l < temp_track.hits.size(); ++l) {
// //         temp_track.hits[l].set_id(index_mapping[temp_track.hits[l].get_id()]);
// //       }
// //     }
// //     tracks.push_back(temp_track);
// //     track_states.push_back(state);
// //     seed_used[seeds[(*cur_seg)[i].seed].index] = true;
// //   }

// //   gettimeofday(&t2, NULL);
// //   time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
// //   time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
// //   KALtime += (time2 - time1);
// // }
