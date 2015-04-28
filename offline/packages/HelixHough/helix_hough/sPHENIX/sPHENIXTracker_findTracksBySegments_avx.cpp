#ifdef AVXHOUGH

#include "vector_math_inline_avx.h"
#include "sPHENIXTracker.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <sys/time.h>


using namespace std;
using namespace Eigen;


static inline void __attribute__((always_inline)) calculateKappaTangents_avx(float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a)
{
  static const __m256 two = {2., 2., 2., 2., 2., 2., 2., 2.};
  
  __m256 x1 = _mm256_load_ps(x1_a);
  __m256 x2 = _mm256_load_ps(x2_a);
  __m256 x3 = _mm256_load_ps(x3_a);
  __m256 y1 = _mm256_load_ps(y1_a);
  __m256 y2 = _mm256_load_ps(y2_a);
  __m256 y3 = _mm256_load_ps(y3_a);
  __m256 z1 = _mm256_load_ps(z1_a);
  __m256 z2 = _mm256_load_ps(z2_a);
  __m256 z3 = _mm256_load_ps(z3_a);
  
  __m256 dx1 = _mm256_load_ps(dx1_a);
  __m256 dx2 = _mm256_load_ps(dx2_a);
  __m256 dx3 = _mm256_load_ps(dx3_a);
  __m256 dy1 = _mm256_load_ps(dy1_a);
  __m256 dy2 = _mm256_load_ps(dy2_a);
  __m256 dy3 = _mm256_load_ps(dy3_a);
  __m256 dz1 = _mm256_load_ps(dz1_a);
  __m256 dz2 = _mm256_load_ps(dz2_a);
  __m256 dz3 = _mm256_load_ps(dz3_a);
  
  
  __m256 D12 = _mm256_sub_ps(x2, x1);
  D12 = _mm256_mul_ps(D12, D12);
  __m256 tmp1 = _mm256_sub_ps(y2, y1);
  tmp1 = _mm256_mul_ps(tmp1, tmp1);
  D12 = _mm256_add_ps(D12, tmp1);
  D12 = _vec256_sqrt_ps(D12);
  
  __m256 D23 = _mm256_sub_ps(x3, x2);
  D23 = _mm256_mul_ps(D23, D23);
  tmp1 = _mm256_sub_ps(y3, y2);
  tmp1 = _mm256_mul_ps(tmp1, tmp1);
  D23 = _mm256_add_ps(D23, tmp1);
  D23 = _vec256_sqrt_ps(D23);
  
  __m256 D31 = _mm256_sub_ps(x1, x3);
  D31 = _mm256_mul_ps(D31, D31);
  tmp1 = _mm256_sub_ps(y1, y3);
  tmp1 = _mm256_mul_ps(tmp1, tmp1);
  D31 = _mm256_add_ps(D31, tmp1);
  D31 = _vec256_sqrt_ps(D31);
  
  __m256 k = _mm256_mul_ps(D12, D23);
  k = _mm256_mul_ps(k, D31);
  k = _vec256_rec_ps(k);
  tmp1 = (D12 + D23 + D31)*(D23 + D31 - D12)*(D12 + D31 - D23)*(D12 + D23 - D31);
  tmp1 = _vec256_sqrt_ps(tmp1);
  k *= tmp1;
  
  __m256 tmp2 = _mm256_cmpgt_ps(tmp1, zero_256);
  tmp1 = _mm256_and_ps(tmp2, k);
  tmp2 = _mm256_andnot_ps(tmp2, zero_256);
  k = _mm256_xor_ps(tmp1, tmp2);
  
  _mm256_store_ps(kappa_a, k);
  __m256 k_inv = _vec256_rec_ps(k);
  
  __m256 D12_inv = _vec256_rec_ps(D12);
  __m256 D23_inv = _vec256_rec_ps(D23);
  __m256 D31_inv = _vec256_rec_ps(D31);
  
  __m256 dr1 = dx1*dx1 + dy1*dy1;
  dr1 = _vec256_sqrt_ps(dr1);
  __m256 dr2 = dx2*dx2 + dy2*dy2;
  dr2 = _vec256_sqrt_ps(dr2);
  __m256 dr3 = dx3*dx3 + dy3*dy3;
  dr3 = _vec256_sqrt_ps(dr3);
  
  __m256 dk1 = (dr1 + dr2)*D12_inv*D12_inv;
  __m256 dk2 = (dr2 + dr3)*D23_inv*D23_inv;
  __m256 dk = dk1 + dk2;
  _mm256_store_ps(dkappa_a, dk);
  
  __m256 ux12 = (x2 - x1)*D12_inv;
  __m256 uy12 = (y2 - y1)*D12_inv;
  __m256 ux23 = (x3 - x2)*D23_inv;
  __m256 uy23 = (y3 - y2)*D23_inv;
  __m256 ux13 = (x3 - x1)*D31_inv;
  __m256 uy13 = (y3 - y1)*D31_inv;
  
  __m256 cosalpha = ux12*ux13 + uy12*uy13;
  __m256 sinalpha = ux13*uy12 - ux12*uy13;
  
  __m256 ux_mid = ux23*cosalpha - uy23*sinalpha;
  __m256 uy_mid = ux23*sinalpha + uy23*cosalpha;
  _mm256_store_ps(ux_mid_a, ux_mid);
  _mm256_store_ps(uy_mid_a, uy_mid);
  
  __m256 ux_end = ux23*cosalpha + uy23*sinalpha;
  __m256 uy_end = uy23*cosalpha - ux23*sinalpha;
  
  _mm256_store_ps(ux_end_a, ux_end);
  _mm256_store_ps(uy_end_a, uy_end);
  
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m256 v = one_256 - sinalpha*sinalpha;
  v = _vec256_sqrt_ps(v);
  v += one_256;
  v = _vec256_rec_ps(v);
  v *= sinalpha;
  __m256 s2 = _vec256_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm256_cmpgt_ps(k, zero_256);
  tmp2 = _mm256_and_ps(tmp1, s2);
  tmp1 = _mm256_andnot_ps(tmp1, D23);
  s2 = _mm256_xor_ps(tmp1, tmp2);
  
  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m256 del_z_2 = z3 - z2;
  __m256 dzdl_2 = s2*s2 + del_z_2*del_z_2;
  dzdl_2 = _vec256_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m256 ddzdl_2 = (dz2 + dz3)*D23_inv;
  _mm256_store_ps(dzdl_2_a, dzdl_2);
  _mm256_store_ps(ddzdl_2_a, ddzdl_2);
  
  sinalpha = ux13*uy23 - ux23*uy13;
  v = one_256 - sinalpha*sinalpha;
  v = _vec256_sqrt_ps(v);
  v += one_256;
  v = _vec256_rec_ps(v);
  v *= sinalpha;
  __m256 s1 = _vec256_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm256_cmpgt_ps(k, zero_256);
  tmp2 = _mm256_and_ps(tmp1, s1);
  tmp1 = _mm256_andnot_ps(tmp1, D12);
  s1 = _mm256_xor_ps(tmp1, tmp2);
  
  __m256 del_z_1 = z2 - z1;
  __m256 dzdl_1 = s1*s1 + del_z_1*del_z_1;
  dzdl_1 = _vec256_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m256 ddzdl_1 = (dz1 + dz2)*D12_inv;
  _mm256_store_ps(dzdl_1_a, dzdl_1);
  _mm256_store_ps(ddzdl_1_a, ddzdl_1);
}


static inline void __attribute__((always_inline)) calculateKappaTangents_avx(float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a, float sinang_cut, float cosang_diff_inv, float* cur_kappa_a, float* cur_dkappa_a, float* cur_ux_a, float* cur_uy_a, float* cur_chi2_a, float* chi2_a)
{
  static const __m256 two = {2., 2., 2., 2., 2., 2., 2., 2.};
  
  __m256 x1 = _mm256_load_ps(x1_a);
  __m256 x2 = _mm256_load_ps(x2_a);
  __m256 x3 = _mm256_load_ps(x3_a);
  __m256 y1 = _mm256_load_ps(y1_a);
  __m256 y2 = _mm256_load_ps(y2_a);
  __m256 y3 = _mm256_load_ps(y3_a);
  __m256 z1 = _mm256_load_ps(z1_a);
  __m256 z2 = _mm256_load_ps(z2_a);
  __m256 z3 = _mm256_load_ps(z3_a);
  
  __m256 dx1 = _mm256_load_ps(dx1_a);
  __m256 dx2 = _mm256_load_ps(dx2_a);
  __m256 dx3 = _mm256_load_ps(dx3_a);
  __m256 dy1 = _mm256_load_ps(dy1_a);
  __m256 dy2 = _mm256_load_ps(dy2_a);
  __m256 dy3 = _mm256_load_ps(dy3_a);
  __m256 dz1 = _mm256_load_ps(dz1_a);
  __m256 dz2 = _mm256_load_ps(dz2_a);
  __m256 dz3 = _mm256_load_ps(dz3_a);
  
  
  __m256 D12 = _mm256_sub_ps(x2, x1);
  D12 = _mm256_mul_ps(D12, D12);
  __m256 tmp1 = _mm256_sub_ps(y2, y1);
  tmp1 = _mm256_mul_ps(tmp1, tmp1);
  D12 = _mm256_add_ps(D12, tmp1);
  D12 = _vec256_sqrt_ps(D12);
  
  __m256 D23 = _mm256_sub_ps(x3, x2);
  D23 = _mm256_mul_ps(D23, D23);
  tmp1 = _mm256_sub_ps(y3, y2);
  tmp1 = _mm256_mul_ps(tmp1, tmp1);
  D23 = _mm256_add_ps(D23, tmp1);
  D23 = _vec256_sqrt_ps(D23);
  
  __m256 D31 = _mm256_sub_ps(x1, x3);
  D31 = _mm256_mul_ps(D31, D31);
  tmp1 = _mm256_sub_ps(y1, y3);
  tmp1 = _mm256_mul_ps(tmp1, tmp1);
  D31 = _mm256_add_ps(D31, tmp1);
  D31 = _vec256_sqrt_ps(D31);
  
  __m256 k = _mm256_mul_ps(D12, D23);
  k = _mm256_mul_ps(k, D31);
  k = _vec256_rec_ps(k);
  tmp1 = (D12 + D23 + D31)*(D23 + D31 - D12)*(D12 + D31 - D23)*(D12 + D23 - D31);
  tmp1 = _vec256_sqrt_ps(tmp1);
  k *= tmp1;
  
  __m256 tmp2 = _mm256_cmpgt_ps(tmp1, zero_256);
  tmp1 = _mm256_and_ps(tmp2, k);
  tmp2 = _mm256_andnot_ps(tmp2, zero_256);
  k = _mm256_xor_ps(tmp1, tmp2);
  
  _mm256_store_ps(kappa_a, k);
  __m256 k_inv = _vec256_rec_ps(k);
  
  __m256 D12_inv = _vec256_rec_ps(D12);
  __m256 D23_inv = _vec256_rec_ps(D23);
  __m256 D31_inv = _vec256_rec_ps(D31);
  
  __m256 dr1 = dx1*dx1 + dy1*dy1;
  dr1 = _vec256_sqrt_ps(dr1);
  __m256 dr2 = dx2*dx2 + dy2*dy2;
  dr2 = _vec256_sqrt_ps(dr2);
  __m256 dr3 = dx3*dx3 + dy3*dy3;
  dr3 = _vec256_sqrt_ps(dr3);
  
  __m256 dk1 = (dr1 + dr2)*D12_inv*D12_inv;
  __m256 dk2 = (dr2 + dr3)*D23_inv*D23_inv;
  __m256 dk = dk1 + dk2;
  _mm256_store_ps(dkappa_a, dk);
  
  __m256 ux12 = (x2 - x1)*D12_inv;
  __m256 uy12 = (y2 - y1)*D12_inv;
  __m256 ux23 = (x3 - x2)*D23_inv;
  __m256 uy23 = (y3 - y2)*D23_inv;
  __m256 ux13 = (x3 - x1)*D31_inv;
  __m256 uy13 = (y3 - y1)*D31_inv;
  
  __m256 cosalpha = ux12*ux13 + uy12*uy13;
  __m256 sinalpha = ux13*uy12 - ux12*uy13;
  
  __m256 ux_mid = ux23*cosalpha - uy23*sinalpha;
  __m256 uy_mid = ux23*sinalpha + uy23*cosalpha;
  _mm256_store_ps(ux_mid_a, ux_mid);
  _mm256_store_ps(uy_mid_a, uy_mid);
  
  __m256 ux_end = ux23*cosalpha + uy23*sinalpha;
  __m256 uy_end = uy23*cosalpha - ux23*sinalpha;
  
  _mm256_store_ps(ux_end_a, ux_end);
  _mm256_store_ps(uy_end_a, uy_end);
  
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m256 v = one_256 - sinalpha*sinalpha;
  v = _vec256_sqrt_ps(v);
  v += one_256;
  v = _vec256_rec_ps(v);
  v *= sinalpha;
  __m256 s2 = _vec256_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm256_cmpgt_ps(k, zero_256);
  tmp2 = _mm256_and_ps(tmp1, s2);
  tmp1 = _mm256_andnot_ps(tmp1, D23);
  s2 = _mm256_xor_ps(tmp1, tmp2);
  
  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m256 del_z_2 = z3 - z2;
  __m256 dzdl_2 = s2*s2 + del_z_2*del_z_2;
  dzdl_2 = _vec256_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m256 ddzdl_2 = (dz2 + dz3)*D23_inv;
  _mm256_store_ps(dzdl_2_a, dzdl_2);
  _mm256_store_ps(ddzdl_2_a, ddzdl_2);
  
  sinalpha = ux13*uy23 - ux23*uy13;
  v = one_256 - sinalpha*sinalpha;
  v = _vec256_sqrt_ps(v);
  v += one_256;
  v = _vec256_rec_ps(v);
  v *= sinalpha;
  __m256 s1 = _vec256_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm256_cmpgt_ps(k, zero_256);
  tmp2 = _mm256_and_ps(tmp1, s1);
  tmp1 = _mm256_andnot_ps(tmp1, D12);
  s1 = _mm256_xor_ps(tmp1, tmp2);
  
  __m256 del_z_1 = z2 - z1;
  __m256 dzdl_1 = s1*s1 + del_z_1*del_z_1;
  dzdl_1 = _vec256_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m256 ddzdl_1 = (dz1 + dz2)*D12_inv;
  _mm256_store_ps(dzdl_1_a, dzdl_1);
  _mm256_store_ps(ddzdl_1_a, ddzdl_1);
  
  __m256 c_dk = _mm256_load_ps(cur_dkappa_a);
  __m256 c_k = _mm256_load_ps(cur_kappa_a);
  __m256 c_ux = _mm256_load_ps(cur_ux_a);
  __m256 c_uy = _mm256_load_ps(cur_uy_a);
  __m256 c_chi2 = _mm256_load_ps(cur_chi2_a);
  __m256 sinang = _mm256_load1_ps(sinang_cut);
  __m256 cosdiff = _mm256_load1_ps(cosang_diff_inv);
  
  __m256 kdiff = c_k - k;
  __m256 n_dk = c_dk + dk + sinang*k;
  __m256 chi2_k = kdiff*kdiff/(n_dk*n_dk);
  __m256 cos_scatter = c_ux*ux_mid + c_uy*uy_mid;
  __m256 chi2_ang = (one_256-cos_scatter)*(one_256-cos_scatter)*cosdiff*cosdiff;
  tmp1 = dzdl_1*sinang;
  _vec256_fabs_ps(tmp1);
  __m256 chi2_dzdl = (dzdl_1 - dzdl_2)/(ddzdl_1 + ddzdl_2 + tmp1);
  chi2_dzdl *= chi2_dzdl;
  chi2_dzdl *= one_o_2_256;
  
  __m256 n_chi2 = c_chi2 + chi2_ang + chi2_k + chi2_dzdl;
  _mm256_store_ps(chi2_a, n_chi2);
}


void sPHENIXTracker::findTracksBySegments_avx(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
  unsigned int allowed_missing = n_layers - req_layers;
  for(unsigned int l=0;l<n_layers;++l)
  {
    layer_sorted_1[findtracks_bin][l].clear();
  }
  for(unsigned int i=0;i<hits.size();++i)
  {
    unsigned int min = (hits[i].layer - allowed_missing);
    if(allowed_missing > hits[i].layer){min = 0;}
    for(unsigned int l=min;l<=hits[i].layer;l+=1){layer_sorted_1[findtracks_bin][l].push_back(hits[i]);}
  }
  findtracks_bin += 1;
  if(findtracks_bin != 4){return;}
  findtracks_bin = 0;
  findTracksBySegments_avx_run(tracks);
}


void sPHENIXTracker::findTracksBySegments_avx_run(vector<SimpleTrack3D>& tracks)
{
  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;
  unsigned int curseg_size = 0;
  unsigned int nextseg_size = 0;
  
  vector<TrackSegment> complete_segments;
  
  unsigned int allowed_missing = n_layers - req_layers;
  
  timeval t1,t2;
  double time1=0.;
  double time2=0.;
  
  gettimeofday(&t1, NULL);
  
  float cosang_diff = 1. - cosang_cut;
  float cosang_diff_inv = 1./cosang_diff;
  float sinang_cut = sqrt(1. - cosang_cut*cosang_cut);
  float easy_chi2_cut = ca_chi2_cut;
  
  vector<float> inv_layer;inv_layer.assign(n_layers, 1.);
  for(unsigned int l=3;l<n_layers;++l)
  {
    inv_layer[l] = 1./(((float)l) - 2.);
  }
  
  unsigned int hit_counter = 0;
  float x1_a[8] __attribute__((aligned(32)));
  float x2_a[8] __attribute__((aligned(32)));
  float x3_a[8] __attribute__((aligned(32)));
  float y1_a[8] __attribute__((aligned(32)));
  float y2_a[8] __attribute__((aligned(32)));
  float y3_a[8] __attribute__((aligned(32)));
  float z1_a[8] __attribute__((aligned(32)));
  float z2_a[8] __attribute__((aligned(32)));
  float z3_a[8] __attribute__((aligned(32)));
  
  float dx1_a[8] __attribute__((aligned(32)));
  float dx2_a[8] __attribute__((aligned(32)));
  float dx3_a[8] __attribute__((aligned(32)));
  float dy1_a[8] __attribute__((aligned(32)));
  float dy2_a[8] __attribute__((aligned(32)));
  float dy3_a[8] __attribute__((aligned(32)));
  float dz1_a[8] __attribute__((aligned(32)));
  float dz2_a[8] __attribute__((aligned(32)));
  float dz3_a[8] __attribute__((aligned(32)));
  
  float kappa_a[8] __attribute__((aligned(32)));
  float dkappa_a[8] __attribute__((aligned(32)));
  
  float ux_mid_a[8] __attribute__((aligned(32)));
  float uy_mid_a[8] __attribute__((aligned(32)));
  float ux_end_a[8] __attribute__((aligned(32)));
  float uy_end_a[8] __attribute__((aligned(32)));
  
  float dzdl_1_a[8] __attribute__((aligned(32)));
  float dzdl_2_a[8] __attribute__((aligned(32)));
  float ddzdl_1_a[8] __attribute__((aligned(32)));
  float ddzdl_2_a[8] __attribute__((aligned(32)));
  
  float cur_kappa_a[8] __attribute__((aligned(32)));
  float cur_dkappa_a[8] __attribute__((aligned(32)));
  float cur_ux_a[8] __attribute__((aligned(32)));
  float cur_uy_a[8] __attribute__((aligned(32)));
  float cur_chi2_a[8] __attribute__((aligned(32)));
  float chi2_a[8] __attribute__((aligned(32)));
  
  unsigned int hit1[8];
  unsigned int hit2[8];
  unsigned int hit3[8];
  unsigned int bins[8];
  
  TrackSegment temp_segment;temp_segment.hits.assign(n_layers, 0);
  // make segments out of first 3 layers
  for(unsigned int b=0;b<4;++b)
  {
    for(unsigned int i=0,sizei=layer_sorted_1[b][0].size();i<sizei;++i)
    {
      for(unsigned int j=0,sizej=layer_sorted_1[b][1].size();j<sizej;++j)
      {
        for(unsigned int k=0,sizek=layer_sorted_1[b][2].size();k<sizek;++k)
        {
          if( (layer_sorted_1[b][0][i].layer >= layer_sorted_1[b][1][j].layer) || (layer_sorted_1[b][1][j].layer >= layer_sorted_1[b][2][k].layer) ){continue;}
          
          x1_a[hit_counter] = layer_sorted_1[b][0][i].x;
          y1_a[hit_counter] = layer_sorted_1[b][0][i].y;
          z1_a[hit_counter] = layer_sorted_1[b][0][i].z;
          dx1_a[hit_counter] = layer_sorted_1[b][0][i].dx;
          dy1_a[hit_counter] = layer_sorted_1[b][0][i].dy;
          dz1_a[hit_counter] = layer_sorted_1[b][0][i].dz;
          
          x2_a[hit_counter] = layer_sorted_1[b][1][j].x;
          y2_a[hit_counter] = layer_sorted_1[b][1][j].y;
          z2_a[hit_counter] = layer_sorted_1[b][1][j].z;
          dx2_a[hit_counter] = layer_sorted_1[b][1][j].dx;
          dy2_a[hit_counter] = layer_sorted_1[b][1][j].dy;
          dz2_a[hit_counter] = layer_sorted_1[b][1][j].dz;
          
          x3_a[hit_counter] = layer_sorted_1[b][2][k].x;
          y3_a[hit_counter] = layer_sorted_1[b][2][k].y;
          z3_a[hit_counter] = layer_sorted_1[b][2][k].z;
          dx3_a[hit_counter] = layer_sorted_1[b][2][k].dx;
          dy3_a[hit_counter] = layer_sorted_1[b][2][k].dy;
          dz3_a[hit_counter] = layer_sorted_1[b][2][k].dz;
          
          hit1[hit_counter] = i;
          hit2[hit_counter] = j;
          hit3[hit_counter] = k;
          bins[hit_counter] = b;
          
          hit_counter += 1;
          
          if(hit_counter == 8)
          {
            calculateKappaTangents_avx(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
            
            for(unsigned int h=0;h<hit_counter;++h)
            {
              temp_segment.chi2 = (dzdl_1_a[h] - dzdl_2_a[h])/(ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h]*sinang_cut));
              temp_segment.chi2 *= temp_segment.chi2;
              if(temp_segment.chi2 > 2.0){continue;}
              temp_segment.ux = ux_end_a[h];
              temp_segment.uy = uy_end_a[h];
              temp_segment.kappa = kappa_a[h];
              if(temp_segment.kappa > top_range.max_k){continue;}
              temp_segment.dkappa = dkappa_a[h];
              temp_segment.hits[0] = hit1[h];
              temp_segment.hits[1] = hit2[h];
              temp_segment.hits[2] = hit3[h];
              temp_segment.bin = bins[h];
              temp_segment.n_hits = 3;
              unsigned int outer_layer = layer_sorted_1[temp_segment.bin][2][temp_segment.hits[2]].layer;
              if( (outer_layer - 2) > allowed_missing ){continue;}
              if( (n_layers - 3) <= allowed_missing ){complete_segments.push_back(temp_segment);}
              if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
              else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
            }
            
            hit_counter=0;
          }
        }
      }
    }
  }
  if(hit_counter != 0)
  {
    calculateKappaTangents_avx(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
    
    for(unsigned int h=0;h<hit_counter;++h)
    {
      temp_segment.chi2 = (dzdl_1_a[h] - dzdl_2_a[h])/(ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h]*sinang_cut));
      temp_segment.chi2 *= temp_segment.chi2;
      if(temp_segment.chi2 > 2.0){continue;}
      temp_segment.ux = ux_end_a[h];
      temp_segment.uy = uy_end_a[h];
      temp_segment.kappa = kappa_a[h];
      if(temp_segment.kappa > top_range.max_k){continue;}
      temp_segment.dkappa = dkappa_a[h];
      temp_segment.hits[0] = hit1[h];
      temp_segment.hits[1] = hit2[h];
      temp_segment.hits[2] = hit3[h];
      temp_segment.bin = bins[h];
      temp_segment.n_hits = 3;
      unsigned int outer_layer = layer_sorted_1[temp_segment.bin][2][temp_segment.hits[2]].layer;
      if( (outer_layer - 2) > allowed_missing ){continue;}
      if( (n_layers - 3) <= allowed_missing ){complete_segments.push_back(temp_segment);}
      if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
      else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
    }
    
    hit_counter=0;
  }
  swap(cur_seg, next_seg);
  swap(curseg_size, nextseg_size);
  
  // add hits to segments layer-by-layer, cutting out bad segments
  unsigned int whichseg[8];
  for(unsigned int l=3;l<n_layers;++l)
  {
    if(l == (n_layers-1)){easy_chi2_cut*=0.25;}
    nextseg_size = 0;
    for(unsigned int i=0,sizei=curseg_size;i<sizei;++i)
    {
      unsigned int b = (*cur_seg)[i].bin;
      for(unsigned int j=0,sizej=layer_sorted_1[b][l].size();j<sizej;++j)
      {
        if( (layer_sorted_1[b][l-1][(*cur_seg)[i].hits[l-1]].layer >= layer_sorted_1[b][l][j].layer) ){continue;}
        
        x1_a[hit_counter] = layer_sorted_1[b][l-2][(*cur_seg)[i].hits[l-2]].x;
        y1_a[hit_counter] = layer_sorted_1[b][l-2][(*cur_seg)[i].hits[l-2]].y;
        z1_a[hit_counter] = layer_sorted_1[b][l-2][(*cur_seg)[i].hits[l-2]].z;
        x2_a[hit_counter] = layer_sorted_1[b][l-1][(*cur_seg)[i].hits[l-1]].x;
        y2_a[hit_counter] = layer_sorted_1[b][l-1][(*cur_seg)[i].hits[l-1]].y;
        z2_a[hit_counter] = layer_sorted_1[b][l-1][(*cur_seg)[i].hits[l-1]].z;
        x3_a[hit_counter] = layer_sorted_1[b][l][j].x;
        y3_a[hit_counter] = layer_sorted_1[b][l][j].y;
        z3_a[hit_counter] = layer_sorted_1[b][l][j].z;
        
        dx1_a[hit_counter] = layer_sorted_1[b][l-2][(*cur_seg)[i].hits[l-2]].dx;
        dy1_a[hit_counter] = layer_sorted_1[b][l-2][(*cur_seg)[i].hits[l-2]].dy;
        dz1_a[hit_counter] = layer_sorted_1[b][l-2][(*cur_seg)[i].hits[l-2]].dz;
        dx2_a[hit_counter] = layer_sorted_1[b][l-1][(*cur_seg)[i].hits[l-1]].dx;
        dy2_a[hit_counter] = layer_sorted_1[b][l-1][(*cur_seg)[i].hits[l-1]].dy;
        dz2_a[hit_counter] = layer_sorted_1[b][l-1][(*cur_seg)[i].hits[l-1]].dz;
        dx3_a[hit_counter] = layer_sorted_1[b][l][j].dx;
        dy3_a[hit_counter] = layer_sorted_1[b][l][j].dy;
        dz3_a[hit_counter] = layer_sorted_1[b][l][j].dz;
        
        cur_kappa_a[hit_counter] = (*cur_seg)[i].kappa;
        cur_dkappa_a[hit_counter] = (*cur_seg)[i].dkappa;
        cur_ux_a[hit_counter] = (*cur_seg)[i].ux;
        cur_uy_a[hit_counter] = (*cur_seg)[i].uy;
        cur_chi2_a[hit_counter] = (*cur_seg)[i].chi2;
        
        whichseg[hit_counter] = i;
        hit1[hit_counter] = j;
        
        hit_counter += 1;
        if(hit_counter == 8)
        {
          calculateKappaTangents_avx(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a, sinang_cut, cosang_diff_inv, cur_kappa_a, cur_dkappa_a, cur_ux_a, cur_uy_a, cur_chi2_a, chi2_a);
          
          for(unsigned int h=0;h<hit_counter;++h)
          {
            if( (chi2_a[h])*inv_layer[l] < easy_chi2_cut )
            {
              temp_segment.chi2 = chi2_a[h];
              temp_segment.ux = ux_end_a[h];
              temp_segment.uy = uy_end_a[h];
              temp_segment.kappa = kappa_a[h];
              if(temp_segment.kappa > top_range.max_k){continue;}
              temp_segment.dkappa = dkappa_a[h];
              temp_segment.bin = (*cur_seg)[whichseg[h]].bin;
              for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
              temp_segment.hits[l] = hit1[h];
              unsigned int outer_layer = layer_sorted_1[temp_segment.bin][l][temp_segment.hits[l]].layer;
              temp_segment.n_hits = l+1;
              if( (n_layers - (l+1)) <= allowed_missing ){complete_segments.push_back(temp_segment);}
              if( (outer_layer - l) > allowed_missing ){continue;}
              if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
              else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
            }
          }
          hit_counter=0;
        }
      }
    }
    if(hit_counter != 0)
    {
      calculateKappaTangents_avx(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a, sinang_cut, cosang_diff_inv, cur_kappa_a, cur_dkappa_a, cur_ux_a, cur_uy_a, cur_chi2_a, chi2_a);
      
      for(unsigned int h=0;h<hit_counter;++h)
      {
        if( (chi2_a[h])*inv_layer[l] < easy_chi2_cut )
        {
          temp_segment.chi2 = chi2_a[h];
          temp_segment.ux = ux_end_a[h];
          temp_segment.uy = uy_end_a[h];
          temp_segment.kappa = kappa_a[h];
          if(temp_segment.kappa > top_range.max_k){continue;}
          temp_segment.dkappa = dkappa_a[h];
          temp_segment.bin = (*cur_seg)[whichseg[h]].bin;
          for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
          temp_segment.hits[l] = hit1[h];
          unsigned int outer_layer = layer_sorted_1[temp_segment.bin][l][temp_segment.hits[l]].layer;
          temp_segment.n_hits = l+1;
          if( (n_layers - (l+1)) <= allowed_missing ){complete_segments.push_back(temp_segment);}
          if( (outer_layer - l) > allowed_missing ){continue;}
          if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
          else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
        }
      }
      hit_counter=0;
    }
    swap(cur_seg, next_seg);
    swap(curseg_size, nextseg_size);
  }
  
  for(unsigned int i=0;i<complete_segments.size();++i)
  {
    if(cur_seg->size() == curseg_size){cur_seg->push_back(complete_segments[i]);curseg_size+=1;}
    else{(*cur_seg)[curseg_size] = complete_segments[i];curseg_size+=1;}
  }
  
  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
  CAtime += (time2 - time1);
  SimpleTrack3D temp_track;
  temp_track.hits.assign(n_layers, SimpleHit3D());
  SimpleTrack3D temp_track_3hits;
  temp_track_3hits.hits.assign(3, SimpleHit3D());
  vector<SimpleHit3D> temp_hits;
  for(unsigned int i=0,sizei=curseg_size;i<sizei;++i)
  {
    temp_track.hits.assign((*cur_seg)[i].n_hits, SimpleHit3D());
    temp_comb.assign((*cur_seg)[i].n_hits, 0);
    for(unsigned int l=0;l<(*cur_seg)[i].n_hits;++l)
    {
      temp_comb[l] = layer_sorted_1[(*cur_seg)[i].bin][l][(*cur_seg)[i].hits[l]].index;
    }
    sort(temp_comb.begin(),temp_comb.end());
    set<vector<unsigned int> >::iterator it = combos.find(temp_comb);
    if(it != combos.end()){continue;}
    combos.insert(temp_comb);
    for(unsigned int l=0;l<(*cur_seg)[i].n_hits;++l)
    {
      temp_track.hits[l] = layer_sorted_1[(*cur_seg)[i].bin][l][(*cur_seg)[i].hits[l]];
    }
    
    gettimeofday(&t1, NULL);
    
    for(unsigned int j=0;j<3;++j){temp_track_3hits.hits[j] = temp_track.hits[j];}
    fitTrack_3(temp_track_3hits);
    HelixKalmanState state;
    state.C = Matrix<float,5,5>::Zero(5,5);
    state.C(0,0) = 1.1;
    state.C(1,1) = 1.1;
    state.C(2,2) = 0.1;
    state.C(3,3) = 1.1;
    state.C(4,4) = 1.1;
    state.phi = temp_track_3hits.phi;
    if(state.phi < 0.){state.phi += 2.*M_PI;}
    state.d = temp_track_3hits.d;
    state.kappa = temp_track_3hits.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track_3hits.z0;
    state.dzdl = temp_track_3hits.dzdl;
    
    bool goodtrack = true;
    for(unsigned int h=0;h<temp_track.hits.size();++h)
    {
      kalman->addHit(temp_track.hits[h], state);nfits+=1;
      if(h >= 3)
      {
        if(state.chi2 != state.chi2)
        {
          goodtrack=false;break;
        }
        if(state.chi2/(2.*((float)(h+1)) - 5.) > chi2_cut)
        {
          goodtrack=false;break;
        }
      }
    }
    
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
    KALtime += (time2 - time1);
    
    if(state.chi2 != state.chi2)
    {
      goodtrack=false;
    }
    if(state.chi2/(2.*((float)(temp_track.hits.size())) - 5.) > chi2_cut)
    {
      goodtrack=false;
    }
    if(goodtrack==false){continue;}
    
    temp_track.phi = state.phi;
    if(temp_track.phi < 0.){temp_track.phi += 2.*M_PI;}
    if(temp_track.phi > 2.*M_PI){temp_track.phi -= 2.*M_PI;}
    temp_track.d = state.d;
    temp_track.kappa = state.kappa;
    temp_track.z0 = state.z0;
    temp_track.dzdl = state.dzdl;
    
    if(temp_track.kappa > top_range.max_k){continue;}
    
    tracks.push_back(temp_track);
    track_states.push_back(state);
    if((remove_hits == true) && (state.chi2 < chi2_removal_cut) && (temp_track.hits.size() >= n_removal_hits) )
    {
      for(unsigned int i=0;i<temp_track.hits.size();++i)
      {
        (*hit_used)[temp_track.hits[i].index] = true;
      }
    }
  }
  for(unsigned int b=0;b<4;++b)
  {
    for(unsigned int l=0;l<n_layers;++l)
    {
      layer_sorted_1[b][l].clear();
    }
  }
}


#endif
