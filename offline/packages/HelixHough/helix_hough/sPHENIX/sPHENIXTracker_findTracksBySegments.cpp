#include "vector_math_inline.h"
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

void sPHENIXTracker::calculateKappaTangents(float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a)
{
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
  tmp1 = (D12 + D23 + D31)*(D23 + D31 - D12)*(D12 + D31 - D23)*(D12 + D23 - D31);
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
  
  __m128 dr1 = dx1*dx1 + dy1*dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2*dx2 + dy2*dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3*dx3 + dy3*dy3;
  dr3 = _vec_sqrt_ps(dr3);
  
  __m128 dk1 = (dr1 + dr2)*D12_inv*D12_inv;
  __m128 dk2 = (dr2 + dr3)*D23_inv*D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);
  
  __m128 ux12 = (x2 - x1)*D12_inv;
  __m128 uy12 = (y2 - y1)*D12_inv;
  __m128 ux23 = (x3 - x2)*D23_inv;
  __m128 uy23 = (y3 - y2)*D23_inv;
  __m128 ux13 = (x3 - x1)*D31_inv;
  __m128 uy13 = (y3 - y1)*D31_inv;
  
  __m128 cosalpha = ux12*ux13 + uy12*uy13;
  __m128 sinalpha = ux13*uy12 - ux12*uy13;
  
  __m128 ux_mid = ux23*cosalpha - uy23*sinalpha;
  __m128 uy_mid = ux23*sinalpha + uy23*cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);
  
  __m128 ux_end = ux23*cosalpha + uy23*sinalpha;
  __m128 uy_end = uy23*cosalpha - ux23*sinalpha;
  
  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);
  
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha*sinalpha;
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
  __m128 dzdl_2 = s2*s2 + del_z_2*del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3)*D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);
  
  sinalpha = ux13*uy23 - ux23*uy13;
  v = one - sinalpha*sinalpha;
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
  __m128 dzdl_1 = s1*s1 + del_z_1*del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2)*D12_inv;
  _mm_store_ps(dzdl_1_a, dzdl_1);
  _mm_store_ps(ddzdl_1_a, ddzdl_1);
}


void sPHENIXTracker::calculateKappaTangents(float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a, float sinang_cut, float cosang_diff_inv, float* cur_kappa_a, float* cur_dkappa_a, float* cur_ux_a, float* cur_uy_a, float* cur_chi2_a, float* chi2_a)
{
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
  tmp1 = (D12 + D23 + D31)*(D23 + D31 - D12)*(D12 + D31 - D23)*(D12 + D23 - D31);
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
  
  __m128 dr1 = dx1*dx1 + dy1*dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2*dx2 + dy2*dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3*dx3 + dy3*dy3;
  dr3 = _vec_sqrt_ps(dr3);
  
  __m128 dk1 = (dr1 + dr2)*D12_inv*D12_inv;
  __m128 dk2 = (dr2 + dr3)*D23_inv*D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);
  
  __m128 ux12 = (x2 - x1)*D12_inv;
  __m128 uy12 = (y2 - y1)*D12_inv;
  __m128 ux23 = (x3 - x2)*D23_inv;
  __m128 uy23 = (y3 - y2)*D23_inv;
  __m128 ux13 = (x3 - x1)*D31_inv;
  __m128 uy13 = (y3 - y1)*D31_inv;
  
  __m128 cosalpha = ux12*ux13 + uy12*uy13;
  __m128 sinalpha = ux13*uy12 - ux12*uy13;
  
  __m128 ux_mid = ux23*cosalpha - uy23*sinalpha;
  __m128 uy_mid = ux23*sinalpha + uy23*cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);
  
  __m128 ux_end = ux23*cosalpha + uy23*sinalpha;
  __m128 uy_end = uy23*cosalpha - ux23*sinalpha;
  
  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);
  
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha*sinalpha;
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
  __m128 dzdl_2 = s2*s2 + del_z_2*del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3)*D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);
  
  sinalpha = ux13*uy23 - ux23*uy13;
  v = one - sinalpha*sinalpha;
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
  __m128 dzdl_1 = s1*s1 + del_z_1*del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2)*D12_inv;
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
  __m128 n_dk = c_dk + dk + sinang*k;
  __m128 chi2_k = kdiff*kdiff/(n_dk*n_dk);
  __m128 cos_scatter = c_ux*ux_mid + c_uy*uy_mid;
  __m128 chi2_ang = (one-cos_scatter)*(one-cos_scatter)*cosdiff*cosdiff;
  tmp1 = dzdl_1*sinang;
  _vec_fabs_ps(tmp1);
  __m128 chi2_dzdl = (dzdl_1 - dzdl_2)/(ddzdl_1 + ddzdl_2 + tmp1);
  chi2_dzdl *= chi2_dzdl;
  chi2_dzdl *= one_o_2;
  
  __m128 n_chi2 = c_chi2 + chi2_ang + chi2_k + chi2_dzdl;
  _mm_store_ps(chi2_a, n_chi2);
}


struct TempComb
{
  TempComb() : ndummies(0) {}
  HelixKalmanState state;
  vector<int> hit_indexes;
  unsigned int ndummies;
};



void sPHENIXTracker::initDummyHits( vector<SimpleHit3D>& dummies, const HelixRange& range, HelixKalmanState& init_state )
{
  SimpleTrack3D dummy_track;dummy_track.hits.push_back( SimpleHit3D() );
  dummy_track.kappa = 0.5*( range.min_k + range.max_k );
  dummy_track.phi = 0.5*( range.min_phi + range.max_phi );
  dummy_track.d = 0.5*( range.min_d + range.max_d );
  dummy_track.dzdl = 0.5*( range.min_dzdl + range.max_dzdl );
  dummy_track.z0 = 0.5*( range.min_z0 + range.max_z0 );
  
  init_state.kappa = dummy_track.kappa;
  init_state.nu = sqrt(dummy_track.kappa);
  init_state.phi = dummy_track.phi;
  init_state.d = dummy_track.d;
  init_state.dzdl = dummy_track.dzdl;
  init_state.z0 = dummy_track.z0;
  
  init_state.C = Matrix<float,5,5>::Zero(5,5);
  init_state.C(0,0) = pow( range.max_phi - range.min_phi , 2.);
  init_state.C(1,1) = pow( range.max_d - range.min_d , 2.);
  init_state.C(2,2) = pow( 10.*sqrt(range.max_k - range.min_k) , 2.);
  init_state.C(3,3) = pow( range.max_z0 - range.min_z0 , 2.);
  init_state.C(4,4) = pow( range.max_dzdl - range.min_dzdl , 2.);
  init_state.chi2 = 0.;
  init_state.position = 0;
  init_state.x_int = 0.;
  init_state.y_int = 0.;
  init_state.z_int = 0.;
  
  for(unsigned int i=0;i<n_layers;++i)
  {
    float x,y,z;
    projectToLayer( dummy_track, i, x, y, z );
    dummies[i].x = x;dummies[i].dx = 5.;
    dummies[i].y = x;dummies[i].dy = 5.;
    dummies[i].z = x;dummies[i].dz = 5.;
    dummies[i].layer = i;
  }
}


static bool next_combo_n( vector<int> const& lsizes, vector<int>& comb_n )
{
  unsigned int n = lsizes.size()/2;
  for(int l=0;l<n;++l)
  {
    if( comb_n[l] == (lsizes[l]-1) )
    {
      comb_n[l] = 0;
    }
    else
    {
      comb_n[l] += 1;
      return true;
    }
  }
  return false;
}

static SimpleHit3D& get_hit( vector<SimpleHit3D>& hits, vector<SimpleHit3D>& dummies, int index )
{
  if(index >= 0)
  {
    return hits[ index ];
  }
  else
  {
    return dummies[ (-index) - 1 ];
  }
}


void sPHENIXTracker::findTracksByCombinatorialKalman(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
  unsigned int half_layers = n_layers/2;
  float err_scale = 2.0;


  vector<SimpleHit3D> dummies(n_layers, SimpleHit3D());
  HelixKalmanState init_state;
  initDummyHits( dummies, range, init_state );
  vector<vector<int> > layer_indexes;layer_indexes.assign( n_layers, vector<int>() );
  for(unsigned int i=0;i<hits.size();++i){layer_indexes[ hits[i].layer ].push_back( i );}
  for(int i =0;i<(int)n_layers;++i){layer_indexes[i].push_back( -(i+1) );}
  
  vector<TempComb> comb1;
  vector<TempComb> comb2;
  
  vector<TempComb>* cur = &comb1;
  vector<TempComb>* prev = &comb2;
  
  unsigned int n_side = 3;
  unsigned int n_req = 3;
  vector<int> lsizes;
  for(int i=0;i<n_side;++i){ lsizes.push_back( layer_indexes[ i ].size() ); }
  for(int i=0;i<n_side;++i){ lsizes.push_back( layer_indexes[ n_layers - n_side + i ].size() ); }
  vector<int> comb_n;comb_n.assign(2*n_side, 0);
  while(true)
  {
    unsigned int ndummies = 0;
    unsigned int n_in = 0;
    unsigned int n_out = 0;
    SimpleTrack3D temp_track;
    TempComb tc;
    for(unsigned int h=0;h<n_side;++h)
    {
      if(layer_indexes[ h ][ comb_n[ h ] ] < 0){ndummies += 1;}
      if(layer_indexes[ h ][ comb_n[ h ] ] >= 0){n_in += 1;}
      SimpleHit3D& hit = get_hit( hits, dummies, layer_indexes[ h ][ comb_n[ h ] ] );
      temp_track.hits.push_back(hit);
      temp_track.hits.back().dx *= 0.5*hit_error_scale[h];
      temp_track.hits.back().dy *= 0.5*hit_error_scale[h];
      temp_track.hits.back().dz *= 0.5*hit_error_scale[h];
      tc.hit_indexes.push_back( layer_indexes[ h ][ comb_n[ h ] ] );
    }
    for(unsigned int h=0;h<n_side;++h)
    {
      if(layer_indexes[ n_layers - n_side + h ][ comb_n[ h+n_side ] ] < 0){ndummies += 1;}
      if(layer_indexes[ n_layers - n_side + h ][ comb_n[ h+n_side ] ] >= 0){n_out += 1;}
      SimpleHit3D& hit = get_hit( hits, dummies, layer_indexes[ n_layers - n_side + h ][ comb_n[ h+n_side ] ] );
      temp_track.hits.push_back(hit);
      temp_track.hits.back().dx *= 0.5*hit_error_scale[n_layers - n_side + h];
      temp_track.hits.back().dy *= 0.5*hit_error_scale[n_layers - n_side + h];
      temp_track.hits.back().dz *= 0.5*hit_error_scale[n_layers - n_side + h];
      tc.hit_indexes.push_back( layer_indexes[ n_layers - n_side + h ][ comb_n[ h+n_side ] ] );
    }
    
    if(layer_indexes[ 0 ][ comb_n[ 0 ] ] < 0){ if( next_combo_n( lsizes, comb_n ) == false ){break;} continue;}
    if(layer_indexes[ 1 ][ comb_n[ 1 ] ] < 0){ if( next_combo_n( lsizes, comb_n ) == false ){break;} continue;}
    
    float init_chi2 = 100.;
    if( (n_in >= n_req) && (n_out >= n_req) ){init_chi2 = fitTrack(temp_track);}
    else
    {
      if( next_combo_n( lsizes, comb_n ) == false ){break;}
      continue;
    }
    if( ( (n_in >= n_req) && (n_out >= n_req) ) )
    {
      cur->push_back( TempComb() );
      TempComb& curcomb = cur->back();
      curcomb.ndummies = ndummies;
      curcomb.state.chi2 = init_chi2/( 2.*( n_in + n_out ) - 5. );
      if( curcomb.state.chi2 > 10. ){cur->pop_back();}
      else{
      curcomb.state.phi = temp_track.phi;
      if(curcomb.state.phi < 0.){curcomb.state.phi += 2.*M_PI;}
      curcomb.state.d = temp_track.d;
      curcomb.state.kappa = temp_track.kappa;
      curcomb.state.nu = sqrt(curcomb.state.kappa);
      curcomb.state.z0 = temp_track.z0;
      curcomb.state.dzdl = temp_track.dzdl;
      curcomb.hit_indexes = tc.hit_indexes;}
    }
    if( next_combo_n( lsizes, comb_n ) == false ){break;}
  }
  if(cur->size() == 0)
  {
    // cout<<"no seeds found"<<endl;
    return;
  }
  
  // // find the best seed
  // int best = 0;
  // float chi2_best = cur->at(0).state.chi2;
  // for( unsigned int i=1;i<cur->size();++i )
  // {
  //   if( cur->at(i).state.chi2 < chi2_best )
  //   {
  //     chi2_best = cur->at(i).state.chi2;
  //     best = i;
  //   }
  // }

  for(ulong sd=0;sd<cur->size();sd++){
  TempComb seed = cur->at(sd);
  
  seed.state.C = Matrix<float,5,5>::Zero(5,5);
  seed.state.C(0,0) = pow(0.01, 2.);
  seed.state.C(1,1) = pow(0.5, 2.);
  seed.state.C(2,2) = pow(0.05*seed.state.nu, 2.);
  seed.state.C(3,3) = pow(0.5, 2.);
  seed.state.C(4,4) = pow(0.05, 2.);
  seed.state.chi2 = 0.;
  seed.state.position = n_layers;
  seed.state.x_int = 0.;
  seed.state.y_int = 0.;
  seed.state.z_int = 0.;
  
  seed.ndummies = 0;
  
  vector<int> in_seed( n_side, -1 );
  for(unsigned int i=0;i<n_side;++i)
  {
    if( seed.hit_indexes[i+n_side] > 0 )
    {
      in_seed[i] = seed.hit_indexes[i+n_side];
    }
  }
  
  seed.hit_indexes.clear();
  seed.hit_indexes.assign(n_layers, 0);
  vector<HelixKalmanState> temp_states;
  vector<int> temp_indexes;
  for(int l=n_layers-1;l>=0;--l)
  {
    temp_states.clear();
    temp_indexes.clear();
    for( unsigned int h=0;h<layer_indexes[l].size();++h )
    {
      if( (n_layers - 1 - l) < n_side )
      {
        if( in_seed[n_side - 1 - (n_layers - 1 - l)] >= 0 )
        {
          if( layer_indexes[l][h] != in_seed[n_side - 1 - (n_layers - 1 - l)] ){continue;}
        }
      }
      
      temp_indexes.push_back( layer_indexes[l][h] );
      temp_states.push_back( seed.state );
      SimpleHit3D hit = get_hit( hits, dummies, temp_indexes.back() );
      hit.dx *= err_scale*hit_error_scale[l];
      hit.dy *= err_scale*hit_error_scale[l];
      hit.dz *= err_scale*hit_error_scale[l];
      kalman->addHit( hit , temp_states.back() );
    }
    // is there a good hit?
    int best = -1;
    float best_chi2 = 99999999.;
    for( unsigned int h=0;h<temp_indexes.size();++h )
    {
      if(temp_indexes[h] < 0){continue;}
      if( (temp_states[h].chi2/(2.*( (float)((n_layers-l)+2 - seed.ndummies) ) ) < (chi2_cut)) && ( (temp_states[h].chi2-seed.state.chi2)<5. ) && (temp_states[h].chi2 < best_chi2) )
      {
        if( (n_layers-l) - (int)(seed.ndummies) > 3 )
        {
          if( temp_states[h].chi2/( 2.*((float)((n_layers-l)-(int)(seed.ndummies))) - 5. ) > chi2_cut )
          {
            continue;
          }
        }
        best = h;
        best_chi2 = temp_states[h].chi2;
      }
    }
    
    if( best >= 0 )
    {
      seed.hit_indexes[l] = ( temp_indexes[best] );
      seed.state = temp_states[ best ];
    }
    else
    {
      seed.hit_indexes[l] = ( -(l+1) );
      seed.ndummies += 1;
    }
    seed.state.position = l+1;
    if( seed.ndummies > (n_layers - req_layers) )
    {
      break;
    }
  }
  if( seed.ndummies > (n_layers - req_layers) ){continue;}

  
  tracks.push_back( SimpleTrack3D() );
  for(unsigned int i=0;i<n_layers;++i)
  {
    if( seed.hit_indexes[i] >= 0 )
    {
      tracks.back().hits.push_back( get_hit( hits, dummies, seed.hit_indexes[i] ) );
    }
  }

  {
    HelixKalmanState state = seed.state;
    
    state.C *= 0.05;
    for(int j=0;j<5;++j)
    {
      state.C(2,j) *= 0.03;
      state.C(j,2) *= 0.03;
    }
    
    state.chi2 = 0.;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;
    state.position = tracks.back().hits.size();
    for(int h=(tracks.back().hits.size() - 1);h>=0;--h)
    {
      SimpleHit3D hit = tracks.back().hits[h];
      float err_scale = 1.;
      int layer = hit.layer;
      if( (layer >= 0) && (layer < (int)(hit_error_scale.size()) ) ){err_scale = hit_error_scale[layer];}
      err_scale *= 0.4;
      hit.dx *= err_scale;hit.dy *= err_scale;hit.dz *= err_scale;
      
      kalman->addHit(hit, state);
    }
    
    if(!(state.kappa == state.kappa)){ tracks.pop_back();continue; }
    seed.state = state;
    
    tracks.back().phi = seed.state.phi;
    tracks.back().d = seed.state.d;
    tracks.back().kappa = seed.state.kappa;
    tracks.back().z0 = seed.state.z0;
    tracks.back().dzdl = seed.state.dzdl;

    if( seed.state.chi2 / (2.*tracks.back().hits.size() - 5.) > chi2_cut )
    {
      tracks.pop_back();continue;
    }
  }

  // cout<<"added track with "<<tracks.back().hits.size()<<" hits"<<endl;
  
  if(seed.state.phi < 0.){seed.state.phi += 2.*M_PI;}
  tracks.back().phi = seed.state.phi;
  tracks.back().d = seed.state.d;
  tracks.back().kappa = seed.state.kappa;
  tracks.back().z0 = seed.state.z0;
  tracks.back().dzdl = seed.state.dzdl;
  track_states.push_back(seed.state);
  if(remove_hits == true)
  {
    for(unsigned int i=0;i<tracks.back().hits.size();++i)
    {
      (*hit_used)[tracks.back().hits[i].index] = true;
    }
  }}
}


void sPHENIXTracker::findTracksBySegments(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
  if(n_layers > 20)
  {
    findTracksByCombinatorialKalman( hits, tracks, range );
    return;
  }
  
  
  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;
  unsigned int curseg_size = 0;
  unsigned int nextseg_size = 0;
  
  vector<TrackSegment> complete_segments;
  
  unsigned int allowed_missing = n_layers - req_layers;
  
  for(unsigned int l=0;l<n_layers;++l)
  {
    layer_sorted[l].clear();
  }
  for(unsigned int i=0;i<hits.size();++i)
  {
    unsigned int min = (hits[i].layer - allowed_missing);
    if(allowed_missing > hits[i].layer){min = 0;}
    for(unsigned int l=min;l<=hits[i].layer;l+=1)
    {
      layer_sorted[l].push_back(hits[i]);
      SimpleHit3D& hit = layer_sorted[l].back();
      float err_scale = 1.;
      int layer = hit.layer;
      if( (layer >= 0) && (layer < (int)(hit_error_scale.size()) ) ){err_scale *= 3.*hit_error_scale[layer];}
      hit.dx *= err_scale;hit.dy *= err_scale;hit.dz *= err_scale;
    }
  }
  for(unsigned int l=0;l<n_layers;++l)
  {
    if(layer_sorted[l].size()==0){return;}
  }
  
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
  
  TrackSegment temp_segment;temp_segment.hits.assign(n_layers, 0);
  // make segments out of first 3 layers
  for(unsigned int i=0,sizei=layer_sorted[0].size();i<sizei;++i)
  {
    for(unsigned int j=0,sizej=layer_sorted[1].size();j<sizej;++j)
    {
      for(unsigned int k=0,sizek=layer_sorted[2].size();k<sizek;++k)
      {
        if( (layer_sorted[0][i].layer >= layer_sorted[1][j].layer) || (layer_sorted[1][j].layer >= layer_sorted[2][k].layer) ){continue;}
        
        x1_a[hit_counter] = layer_sorted[0][i].x;
        y1_a[hit_counter] = layer_sorted[0][i].y;
        z1_a[hit_counter] = layer_sorted[0][i].z;
        dx1_a[hit_counter] = layer_sorted[0][i].dx;
        dy1_a[hit_counter] = layer_sorted[0][i].dy;
        dz1_a[hit_counter] = layer_sorted[0][i].dz;
        
        x2_a[hit_counter] = layer_sorted[1][j].x;
        y2_a[hit_counter] = layer_sorted[1][j].y;
        z2_a[hit_counter] = layer_sorted[1][j].z;
        dx2_a[hit_counter] = layer_sorted[1][j].dx;
        dy2_a[hit_counter] = layer_sorted[1][j].dy;
        dz2_a[hit_counter] = layer_sorted[1][j].dz;
        
        x3_a[hit_counter] = layer_sorted[2][k].x;
        y3_a[hit_counter] = layer_sorted[2][k].y;
        z3_a[hit_counter] = layer_sorted[2][k].z;
        dx3_a[hit_counter] = layer_sorted[2][k].dx;
        dy3_a[hit_counter] = layer_sorted[2][k].dy;
        dz3_a[hit_counter] = layer_sorted[2][k].dz;
        
        hit1[hit_counter] = i;
        hit2[hit_counter] = j;
        hit3[hit_counter] = k;
        
        hit_counter += 1;
        
        if(hit_counter == 4)
        {
          calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
          
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
            temp_segment.n_hits = 3;
            unsigned int outer_layer = layer_sorted[2][temp_segment.hits[2]].layer;
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
  if(hit_counter != 0)
  {
    calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
    
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
      temp_segment.n_hits = 3;
      unsigned int outer_layer = layer_sorted[2][temp_segment.hits[2]].layer;
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
  unsigned int whichseg[4];
  for(unsigned int l=3;l<n_layers;++l)
  {
    if(l == (n_layers-1)){easy_chi2_cut*=0.25;}
    nextseg_size = 0;
    for(unsigned int i=0,sizei=curseg_size;i<sizei;++i)
    {
      for(unsigned int j=0,sizej=layer_sorted[l].size();j<sizej;++j)
      {
        if( (layer_sorted[l-1][(*cur_seg)[i].hits[l-1]].layer >= layer_sorted[l][j].layer) ){continue;}
        
        x1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[i].hits[l-2]].x;
        y1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[i].hits[l-2]].y;
        z1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[i].hits[l-2]].z;
        x2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[i].hits[l-1]].x;
        y2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[i].hits[l-1]].y;
        z2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[i].hits[l-1]].z;
        x3_a[hit_counter] = layer_sorted[l][j].x;
        y3_a[hit_counter] = layer_sorted[l][j].y;
        z3_a[hit_counter] = layer_sorted[l][j].z;
        
        dx1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[i].hits[l-2]].dx;
        dy1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[i].hits[l-2]].dy;
        dz1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[i].hits[l-2]].dz;
        dx2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[i].hits[l-1]].dx;
        dy2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[i].hits[l-1]].dy;
        dz2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[i].hits[l-1]].dz;
        dx3_a[hit_counter] = layer_sorted[l][j].dx;
        dy3_a[hit_counter] = layer_sorted[l][j].dy;
        dz3_a[hit_counter] = layer_sorted[l][j].dz;
        
        cur_kappa_a[hit_counter] = (*cur_seg)[i].kappa;
        cur_dkappa_a[hit_counter] = (*cur_seg)[i].dkappa;
        cur_ux_a[hit_counter] = (*cur_seg)[i].ux;
        cur_uy_a[hit_counter] = (*cur_seg)[i].uy;
        cur_chi2_a[hit_counter] = (*cur_seg)[i].chi2;
        
        whichseg[hit_counter] = i;
        hit1[hit_counter] = j;
        
        hit_counter += 1;
        if(hit_counter == 4)
        {
          calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a, sinang_cut, cosang_diff_inv, cur_kappa_a, cur_dkappa_a, cur_ux_a, cur_uy_a, cur_chi2_a, chi2_a);
          
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
              for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
              temp_segment.hits[l] = hit1[h];
              unsigned int outer_layer = layer_sorted[l][temp_segment.hits[l]].layer;
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
      calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a, sinang_cut, cosang_diff_inv, cur_kappa_a, cur_dkappa_a, cur_ux_a, cur_uy_a, cur_chi2_a, chi2_a);
      
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
          for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
          temp_segment.hits[l] = hit1[h];
          unsigned int outer_layer = layer_sorted[l][temp_segment.hits[l]].layer;
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
  vector<SimpleHit3D> temp_hits;
  for(unsigned int i=0,sizei=curseg_size;i<sizei;++i)
  {
    temp_track.hits.assign((*cur_seg)[i].n_hits, SimpleHit3D());
    
    temp_comb.assign((*cur_seg)[i].n_hits, 0);
    for(unsigned int l=0;l<(*cur_seg)[i].n_hits;++l)
    {
      temp_comb[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].index;
    }
    sort(temp_comb.begin(),temp_comb.end());
    set<vector<unsigned int> >::iterator it = combos.find(temp_comb);
    if(it != combos.end()){continue;}
    if(combos.size() > 10000){combos.clear();}
    combos.insert(temp_comb);
    
    for(unsigned int l=0;l<(*cur_seg)[i].n_hits;++l)
    {
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
      if (init_chi2 > fast_chi2_cut_par0+fast_chi2_cut_par1/kappaToPt(temp_track.kappa)) {
  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
  KALtime += (time2 - time1);
  continue;
      }
    }
    HelixKalmanState state;
    state.phi = temp_track.phi;
    if(state.phi < 0.){state.phi += 2.*M_PI;}
    state.d = temp_track.d;
    state.kappa = temp_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track.z0;
    state.dzdl = temp_track.dzdl;
    state.C = Matrix<float,5,5>::Zero(5,5);
    state.C(0,0) = pow(0.01, 2.);
    state.C(1,1) = pow(0.01, 2.);
    state.C(2,2) = pow(0.01*state.nu, 2.);
    state.C(3,3) = pow(0.05, 2.);
    state.C(4,4) = pow(0.05, 2.);
    state.chi2 = 0.;
    state.position = 0;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;
        
    for(unsigned int h=0;h<temp_track.hits.size();++h)
    {
      kalman->addHit(temp_track.hits[h], state);nfits+=1;
    }
    
    // fudge factor for non-gaussian hit sizes
    state.C *= 3.;
    state.chi2 *= 6.;
    
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
    KALtime += (time2 - time1);
    
    if( !(temp_track.kappa == temp_track.kappa) ){continue;}
    if(temp_track.kappa > top_range.max_k){continue;}
    if( !( state.chi2 == state.chi2 ) ){continue;}
    if(state.chi2/(2.*((float)(temp_track.hits.size())) - 5.) > chi2_cut){continue;}
    
    if( cut_on_dca == true )
    {
      if( fabs(temp_track.d) > dca_cut ){continue;}
      if( fabs(temp_track.z0) > dca_cut ){continue;}
    }

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
}


void sPHENIXTracker::findSeededTracksbySegments(vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
  unsigned int allowed_missing = n_layers - seed_layer - req_layers;
  
  for(unsigned int l=0;l<n_layers;++l)
  {
    layer_sorted[l].clear();
  }
  for(unsigned int i=0;i<hits.size();++i)
  {
    unsigned int min = (hits[i].layer - allowed_missing);
    if(allowed_missing > hits[i].layer){min = 0;}
    for(unsigned int l=min;l<=hits[i].layer;l+=1){layer_sorted[l].push_back(hits[i]);}
  }
  
  findSeededTracksbySegments_run(seeds, hits, tracks);
}


void sPHENIXTracker::findSeededTracksbySegments_run(vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks)
{
  if(seeds.size() == 0){return;}
  
  timeval t1,t2;
  double time1=0.;
  double time2=0.;
  
  gettimeofday(&t1, NULL);
  
  unsigned int first_new_layer = seed_layer;
  
  unsigned int allowed_missing = n_layers - seed_layer - req_layers;
  
  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;
  unsigned int curseg_size = 0;
  unsigned int nextseg_size = 0;
  
  vector<TrackSegment> complete_segments;
  
  float cosang_diff = 1. - cosang_cut;
  float cosang_diff_inv = 1./cosang_diff;
  float sinang_cut = sqrt(1. - cosang_cut*cosang_cut);
  float easy_chi2_cut = ca_chi2_cut;
  
  unsigned int hit_counter = 0;
  float x1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float x2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float x3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float y1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float y2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float y3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float z1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float z2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float z3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float dx1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dx2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dx3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dy1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dy2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dy3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dz1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dz2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dz3_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float kappa_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dkappa_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float ux_mid_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float uy_mid_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float ux_end_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float uy_end_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  float dzdl_1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float dzdl_2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float ddzdl_1_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  float ddzdl_2_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  TrackSegment temp_segment;temp_segment.hits.assign(n_layers, 0);
  unsigned int whichhit[4];
  unsigned int whichseed[4];
  for(unsigned int seed=0,seedsize=seeds.size();seed<seedsize;++seed)
  {
    unsigned int firsthit = seeds[seed].hits.size();
    
    x1_a[hit_counter] = seeds[seed].hits[firsthit-3].x;
    y1_a[hit_counter] = seeds[seed].hits[firsthit-3].y;
    z1_a[hit_counter] = seeds[seed].hits[firsthit-3].z;
    x2_a[hit_counter] = seeds[seed].hits[firsthit-2].x;
    y2_a[hit_counter] = seeds[seed].hits[firsthit-2].y;
    z2_a[hit_counter] = seeds[seed].hits[firsthit-2].z;
    x3_a[hit_counter] = seeds[seed].hits[firsthit-1].x;
    y3_a[hit_counter] = seeds[seed].hits[firsthit-1].y;
    z3_a[hit_counter] = seeds[seed].hits[firsthit-1].z;
    
    dx1_a[hit_counter] = seeds[seed].hits[firsthit-3].dx;
    dy1_a[hit_counter] = seeds[seed].hits[firsthit-3].dy;
    dz1_a[hit_counter] = seeds[seed].hits[firsthit-3].dz;
    dx2_a[hit_counter] = seeds[seed].hits[firsthit-2].dx;
    dy2_a[hit_counter] = seeds[seed].hits[firsthit-2].dy;
    dz2_a[hit_counter] = seeds[seed].hits[firsthit-2].dz;
    dx3_a[hit_counter] = seeds[seed].hits[firsthit-1].dx;
    dy3_a[hit_counter] = seeds[seed].hits[firsthit-1].dy;
    dz3_a[hit_counter] = seeds[seed].hits[firsthit-1].dz;
    
    whichseed[hit_counter] = seed;
    hit_counter += 1;
    if(hit_counter == 4)
    {
      calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
      
      for(unsigned int h=0;h<hit_counter;++h)
      {
        temp_segment.chi2 = 0.;
        temp_segment.ux = ux_end_a[h];
        temp_segment.uy = uy_end_a[h];
        temp_segment.kappa = kappa_a[h];
        temp_segment.dkappa = dkappa_a[h];
        temp_segment.seed = whichseed[h];
        if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
        else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
      }
      
      hit_counter = 0;
    }
  }
  if(hit_counter != 0)
  {
    calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
    
    for(unsigned int h=0;h<hit_counter;++h)
    {
      temp_segment.chi2 = 0.;
      temp_segment.ux = ux_end_a[h];
      temp_segment.uy = uy_end_a[h];
      temp_segment.kappa = kappa_a[h];
      temp_segment.dkappa = dkappa_a[h];
      temp_segment.seed = whichseed[h];
      if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
      else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
    }
    
    hit_counter = 0;
  }
  swap(cur_seg, next_seg);
  swap(curseg_size, nextseg_size);
  unsigned int whichseg[4];
  for(unsigned int l=first_new_layer;l<n_layers;++l)
  {
    nextseg_size = 0;
    for(unsigned int j=0,sizej=curseg_size;j<sizej;++j)
    {
      unsigned int firsthit = seeds[(*cur_seg)[j].seed].hits.size();
      for(unsigned int i=0,sizei=layer_sorted[l].size();i<sizei;++i)
      {
        if((l-2)<=(first_new_layer-1))
        {
          x1_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-2-(first_new_layer-firsthit)].x;
          y1_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-2-(first_new_layer-firsthit)].y;
          z1_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-2-(first_new_layer-firsthit)].z;
          dx1_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-2-(first_new_layer-firsthit)].dx;
          dy1_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-2-(first_new_layer-firsthit)].dy;
          dz1_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-2-(first_new_layer-firsthit)].dz;
        }
        else
        {
          x1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[j].hits[l-2]].x;
          y1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[j].hits[l-2]].y;
          z1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[j].hits[l-2]].z;
          dx1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[j].hits[l-2]].dx;
          dy1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[j].hits[l-2]].dy;
          dz1_a[hit_counter] = layer_sorted[l-2][(*cur_seg)[j].hits[l-2]].dz;
        }
        if((l-1)<=(first_new_layer-1))
        {
          x2_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-1-(first_new_layer-firsthit)].x;
          y2_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-1-(first_new_layer-firsthit)].y;
          z2_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-1-(first_new_layer-firsthit)].z;
          dx2_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-1-(first_new_layer-firsthit)].dx;
          dy2_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-1-(first_new_layer-firsthit)].dy;
          dz2_a[hit_counter] = seeds[(*cur_seg)[j].seed].hits[l-1-(first_new_layer-firsthit)].dz;
        }
        else
        {
          x2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[j].hits[l-1]].x;
          y2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[j].hits[l-1]].y;
          z2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[j].hits[l-1]].z;
          dx2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[j].hits[l-1]].dx;
          dy2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[j].hits[l-1]].dy;
          dz2_a[hit_counter] = layer_sorted[l-1][(*cur_seg)[j].hits[l-1]].dz;
        }
        x3_a[hit_counter] = layer_sorted[l][i].x;
        y3_a[hit_counter] = layer_sorted[l][i].y;
        z3_a[hit_counter] = layer_sorted[l][i].z;
        dx3_a[hit_counter] = layer_sorted[l][i].dx;
        dy3_a[hit_counter] = layer_sorted[l][i].dy;
        dz3_a[hit_counter] = layer_sorted[l][i].dz;
        
        whichhit[hit_counter] = i;
        whichseg[hit_counter] = j;
        hit_counter += 1;
        if(hit_counter == 4)
        {
          calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
          
          for(unsigned int h=0;h<hit_counter;++h)
          {
            float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
            float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
            dk += sinang_cut*kappa_a[h];
            float chi2_k = kdiff*kdiff/(dk*dk);
            float cos_scatter = (*cur_seg)[whichseg[h]].ux*ux_mid_a[h] + (*cur_seg)[whichseg[h]].uy*uy_mid_a[h];
            float chi2_ang = (1.-cos_scatter)*(1.-cos_scatter)*cosang_diff_inv*cosang_diff_inv;
            float chi2_dzdl = (dzdl_1_a[h] - dzdl_2_a[h])/(ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h]*sinang_cut));
            chi2_dzdl *= chi2_dzdl;
            if( ((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl)/((float)l - 2.) < easy_chi2_cut )
            {
              temp_segment.chi2 = (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
              temp_segment.ux = ux_end_a[h];
              temp_segment.uy = uy_end_a[h];
              temp_segment.kappa = kappa_a[h];
              temp_segment.dkappa = dkappa_a[h];
              for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
              temp_segment.hits[l] = whichhit[h];
              temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
              unsigned int outer_layer = layer_sorted[l][temp_segment.hits[l]].layer;
              temp_segment.n_hits = l+1;
              if( (n_layers - (l+1)) <= allowed_missing ){complete_segments.push_back(temp_segment);}
              if( (outer_layer - l) > allowed_missing ){continue;}
              if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
              else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
            }
          }
          hit_counter = 0;
        }
      }
    }
    if(hit_counter != 0)
    {
      calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);
      
      for(unsigned int h=0;h<hit_counter;++h)
      {
        float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
        float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
        dk += sinang_cut*kappa_a[h];
        float chi2_k = kdiff*kdiff/(dk*dk);
        float cos_scatter = (*cur_seg)[whichseg[h]].ux*ux_mid_a[h] + (*cur_seg)[whichseg[h]].uy*uy_mid_a[h];
        float chi2_ang = (1.-cos_scatter)*(1.-cos_scatter)*cosang_diff_inv*cosang_diff_inv;
        float chi2_dzdl = (dzdl_1_a[h] - dzdl_2_a[h])/(ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h]*sinang_cut));
        chi2_dzdl *= chi2_dzdl;
        if( ((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl)/((float)l - 2.) < easy_chi2_cut )
        {
          temp_segment.chi2 = (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
          temp_segment.ux = ux_end_a[h];
          temp_segment.uy = uy_end_a[h];
          temp_segment.kappa = kappa_a[h];
          temp_segment.dkappa = dkappa_a[h];
          for(unsigned int ll=0;ll<l;++ll){temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];}
          temp_segment.hits[l] = whichhit[h];
          temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
          unsigned int outer_layer = layer_sorted[l][temp_segment.hits[l]].layer;
          temp_segment.n_hits = l+1;
          if( (n_layers - (l+1)) <= allowed_missing ){complete_segments.push_back(temp_segment);}
          if( (outer_layer - l) > allowed_missing ){continue;}
          if(next_seg->size() == nextseg_size){next_seg->push_back(temp_segment);nextseg_size+=1;}
          else{(*next_seg)[nextseg_size] = temp_segment;nextseg_size+=1;}
        }
      }
      hit_counter = 0;
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
  
  
  gettimeofday(&t1, NULL);
  
  SimpleTrack3D temp_track;
  temp_track.hits.assign(n_layers, SimpleHit3D());
  for(unsigned int i=0,sizei=curseg_size;i<sizei;++i)
  {
    temp_track.hits.clear();
    for(unsigned int l=0;l<seeds[(*cur_seg)[i].seed].hits.size();++l)
    {
      temp_track.hits.push_back(seeds[(*cur_seg)[i].seed].hits[l]);
    }
    for(unsigned int l=seed_layer;l<(*cur_seg)[i].n_hits;++l)
    {
      temp_track.hits.push_back(layer_sorted[l][(*cur_seg)[i].hits[l]]);
    }
    
    
    HelixKalmanState state;
    
    fitTrack(temp_track);
    state.C = Matrix<float,5,5>::Zero(5,5);
    state.phi = temp_track.phi;
    if(state.phi < 0.){state.phi += 2.*M_PI;}
    state.d = temp_track.d;
    state.kappa = temp_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track.z0;
    state.dzdl = temp_track.dzdl;
    state.C(0,0) = pow(0.1*state.phi, 2.);
    state.C(1,1) = pow(0.1*state.d, 2.);
    state.C(2,2) = pow(0.1*state.nu, 2.);
    state.C(3,3) = pow(0.1*state.z0, 2.);
    state.C(4,4) = pow(0.1*state.dzdl, 2.);
    
    state.chi2 = 0.;
    state.position = 0;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;
    for(unsigned int l=0;l<temp_track.hits.size();++l)
    {
      kalman->addHit(temp_track.hits[l], state);nfits+=1;
    }

    
    
    if(state.chi2/(2.*((float)(temp_track.hits.size()+1)) - 5.) > chi2_cut)
    {
      continue;
    }
    if(state.chi2 != state.chi2)
    {
      continue;
    }
    
    temp_track.phi = state.phi;
    if(temp_track.phi < 0.){temp_track.phi += 2.*M_PI;}
    if(temp_track.phi > 2.*M_PI){temp_track.phi -= 2.*M_PI;}
    temp_track.d = state.d;
    temp_track.kappa = state.kappa;
    temp_track.z0 = state.z0;
    temp_track.dzdl = state.dzdl;
    
    if((remove_hits == true) && (state.chi2 < chi2_removal_cut) && (temp_track.hits.size() >= n_removal_hits))
    {
      for(unsigned int l=(seeds[(*cur_seg)[i].seed].hits.size());l<temp_track.hits.size();++l)
      {
        (*hit_used)[temp_track.hits[l].index] = true;
        temp_track.hits[l].index = index_mapping[temp_track.hits[l].index];
      }
    }
    else
    {
      for(unsigned int l=(seeds[(*cur_seg)[i].seed].hits.size());l<temp_track.hits.size();++l)
      {
        temp_track.hits[l].index = index_mapping[temp_track.hits[l].index];
      }
    }
    tracks.push_back(temp_track);
    track_states.push_back(state);
    seed_used[seeds[(*cur_seg)[i].seed].index] = true;
  }
  
  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
  KALtime += (time2 - time1);
}





