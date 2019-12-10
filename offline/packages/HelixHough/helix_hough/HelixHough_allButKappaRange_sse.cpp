#include "HelixHough.h"

#include "vector_math_inline.h"

#include <emmintrin.h>
#include <xmmintrin.h>

using namespace std;


static const __m128 one_o_4 = {0.25,0.25,0.25,0.25};
static const __m128 two = {2., 2., 2., 2.};
static const __m128 one_o_100 = {0.01,0.01,0.01,0.01};
static const __m128 close_one = {0.999,0.999,0.999,0.999};
static const __m128 four = {4., 4., 4., 4.};
static const __m128 one_o_3 = {0.3333333333333333333, 0.3333333333333333333, 0.3333333333333333333, 0.3333333333333333333};
static const __m128 _3_o_20 = {0.15, 0.15, 0.15, 0.15};
static const __m128 _5_o_56 = {8.92857142857142877e-02, 8.92857142857142877e-02, 8.92857142857142877e-02, 8.92857142857142877e-02};
static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
static const __m128 three_pi_over_two = {3.*0x1.921fb54442d1846ap0f, 3.*0x1.921fb54442d1846ap0f, 3.*0x1.921fb54442d1846ap0f, 3.*0x1.921fb54442d1846ap0f};


static inline void __attribute__((always_inline)) calculate_phi_d(__m128& k, __m128& Delta, __m128& Delta2, __m128& Delta_inv, __m128& ux, __m128& uy, __m128& x3, __m128& y3, __m128& phi_1, __m128& phi_2, __m128& d_1, __m128& d_2)
{
  __m128 k_inv = _vec_rec_ps(k);
  __m128 k2 = _mm_mul_ps(k, k);
  
  __m128 uscale2 = _mm_mul_ps(Delta2, k2);
  uscale2 = _mm_mul_ps(one_o_4, uscale2);
  uscale2 = _mm_sub_ps(one, uscale2);
  __m128 uscale = _vec_sqrt_ps(uscale2);
  
  __m128 tmp1 = _mm_mul_ps(k, y3);
  __m128 tmp2 = _mm_mul_ps(uscale, uy);
  __m128 tmp3 = _mm_mul_ps(k, x3);
  __m128 tmp4 = _mm_mul_ps(uscale, ux);
  
  __m128 tmp5 = _mm_add_ps(tmp1, tmp2);
  __m128 tmp6 = _mm_add_ps(tmp3, tmp4);
  phi_1 = _vec_atan2_ps(tmp5, tmp6);
  tmp5 = _mm_sub_ps(tmp1, tmp2);
  tmp6 = _mm_sub_ps(tmp3, tmp4);
  phi_2 = _vec_atan2_ps(tmp5, tmp6);
  
  tmp1 = _mm_cmplt_ps(phi_1, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_1 = _mm_add_ps(phi_1, tmp1);
  
  tmp1 = _mm_cmplt_ps(phi_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_2 = _mm_add_ps(phi_2, tmp1);
  
  tmp1 = _mm_mul_ps(x3, x3);
  tmp2 = _mm_mul_ps(y3, y3);
  tmp1 = _mm_add_ps(tmp1, tmp2);
  __m128 kd1 = _mm_mul_ps(k2, tmp1);
  kd1 = _mm_add_ps(kd1, uscale2);
  tmp1 = _mm_mul_ps(two, k);
  tmp1 = _mm_mul_ps(tmp1, uscale);
  __m128 k0val = _mm_mul_ps(x3, ux); // value of d when k = 0
  tmp3 = _mm_mul_ps(y3, uy);
  k0val = _mm_add_ps(k0val, tmp3);
  tmp1 = _mm_mul_ps(tmp1, k0val);
  __m128 kd2 = _mm_sub_ps(kd1, tmp1);
  kd1 = _mm_add_ps(kd1, tmp1);
  kd1 = _vec_sqrt_ps(kd1);
  kd2 = _vec_sqrt_ps(kd2);
  
  //helicity 1, k 1
  d_1 = kd1;
  d_1 = _mm_sub_ps(d_1, one);
  d_1 = _mm_mul_ps(d_1, k_inv);
  //helicity 2, k 1
  d_2 = kd2;
  d_2 = _mm_sub_ps(d_2, one);
  d_2 = _mm_mul_ps(d_2, k_inv);
  // if k=0 , set d to k0val
  tmp1 = _mm_cmpeq_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, k0val);
  tmp3 = _mm_andnot_ps(tmp1, d_1);
  d_1 = _mm_xor_ps(tmp2, tmp3);
  tmp3 = _mm_andnot_ps(tmp1, d_2);
  d_2 = _mm_xor_ps(tmp2, tmp3);
}


static inline __m128 __attribute__((always_inline)) calculate_dzdl(__m128& Delta, __m128& z1, __m128& z2, __m128& k)
{
  __m128 v = _mm_mul_ps(one_o_2, k);
  v = _mm_mul_ps(v, Delta);
  //if(v > 0.999){v = 0.999;}
  __m128 tmp1 = _mm_cmpgt_ps(v, close_one);
  __m128 tmp2 = _mm_and_ps(tmp1, close_one);
  __m128 tmp3 = _mm_andnot_ps(tmp1, v);
  v = _mm_xor_ps(tmp2, tmp3);
  __m128 one_o_v = _vec_rec_ps(v);
  //power series assuming v_v is small
  __m128 s = zero;
  __m128 temp1 = _mm_mul_ps(v, v);
  __m128 temp2 = _mm_mul_ps(one_o_2, Delta);
  tmp1 = _mm_mul_ps(two, temp2);
  s = _mm_add_ps(s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, one_o_3);
  s = _mm_add_ps(s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, _3_o_20);
  s = _mm_add_ps(s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, _5_o_56);
  s = _mm_add_ps(s, tmp1);
  ////////////////////////////////////
  //otherwise we calculate an arcsin
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  //s = 2*asin(v)/k
  tmp1 = _mm_mul_ps(v, v);
  tmp1 = _mm_sub_ps(one, tmp1);
  tmp1 = _vec_sqrt_ps(tmp1);
  tmp1 = _mm_add_ps(one, tmp1);
  tmp1 = _mm_mul_ps(tmp1, one_o_v);
  tmp2 = _vec_atan_ps(tmp1);
  tmp2 = _mm_sub_ps(pi_over_two, tmp2);
  tmp2 = _mm_mul_ps(four, tmp2);
  tmp2 = _mm_mul_ps(tmp2, one_o_v);
  tmp2 = _mm_mul_ps(tmp2, Delta);
  tmp2 = _mm_mul_ps(tmp2, one_o_2);
  ////////////////////////////////////
  //choose between the two methods to calculate s
  tmp1 = _mm_cmpgt_ps(v, one_o_100);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  tmp2 = _mm_andnot_ps(tmp1, s);
  __m128 s1 = _mm_xor_ps(tmp3, tmp2);
  
  __m128 dz2 = _mm_sub_ps(z2, z1);
  dz2 = _mm_mul_ps(dz2, dz2);
  tmp1 = _mm_mul_ps(s1, s1);
  tmp1 = _mm_add_ps(tmp1, dz2);
  __m128 dzdl = _mm_div_ps(dz2, tmp1);
  dzdl = _vec_sqrt_ps(dzdl);
  //if z2 < z1, dzdl = -dzdl
  tmp1 = _mm_cmplt_ps(z2, z1);
  tmp2 = _mm_xor_ps(dzdl, SIGNMASK);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  dzdl = _mm_andnot_ps(tmp1, dzdl);
  dzdl = _mm_xor_ps(dzdl, tmp3);
  
  return dzdl;
}


static inline __m128 __attribute__((always_inline)) calculate_z0(__m128& x, __m128& y, __m128& z, __m128& k, __m128& d, __m128& phi, __m128& dzdl)
{
  __m128 cs, sn;
  _vec_sin_cos_ps(phi, sn, cs);
  __m128 dx = d*cs;
  __m128 dy = d*sn;
  
  __m128 Dx = dx - x;
  __m128 Dy = dy - y;
  __m128 Delta = _vec_sqrt_ps(Dx*Dx + Dy*Dy);
  
  __m128 v = _mm_mul_ps(one_o_2, k);
  v = _mm_mul_ps(v, Delta);
  //if(v > 0.999){v = 0.999;}
  __m128 tmp1 = _mm_cmpgt_ps(v, close_one);
  __m128 tmp2 = _mm_and_ps(tmp1, close_one);
  __m128 tmp3 = _mm_andnot_ps(tmp1, v);
  v = _mm_xor_ps(tmp2, tmp3);
  __m128 one_o_v = _vec_rec_ps(v);
  //power series assuming v_v is small
  __m128 s = zero;
  __m128 temp1 = _mm_mul_ps(v, v);
  __m128 temp2 = _mm_mul_ps(one_o_2, Delta);
  tmp1 = _mm_mul_ps(two, temp2);
  s = _mm_add_ps(s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, one_o_3);
  s = _mm_add_ps(s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, _3_o_20);
  s = _mm_add_ps(s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, _5_o_56);
  s = _mm_add_ps(s, tmp1);
  ////////////////////////////////////
  //otherwise we calculate an arcsin
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  //s = 2*asin(v)/k
  tmp1 = _mm_mul_ps(v, v);
  tmp1 = _mm_sub_ps(one, tmp1);
  tmp1 = _vec_sqrt_ps(tmp1);
  tmp1 = _mm_add_ps(one, tmp1);
  tmp1 = _mm_mul_ps(tmp1, one_o_v);
  tmp2 = _vec_atan_ps(tmp1);
  tmp2 = _mm_sub_ps(pi_over_two, tmp2);
  tmp2 = _mm_mul_ps(four, tmp2);
  tmp2 = _mm_mul_ps(tmp2, one_o_v);
  tmp2 = _mm_mul_ps(tmp2, Delta);
  tmp2 = _mm_mul_ps(tmp2, one_o_2);
  ////////////////////////////////////
  //choose between the two methods to calculate s
  tmp1 = _mm_cmpgt_ps(v, one_o_100);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  tmp2 = _mm_andnot_ps(tmp1, s);
  __m128 s1 = _mm_xor_ps(tmp3, tmp2);
  
  // dz/ds = dzdl/(1 - dzdl^2)
  __m128 dzds = dzdl*_vec_rsqrt_ps(one - dzdl*dzdl);
  
  return (z - dzds*s1);
}


static inline void __attribute__((always_inline)) find_min_max(__m128 val1, __m128 val2, __m128& min, __m128& max)
{
  __m128 tmp1 = _mm_cmplt_ps(val1, val2);
  __m128 tmp2 = _mm_and_ps(tmp1, val1);
  __m128 tmp3 = _mm_andnot_ps(tmp1, val2);
  min = _mm_xor_ps(tmp2, tmp3);
  tmp2 = _mm_and_ps(tmp1, val2);
  tmp3 = _mm_andnot_ps(tmp1, val1);
  max = _mm_xor_ps(tmp2, tmp3);
}


void HelixHough::allButKappaRange_sse(float* x1_a,float* x2_a,float* y1_a,float* y2_a,float* z1_a,float* z2_a, float* min_k_a,float* max_k_a, float* min_phi_1_a,float* max_phi_1_a,float* min_phi_2_a,float* max_phi_2_a, float* min_d_1_a,float* max_d_1_a,float* min_d_2_a,float* max_d_2_a, float* min_dzdl_a,float* max_dzdl_a, float* min_z0_1_a,float* max_z0_1_a,float* min_z0_2_a,float* max_z0_2_a)
{
  __m128 x1 = _mm_load_ps(x1_a);
  __m128 y1 = _mm_load_ps(y1_a);
  __m128 z1 = _mm_load_ps(z1_a);
  
  __m128 x2 = _mm_load_ps(x2_a);
  __m128 y2 = _mm_load_ps(y2_a);
  __m128 z2 = _mm_load_ps(z2_a);
  
  __m128 x3 = _mm_add_ps(x1, x2);
  x3 = _mm_mul_ps(one_o_2, x3);
  __m128 y3 = _mm_add_ps(y1, y2);
  y3 = _mm_mul_ps(one_o_2, y3);
  __m128 Delta2 = _mm_sub_ps(x2, x1);
  Delta2 = _mm_mul_ps(Delta2, Delta2);
  __m128 tmp1 = _mm_sub_ps(y2, y1);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  Delta2 = _mm_add_ps(Delta2, tmp1);
  __m128 Delta = _vec_sqrt_ps(Delta2);
  __m128 Delta_inv = _vec_rec_ps(Delta);
  __m128 ux = _mm_sub_ps(y2, y1);
  ux = _mm_mul_ps(ux, Delta_inv);
  __m128 uy = _mm_sub_ps(x1, x2);
  uy = _mm_mul_ps(uy, Delta_inv);
  
  
  __m128 k = _mm_load_ps(min_k_a);
  __m128 phi_1_1, phi_2_1;
  __m128 d_1_1, d_2_1;
  calculate_phi_d(k, Delta, Delta2, Delta_inv, ux, uy, x3, y3, phi_1_1, phi_2_1, d_1_1, d_2_1);
  
  __m128 dzdl_1 = calculate_dzdl(Delta, z1, z2, k);
  __m128 z0_1_1 = calculate_z0(x1, y1, z1, k, d_1_1, phi_1_1, dzdl_1);
  __m128 z0_2_1 = calculate_z0(x1, y1, z1, k, d_2_1, phi_2_1, dzdl_1);
  
  
  k = _mm_load_ps(max_k_a);
  __m128 phi_1_2, phi_2_2;
  __m128 d_1_2, d_2_2;
  calculate_phi_d(k, Delta, Delta2, Delta_inv, ux, uy, x3, y3, phi_1_2, phi_2_2, d_1_2, d_2_2);
  
  __m128 dzdl_2 = calculate_dzdl(Delta, z1, z2, k);
  __m128 z0_1_2 = calculate_z0(x1, y1, z1, k, d_1_2, phi_1_2, dzdl_2);
  __m128 z0_2_2 = calculate_z0(x1, y1, z1, k, d_2_2, phi_2_2, dzdl_2);
  
  // choose the min and max for each helicity
  // start with helicity 1 :
  // check if phi overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_1_1, pi_over_two);
  __m128 tmp2 = _mm_cmplt_ps(phi_1_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmpgt_ps(phi_1_1, three_pi_over_two);
  __m128 tmp3 = _mm_cmpgt_ps(phi_1_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp1 = _mm_and_ps(tmp1, tmp2);
  // tmp1 is now all ones if phi_r overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);
  __m128 tmp4 = _mm_cmpgt_ps(phi_1_1, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  __m128 tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_1_1 = _mm_sub_ps(phi_1_1, tmp3);
  tmp4 = _mm_cmpgt_ps(phi_1_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_1_2 = _mm_sub_ps(phi_1_2, tmp3);
  
  __m128 phi_1_min, phi_1_max;
  find_min_max(phi_1_1, phi_1_2, phi_1_min, phi_1_max);
  
  __m128 d_1_min, d_1_max;
  find_min_max(d_1_1, d_1_2, d_1_min, d_1_max);
  
  __m128 z0_1_min, z0_1_max;
  find_min_max(z0_1_1, z0_1_2, z0_1_min, z0_1_max);
  
  _mm_store_ps(min_phi_1_a, phi_1_min);
  _mm_store_ps(max_phi_1_a, phi_1_max);
  _mm_store_ps(min_d_1_a, d_1_min);
  _mm_store_ps(max_d_1_a, d_1_max);
  _mm_store_ps(min_z0_1_a, z0_1_min);
  _mm_store_ps(max_z0_1_a, z0_1_max);
  
  
  // choose the min and max for each helicity
  // now for helicity 2 :
  // check if phi overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_2_1, pi_over_two);
  tmp2 = _mm_cmplt_ps(phi_2_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmpgt_ps(phi_2_1, three_pi_over_two);
  tmp3 = _mm_cmpgt_ps(phi_2_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp1 = _mm_and_ps(tmp1, tmp2);
  // tmp1 is now all ones if phi_r overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);
  tmp4 = _mm_cmpgt_ps(phi_2_1, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_2_1 = _mm_sub_ps(phi_2_1, tmp3);
  tmp4 = _mm_cmpgt_ps(phi_2_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_2_2 = _mm_sub_ps(phi_2_2, tmp3);
  
  __m128 phi_2_min, phi_2_max;
  find_min_max(phi_2_1, phi_2_2, phi_2_min, phi_2_max);
  
  __m128 d_2_min, d_2_max;
  find_min_max(d_2_1, d_2_2, d_2_min, d_2_max);
  
  __m128 z0_2_min, z0_2_max;
  find_min_max(z0_2_1, z0_2_2, z0_2_min, z0_2_max);
  
  _mm_store_ps(min_phi_2_a, phi_2_min);
  _mm_store_ps(max_phi_2_a, phi_2_max);
  _mm_store_ps(min_d_2_a, d_2_min);
  _mm_store_ps(max_d_2_a, d_2_max);
  _mm_store_ps(min_z0_2_a, z0_2_min);
  _mm_store_ps(max_z0_2_a, z0_2_max);
  
  __m128 dzdl_min, dzdl_max;
  find_min_max(dzdl_1, dzdl_2, dzdl_min, dzdl_max);
  
  _mm_store_ps(min_dzdl_a, dzdl_min);
  _mm_store_ps(max_dzdl_a, dzdl_max);
}


