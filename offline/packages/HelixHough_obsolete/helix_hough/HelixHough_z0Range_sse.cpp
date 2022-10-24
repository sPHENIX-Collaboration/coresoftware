#include "vector_math_inline.h"
#include "HelixHough.h"
#include <math.h>
#include <iostream>

using namespace std;


static const __m128 one_o_10 = {0.1,0.1,0.1,0.1};
static const __m128 one_o_2 = {0.5,0.5,0.5,0.5};
static const __m128 close_one = {0.999,0.999,0.999,0.999};
static const __m128 two = {2., 2., 2., 2.};
static const __m128 four = {4., 4., 4., 4.};
static const __m128 one_o_3 = {0.3333333333333333333, 0.3333333333333333333, 0.3333333333333333333, 0.3333333333333333333};
static const __m128 _3_o_20 = {0.15, 0.15, 0.15, 0.15};
static const __m128 _5_o_56 = {8.92857142857142877e-02, 8.92857142857142877e-02, 8.92857142857142877e-02, 8.92857142857142877e-02};
static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));

void HelixHough::z0Range_sse(const SimpleHit3D& hit, float cosphi1, float sinphi1, float cosphi2, float sinphi2, float min_k, float max_k, float min_phi, float max_phi, float min_d, float max_d, float min_theta, float max_theta, float& min_z0, float& max_z0)
{
  float z0_1[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  
  __m128 v_z0_1 = _mm_load_ps(z0_1);
  
  float x_arr[4] __attribute__((aligned(16))) = {hit.x,hit.x,hit.x,hit.x};
  float y_arr[4] __attribute__((aligned(16))) = {hit.y,hit.y,hit.y,hit.y};
  float z_arr[4] __attribute__((aligned(16))) = {hit.z,hit.z,hit.z,hit.z};
  
  __m128 v_x = _mm_load_ps(x_arr);
  __m128 v_y = _mm_load_ps(y_arr);
  __m128 v_z = _mm_load_ps(z_arr);
  
  
  //   p0,p0,p1,p1  d0,d0,d1,d1  k0,k0,k1,k1  t0,t1,t0,t1
  
  float cosphi_1[4] __attribute__((aligned(16))) = {cosphi1, cosphi1, cosphi2, cosphi2}; __m128 v_cosphi_1 = _mm_load_ps(cosphi_1);
  float sinphi_1[4] __attribute__((aligned(16))) = {sinphi1, sinphi1, sinphi2, sinphi2}; __m128 v_sinphi_1 = _mm_load_ps(sinphi_1);
  float d_1[4] __attribute__((aligned(16))) = {min_d, min_d, max_d, max_d}; __m128 v_d_1 = _mm_load_ps(d_1);
  float k_1[4] __attribute__((aligned(16))) = {min_k, min_k, max_k, max_k}; __m128 v_k_1 = _mm_load_ps(k_1);
  float t_1[4] __attribute__((aligned(16))) = {min_theta, max_theta, min_theta, max_theta}; __m128 v_t_1 = _mm_load_ps(t_1);
  
  
  //first 4
  __m128 v_dx = _mm_mul_ps(v_cosphi_1, v_d_1);
  __m128 v_dy = _mm_mul_ps(v_sinphi_1, v_d_1);
  __m128 tmp1 = _mm_sub_ps(v_x, v_dx);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  __m128 v_D = _mm_sub_ps(v_y, v_dy);
  v_D = _mm_mul_ps(v_D, v_D);
  v_D = _mm_add_ps(v_D, tmp1);
  v_D = _mm_sqrt_ps(v_D);
  __m128 v_v = _mm_mul_ps(one_o_2, v_k_1);
  v_v = _mm_mul_ps(v_v, v_D);
  //if(v > 0.999){v = 0.999;}
  tmp1 = _mm_cmpgt_ps(v_v, close_one);
  __m128 tmp2 = _mm_and_ps(tmp1, close_one);
  __m128 tmp3 = _mm_andnot_ps(tmp1, v_v);
  v_v = _mm_xor_ps(tmp2, tmp3);
  //power series assuming v_v is small
  __m128 v_s = zero;
  __m128 temp1 = _mm_mul_ps(v_v, v_v);
  __m128 temp2 = _mm_mul_ps(one_o_2, v_D);
  tmp1 = _mm_mul_ps(two, temp2);
  v_s = _mm_add_ps(v_s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, one_o_3);
  v_s = _mm_add_ps(v_s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, _3_o_20);
  v_s = _mm_add_ps(v_s, tmp1);
  temp2 = _mm_mul_ps(temp2, temp1);
  tmp1 = _mm_mul_ps(temp2, _5_o_56);
  v_s = _mm_add_ps(v_s, tmp1);
  ////////////////////////////////////
  //otherwise we calculate an arcsin
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  //s = 2*asin(v)/k
  tmp1 = _mm_mul_ps(v_v, v_v);
  tmp1 = _mm_sub_ps(one, tmp1);
  tmp1 = _mm_sqrt_ps(tmp1);
  tmp1 = _mm_add_ps(one, tmp1);
  tmp1 = _mm_div_ps(v_v, tmp1);
  tmp2 = _vec_atan_ps(tmp1);
  tmp2 = _mm_mul_ps(four, tmp2);
  tmp2 = _mm_div_ps(tmp2, v_k_1);
  ////////////////////////////////////
  //choose between the two methods to calculate s
  tmp1 = _mm_cmpgt_ps(v_v, one_o_10);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  tmp2 = _mm_andnot_ps(tmp1, v_s);
  v_s = _mm_xor_ps(tmp3, tmp2);
  ////////////////////////////////////
  tmp1 = _mm_mul_ps(v_s, v_s);
  tmp2 = _mm_mul_ps(v_t_1, v_t_1);
  tmp3 = _mm_mul_ps(tmp1, tmp2);
  __m128 tmp4 = _mm_sub_ps(one, tmp2);
  tmp4 = _mm_div_ps(tmp3, tmp4);
  tmp4 = _mm_sqrt_ps(tmp4);
  //if(t > 0){dz = -tmp4;}
  //else{dz = tmp4;}
  tmp3 = _mm_xor_ps(tmp4, SIGNMASK);
  tmp1 = _mm_cmpgt_ps(v_t_1, zero);
  tmp2 = _mm_and_ps(tmp1, tmp3);
  __m128 v_dz = _mm_andnot_ps(tmp1, tmp4);
  v_dz = _mm_xor_ps(tmp2, v_dz);
  v_z0_1 = _mm_add_ps(v_dz, v_z);
  
  
  float z0_a[4] __attribute__((aligned(16))) = {0.,0.,0.,0.};
  _mm_store_ps(z0_a, v_z0_1);
  
  min_z0 = z0_a[0];
  if (z0_a[1] < min_z0)  min_z0 = z0_a[1] ;
  if (z0_a[2] < min_z0)  min_z0 = z0_a[2] ;
  if (z0_a[3] < min_z0)  min_z0 = z0_a[3] ;
  
  max_z0 = z0_a[0];
  if (z0_a[1] > max_z0)  max_z0 = z0_a[1] ;
  if (z0_a[2] > max_z0)  max_z0 = z0_a[2] ;
  if (z0_a[3] > max_z0)  max_z0 = z0_a[3] ;
  
  
  min_z0 -= hit.dz;
  max_z0 += hit.dz;
}
