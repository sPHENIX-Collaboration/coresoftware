#include "HelixHough.h"

#include "vector_math_inline.h"

#include <emmintrin.h>
#include <xmmintrin.h>

using namespace std;


static const __m128 one_o_100 = {0.01,0.01,0.01,0.01};
static const __m128 close_one = {0.999,0.999,0.999,0.999};
static const __m128 two = {2., 2., 2., 2.};
static const __m128 four = {4., 4., 4., 4.};
static const __m128 one_o_3 = {0.3333333333333333333, 0.3333333333333333333, 0.3333333333333333333, 0.3333333333333333333};
static const __m128 _3_o_20 = {0.15, 0.15, 0.15, 0.15};
static const __m128 _5_o_56 = {8.92857142857142877e-02, 8.92857142857142877e-02, 8.92857142857142877e-02, 8.92857142857142877e-02};
static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));


void HelixHough::dzdlRange_sse(float* x_a, float* y_a, float* z_a, float cosphi1, float sinphi1, float cosphi2, float sinphi2, float min_k_val, float max_k_val, float min_d_val, float max_d_val, float* min_z0_val, float* max_z0_val, float* min_dzdl_a, float* max_dzdl_a)
{
  __m128 x = _mm_load_ps(x_a);
  __m128 y = _mm_load_ps(y_a);
  __m128 z = _mm_load_ps(z_a);
  
  __m128 cosphi_min = _mm_load1_ps(&(cosphi1));
  __m128 cosphi_max = _mm_load1_ps(&(cosphi2));
  __m128 sinphi_min = _mm_load1_ps(&(sinphi1));
  __m128 sinphi_max = _mm_load1_ps(&(sinphi2));
  
  __m128 min_k = _mm_load1_ps(&(min_k_val));
  __m128 max_k = _mm_load1_ps(&(max_k_val));
  __m128 min_d = _mm_load1_ps(&(min_d_val));
  __m128 max_d = _mm_load1_ps(&(max_d_val));
  __m128 min_z0 = _mm_load_ps((min_z0_val));
  __m128 max_z0 = _mm_load_ps((max_z0_val));
  
  __m128 d = min_d;                                            __m128 d_2 = max_d;
  __m128 k = min_k;                                            __m128 k_2 = max_k;
  __m128 dx = _mm_mul_ps(cosphi_min, d);                       __m128 dx_2 = _mm_mul_ps(cosphi_max,d_2);
  __m128 dy = _mm_mul_ps(sinphi_min,d);                        __m128 dy_2 = _mm_mul_ps(sinphi_max,d_2);
  __m128 D = _mm_sub_ps(x, dx);                                __m128 D_2 = _mm_sub_ps(x, dx_2);
  D = _mm_mul_ps(D, D);                                        D_2 = _mm_mul_ps(D_2, D_2);
  __m128 tmp1 = _mm_sub_ps(y, dy);                             __m128 tmp1_2 = _mm_sub_ps(y, dy_2);
  tmp1 = _mm_mul_ps(tmp1, tmp1);                               tmp1_2 = _mm_mul_ps(tmp1_2, tmp1_2);
  D = _mm_add_ps(D, tmp1);                                     D_2 = _mm_add_ps(D_2, tmp1_2);
  D = _vec_sqrt_ps(D);                                         D_2 = _vec_sqrt_ps(D_2);
  __m128 v = _mm_mul_ps(one_o_2, k);                           __m128 v_2 = _mm_mul_ps(one_o_2, k_2);
  v = _mm_mul_ps(v, D);                                        v_2 = _mm_mul_ps(v_2, D_2);
  //if(v > 0.999){v = 0.999;}                                  //if(v > 0.999){v = 0.999;}
  tmp1 = _mm_cmpgt_ps(v, close_one);                           tmp1_2 = _mm_cmpgt_ps(v_2, close_one);
  __m128 tmp2 = _mm_and_ps(tmp1, close_one);                   __m128 tmp2_2 = _mm_and_ps(tmp1_2, close_one);
  __m128 tmp3 = _mm_andnot_ps(tmp1, v);                        __m128 tmp3_2 = _mm_andnot_ps(tmp1_2, v_2);
  v = _mm_xor_ps(tmp2, tmp3);                                  v_2 = _mm_xor_ps(tmp2_2, tmp3_2);
  __m128 one_o_v = _vec_rec_ps(v);                          __m128 one_o_v_2 = _vec_rec_ps(v_2);
  //power series assuming v_v is small                         //power series assuming v_v is small
  __m128 s = zero;                                             __m128 s_2 = zero;
  __m128 temp1 = _mm_mul_ps(v, v);                             __m128 temp1_2 = _mm_mul_ps(v_2, v_2);
  __m128 temp2 = _mm_mul_ps(one_o_2, D);                       __m128 temp2_2 = _mm_mul_ps(one_o_2, D_2);
  tmp1 = _mm_mul_ps(two, temp2);                               tmp1_2 = _mm_mul_ps(two, temp2_2);
  s = _mm_add_ps(s, tmp1);                                     s_2 = _mm_add_ps(s_2, tmp1_2);
  temp2 = _mm_mul_ps(temp2, temp1);                            temp2_2 = _mm_mul_ps(temp2_2, temp1_2);
  tmp1 = _mm_mul_ps(temp2, one_o_3);                           tmp1_2 = _mm_mul_ps(temp2_2, one_o_3);
  s = _mm_add_ps(s, tmp1);                                     s_2 = _mm_add_ps(s_2, tmp1_2);
  temp2 = _mm_mul_ps(temp2, temp1);                            temp2_2 = _mm_mul_ps(temp2_2, temp1_2);
  tmp1 = _mm_mul_ps(temp2, _3_o_20);                           tmp1_2 = _mm_mul_ps(temp2_2, _3_o_20);
  s = _mm_add_ps(s, tmp1);                                     s_2 = _mm_add_ps(s_2, tmp1_2);
  temp2 = _mm_mul_ps(temp2, temp1);                            temp2_2 = _mm_mul_ps(temp2_2, temp1_2);
  tmp1 = _mm_mul_ps(temp2, _5_o_56);                           tmp1_2 = _mm_mul_ps(temp2_2, _5_o_56);
  s = _mm_add_ps(s, tmp1);                                     s_2 = _mm_add_ps(s_2, tmp1_2);
  ////////////////////////////////////                         ////////////////////////////////////
  //otherwise we calculate an arcsin                           //otherwise we calculate an arcsin
  //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )              //asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  //s = 2*asin(v)/k                                            //s = 2*asin(v)/k
  tmp1 = _mm_mul_ps(v, v);                                     tmp1_2 = _mm_mul_ps(v_2, v_2);
  tmp1 = _mm_sub_ps(one, tmp1);                                tmp1_2 = _mm_sub_ps(one, tmp1_2);
  tmp1 = _vec_sqrt_ps(tmp1);                                   tmp1_2 = _vec_sqrt_ps(tmp1_2);
  tmp1 = _mm_add_ps(one, tmp1);                                tmp1_2 = _mm_add_ps(one, tmp1_2);
  tmp1 = _mm_mul_ps(tmp1, one_o_v);                            tmp1_2 = _mm_mul_ps(tmp1_2, one_o_v_2);
  tmp2 = _vec_atan_ps(tmp1);                                   tmp2_2 = _vec_atan_ps(tmp1_2);
  tmp2 = _mm_sub_ps(pi_over_two, tmp2);                        tmp2_2 = _mm_sub_ps(pi_over_two, tmp2_2);
  tmp2 = _mm_mul_ps(four, tmp2);                               tmp2_2 = _mm_mul_ps(four, tmp2_2);
  tmp2 = _mm_mul_ps(tmp2, one_o_v);                            tmp2_2 = _mm_mul_ps(tmp2_2, one_o_v_2);
  tmp2 = _mm_mul_ps(tmp2, D);                                  tmp2_2 = _mm_mul_ps(tmp2_2, D_2);
  tmp2 = _mm_mul_ps(tmp2, one_o_2);                            tmp2_2 = _mm_mul_ps(tmp2_2, one_o_2);
  ////////////////////////////////////                         ////////////////////////////////////
  //choose between the two methods to calculate s              //choose between the two methods to calculate s
  tmp1 = _mm_cmpgt_ps(v, one_o_100);                           tmp1_2 = _mm_cmpgt_ps(v_2, one_o_100);
  tmp3 = _mm_and_ps(tmp1, tmp2);                               tmp3_2 = _mm_and_ps(tmp1_2, tmp2_2);
  tmp2 = _mm_andnot_ps(tmp1, s);                               tmp2_2 = _mm_andnot_ps(tmp1_2, s_2);
  __m128 s1 = _mm_xor_ps(tmp3, tmp2);                          __m128 s2 = _mm_xor_ps(tmp3_2, tmp2_2);
  
  
  
  
  ////////////////////////////////////
  __m128 dz2 = _mm_sub_ps(z, max_z0);
  dz2 = _mm_mul_ps(dz2, dz2);
  tmp1 = _mm_mul_ps(s1, s1);
  tmp1 = _mm_add_ps(tmp1, dz2);
  __m128 dzdl_1 = _mm_div_ps(dz2, tmp1);
  dzdl_1 = _vec_sqrt_ps(dzdl_1);
  //if z < max_z0, dzdl = -dzdl_1
  tmp1 = _mm_cmplt_ps(z, max_z0);
  tmp2 = _mm_xor_ps(dzdl_1, SIGNMASK);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  dzdl_1 = _mm_andnot_ps(tmp1, dzdl_1);
  dzdl_1 = _mm_xor_ps(dzdl_1, tmp3);
  
  ////////////////////////////////////
  dz2 = _mm_sub_ps(z, min_z0);
  dz2 = _mm_mul_ps(dz2, dz2);
  tmp1 = _mm_mul_ps(s1, s1);
  tmp1 = _mm_add_ps(tmp1, dz2);
  __m128 dzdl_2 = _mm_div_ps(dz2, tmp1);
  dzdl_2 = _vec_sqrt_ps(dzdl_2);
  //if z < min_z0, dzdl = -dzdl_2
  tmp1 = _mm_cmplt_ps(z, min_z0);
  tmp2 = _mm_xor_ps(dzdl_2, SIGNMASK);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  dzdl_2 = _mm_andnot_ps(tmp1, dzdl_2);
  dzdl_2 = _mm_xor_ps(dzdl_2, tmp3);
  
  ////////////////////////////////////
  dz2 = _mm_sub_ps(z, max_z0);
  dz2 = _mm_mul_ps(dz2, dz2);
  tmp1 = _mm_mul_ps(s2, s2);
  tmp1 = _mm_add_ps(tmp1, dz2);
  __m128 dzdl_3 = _mm_div_ps(dz2, tmp1);
  dzdl_3 = _vec_sqrt_ps(dzdl_3);
  //if z < max_z0, dzdl = -dzdl_3
  tmp1 = _mm_cmplt_ps(z, max_z0);
  tmp2 = _mm_xor_ps(dzdl_3, SIGNMASK);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  dzdl_3 = _mm_andnot_ps(tmp1, dzdl_3);
  dzdl_3 = _mm_xor_ps(dzdl_3, tmp3);
  
  ////////////////////////////////////
  dz2 = _mm_sub_ps(z, min_z0);
  dz2 = _mm_mul_ps(dz2, dz2);
  tmp1 = _mm_mul_ps(s2, s2);
  tmp1 = _mm_add_ps(tmp1, dz2);
  __m128 dzdl_4 = _mm_div_ps(dz2, tmp1);
  dzdl_4 = _vec_sqrt_ps(dzdl_4);
  //if z < min_z0, dzdl = -dzdl_4
  tmp1 = _mm_cmplt_ps(z, min_z0);
  tmp2 = _mm_xor_ps(dzdl_4, SIGNMASK);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  dzdl_4 = _mm_andnot_ps(tmp1, dzdl_4);
  dzdl_4 = _mm_xor_ps(dzdl_4, tmp3);
  
  
  
  
  __m128 dzdl_max = dzdl_1;
  tmp1 = _mm_cmpgt_ps(dzdl_2, dzdl_max);
  tmp2 = _mm_and_ps(tmp1, dzdl_2);
  tmp3 = _mm_andnot_ps(tmp1, dzdl_max);
  dzdl_max = _mm_xor_ps(tmp2, tmp3);
  tmp1 = _mm_cmpgt_ps(dzdl_3, dzdl_max);
  tmp2 = _mm_and_ps(tmp1, dzdl_3);
  tmp3 = _mm_andnot_ps(tmp1, dzdl_max);
  dzdl_max = _mm_xor_ps(tmp2, tmp3);
  tmp1 = _mm_cmpgt_ps(dzdl_4, dzdl_max);
  tmp2 = _mm_and_ps(tmp1, dzdl_4);
  tmp3 = _mm_andnot_ps(tmp1, dzdl_max);
  dzdl_max = _mm_xor_ps(tmp2, tmp3);
  
  __m128 dzdl_min = dzdl_1;
  tmp1 = _mm_cmplt_ps(dzdl_2, dzdl_min);
  tmp2 = _mm_and_ps(tmp1, dzdl_2);
  tmp3 = _mm_andnot_ps(tmp1, dzdl_min);
  dzdl_min = _mm_xor_ps(tmp2, tmp3);
  tmp1 = _mm_cmplt_ps(dzdl_3, dzdl_min);
  tmp2 = _mm_and_ps(tmp1, dzdl_3);
  tmp3 = _mm_andnot_ps(tmp1, dzdl_min);
  dzdl_min = _mm_xor_ps(tmp2, tmp3);
  tmp1 = _mm_cmplt_ps(dzdl_4, dzdl_min);
  tmp2 = _mm_and_ps(tmp1, dzdl_4);
  tmp3 = _mm_andnot_ps(tmp1, dzdl_min);
  dzdl_min = _mm_xor_ps(tmp2, tmp3);
  
  
  
  _mm_store_ps(min_dzdl_a, dzdl_min);
  _mm_store_ps(max_dzdl_a, dzdl_max);
}
