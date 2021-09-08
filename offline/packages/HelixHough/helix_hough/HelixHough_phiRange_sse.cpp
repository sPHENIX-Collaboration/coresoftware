#include "HelixHough.h"
#include "vector_math_inline.h"

#include <emmintrin.h>
#include <xmmintrin.h>

using namespace std;


static const __m128 two = {2., 2., 2., 2.};
static const __m128 three_pi_over_two = {3.*0x1.921fb54442d1846ap0f, 3.*0x1.921fb54442d1846ap0f, 3.*0x1.921fb54442d1846ap0f, 3.*0x1.921fb54442d1846ap0f};
static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));

void HelixHough::phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi_1, float* max_phi_1, float* min_phi_2, float* max_phi_2)
{
  __m128 x = _mm_load_ps(hit_x);
  __m128 y = _mm_load_ps(hit_y);
  
  __m128 hit_phi = _vec_atan2_ps(y,x);
  //if phi < 0, phi += 2*pi
  __m128 tmp1 = _mm_cmplt_ps(hit_phi, zero);
  __m128 tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  hit_phi = _mm_add_ps(hit_phi, tmp1);
  
  // min_d, min_k
  __m128 d = _mm_load_ps(min_d);
  __m128 k = _mm_load_ps(min_k);
  __m128 D = _mm_mul_ps(x,x);
  tmp1 = _mm_mul_ps(y,y);
  D = _mm_add_ps(D,tmp1);
  D = _vec_sqrt_ps(D);
  __m128 D_inv = _vec_rec_ps(D);
  __m128 ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  __m128 hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  __m128 neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  __m128 xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  __m128 xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  __m128 yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  __m128 yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  __m128 phi_r_1 = _vec_atan2_ps(yk1, xk1);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_r_1, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_r_1 = _mm_add_ps(phi_r_1, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_r_1 = _mm_andnot_ps(neg, phi_r_1);
  phi_r_1 = _mm_xor_ps(tmp1, phi_r_1);
  
  __m128 phi_l_1 = _vec_atan2_ps(yk2, xk2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_l_1, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_l_1 = _mm_add_ps(phi_l_1, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_l_1 = _mm_andnot_ps(neg, phi_l_1);
  phi_l_1 = _mm_xor_ps(tmp1, phi_l_1);
  
  
  // min_d, max_k
  d = _mm_load_ps(min_d);
  k = _mm_load_ps(max_k);
  D = _mm_mul_ps(x,x);
  tmp1 = _mm_mul_ps(y,y);
  D = _mm_add_ps(D,tmp1);
  D = _vec_sqrt_ps(D);
  D_inv = _vec_rec_ps(D);
  ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  __m128 phi_r_2 = _vec_atan2_ps(yk1, xk1);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_r_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_r_2 = _mm_add_ps(phi_r_2, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_r_2 = _mm_andnot_ps(neg, phi_r_2);
  phi_r_2 = _mm_xor_ps(tmp1, phi_r_2);
  
  __m128 phi_l_2 = _vec_atan2_ps(yk2, xk2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_l_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_l_2 = _mm_add_ps(phi_l_2, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_l_2 = _mm_andnot_ps(neg, phi_l_2);
  phi_l_2 = _mm_xor_ps(tmp1, phi_l_2);
  
  // max_d, min_k
  d = _mm_load_ps(max_d);
  k = _mm_load_ps(max_k);
  D = _mm_mul_ps(x,x);
  tmp1 = _mm_mul_ps(y,y);
  D = _mm_add_ps(D,tmp1);
  D = _vec_sqrt_ps(D);
  D_inv = _vec_rec_ps(D);
  ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  __m128 phi_r_3 = _vec_atan2_ps(yk1, xk1);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_r_3, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_r_3 = _mm_add_ps(phi_r_3, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_r_3 = _mm_andnot_ps(neg, phi_r_3);
  phi_r_3 = _mm_xor_ps(tmp1, phi_r_3);
  
  __m128 phi_l_3 = _vec_atan2_ps(yk2, xk2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_l_3, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_l_3 = _mm_add_ps(phi_l_3, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_l_3 = _mm_andnot_ps(neg, phi_l_3);
  phi_l_3 = _mm_xor_ps(tmp1, phi_l_3);
  
  // max_d, max_k
  d = _mm_load_ps(max_d);
  k = _mm_load_ps(max_k);
  D = _mm_mul_ps(x,x);
  tmp1 = _mm_mul_ps(y,y);
  D = _mm_add_ps(D,tmp1);
  D = _vec_sqrt_ps(D);
  D_inv = _vec_rec_ps(D);
  ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  __m128 phi_r_4 = _vec_atan2_ps(yk1, xk1);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_r_4, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_r_4 = _mm_add_ps(phi_r_4, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_r_4 = _mm_andnot_ps(neg, phi_r_4);
  phi_r_4 = _mm_xor_ps(tmp1, phi_r_4);
  
  __m128 phi_l_4 = _vec_atan2_ps(yk2, xk2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_l_4, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_l_4 = _mm_add_ps(phi_l_4, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_l_4 = _mm_andnot_ps(neg, phi_l_4);
  phi_l_4 = _mm_xor_ps(tmp1, phi_l_4);
  
  ////////////////////////////////////////////////////////////////
  
  // check if phi_r overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_r_1, pi_over_two);
  tmp2 = _mm_cmplt_ps(phi_r_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_r_3, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_r_4, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  
  tmp2 = _mm_cmpgt_ps(phi_r_1, three_pi_over_two);
  __m128 tmp3 = _mm_cmpgt_ps(phi_r_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_r_3, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_r_4, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  
  tmp1 = _mm_and_ps(tmp1, tmp2);
  
  // tmp1 is now all ones if phi_r overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi_r values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);
  
  __m128 tmp4 = _mm_cmpgt_ps(phi_r_1, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  __m128 tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_r_1 = _mm_sub_ps(phi_r_1, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_r_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_r_2 = _mm_sub_ps(phi_r_2, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_r_3, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_r_3 = _mm_sub_ps(phi_r_3, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_r_4, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_r_4 = _mm_sub_ps(phi_r_4, tmp3);
  
  
  // find the minimum phi_r
  __m128 phi_r_min = phi_r_1;
  tmp2 = _mm_cmplt_ps(phi_r_2, phi_r_min);
  tmp3 = _mm_and_ps(tmp2, phi_r_2);
  phi_r_min = _mm_andnot_ps(tmp2, phi_r_min);
  phi_r_min = _mm_xor_ps(phi_r_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_r_3, phi_r_min);
  tmp3 = _mm_and_ps(tmp2, phi_r_3);
  phi_r_min = _mm_andnot_ps(tmp2, phi_r_min);
  phi_r_min = _mm_xor_ps(phi_r_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_r_4, phi_r_min);
  tmp3 = _mm_and_ps(tmp2, phi_r_4);
  phi_r_min = _mm_andnot_ps(tmp2, phi_r_min);
  phi_r_min = _mm_xor_ps(phi_r_min, tmp3);
  
  // find the maximum phi_r
  __m128 phi_r_max = phi_r_1;
  tmp2 = _mm_cmpgt_ps(phi_r_2, phi_r_max);
  tmp3 = _mm_and_ps(tmp2, phi_r_2);
  phi_r_max = _mm_andnot_ps(tmp2, phi_r_max);
  phi_r_max = _mm_xor_ps(phi_r_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_r_3, phi_r_max);
  tmp3 = _mm_and_ps(tmp2, phi_r_3);
  phi_r_max = _mm_andnot_ps(tmp2, phi_r_max);
  phi_r_max = _mm_xor_ps(phi_r_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_r_4, phi_r_max);
  tmp3 = _mm_and_ps(tmp2, phi_r_4);
  phi_r_max = _mm_andnot_ps(tmp2, phi_r_max);
  phi_r_max = _mm_xor_ps(phi_r_max, tmp3);
  
  _mm_store_ps(min_phi_1, phi_r_min);
  _mm_store_ps(max_phi_1, phi_r_max);
  
  
  
  // check if phi_l overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_l_1, pi_over_two);
  tmp2 = _mm_cmplt_ps(phi_l_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_l_3, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_l_4, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  
  tmp2 = _mm_cmpgt_ps(phi_l_1, three_pi_over_two);
  tmp3 = _mm_cmpgt_ps(phi_l_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_l_3, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_l_4, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  
  tmp1 = _mm_and_ps(tmp1, tmp2);
  
  // tmp1 is now all ones if phi_l overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi_l values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_l_1, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_l_1 = _mm_sub_ps(phi_l_1, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_l_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_l_2 = _mm_sub_ps(phi_l_2, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_l_3, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_l_3 = _mm_sub_ps(phi_l_3, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_l_4, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_l_4 = _mm_sub_ps(phi_l_4, tmp3);
  
  
  // find the minimum phi_l
  __m128 phi_l_min = phi_l_1;
  tmp2 = _mm_cmplt_ps(phi_l_2, phi_l_min);
  tmp3 = _mm_and_ps(tmp2, phi_l_2);
  phi_l_min = _mm_andnot_ps(tmp2, phi_l_min);
  phi_l_min = _mm_xor_ps(phi_l_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_l_3, phi_l_min);
  tmp3 = _mm_and_ps(tmp2, phi_l_3);
  phi_l_min = _mm_andnot_ps(tmp2, phi_l_min);
  phi_l_min = _mm_xor_ps(phi_l_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_l_4, phi_l_min);
  tmp3 = _mm_and_ps(tmp2, phi_l_4);
  phi_l_min = _mm_andnot_ps(tmp2, phi_l_min);
  phi_l_min = _mm_xor_ps(phi_l_min, tmp3);
  
  // find the maximum phi_l
  __m128 phi_l_max = phi_l_1;
  tmp2 = _mm_cmpgt_ps(phi_l_2, phi_l_max);
  tmp3 = _mm_and_ps(tmp2, phi_l_2);
  phi_l_max = _mm_andnot_ps(tmp2, phi_l_max);
  phi_l_max = _mm_xor_ps(phi_l_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_l_3, phi_l_max);
  tmp3 = _mm_and_ps(tmp2, phi_l_3);
  phi_l_max = _mm_andnot_ps(tmp2, phi_l_max);
  phi_l_max = _mm_xor_ps(phi_l_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_l_4, phi_l_max);
  tmp3 = _mm_and_ps(tmp2, phi_l_4);
  phi_l_max = _mm_andnot_ps(tmp2, phi_l_max);
  phi_l_max = _mm_xor_ps(phi_l_max, tmp3);
  
  _mm_store_ps(min_phi_2, phi_l_min);
  _mm_store_ps(max_phi_2, phi_l_max);
  
}


static inline __m128 compare_sign(__m128 a, __m128 b)
{
  const __m128i MASK = _mm_set1_epi32(0xffffffff);
  
  __m128  f = _mm_xor_ps(a,b);
  __m128i i = _mm_castps_si128(f);
  
  i = _mm_srai_epi32(i,31);
  i = _mm_xor_si128(i,MASK);
  
  f = _mm_castsi128_ps(i);
  
  return f;
}


void HelixHough::phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi, float* max_phi, float hel, __m128& phi_3_out, __m128& phi_4_out)
{
  __m128 helicity_vec = _mm_load1_ps(&(hel));
  
  __m128 x = _mm_load_ps(hit_x);
  __m128 y = _mm_load_ps(hit_y);
  
  __m128 d_min = _mm_load_ps(min_d);
  __m128 d_max = _mm_load_ps(max_d);
  __m128 k_min = _mm_load_ps(min_k);
  __m128 k_max = _mm_load_ps(max_k);
  
  __m128 hit_phi = _vec_atan2_ps(y,x);
  //if phi < 0, phi += 2*pi
  __m128 tmp1 = _mm_cmplt_ps(hit_phi, zero);
  __m128 tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  hit_phi = _mm_add_ps(hit_phi, tmp1);
  
  // min_d, max_k
  __m128 d = d_min;
  __m128 k = k_max;
  __m128 D = _mm_mul_ps(x,x);
  tmp1 = _mm_mul_ps(y,y);
  D = _mm_add_ps(D,tmp1);
  D = _vec_sqrt_ps(D);
  __m128 D_inv = _vec_rec_ps(D);
  __m128 ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  __m128 hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  __m128 neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  __m128 xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  __m128 xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  __m128 yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  __m128 yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  __m128 crossproduct = _mm_mul_ps(x, yk1);
  __m128 crosstemp = _mm_mul_ps(y, xk1);
  crossproduct = _mm_sub_ps(crossproduct, crosstemp);
  __m128 correct_helicity = compare_sign(crossproduct, helicity_vec);
  
  __m128 xk = _mm_and_ps(correct_helicity, xk1);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);
  xk = _mm_xor_ps(xk, tmp1);
  __m128 yk = _mm_and_ps(correct_helicity, yk1);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);
  yk = _mm_xor_ps(yk, tmp1);
  
  
  __m128 phi_1 = _vec_atan2_ps(yk, xk);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_1, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_1 = _mm_add_ps(phi_1, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_1 = _mm_andnot_ps(neg, phi_1);
  phi_1 = _mm_xor_ps(tmp1, phi_1);
  phi_3_out = phi_1;
  
  // max_d, max_k
  d = d_max;
  k = k_max;
  ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  xk = _mm_and_ps(correct_helicity, xk1);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);
  xk = _mm_xor_ps(xk, tmp1);
  yk = _mm_and_ps(correct_helicity, yk1);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);
  yk = _mm_xor_ps(yk, tmp1);
  
  __m128 phi_2 = _vec_atan2_ps(yk, xk);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_2 = _mm_add_ps(phi_2, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_2 = _mm_andnot_ps(neg, phi_2);
  phi_2 = _mm_xor_ps(tmp1, phi_2);
  phi_4_out = phi_2;
  
  
  // min_d, min_k
  d = d_min;
  k = k_min;
  ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  xk = _mm_and_ps(correct_helicity, xk1);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);
  xk = _mm_xor_ps(xk, tmp1);
  yk = _mm_and_ps(correct_helicity, yk1);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);
  yk = _mm_xor_ps(yk, tmp1);
  
  __m128 phi_3 = _vec_atan2_ps(yk, xk);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_3, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_3 = _mm_add_ps(phi_3, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_3 = _mm_andnot_ps(neg, phi_3);
  phi_3 = _mm_xor_ps(tmp1, phi_3);
  
  
  // max_d, min_k
  d = d_max;
  k = k_min;
  ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  xk = _mm_and_ps(correct_helicity, xk1);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);
  xk = _mm_xor_ps(xk, tmp1);
  yk = _mm_and_ps(correct_helicity, yk1);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);
  yk = _mm_xor_ps(yk, tmp1);
  
  __m128 phi_4 = _vec_atan2_ps(yk, xk);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_4, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_4 = _mm_add_ps(phi_4, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_4 = _mm_andnot_ps(neg, phi_4);
  phi_4 = _mm_xor_ps(tmp1, phi_4);
  
  ////////////////////////////////////////////////////////////////
  
  // check if phi overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_1, pi_over_two);
  tmp2 = _mm_cmplt_ps(phi_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_3, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_4, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  
  tmp2 = _mm_cmpgt_ps(phi_1, three_pi_over_two);
  __m128 tmp3 = _mm_cmpgt_ps(phi_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_3, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_4, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  
  tmp1 = _mm_and_ps(tmp1, tmp2);
  
  // tmp1 is now all ones if phi overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);
  
  __m128 tmp4 = _mm_cmpgt_ps(phi_1, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  __m128 tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_1 = _mm_sub_ps(phi_1, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_2 = _mm_sub_ps(phi_2, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_3, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_3 = _mm_sub_ps(phi_3, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_4, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_4 = _mm_sub_ps(phi_4, tmp3);
  
  
  // find the minimum phi
  __m128 phi_min = phi_1;
  tmp2 = _mm_cmplt_ps(phi_2, phi_min);
  tmp3 = _mm_and_ps(tmp2, phi_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);
  phi_min = _mm_xor_ps(phi_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_3, phi_min);
  tmp3 = _mm_and_ps(tmp2, phi_3);
  phi_min = _mm_andnot_ps(tmp2, phi_min);
  phi_min = _mm_xor_ps(phi_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_4, phi_min);
  tmp3 = _mm_and_ps(tmp2, phi_4);
  phi_min = _mm_andnot_ps(tmp2, phi_min);
  phi_min = _mm_xor_ps(phi_min, tmp3);
  
  // find the maximum phi
  __m128 phi_max = phi_1;
  tmp2 = _mm_cmpgt_ps(phi_2, phi_max);
  tmp3 = _mm_and_ps(tmp2, phi_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);
  phi_max = _mm_xor_ps(phi_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_3, phi_max);
  tmp3 = _mm_and_ps(tmp2, phi_3);
  phi_max = _mm_andnot_ps(tmp2, phi_max);
  phi_max = _mm_xor_ps(phi_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_4, phi_max);
  tmp3 = _mm_and_ps(tmp2, phi_4);
  phi_max = _mm_andnot_ps(tmp2, phi_max);
  phi_max = _mm_xor_ps(phi_max, tmp3);
  
  
  _mm_store_ps(min_phi, phi_min);
  _mm_store_ps(max_phi, phi_max);
}


void HelixHough::phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* max_k, float* min_phi, float* max_phi, float hel, __m128& phi_3, __m128& phi_4, __m128& phi_3_out, __m128& phi_4_out)
{
  __m128 helicity_vec = _mm_load1_ps(&(hel));
  
  __m128 x = _mm_load_ps(hit_x);
  __m128 y = _mm_load_ps(hit_y);
  
  __m128 d_min = _mm_load_ps(min_d);
  __m128 d_max = _mm_load_ps(max_d);
  __m128 k_max = _mm_load_ps(max_k);
  
  __m128 hit_phi = _vec_atan2_ps(y,x);
  //if phi < 0, phi += 2*pi
  __m128 tmp1 = _mm_cmplt_ps(hit_phi, zero);
  __m128 tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  hit_phi = _mm_add_ps(hit_phi, tmp1);
  
  // min_d, max_k
  __m128 d = d_min;
  __m128 k = k_max;
  __m128 D = _mm_mul_ps(x,x);
  tmp1 = _mm_mul_ps(y,y);
  D = _mm_add_ps(D,tmp1);
  D = _vec_sqrt_ps(D);
  __m128 D_inv = _vec_rec_ps(D);
  __m128 ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  __m128 hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  __m128 neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  __m128 xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  __m128 xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  __m128 yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  __m128 yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  __m128 crossproduct = _mm_mul_ps(x, yk1);
  __m128 crosstemp = _mm_mul_ps(y, xk1);
  crossproduct = _mm_sub_ps(crossproduct, crosstemp);
  __m128 correct_helicity = compare_sign(crossproduct, helicity_vec);
  
  __m128 xk = _mm_and_ps(correct_helicity, xk1);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);
  xk = _mm_xor_ps(xk, tmp1);
  __m128 yk = _mm_and_ps(correct_helicity, yk1);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);
  yk = _mm_xor_ps(yk, tmp1);
  
  
  __m128 phi_1 = _vec_atan2_ps(yk, xk);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_1, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_1 = _mm_add_ps(phi_1, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_1 = _mm_andnot_ps(neg, phi_1);
  phi_1 = _mm_xor_ps(tmp1, phi_1);
  phi_3_out = phi_1;
  
  // max_d, max_k
  d = d_max;
  k = k_max;
  ak = d;
  ak = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  tmp1 = _mm_mul_ps(D,D);
  tmp1 = _mm_mul_ps(tmp1, k);
  ak = _mm_add_ps(ak, tmp1);
  ak = _mm_mul_ps(ak, D_inv);
  ak = _mm_mul_ps(ak, one_o_2);
  hk = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);
  hk = _mm_mul_ps(hk,hk);
  tmp1 = _mm_mul_ps(ak,ak);
  hk = _mm_sub_ps(hk, tmp1);
  neg = _mm_cmple_ps(hk, zero);
  hk = _vec_sqrt_ps(hk);
  
  xk1 = _mm_mul_ps(ak, x);
  tmp1 = _mm_mul_ps(hk,y);
  xk2 = _mm_sub_ps(xk1, tmp1);
  xk1 = _mm_add_ps(xk1, tmp1);
  xk1 = _mm_mul_ps(xk1, D_inv);
  xk2 = _mm_mul_ps(xk2, D_inv);
  
  yk1 = _mm_mul_ps(ak, y);
  tmp1 = _mm_mul_ps(hk,x);
  yk2 = _mm_add_ps(yk1, tmp1);
  yk1 = _mm_sub_ps(yk1, tmp1);
  yk1 = _mm_mul_ps(yk1, D_inv);
  yk2 = _mm_mul_ps(yk2, D_inv);
  
  xk = _mm_and_ps(correct_helicity, xk1);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);
  xk = _mm_xor_ps(xk, tmp1);
  yk = _mm_and_ps(correct_helicity, yk1);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);
  yk = _mm_xor_ps(yk, tmp1);
  
  __m128 phi_2 = _vec_atan2_ps(yk, xk);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);
  phi_2 = _mm_add_ps(phi_2, tmp1);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);
  phi_2 = _mm_andnot_ps(neg, phi_2);
  phi_2 = _mm_xor_ps(tmp1, phi_2);
  phi_4_out = phi_2;
  
  // check if phi overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_1, pi_over_two);
  tmp2 = _mm_cmplt_ps(phi_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_3, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  tmp2 = _mm_cmplt_ps(phi_4, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);
  
  tmp2 = _mm_cmpgt_ps(phi_1, three_pi_over_two);
  __m128 tmp3 = _mm_cmpgt_ps(phi_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_3, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  tmp3 = _mm_cmpgt_ps(phi_4, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);
  
  tmp1 = _mm_and_ps(tmp1, tmp2);
  
  // tmp1 is now all ones if phi overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);
  
  __m128 tmp4 = _mm_cmpgt_ps(phi_1, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  __m128 tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_1 = _mm_sub_ps(phi_1, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_2 = _mm_sub_ps(phi_2, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_3, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_3 = _mm_sub_ps(phi_3, tmp3);
  
  tmp4 = _mm_cmpgt_ps(phi_4, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);
  tmp5 = _mm_andnot_ps(tmp4, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);
  phi_4 = _mm_sub_ps(phi_4, tmp3);
  
  
  // find the minimum phi
  __m128 phi_min = phi_1;
  tmp2 = _mm_cmplt_ps(phi_2, phi_min);
  tmp3 = _mm_and_ps(tmp2, phi_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);
  phi_min = _mm_xor_ps(phi_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_3, phi_min);
  tmp3 = _mm_and_ps(tmp2, phi_3);
  phi_min = _mm_andnot_ps(tmp2, phi_min);
  phi_min = _mm_xor_ps(phi_min, tmp3);
  tmp2 = _mm_cmplt_ps(phi_4, phi_min);
  tmp3 = _mm_and_ps(tmp2, phi_4);
  phi_min = _mm_andnot_ps(tmp2, phi_min);
  phi_min = _mm_xor_ps(phi_min, tmp3);
  
  // find the maximum phi
  __m128 phi_max = phi_1;
  tmp2 = _mm_cmpgt_ps(phi_2, phi_max);
  tmp3 = _mm_and_ps(tmp2, phi_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);
  phi_max = _mm_xor_ps(phi_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_3, phi_max);
  tmp3 = _mm_and_ps(tmp2, phi_3);
  phi_max = _mm_andnot_ps(tmp2, phi_max);
  phi_max = _mm_xor_ps(phi_max, tmp3);
  tmp2 = _mm_cmpgt_ps(phi_4, phi_max);
  tmp3 = _mm_and_ps(tmp2, phi_4);
  phi_max = _mm_andnot_ps(tmp2, phi_max);
  phi_max = _mm_xor_ps(phi_max, tmp3);
  
  
  _mm_store_ps(min_phi, phi_min);
  _mm_store_ps(max_phi, phi_max);
}


void HelixHough::phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi, float* max_phi, float* min_phi_2, float* max_phi_2, float hel, __m128& phi_3_out, __m128& phi_4_out, float* hit_x_2, float* hit_y_2, __m128& phi_3_out_2, __m128& phi_4_out_2)
{
  __m128 helicity_vec = _mm_load1_ps(&(hel));
  
  __m128 x = _mm_load_ps(hit_x);                                    __m128 x_2 = _mm_load_ps(hit_x_2);
  __m128 y = _mm_load_ps(hit_y);                                    __m128 y_2 = _mm_load_ps(hit_y_2);
  
  __m128 d_min = _mm_load_ps(min_d);
  __m128 d_max = _mm_load_ps(max_d);
  __m128 k_min = _mm_load_ps(min_k);
  __m128 k_max = _mm_load_ps(max_k);
  
  __m128 hit_phi = _vec_atan2_ps(y,x);                               __m128 hit_phi_2 = _vec_atan2_ps(y_2,x_2);
  //if phi < 0, phi += 2*pi
  __m128 tmp1 = _mm_cmplt_ps(hit_phi, zero);                         __m128 tmp1_2 = _mm_cmplt_ps(hit_phi_2, zero);
  __m128 tmp2 = _mm_and_ps(tmp1, twopi);                             __m128 tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  hit_phi = _mm_add_ps(hit_phi, tmp1);                               hit_phi_2 = _mm_add_ps(hit_phi_2, tmp1_2);
  
  // min_d, max_k
  __m128 d = d_min;                                                  
  __m128 k = k_max;                                                  
  __m128 D = _mm_mul_ps(x,x);                                        __m128 D_2 = _mm_mul_ps(x_2,x_2);
  tmp1 = _mm_mul_ps(y,y);                                            tmp1_2 = _mm_mul_ps(y_2,y_2);
  D = _mm_add_ps(D,tmp1);                                            D_2 = _mm_add_ps(D_2,tmp1_2);
  D = _vec_sqrt_ps(D);                                                D_2 = _vec_sqrt_ps(D_2);
  __m128 D_inv = _vec_rec_ps(D);                                     __m128 D_inv_2 = _vec_rec_ps(D_2);
  __m128 ak = d;                                                     __m128 ak_2 = d;
  ak = _mm_mul_ps(d, two);                                           ak_2 = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);                                            tmp1_2 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  tmp1 = _mm_mul_ps(D,D);                                            tmp1_2 = _mm_mul_ps(D_2,D_2);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  ak = _mm_mul_ps(ak, D_inv);                                        ak_2 = _mm_mul_ps(ak_2, D_inv_2);
  ak = _mm_mul_ps(ak, one_o_2);                                      ak_2 = _mm_mul_ps(ak_2, one_o_2);
  __m128 hk = _mm_mul_ps(d,k);                                       __m128 hk_2 = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);                                          hk_2 = _mm_add_ps(hk_2, one);
  hk = _mm_mul_ps(hk,hk);                                            hk_2 = _mm_mul_ps(hk_2,hk_2);
  tmp1 = _mm_mul_ps(ak,ak);                                          tmp1_2 = _mm_mul_ps(ak_2,ak_2);
  hk = _mm_sub_ps(hk, tmp1);                                         hk_2 = _mm_sub_ps(hk_2, tmp1_2);
  __m128 neg = _mm_cmple_ps(hk, zero);                               __m128 neg_2 = _mm_cmple_ps(hk_2, zero);
  hk = _vec_sqrt_ps(hk);                                              hk_2 = _vec_sqrt_ps(hk_2);
  
  __m128 xk1 = _mm_mul_ps(ak, x);                                    __m128 xk1_2 = _mm_mul_ps(ak_2, x_2);
  tmp1 = _mm_mul_ps(hk,y);                                           tmp1_2 = _mm_mul_ps(hk_2,y_2);
  __m128 xk2 = _mm_sub_ps(xk1, tmp1);                                __m128 xk2_2 = _mm_sub_ps(xk1_2, tmp1_2);
  xk1 = _mm_add_ps(xk1, tmp1);                                       xk1_2 = _mm_add_ps(xk1_2, tmp1_2);
  xk1 = _mm_mul_ps(xk1, D_inv);                                      xk1_2 = _mm_mul_ps(xk1_2, D_inv_2);
  xk2 = _mm_mul_ps(xk2, D_inv);                                      xk2_2 = _mm_mul_ps(xk2_2, D_inv_2);
  
  __m128 yk1 = _mm_mul_ps(ak, y);                                    __m128 yk1_2 = _mm_mul_ps(ak_2, y_2);
  tmp1 = _mm_mul_ps(hk,x);                                           tmp1_2 = _mm_mul_ps(hk_2,x_2);
  __m128 yk2 = _mm_add_ps(yk1, tmp1);                                __m128 yk2_2 = _mm_add_ps(yk1_2, tmp1_2);
  yk1 = _mm_sub_ps(yk1, tmp1);                                       yk1_2 = _mm_sub_ps(yk1_2, tmp1_2);
  yk1 = _mm_mul_ps(yk1, D_inv);                                      yk1_2 = _mm_mul_ps(yk1_2, D_inv_2);
  yk2 = _mm_mul_ps(yk2, D_inv);                                      yk2_2 = _mm_mul_ps(yk2_2, D_inv_2);
  
  __m128 crossproduct = _mm_mul_ps(x, yk1);                               __m128 crossproduct_2 = _mm_mul_ps(x_2, yk1_2);
  __m128 crosstemp = _mm_mul_ps(y, xk1);                                  __m128 crosstemp_2 = _mm_mul_ps(y_2, xk1_2);
  crossproduct = _mm_sub_ps(crossproduct, crosstemp);                     crossproduct_2 = _mm_sub_ps(crossproduct_2, crosstemp_2);
  __m128 correct_helicity = compare_sign(crossproduct, helicity_vec);     __m128 correct_helicity_2 = compare_sign(crossproduct_2, helicity_vec);
  
  __m128 xk = _mm_and_ps(correct_helicity, xk1);                     __m128 xk_2 = _mm_and_ps(correct_helicity_2, xk1_2);  
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, xk2_2);
  xk = _mm_xor_ps(xk, tmp1);                                         xk_2 = _mm_xor_ps(xk_2, tmp1_2);
  __m128 yk = _mm_and_ps(correct_helicity, yk1);                     __m128 yk_2 = _mm_and_ps(correct_helicity_2, yk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, yk2_2);
  yk = _mm_xor_ps(yk, tmp1);                                         yk_2 = _mm_xor_ps(yk_2, tmp1_2);
  
  
  __m128 phi_1 = _vec_atan2_ps(yk, xk);                              __m128 phi_1_2 = _vec_atan2_ps(yk_2, xk_2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_1, zero);                                  tmp1_2 = _mm_cmplt_ps(phi_1_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  phi_1 = _mm_add_ps(phi_1, tmp1);                                   phi_1_2 = _mm_add_ps(phi_1_2, tmp1_2);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);                                   tmp1_2 = _mm_and_ps(neg_2, hit_phi_2);
  phi_1 = _mm_andnot_ps(neg, phi_1);                                 phi_1_2 = _mm_andnot_ps(neg_2, phi_1_2);
  phi_1 = _mm_xor_ps(tmp1, phi_1);                                   phi_1_2 = _mm_xor_ps(tmp1_2, phi_1_2);
  phi_3_out = phi_1;                                                 phi_3_out_2 = phi_1_2;
  
  // max_d, max_k
  d = d_max;                                                         
  k = k_max;                                                         
  ak = d;                                                            ak_2 = d;
  ak = _mm_mul_ps(d, two);                                           ak_2 = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);                                            tmp1_2 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  tmp1 = _mm_mul_ps(D,D);                                            tmp1_2 = _mm_mul_ps(D_2,D_2);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  ak = _mm_mul_ps(ak, D_inv);                                        ak_2 = _mm_mul_ps(ak_2, D_inv_2);
  ak = _mm_mul_ps(ak, one_o_2);                                      ak_2 = _mm_mul_ps(ak_2, one_o_2);
  hk = _mm_mul_ps(d,k);                                              hk_2 = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);                                          hk_2 = _mm_add_ps(hk_2, one);
  hk = _mm_mul_ps(hk,hk);                                            hk_2 = _mm_mul_ps(hk_2,hk_2);
  tmp1 = _mm_mul_ps(ak,ak);                                          tmp1_2 = _mm_mul_ps(ak_2,ak_2);
  hk = _mm_sub_ps(hk, tmp1);                                         hk_2 = _mm_sub_ps(hk_2, tmp1_2);
  neg = _mm_cmple_ps(hk, zero);                                      neg_2 = _mm_cmple_ps(hk_2, zero);
  hk = _vec_sqrt_ps(hk);                                              hk_2 = _vec_sqrt_ps(hk_2);
  
  xk1 = _mm_mul_ps(ak, x);                                           xk1_2 = _mm_mul_ps(ak_2, x_2);
  tmp1 = _mm_mul_ps(hk,y);                                           tmp1_2 = _mm_mul_ps(hk_2,y_2);
  xk2 = _mm_sub_ps(xk1, tmp1);                                       xk2_2 = _mm_sub_ps(xk1_2, tmp1_2);
  xk1 = _mm_add_ps(xk1, tmp1);                                       xk1_2 = _mm_add_ps(xk1_2, tmp1_2);
  xk1 = _mm_mul_ps(xk1, D_inv);                                      xk1_2 = _mm_mul_ps(xk1_2, D_inv_2);
  xk2 = _mm_mul_ps(xk2, D_inv);                                      xk2_2 = _mm_mul_ps(xk2_2, D_inv_2);
  
  yk1 = _mm_mul_ps(ak, y);                                           yk1_2 = _mm_mul_ps(ak_2, y_2);
  tmp1 = _mm_mul_ps(hk,x);                                           tmp1_2 = _mm_mul_ps(hk_2, x_2);
  yk2 = _mm_add_ps(yk1, tmp1);                                       yk2_2 = _mm_add_ps(yk1_2, tmp1_2);
  yk1 = _mm_sub_ps(yk1, tmp1);                                       yk1_2 = _mm_sub_ps(yk1_2, tmp1_2);
  yk1 = _mm_mul_ps(yk1, D_inv);                                      yk1_2 = _mm_mul_ps(yk1_2, D_inv_2);
  yk2 = _mm_mul_ps(yk2, D_inv);                                      yk2_2 = _mm_mul_ps(yk2_2, D_inv_2);
  
  xk = _mm_and_ps(correct_helicity, xk1);                            xk_2 = _mm_and_ps(correct_helicity_2, xk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, xk2_2);
  xk = _mm_xor_ps(xk, tmp1);                                         xk_2 = _mm_xor_ps(xk_2, tmp1_2);
  yk = _mm_and_ps(correct_helicity, yk1);                            yk_2 = _mm_and_ps(correct_helicity_2, yk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, yk2_2);
  yk = _mm_xor_ps(yk, tmp1);                                         yk_2 = _mm_xor_ps(yk_2, tmp1_2);
  
  __m128 phi_2 = _vec_atan2_ps(yk, xk);                              __m128 phi_2_2 = _vec_atan2_ps(yk_2, xk_2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_2, zero);                                  tmp1_2 = _mm_cmplt_ps(phi_2_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  phi_2 = _mm_add_ps(phi_2, tmp1);                                   phi_2_2 = _mm_add_ps(phi_2_2, tmp1_2);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);                                   tmp1_2 = _mm_and_ps(neg_2, hit_phi_2);
  phi_2 = _mm_andnot_ps(neg, phi_2);                                 phi_2_2 = _mm_andnot_ps(neg_2, phi_2_2);
  phi_2 = _mm_xor_ps(tmp1, phi_2);                                   phi_2_2 = _mm_xor_ps(tmp1_2, phi_2_2);
  phi_4_out = phi_2;                                                 phi_4_out_2 = phi_2_2;
  
  
  // min_d, min_k
  d = d_min;
  k = k_min;
  ak = d;                                                            ak_2 = d;
  ak = _mm_mul_ps(d, two);                                           ak_2 = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);                                            tmp1_2 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  tmp1 = _mm_mul_ps(D,D);                                            tmp1_2 = _mm_mul_ps(D_2,D_2);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  ak = _mm_mul_ps(ak, D_inv);                                        ak_2 = _mm_mul_ps(ak_2, D_inv_2);
  ak = _mm_mul_ps(ak, one_o_2);                                      ak_2 = _mm_mul_ps(ak_2, one_o_2);
  hk = _mm_mul_ps(d,k);                                              hk_2 = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);                                          hk_2 = _mm_add_ps(hk_2, one);
  hk = _mm_mul_ps(hk,hk);                                            hk_2 = _mm_mul_ps(hk_2,hk_2);
  tmp1 = _mm_mul_ps(ak,ak);                                          tmp1_2 = _mm_mul_ps(ak_2,ak_2);
  hk = _mm_sub_ps(hk, tmp1);                                         hk_2 = _mm_sub_ps(hk_2, tmp1_2);
  neg = _mm_cmple_ps(hk, zero);                                      neg_2 = _mm_cmple_ps(hk_2, zero);
  hk = _vec_sqrt_ps(hk);                                              hk_2 = _vec_sqrt_ps(hk_2);
  
  xk1 = _mm_mul_ps(ak, x);                                           xk1_2 = _mm_mul_ps(ak_2, x_2);
  tmp1 = _mm_mul_ps(hk, y);                                           tmp1_2 = _mm_mul_ps(hk_2, y_2);
  xk2 = _mm_sub_ps(xk1, tmp1);                                       xk2_2 = _mm_sub_ps(xk1_2, tmp1_2);
  xk1 = _mm_add_ps(xk1, tmp1);                                       xk1_2 = _mm_add_ps(xk1_2, tmp1_2);
  xk1 = _mm_mul_ps(xk1, D_inv);                                      xk1_2 = _mm_mul_ps(xk1_2, D_inv_2);
  xk2 = _mm_mul_ps(xk2, D_inv);                                      xk2_2 = _mm_mul_ps(xk2_2, D_inv_2);
  
  yk1 = _mm_mul_ps(ak, y);                                           yk1_2 = _mm_mul_ps(ak_2, y_2);
  tmp1 = _mm_mul_ps(hk,x);                                           tmp1_2 = _mm_mul_ps(hk_2, x_2);
  yk2 = _mm_add_ps(yk1, tmp1);                                       yk2_2 = _mm_add_ps(yk1_2, tmp1_2);
  yk1 = _mm_sub_ps(yk1, tmp1);                                       yk1_2 = _mm_sub_ps(yk1_2, tmp1_2);
  yk1 = _mm_mul_ps(yk1, D_inv);                                      yk1_2 = _mm_mul_ps(yk1_2, D_inv_2);
  yk2 = _mm_mul_ps(yk2, D_inv);                                      yk2_2 = _mm_mul_ps(yk2_2, D_inv_2);
  
  xk = _mm_and_ps(correct_helicity, xk1);                            xk_2 = _mm_and_ps(correct_helicity_2, xk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, xk2_2);
  xk = _mm_xor_ps(xk, tmp1);                                         xk_2 = _mm_xor_ps(xk_2, tmp1_2);
  yk = _mm_and_ps(correct_helicity, yk1);                            yk_2 = _mm_and_ps(correct_helicity_2, yk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, yk2_2);
  yk = _mm_xor_ps(yk, tmp1);                                         yk_2 = _mm_xor_ps(yk_2, tmp1_2);
  
  __m128 phi_3 = _vec_atan2_ps(yk, xk);                              __m128 phi_3_2 = _vec_atan2_ps(yk_2, xk_2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_3, zero);                                  tmp1_2 = _mm_cmplt_ps(phi_3_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  phi_3 = _mm_add_ps(phi_3, tmp1);                                   phi_3_2 = _mm_add_ps(phi_3_2, tmp1_2);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);                                   tmp1_2 = _mm_and_ps(neg_2, hit_phi_2);
  phi_3 = _mm_andnot_ps(neg, phi_3);                                 phi_3_2 = _mm_andnot_ps(neg_2, phi_3_2);
  phi_3 = _mm_xor_ps(tmp1, phi_3);                                   phi_3_2 = _mm_xor_ps(tmp1_2, phi_3_2);
  
  
  // max_d, min_k
  d = d_max;
  k = k_min;
  ak = d;                                                            ak_2 = d;
  ak = _mm_mul_ps(d, two);                                           ak_2 = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);                                            tmp1_2 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  tmp1 = _mm_mul_ps(D,D);                                            tmp1_2 = _mm_mul_ps(D_2,D_2);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  ak = _mm_mul_ps(ak, D_inv);                                        ak_2 = _mm_mul_ps(ak_2, D_inv_2);
  ak = _mm_mul_ps(ak, one_o_2);                                      ak_2 = _mm_mul_ps(ak_2, one_o_2);
  hk = _mm_mul_ps(d,k);                                              hk_2 = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);                                          hk_2 = _mm_add_ps(hk_2, one);
  hk = _mm_mul_ps(hk,hk);                                            hk_2 = _mm_mul_ps(hk_2,hk_2);
  tmp1 = _mm_mul_ps(ak,ak);                                          tmp1_2 = _mm_mul_ps(ak_2,ak_2);
  hk = _mm_sub_ps(hk, tmp1);                                         hk_2 = _mm_sub_ps(hk_2, tmp1_2);
  neg = _mm_cmple_ps(hk, zero);                                      neg_2 = _mm_cmple_ps(hk_2, zero);
  hk = _vec_sqrt_ps(hk);                                              hk_2 = _vec_sqrt_ps(hk_2);
  
  xk1 = _mm_mul_ps(ak, x);                                           xk1_2 = _mm_mul_ps(ak_2, x_2);
  tmp1 = _mm_mul_ps(hk, y);                                           tmp1_2 = _mm_mul_ps(hk_2, y_2);
  xk2 = _mm_sub_ps(xk1, tmp1);                                       xk2_2 = _mm_sub_ps(xk1_2, tmp1_2);
  xk1 = _mm_add_ps(xk1, tmp1);                                       xk1_2 = _mm_add_ps(xk1_2, tmp1_2);
  xk1 = _mm_mul_ps(xk1, D_inv);                                      xk1_2 = _mm_mul_ps(xk1_2, D_inv_2);
  xk2 = _mm_mul_ps(xk2, D_inv);                                      xk2_2 = _mm_mul_ps(xk2_2, D_inv_2);
  
  yk1 = _mm_mul_ps(ak, y);                                           yk1_2 = _mm_mul_ps(ak_2, y_2);
  tmp1 = _mm_mul_ps(hk,x);                                           tmp1_2 = _mm_mul_ps(hk_2, x_2);
  yk2 = _mm_add_ps(yk1, tmp1);                                       yk2_2 = _mm_add_ps(yk1_2, tmp1_2);
  yk1 = _mm_sub_ps(yk1, tmp1);                                       yk1_2 = _mm_sub_ps(yk1_2, tmp1_2);
  yk1 = _mm_mul_ps(yk1, D_inv);                                      yk1_2 = _mm_mul_ps(yk1_2, D_inv_2);
  yk2 = _mm_mul_ps(yk2, D_inv);                                      yk2_2 = _mm_mul_ps(yk2_2, D_inv_2);
  
  xk = _mm_and_ps(correct_helicity, xk1);                            xk_2 = _mm_and_ps(correct_helicity_2, xk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, xk2_2);
  xk = _mm_xor_ps(xk, tmp1);                                         xk_2 = _mm_xor_ps(xk_2, tmp1_2);
  yk = _mm_and_ps(correct_helicity, yk1);                            yk_2 = _mm_and_ps(correct_helicity_2, yk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, yk2_2);
  yk = _mm_xor_ps(yk, tmp1);                                         yk_2 = _mm_xor_ps(yk_2, tmp1_2);
  
  __m128 phi_4 = _vec_atan2_ps(yk, xk);                              __m128 phi_4_2 = _vec_atan2_ps(yk_2, xk_2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_4, zero);                                  tmp1_2 = _mm_cmplt_ps(phi_4_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  phi_4 = _mm_add_ps(phi_4, tmp1);                                   phi_4_2 = _mm_add_ps(phi_4_2, tmp1_2);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);                                   tmp1_2 = _mm_and_ps(neg_2, hit_phi_2);
  phi_4 = _mm_andnot_ps(neg, phi_4);                                 phi_4_2 = _mm_andnot_ps(neg_2, phi_4_2);
  phi_4 = _mm_xor_ps(tmp1, phi_4);                                   phi_4_2 = _mm_xor_ps(tmp1_2, phi_4_2);
  
  ////////////////////////////////////////////////////////////////
  
  // check if phi overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_1, pi_over_two);                           tmp1_2 = _mm_cmplt_ps(phi_1_2, pi_over_two);
  tmp2 = _mm_cmplt_ps(phi_2, pi_over_two);                           tmp2_2 = _mm_cmplt_ps(phi_2_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);                                      tmp1_2 = _mm_or_ps(tmp1_2, tmp2_2);
  tmp2 = _mm_cmplt_ps(phi_3, pi_over_two);                           tmp2_2 = _mm_cmplt_ps(phi_3_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);                                      tmp1_2 = _mm_or_ps(tmp1_2, tmp2_2);
  tmp2 = _mm_cmplt_ps(phi_4, pi_over_two);                           tmp2_2 = _mm_cmplt_ps(phi_4_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);                                      tmp1_2 = _mm_or_ps(tmp1_2, tmp2_2);
  
  tmp2 = _mm_cmpgt_ps(phi_1, three_pi_over_two);                     tmp2_2 = _mm_cmpgt_ps(phi_1_2, three_pi_over_two);
  __m128 tmp3 = _mm_cmpgt_ps(phi_2, three_pi_over_two);              __m128 tmp3_2 = _mm_cmpgt_ps(phi_2_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);                                      tmp2_2 = _mm_or_ps(tmp2_2, tmp3_2);
  tmp3 = _mm_cmpgt_ps(phi_3, three_pi_over_two);                     tmp3_2 = _mm_cmpgt_ps(phi_3_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);                                      tmp2_2 = _mm_or_ps(tmp2_2, tmp3_2);
  tmp3 = _mm_cmpgt_ps(phi_4, three_pi_over_two);                     tmp3_2 = _mm_cmpgt_ps(phi_4_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);                                      tmp2_2 = _mm_or_ps(tmp2_2, tmp3_2);
  
  tmp1 = _mm_and_ps(tmp1, tmp2);                                     tmp1_2 = _mm_and_ps(tmp1_2, tmp2_2);
  
  // tmp1 is now all ones if phi overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);                                  tmp3_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);                                     tmp2_2 = _mm_xor_ps(tmp2_2, tmp3_2);
  
  __m128 tmp4 = _mm_cmpgt_ps(phi_1, three_pi_over_two);              __m128 tmp4_2 = _mm_cmpgt_ps(phi_1_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  __m128 tmp5 = _mm_andnot_ps(tmp4, zero);                           __m128 tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_1 = _mm_sub_ps(phi_1, tmp3);                                   phi_1_2 = _mm_sub_ps(phi_1_2, tmp3_2);
  
  tmp4 = _mm_cmpgt_ps(phi_2, three_pi_over_two);                     tmp4_2 = _mm_cmpgt_ps(phi_2_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  tmp5 = _mm_andnot_ps(tmp4, zero);                                  tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_2 = _mm_sub_ps(phi_2, tmp3);                                   phi_2_2 = _mm_sub_ps(phi_2_2, tmp3_2);
  
  tmp4 = _mm_cmpgt_ps(phi_3, three_pi_over_two);                     tmp4_2 = _mm_cmpgt_ps(phi_3_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  tmp5 = _mm_andnot_ps(tmp4, zero);                                  tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_3 = _mm_sub_ps(phi_3, tmp3);                                   phi_3_2 = _mm_sub_ps(phi_3_2, tmp3_2);
  
  tmp4 = _mm_cmpgt_ps(phi_4, three_pi_over_two);                     tmp4_2 = _mm_cmpgt_ps(phi_4_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  tmp5 = _mm_andnot_ps(tmp4, zero);                                  tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_4 = _mm_sub_ps(phi_4, tmp3);                                   phi_4_2 = _mm_sub_ps(phi_4_2, tmp3_2);
  
  
  // find the minimum phi
  __m128 phi_min = phi_1;                                            __m128 phi_min_2 = phi_1_2;
  tmp2 = _mm_cmplt_ps(phi_2, phi_min);                               tmp2_2 = _mm_cmplt_ps(phi_2_2, phi_min_2);
  tmp3 = _mm_and_ps(tmp2, phi_2);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_2_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);                            phi_min_2 = _mm_andnot_ps(tmp2_2, phi_min_2);
  phi_min = _mm_xor_ps(phi_min, tmp3);                               phi_min_2 = _mm_xor_ps(phi_min_2, tmp3_2);
  tmp2 = _mm_cmplt_ps(phi_3, phi_min);                               tmp2_2 = _mm_cmplt_ps(phi_3_2, phi_min_2);
  tmp3 = _mm_and_ps(tmp2, phi_3);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_3_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);                            phi_min_2 = _mm_andnot_ps(tmp2_2, phi_min_2);
  phi_min = _mm_xor_ps(phi_min, tmp3);                               phi_min_2 = _mm_xor_ps(phi_min_2, tmp3_2);
  tmp2 = _mm_cmplt_ps(phi_4, phi_min);                               tmp2_2 = _mm_cmplt_ps(phi_4_2, phi_min_2);
  tmp3 = _mm_and_ps(tmp2, phi_4);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_4_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);                            phi_min_2 = _mm_andnot_ps(tmp2_2, phi_min_2);
  phi_min = _mm_xor_ps(phi_min, tmp3);                               phi_min_2 = _mm_xor_ps(phi_min_2, tmp3_2);
  
  // find the maximum phi
  __m128 phi_max = phi_1;                                            __m128 phi_max_2 = phi_1_2;
  tmp2 = _mm_cmpgt_ps(phi_2, phi_max);                               tmp2_2 = _mm_cmpgt_ps(phi_2_2, phi_max_2);
  tmp3 = _mm_and_ps(tmp2, phi_2);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_2_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);                            phi_max_2 = _mm_andnot_ps(tmp2_2, phi_max_2);
  phi_max = _mm_xor_ps(phi_max, tmp3);                               phi_max_2 = _mm_xor_ps(phi_max_2, tmp3_2);
  tmp2 = _mm_cmpgt_ps(phi_3, phi_max);                               tmp2_2 = _mm_cmpgt_ps(phi_3_2, phi_max_2);
  tmp3 = _mm_and_ps(tmp2, phi_3);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_3_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);                            phi_max_2 = _mm_andnot_ps(tmp2_2, phi_max_2);
  phi_max = _mm_xor_ps(phi_max, tmp3);                               phi_max_2 = _mm_xor_ps(phi_max_2, tmp3_2);
  tmp2 = _mm_cmpgt_ps(phi_4, phi_max);                               tmp2_2 = _mm_cmpgt_ps(phi_4_2, phi_max_2);
  tmp3 = _mm_and_ps(tmp2, phi_4);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_4_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);                            phi_max_2 = _mm_andnot_ps(tmp2_2, phi_max_2);
  phi_max = _mm_xor_ps(phi_max, tmp3);                               phi_max_2 = _mm_xor_ps(phi_max_2, tmp3_2);
  
  
  _mm_store_ps(min_phi, phi_min);                                    _mm_store_ps(min_phi_2, phi_min_2);
  _mm_store_ps(max_phi, phi_max);                                    _mm_store_ps(max_phi_2, phi_max_2);
}


void HelixHough::phiRange_sse(float* hit_x, float* hit_y, float* min_d, float* max_d, float* min_k, float* max_k, float* min_phi, float* max_phi, float* min_phi_2, float* max_phi_2, float hel, __m128& phi_3, __m128& phi_4, __m128& phi_3_out, __m128& phi_4_out, float* hit_x_2, float* hit_y_2, __m128& phi_3_2, __m128& phi_4_2, __m128& phi_3_out_2, __m128& phi_4_out_2)
{
  __m128 helicity_vec = _mm_load1_ps(&(hel));
  
  __m128 x = _mm_load_ps(hit_x);                                    __m128 x_2 = _mm_load_ps(hit_x_2);
  __m128 y = _mm_load_ps(hit_y);                                    __m128 y_2 = _mm_load_ps(hit_y_2);
  
  __m128 d_min = _mm_load_ps(min_d);
  __m128 d_max = _mm_load_ps(max_d);
  __m128 k_max = _mm_load_ps(max_k);
  
  __m128 hit_phi = _vec_atan2_ps(y,x);                               __m128 hit_phi_2 = _vec_atan2_ps(y_2,x_2);
  //if phi < 0, phi += 2*pi
  __m128 tmp1 = _mm_cmplt_ps(hit_phi, zero);                         __m128 tmp1_2 = _mm_cmplt_ps(hit_phi_2, zero);
  __m128 tmp2 = _mm_and_ps(tmp1, twopi);                             __m128 tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  hit_phi = _mm_add_ps(hit_phi, tmp1);                               hit_phi_2 = _mm_add_ps(hit_phi_2, tmp1_2);
  
  // min_d, max_k
  __m128 d = d_min;                                                  
  __m128 k = k_max;                                                  
  __m128 D = _mm_mul_ps(x,x);                                        __m128 D_2 = _mm_mul_ps(x_2,x_2);
  tmp1 = _mm_mul_ps(y,y);                                            tmp1_2 = _mm_mul_ps(y_2,y_2);
  D = _mm_add_ps(D,tmp1);                                            D_2 = _mm_add_ps(D_2,tmp1_2);
  D = _vec_sqrt_ps(D);                                                D_2 = _vec_sqrt_ps(D_2);
  __m128 D_inv = _vec_rec_ps(D);                                     __m128 D_inv_2 = _vec_rec_ps(D_2);
  __m128 ak = d;                                                     __m128 ak_2 = d;
  ak = _mm_mul_ps(d, two);                                           ak_2 = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);                                            tmp1_2 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  tmp1 = _mm_mul_ps(D,D);                                            tmp1_2 = _mm_mul_ps(D_2,D_2);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  ak = _mm_mul_ps(ak, D_inv);                                        ak_2 = _mm_mul_ps(ak_2, D_inv_2);
  ak = _mm_mul_ps(ak, one_o_2);                                      ak_2 = _mm_mul_ps(ak_2, one_o_2);
  __m128 hk = _mm_mul_ps(d,k);                                       __m128 hk_2 = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);                                          hk_2 = _mm_add_ps(hk_2, one);
  hk = _mm_mul_ps(hk,hk);                                            hk_2 = _mm_mul_ps(hk_2,hk_2);
  tmp1 = _mm_mul_ps(ak,ak);                                          tmp1_2 = _mm_mul_ps(ak_2,ak_2);
  hk = _mm_sub_ps(hk, tmp1);                                         hk_2 = _mm_sub_ps(hk_2, tmp1_2);
  __m128 neg = _mm_cmple_ps(hk, zero);                               __m128 neg_2 = _mm_cmple_ps(hk_2, zero);
  hk = _vec_sqrt_ps(hk);                                              hk_2 = _vec_sqrt_ps(hk_2);
  
  __m128 xk1 = _mm_mul_ps(ak, x);                                    __m128 xk1_2 = _mm_mul_ps(ak_2, x_2);
  tmp1 = _mm_mul_ps(hk,y);                                           tmp1_2 = _mm_mul_ps(hk_2,y_2);
  __m128 xk2 = _mm_sub_ps(xk1, tmp1);                                __m128 xk2_2 = _mm_sub_ps(xk1_2, tmp1_2);
  xk1 = _mm_add_ps(xk1, tmp1);                                       xk1_2 = _mm_add_ps(xk1_2, tmp1_2);
  xk1 = _mm_mul_ps(xk1, D_inv);                                      xk1_2 = _mm_mul_ps(xk1_2, D_inv_2);
  xk2 = _mm_mul_ps(xk2, D_inv);                                      xk2_2 = _mm_mul_ps(xk2_2, D_inv_2);
  
  __m128 yk1 = _mm_mul_ps(ak, y);                                    __m128 yk1_2 = _mm_mul_ps(ak_2, y_2);
  tmp1 = _mm_mul_ps(hk,x);                                           tmp1_2 = _mm_mul_ps(hk_2,x_2);
  __m128 yk2 = _mm_add_ps(yk1, tmp1);                                __m128 yk2_2 = _mm_add_ps(yk1_2, tmp1_2);
  yk1 = _mm_sub_ps(yk1, tmp1);                                       yk1_2 = _mm_sub_ps(yk1_2, tmp1_2);
  yk1 = _mm_mul_ps(yk1, D_inv);                                      yk1_2 = _mm_mul_ps(yk1_2, D_inv_2);
  yk2 = _mm_mul_ps(yk2, D_inv);                                      yk2_2 = _mm_mul_ps(yk2_2, D_inv_2);
  
  __m128 crossproduct = _mm_mul_ps(x, yk1);                               __m128 crossproduct_2 = _mm_mul_ps(x_2, yk1_2);
  __m128 crosstemp = _mm_mul_ps(y, xk1);                                  __m128 crosstemp_2 = _mm_mul_ps(y_2, xk1_2);
  crossproduct = _mm_sub_ps(crossproduct, crosstemp);                     crossproduct_2 = _mm_sub_ps(crossproduct_2, crosstemp_2);
  __m128 correct_helicity = compare_sign(crossproduct, helicity_vec);     __m128 correct_helicity_2 = compare_sign(crossproduct_2, helicity_vec);
  
  __m128 xk = _mm_and_ps(correct_helicity, xk1);                     __m128 xk_2 = _mm_and_ps(correct_helicity_2, xk1_2);  
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, xk2_2);
  xk = _mm_xor_ps(xk, tmp1);                                         xk_2 = _mm_xor_ps(xk_2, tmp1_2);
  __m128 yk = _mm_and_ps(correct_helicity, yk1);                     __m128 yk_2 = _mm_and_ps(correct_helicity_2, yk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, yk2_2);
  yk = _mm_xor_ps(yk, tmp1);                                         yk_2 = _mm_xor_ps(yk_2, tmp1_2);
  
  
  __m128 phi_1 = _vec_atan2_ps(yk, xk);                              __m128 phi_1_2 = _vec_atan2_ps(yk_2, xk_2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_1, zero);                                  tmp1_2 = _mm_cmplt_ps(phi_1_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  phi_1 = _mm_add_ps(phi_1, tmp1);                                   phi_1_2 = _mm_add_ps(phi_1_2, tmp1_2);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);                                   tmp1_2 = _mm_and_ps(neg_2, hit_phi_2);
  phi_1 = _mm_andnot_ps(neg, phi_1);                                 phi_1_2 = _mm_andnot_ps(neg_2, phi_1_2);
  phi_1 = _mm_xor_ps(tmp1, phi_1);                                   phi_1_2 = _mm_xor_ps(tmp1_2, phi_1_2);
  phi_3_out = phi_1;                                                 phi_3_out_2 = phi_1_2;
  
  // max_d, max_k
  d = d_max;                                                         
  k = k_max;                                                         
  ak = d;                                                            ak_2 = d;
  ak = _mm_mul_ps(d, two);                                           ak_2 = _mm_mul_ps(d, two);
  tmp1 = _mm_mul_ps(d,d);                                            tmp1_2 = _mm_mul_ps(d,d);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  tmp1 = _mm_mul_ps(D,D);                                            tmp1_2 = _mm_mul_ps(D_2,D_2);
  tmp1 = _mm_mul_ps(tmp1, k);                                        tmp1_2 = _mm_mul_ps(tmp1_2, k);
  ak = _mm_add_ps(ak, tmp1);                                         ak_2 = _mm_add_ps(ak_2, tmp1_2);
  ak = _mm_mul_ps(ak, D_inv);                                        ak_2 = _mm_mul_ps(ak_2, D_inv_2);
  ak = _mm_mul_ps(ak, one_o_2);                                      ak_2 = _mm_mul_ps(ak_2, one_o_2);
  hk = _mm_mul_ps(d,k);                                              hk_2 = _mm_mul_ps(d,k);
  hk = _mm_add_ps(hk, one);                                          hk_2 = _mm_add_ps(hk_2, one);
  hk = _mm_mul_ps(hk,hk);                                            hk_2 = _mm_mul_ps(hk_2,hk_2);
  tmp1 = _mm_mul_ps(ak,ak);                                          tmp1_2 = _mm_mul_ps(ak_2,ak_2);
  hk = _mm_sub_ps(hk, tmp1);                                         hk_2 = _mm_sub_ps(hk_2, tmp1_2);
  neg = _mm_cmple_ps(hk, zero);                                      neg_2 = _mm_cmple_ps(hk_2, zero);
  hk = _vec_sqrt_ps(hk);                                              hk_2 = _vec_sqrt_ps(hk_2);
  
  xk1 = _mm_mul_ps(ak, x);                                           xk1_2 = _mm_mul_ps(ak_2, x_2);
  tmp1 = _mm_mul_ps(hk, y);                                           tmp1_2 = _mm_mul_ps(hk_2, y_2);
  xk2 = _mm_sub_ps(xk1, tmp1);                                       xk2_2 = _mm_sub_ps(xk1_2, tmp1_2);
  xk1 = _mm_add_ps(xk1, tmp1);                                       xk1_2 = _mm_add_ps(xk1_2, tmp1_2);
  xk1 = _mm_mul_ps(xk1, D_inv);                                      xk1_2 = _mm_mul_ps(xk1_2, D_inv_2);
  xk2 = _mm_mul_ps(xk2, D_inv);                                      xk2_2 = _mm_mul_ps(xk2_2, D_inv_2);
  
  yk1 = _mm_mul_ps(ak, y);                                           yk1_2 = _mm_mul_ps(ak_2, y_2);
  tmp1 = _mm_mul_ps(hk,x);                                           tmp1_2 = _mm_mul_ps(hk_2, x_2);
  yk2 = _mm_add_ps(yk1, tmp1);                                       yk2_2 = _mm_add_ps(yk1_2, tmp1_2);
  yk1 = _mm_sub_ps(yk1, tmp1);                                       yk1_2 = _mm_sub_ps(yk1_2, tmp1_2);
  yk1 = _mm_mul_ps(yk1, D_inv);                                      yk1_2 = _mm_mul_ps(yk1_2, D_inv_2);
  yk2 = _mm_mul_ps(yk2, D_inv);                                      yk2_2 = _mm_mul_ps(yk2_2, D_inv_2);
  
  xk = _mm_and_ps(correct_helicity, xk1);                            xk_2 = _mm_and_ps(correct_helicity_2, xk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, xk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, xk2_2);
  xk = _mm_xor_ps(xk, tmp1);                                         xk_2 = _mm_xor_ps(xk_2, tmp1_2);
  yk = _mm_and_ps(correct_helicity, yk1);                            yk_2 = _mm_and_ps(correct_helicity_2, yk1_2);
  tmp1 = _mm_andnot_ps(correct_helicity, yk2);                       tmp1_2 = _mm_andnot_ps(correct_helicity_2, yk2_2);
  yk = _mm_xor_ps(yk, tmp1);                                         yk_2 = _mm_xor_ps(yk_2, tmp1_2);
  
  __m128 phi_2 = _vec_atan2_ps(yk, xk);                              __m128 phi_2_2 = _vec_atan2_ps(yk_2, xk_2);
  //if phi < 0, phi += 2*pi
  tmp1 = _mm_cmplt_ps(phi_2, zero);                                  tmp1_2 = _mm_cmplt_ps(phi_2_2, zero);
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp1 = _mm_andnot_ps(tmp1, zero);                                  tmp1_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp1 = _mm_xor_ps(tmp1, tmp2);                                     tmp1_2 = _mm_xor_ps(tmp1_2, tmp2_2);
  phi_2 = _mm_add_ps(phi_2, tmp1);                                   phi_2_2 = _mm_add_ps(phi_2_2, tmp1_2);
  // if neg==true, phi = hit_phi
  tmp1 = _mm_and_ps(neg, hit_phi);                                   tmp1_2 = _mm_and_ps(neg_2, hit_phi_2);
  phi_2 = _mm_andnot_ps(neg, phi_2);                                 phi_2_2 = _mm_andnot_ps(neg_2, phi_2_2);
  phi_2 = _mm_xor_ps(tmp1, phi_2);                                   phi_2_2 = _mm_xor_ps(tmp1_2, phi_2_2);
  phi_4_out = phi_2;                                                 phi_4_out_2 = phi_2_2;
  
  
  ////////////////////////////////////////////////////////////////
  
  // check if phi overlaps the 0,2pi jump
  tmp1 = _mm_cmplt_ps(phi_1, pi_over_two);                           tmp1_2 = _mm_cmplt_ps(phi_1_2, pi_over_two);
  tmp2 = _mm_cmplt_ps(phi_2, pi_over_two);                           tmp2_2 = _mm_cmplt_ps(phi_2_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);                                      tmp1_2 = _mm_or_ps(tmp1_2, tmp2_2);
  tmp2 = _mm_cmplt_ps(phi_3, pi_over_two);                           tmp2_2 = _mm_cmplt_ps(phi_3_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);                                      tmp1_2 = _mm_or_ps(tmp1_2, tmp2_2);
  tmp2 = _mm_cmplt_ps(phi_4, pi_over_two);                           tmp2_2 = _mm_cmplt_ps(phi_4_2, pi_over_two);
  tmp1 = _mm_or_ps(tmp1, tmp2);                                      tmp1_2 = _mm_or_ps(tmp1_2, tmp2_2);
  
  tmp2 = _mm_cmpgt_ps(phi_1, three_pi_over_two);                     tmp2_2 = _mm_cmpgt_ps(phi_1_2, three_pi_over_two);
  __m128 tmp3 = _mm_cmpgt_ps(phi_2, three_pi_over_two);              __m128 tmp3_2 = _mm_cmpgt_ps(phi_2_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);                                      tmp2_2 = _mm_or_ps(tmp2_2, tmp3_2);
  tmp3 = _mm_cmpgt_ps(phi_3, three_pi_over_two);                     tmp3_2 = _mm_cmpgt_ps(phi_3_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);                                      tmp2_2 = _mm_or_ps(tmp2_2, tmp3_2);
  tmp3 = _mm_cmpgt_ps(phi_4, three_pi_over_two);                     tmp3_2 = _mm_cmpgt_ps(phi_4_2, three_pi_over_two);
  tmp2 = _mm_or_ps(tmp2, tmp3);                                      tmp2_2 = _mm_or_ps(tmp2_2, tmp3_2);
  
  tmp1 = _mm_and_ps(tmp1, tmp2);                                     tmp1_2 = _mm_and_ps(tmp1_2, tmp2_2);
  
  // tmp1 is now all ones if phi overlaps the jump, all zeros otherwise
  // if tmp1 is true, then subtract 2*pi from all of the phi values > 3*pi/2
  tmp2 = _mm_and_ps(tmp1, twopi);                                    tmp2_2 = _mm_and_ps(tmp1_2, twopi);
  tmp3 = _mm_andnot_ps(tmp1, zero);                                  tmp3_2 = _mm_andnot_ps(tmp1_2, zero);
  tmp2 = _mm_xor_ps(tmp2, tmp3);                                     tmp2_2 = _mm_xor_ps(tmp2_2, tmp3_2);
  
  __m128 tmp4 = _mm_cmpgt_ps(phi_1, three_pi_over_two);              __m128 tmp4_2 = _mm_cmpgt_ps(phi_1_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  __m128 tmp5 = _mm_andnot_ps(tmp4, zero);                           __m128 tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_1 = _mm_sub_ps(phi_1, tmp3);                                   phi_1_2 = _mm_sub_ps(phi_1_2, tmp3_2);
  
  tmp4 = _mm_cmpgt_ps(phi_2, three_pi_over_two);                     tmp4_2 = _mm_cmpgt_ps(phi_2_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  tmp5 = _mm_andnot_ps(tmp4, zero);                                  tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_2 = _mm_sub_ps(phi_2, tmp3);                                   phi_2_2 = _mm_sub_ps(phi_2_2, tmp3_2);
  
  tmp4 = _mm_cmpgt_ps(phi_3, three_pi_over_two);                     tmp4_2 = _mm_cmpgt_ps(phi_3_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  tmp5 = _mm_andnot_ps(tmp4, zero);                                  tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_3 = _mm_sub_ps(phi_3, tmp3);                                   phi_3_2 = _mm_sub_ps(phi_3_2, tmp3_2);
  
  tmp4 = _mm_cmpgt_ps(phi_4, three_pi_over_two);                     tmp4_2 = _mm_cmpgt_ps(phi_4_2, three_pi_over_two);
  tmp3 = _mm_and_ps(tmp4, tmp2);                                     tmp3_2 = _mm_and_ps(tmp4_2, tmp2_2);
  tmp5 = _mm_andnot_ps(tmp4, zero);                                  tmp5_2 = _mm_andnot_ps(tmp4_2, zero);
  tmp3 = _mm_xor_ps(tmp3, tmp5);                                     tmp3_2 = _mm_xor_ps(tmp3_2, tmp5_2);
  phi_4 = _mm_sub_ps(phi_4, tmp3);                                   phi_4_2 = _mm_sub_ps(phi_4_2, tmp3_2);
  
  
  // find the minimum phi
  __m128 phi_min = phi_1;                                            __m128 phi_min_2 = phi_1_2;
  tmp2 = _mm_cmplt_ps(phi_2, phi_min);                               tmp2_2 = _mm_cmplt_ps(phi_2_2, phi_min_2);
  tmp3 = _mm_and_ps(tmp2, phi_2);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_2_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);                            phi_min_2 = _mm_andnot_ps(tmp2_2, phi_min_2);
  phi_min = _mm_xor_ps(phi_min, tmp3);                               phi_min_2 = _mm_xor_ps(phi_min_2, tmp3_2);
  tmp2 = _mm_cmplt_ps(phi_3, phi_min);                               tmp2_2 = _mm_cmplt_ps(phi_3_2, phi_min_2);
  tmp3 = _mm_and_ps(tmp2, phi_3);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_3_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);                            phi_min_2 = _mm_andnot_ps(tmp2_2, phi_min_2);
  phi_min = _mm_xor_ps(phi_min, tmp3);                               phi_min_2 = _mm_xor_ps(phi_min_2, tmp3_2);
  tmp2 = _mm_cmplt_ps(phi_4, phi_min);                               tmp2_2 = _mm_cmplt_ps(phi_4_2, phi_min_2);
  tmp3 = _mm_and_ps(tmp2, phi_4);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_4_2);
  phi_min = _mm_andnot_ps(tmp2, phi_min);                            phi_min_2 = _mm_andnot_ps(tmp2_2, phi_min_2);
  phi_min = _mm_xor_ps(phi_min, tmp3);                               phi_min_2 = _mm_xor_ps(phi_min_2, tmp3_2);
  
  // find the maximum phi
  __m128 phi_max = phi_1;                                            __m128 phi_max_2 = phi_1_2;
  tmp2 = _mm_cmpgt_ps(phi_2, phi_max);                               tmp2_2 = _mm_cmpgt_ps(phi_2_2, phi_max_2);
  tmp3 = _mm_and_ps(tmp2, phi_2);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_2_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);                            phi_max_2 = _mm_andnot_ps(tmp2_2, phi_max_2);
  phi_max = _mm_xor_ps(phi_max, tmp3);                               phi_max_2 = _mm_xor_ps(phi_max_2, tmp3_2);
  tmp2 = _mm_cmpgt_ps(phi_3, phi_max);                               tmp2_2 = _mm_cmpgt_ps(phi_3_2, phi_max_2);
  tmp3 = _mm_and_ps(tmp2, phi_3);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_3_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);                            phi_max_2 = _mm_andnot_ps(tmp2_2, phi_max_2);
  phi_max = _mm_xor_ps(phi_max, tmp3);                               phi_max_2 = _mm_xor_ps(phi_max_2, tmp3_2);
  tmp2 = _mm_cmpgt_ps(phi_4, phi_max);                               tmp2_2 = _mm_cmpgt_ps(phi_4_2, phi_max_2);
  tmp3 = _mm_and_ps(tmp2, phi_4);                                    tmp3_2 = _mm_and_ps(tmp2_2, phi_4_2);
  phi_max = _mm_andnot_ps(tmp2, phi_max);                            phi_max_2 = _mm_andnot_ps(tmp2_2, phi_max_2);
  phi_max = _mm_xor_ps(phi_max, tmp3);                               phi_max_2 = _mm_xor_ps(phi_max_2, tmp3_2);
  
  
  _mm_store_ps(min_phi, phi_min);                                    _mm_store_ps(min_phi_2, phi_min_2);
  _mm_store_ps(max_phi, phi_max);                                    _mm_store_ps(max_phi_2, phi_max_2);
}


