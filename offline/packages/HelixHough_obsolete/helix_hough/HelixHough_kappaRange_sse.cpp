#include "vector_math_inline.h"
#include "HelixHough.h"
#include <cmath>
#include <iostream>

using namespace std;


static const unsigned int allones[4] __attribute__((aligned(16))) = {0xffffffff,0xffffffff,0xffffffff,0xffffffff};
static const unsigned int allzeroes[4] __attribute__((aligned(16))) = {0,0,0,0};
static const __m128 one_o_2 = {0.5,0.5,0.5,0.5};


static inline void __attribute__((always_inline)) _vec_select_ux_uy(__m128 v_phi, __m128& ux, __m128& uy, __m128& tanphi, __m128& phisel)
{
  //extract the sign of phi and set pos_phi = |phi|
  __m128i vec_sgn = _mm_load_si128((__m128i*)sign_int);
  __m128 phi_sign = _mm_and_ps((__m128)vec_sgn, v_phi);
  __m128 pos_phi = _mm_xor_ps(phi_sign, v_phi);
  
  //|phi - pi|
  __m128 v_phi_m_pi = _mm_sub_ps(v_phi, pi);
  phi_sign = _mm_and_ps((__m128)vec_sgn, v_phi_m_pi);
  __m128 pos_phi_m_pi = _mm_xor_ps(phi_sign, v_phi_m_pi);
  
  //|twopi - phi|
  __m128 v_twopi_m_phi = _mm_sub_ps(twopi, v_phi);
  phi_sign = _mm_and_ps((__m128)vec_sgn, v_twopi_m_phi);
  __m128 pos_twopi_m_phi = _mm_xor_ps(phi_sign, v_twopi_m_phi);
  
  phisel = _mm_cmplt_ps(pos_phi, pi_over_four);
  __m128 tmp2 = _mm_cmplt_ps(pos_phi_m_pi, pi_over_four);
  phisel = _mm_or_ps(phisel, tmp2);
  tmp2 = _mm_cmplt_ps(pos_twopi_m_phi, pi_over_four);
  phisel = _mm_or_ps(phisel, tmp2);
  
  __m128 sinval = {0.,0.,0.,0.};
  __m128 cosval = {0.,0.,0.,0.};
  _vec_sin_cos_ps(v_phi, sinval, cosval);
  
  __m128 num = _mm_and_ps(phisel, sinval);
  tmp2 = _mm_andnot_ps(phisel, cosval);
  num = _mm_xor_ps(num, tmp2);
  __m128 denom = _mm_and_ps(phisel, cosval);
  tmp2 = _mm_andnot_ps(phisel, sinval);
  denom = _mm_xor_ps(denom, tmp2);
  tanphi = _mm_div_ps(num, denom);
  __m128 v_phi_p_pi02 = _mm_add_ps(v_phi, pi_over_two);
  __m128 phi_compare = _mm_and_ps(phisel, v_phi_p_pi02);
  tmp2 = _mm_andnot_ps(phisel, v_phi);
  phi_compare = _mm_xor_ps(phi_compare, tmp2);
  __m128 tmp3 = _mm_mul_ps(tanphi, tanphi);
  tmp3 = _mm_add_ps(one, tmp3);
  tmp3 = _mm_div_ps(one, tmp3);
  tmp3 = _mm_sqrt_ps(tmp3);
  tmp2 = _mm_cmpgt_ps(phi_compare, pi);
  __m128 tmp4 = _mm_cmplt_ps(phi_compare, twopi);
  tmp2 = _mm_and_ps(tmp2, tmp4);
  tmp2 = _mm_and_ps((__m128)vec_sgn, tmp2);
  tmp3 = _mm_xor_ps(tmp3, tmp2);
  tmp4 = _mm_mul_ps(tanphi, tmp3);
  ux = _mm_and_ps(phisel, tmp3);
  tmp2 = _mm_andnot_ps(phisel, tmp4);
  ux = _mm_xor_ps(ux, tmp2);
  uy = _mm_and_ps(phisel, tmp4);
  tmp2 = _mm_andnot_ps(phisel, tmp3);
  uy = _mm_xor_ps(uy, tmp2);
}
      
      
static inline void __attribute__((always_inline)) _vec_invRRangeOnePhi(__m128 x, __m128 y, __m128 v_phi, __m128 v_mind, __m128 v_maxd, __m128& v_mink, __m128& v_maxk, __m128& reflection_min, __m128& reflection_max, __m128& incorrect_quadrant, __m128& correct_quadrant)
{
  __m128 dx = {0.,0.,0.,0.};
  __m128 dy = {0.,0.,0.,0.};
  __m128 ux = {0.,0.,0.,0.};
  __m128 uy = {0.,0.,0.,0.};
  __m128 tanphi = {0.,0.,0.,0.};
  __m128 phisel = {0.,0.,0.,0.};
  __m128 t = {0.,0.,0.,0.};
  __m128 tP = {0.,0.,0.,0.};
  __m128 mx = {0.,0.,0.,0.};
  __m128 my = {0.,0.,0.,0.};
  __m128 Ix = {0.,0.,0.,0.};
  __m128 Iy = {0.,0.,0.,0.};
  __m128 r1 = {0.,0.,0.,0.};
  __m128 r2 = {0.,0.,0.,0.};
  
  __m128i vec_ones = _mm_load_si128((__m128i*)allones);
  __m128i vec_sgn = _mm_load_si128((__m128i*)sign_int);
  
  _vec_select_ux_uy(v_phi, ux, uy, tanphi, phisel);
  
  tP = _mm_mul_ps(ux, x);
  __m128 tmp1 = _mm_mul_ps(uy, y);
  tP = _mm_add_ps(tP, tmp1);
  dx = _mm_mul_ps(ux, v_mind);
  dy = _mm_mul_ps(uy, v_mind);
  
  tmp1 = _mm_cmplt_ps(v_mind, tP);
  reflection_min = _mm_or_ps(reflection_min, tmp1);
  tmp1 = _mm_xor_ps(tmp1, (__m128)vec_ones);
  reflection_max = _mm_or_ps(reflection_max, tmp1);
  mx = _mm_add_ps(dx, x);
  mx = _mm_mul_ps(mx, one_o_2);
  my = _mm_add_ps(dy, y);
  my = _mm_mul_ps(my, one_o_2);
  
  __m128 tmp2 = _mm_mul_ps(mx, tanphi);
  tmp2 = _mm_sub_ps(my, tmp2);
  tmp1 = _mm_sub_ps(y, dy);
  __m128 tmp3 = _mm_mul_ps(tanphi, tmp1);
  tmp1 = _mm_sub_ps(x, dx);
  tmp1 = _mm_add_ps(tmp1, tmp3);
  tmp1 = _mm_div_ps(tmp2, tmp1);
  tmp2 = _mm_mul_ps(tanphi, my);
  tmp2 = _mm_sub_ps(tmp2, mx);
  tmp3 = _mm_sub_ps(x, dx);
  tmp3 = _mm_mul_ps(tmp3, tanphi);
  __m128 tmp4 = _mm_sub_ps(y, dy);
  tmp3 = _mm_add_ps(tmp3, tmp4);
  tmp2 = _mm_div_ps(tmp2, tmp3);
  t = _mm_and_ps(phisel, tmp1);
  tmp3 = _mm_andnot_ps(phisel, tmp2);
  t = _mm_xor_ps(t, tmp3);
  
  tmp1 = _mm_sub_ps(y,dy);
  tmp1 = _mm_mul_ps(tmp1, t);
  Ix = _mm_add_ps(mx, tmp1);
  tmp2 = _mm_sub_ps(dx,x);
  tmp2 = _mm_mul_ps(tmp2, t);
  Iy = _mm_add_ps(my, tmp2);
  
  //extract sign of Ix
  __m128 Ix_sign = _mm_and_ps((__m128)vec_sgn, Ix);
  //sign of ux
  __m128 ux_sign = _mm_and_ps((__m128)vec_sgn, ux);
  //compare signs of Ix and ux
  tmp1 = _mm_cmpeq_ps(Ix_sign, ux_sign);
  //now for Iy,uy
  __m128 Iy_sign = _mm_and_ps((__m128)vec_sgn, Iy);
  __m128 uy_sign = _mm_and_ps((__m128)vec_sgn, uy);
  tmp2 = _mm_cmpeq_ps(Iy_sign, uy_sign);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  correct_quadrant = _mm_or_ps(correct_quadrant, tmp3);
  tmp3 = _mm_xor_ps(tmp3, (__m128)vec_ones);
  incorrect_quadrant = _mm_or_ps(incorrect_quadrant, tmp3);
  
  tmp1 = _mm_sub_ps(Ix, dx);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  tmp2 = _mm_sub_ps(Iy, dy);
  tmp2 = _mm_mul_ps(tmp2, tmp2);
  r1 = _mm_add_ps(tmp1, tmp2);
  r1 = _mm_sqrt_ps(r1);
  
  //now repeat for maxd
  dx = _mm_mul_ps(ux, v_maxd);
  dy = _mm_mul_ps(uy, v_maxd);
  
  tmp1 = _mm_cmplt_ps(v_maxd, tP);
  reflection_min = _mm_or_ps(reflection_min, tmp1);
  tmp1 = _mm_xor_ps(tmp1, (__m128)vec_ones);
  reflection_max = _mm_or_ps(reflection_max, tmp1);
  mx = _mm_add_ps(dx, x);
  mx = _mm_mul_ps(mx, one_o_2);
  my = _mm_add_ps(dy, y);
  my = _mm_mul_ps(my, one_o_2);
  
  tmp2 = _mm_mul_ps(mx, tanphi);
  tmp2 = _mm_sub_ps(my, tmp2);
  tmp1 = _mm_sub_ps(y, dy);
  tmp3 = _mm_mul_ps(tanphi, tmp1);
  tmp1 = _mm_sub_ps(x, dx);
  tmp1 = _mm_add_ps(tmp1, tmp3);
  tmp1 = _mm_div_ps(tmp2, tmp1);
  tmp2 = _mm_mul_ps(tanphi, my);
  tmp2 = _mm_sub_ps(tmp2, mx);
  tmp3 = _mm_sub_ps(x, dx);
  tmp3 = _mm_mul_ps(tmp3, tanphi);
  tmp4 = _mm_sub_ps(y, dy);
  tmp3 = _mm_add_ps(tmp3, tmp4);
  tmp2 = _mm_div_ps(tmp2, tmp3);
  t = _mm_and_ps(phisel, tmp1);
  tmp3 = _mm_andnot_ps(phisel, tmp2);
  t = _mm_xor_ps(t, tmp3);
  
  tmp1 = _mm_sub_ps(y,dy);
  tmp1 = _mm_mul_ps(tmp1, t);
  Ix = _mm_add_ps(mx, tmp1);
  tmp2 = _mm_sub_ps(dx,x);
  tmp2 = _mm_mul_ps(tmp2, t);
  Iy = _mm_add_ps(my, tmp2);
  
  //extract sign of Ix
  Ix_sign = _mm_and_ps((__m128)vec_sgn, Ix);
  //sign of ux
  ux_sign = _mm_and_ps((__m128)vec_sgn, ux);
  //compare signs of Ix and ux
  tmp1 = _mm_cmpeq_ps(Ix_sign, ux_sign);
  //now for Iy,uy
  Iy_sign = _mm_and_ps((__m128)vec_sgn, Iy);
  uy_sign = _mm_and_ps((__m128)vec_sgn, uy);
  tmp2 = _mm_cmpeq_ps(Iy_sign, uy_sign);
  tmp3 = _mm_and_ps(tmp1, tmp2);
  correct_quadrant = _mm_or_ps(correct_quadrant, tmp3);
  tmp3 = _mm_xor_ps(tmp3, (__m128)vec_ones);
  incorrect_quadrant = _mm_or_ps(incorrect_quadrant, tmp3);
  
  tmp1 = _mm_sub_ps(Ix, dx);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  tmp2 = _mm_sub_ps(Iy, dy);
  tmp2 = _mm_mul_ps(tmp2, tmp2);
  r2 = _mm_add_ps(tmp1, tmp2);
  r2 = _mm_sqrt_ps(r2);
  
  __m128 invr1 = _mm_div_ps(one, r1);
  __m128 invr2 = _mm_div_ps(one, r2);
  tmp1 = _mm_cmplt_ps(r1, r2);
  v_mink = _mm_and_ps(tmp1, invr2);
  tmp2 = _mm_andnot_ps(tmp1, invr1);
  v_mink = _mm_xor_ps(v_mink, tmp2);
  v_maxk = _mm_and_ps(tmp1, invr1);
  tmp2 = _mm_andnot_ps(tmp1, invr2);
  v_maxk = _mm_xor_ps(v_maxk, tmp2);
}
      
      
void HelixHough::kappaRange_sse(const SimpleHit3D& hit, float* minphi, float* maxphi, float* mind, float* maxd, float* min_invR, float* max_invR)
{
  __m128 v_minphi = _mm_load_ps(minphi);
  __m128 v_maxphi = _mm_load_ps(maxphi);
  __m128 v_mind = _mm_load_ps(mind);
  __m128 v_maxd = _mm_load_ps(maxd);
  __m128 v_mink = _mm_load_ps(min_invR);
  __m128 v_maxk = _mm_load_ps(max_invR);
  
  float x_arr[4] __attribute__((aligned(16))) = {hit.x,hit.x,hit.x,hit.x};
  float y_arr[4] __attribute__((aligned(16))) = {hit.y,hit.y,hit.y,hit.y};
  
  //   __m128 x = _mm_load1_ps(&(hit.x));
  //   __m128 y = _mm_load1_ps(&(hit.y));
  __m128 x = _mm_load_ps(x_arr);
  __m128 y = _mm_load_ps(y_arr);
  
  __m128i vec_zeroes = _mm_load_si128((__m128i*)allzeroes);
  
  __m128 reflection_min = (__m128)(vec_zeroes);
  __m128 reflection_max = (__m128)(vec_zeroes);
  __m128 correct_quadrant = (__m128)(vec_zeroes);
  __m128 incorrect_quadrant = (__m128)(vec_zeroes);
  
  __m128 v_mink1 = {0.,0.,0.,0.};
  __m128 v_maxk1 = {0.,0.,0.,0.};
  _vec_invRRangeOnePhi(x, y, v_minphi, v_mind, v_maxd, v_mink1, v_maxk1, reflection_min, reflection_max, incorrect_quadrant, correct_quadrant);
  
  __m128 v_mink2 = {0.,0.,0.,0.};
  __m128 v_maxk2 = {0.,0.,0.,0.};
  _vec_invRRangeOnePhi(x, y, v_maxphi, v_mind, v_maxd, v_mink2, v_maxk2, reflection_min, reflection_max, incorrect_quadrant, correct_quadrant);
  
  __m128 tmp1 = _mm_cmplt_ps(v_mink1, v_mink2);
  v_mink = _mm_and_ps(tmp1, v_mink1);
  __m128 tmp2 = _mm_andnot_ps(tmp1, v_mink2);
  v_mink = _mm_xor_ps(v_mink, tmp2);
  
  tmp1 = _mm_cmpgt_ps(v_maxk1, v_maxk2);
  v_maxk = _mm_and_ps(tmp1, v_maxk1);
  tmp2 = _mm_andnot_ps(tmp1, v_maxk2);
  v_maxk = _mm_xor_ps(v_maxk, tmp2);
  
  tmp1 = _mm_and_ps(reflection_min, reflection_max);
  v_mink = _mm_andnot_ps(tmp1, v_mink);
  tmp2 = _mm_and_ps(tmp1, zero);
  v_mink = _mm_xor_ps(v_mink, tmp2);
  
  tmp1 = _mm_and_ps(correct_quadrant, incorrect_quadrant);
  v_mink = _mm_andnot_ps(tmp1, v_mink);
  tmp2 = _mm_and_ps(tmp1, zero);
  v_mink = _mm_xor_ps(v_mink, tmp2);
  
  _mm_store_ps(min_invR, v_mink);
  _mm_store_ps(max_invR, v_maxk);
}


