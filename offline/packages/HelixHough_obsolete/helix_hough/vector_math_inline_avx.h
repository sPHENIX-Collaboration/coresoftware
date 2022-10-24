#ifdef AVXHOUGH

#include <immintrin.h>


static const __m256 pi_over_two_256 = {0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f};
static const __m256 pi_over_four_256 = {0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f};
static const __m256 three_256 = {3., 3., 3., 3., 3., 3., 3., 3.};
static const __m256 one_o_2_256 = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
static const __m256 zero_256 = {0., 0., 0., 0., 0., 0., 0., 0.};
static const __m256 one_256 = {1., 1., 1., 1., 1., 1., 1., 1.};
static const __m256 negone_256 = {-1., -1., -1., -1., -1., -1., -1., -1.};
static const __m256 atanC0_256 = {0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f};
static const __m256 atanC1_256 = {-0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f};
static const __m256 atanC2_256 = {0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f};
static const __m256 atanC3_256 = {-0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f};
static const __m256 sqrt2_minus_1_256 = {0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f};
static const __m256 sqrt2_plus_1_256 = {0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f};

static const unsigned int sign_int_256[8] __attribute__((aligned(32))) = {0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000};


static inline __m256  __attribute__((always_inline)) _mm256_cmplt_ps(__m256 a, __m256 b){ return _mm256_cmp_ps(a, b, _CMP_LT_OQ); }
static inline __m256  __attribute__((always_inline)) _mm256_cmpgt_ps(__m256 a, __m256 b){ return _mm256_cmp_ps(a, b, _CMP_GT_OQ); }
static inline __m256  __attribute__((always_inline)) _mm256_cmpeq_ps(__m256 a, __m256 b){ return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }


static inline __m256  __attribute__((always_inline)) _mm256_load1_ps(float x){ return _mm256_set_ps(x, x, x, x, x, x, x, x); }


inline __m256 __attribute__((always_inline)) _vec256_rsqrt_ps(__m256 x)
{
  __m256 x0 = _mm256_rsqrt_ps(x);
  return _mm256_mul_ps( one_o_2_256, _mm256_mul_ps( x0, _mm256_sub_ps( three_256, _mm256_mul_ps( x0, _mm256_mul_ps(x0, x) ) ) ) );
}


inline __m256 __attribute__((always_inline)) _vec256_sqrt_ps(__m256 x)
{
  __m256 x0 = _mm256_rsqrt_ps(x);
  return _mm256_mul_ps(_mm256_mul_ps( one_o_2_256, _mm256_mul_ps( x0, _mm256_sub_ps( three_256, _mm256_mul_ps( x0, _mm256_mul_ps(x0, x) ) ) ) ), x);
}


inline __m256 __attribute__((always_inline)) _vec256_rec_ps(__m256 x)
{
  __m256 a = _vec256_rsqrt_ps(x);
  return _mm256_mul_ps(a, a);
}


inline void __attribute__((always_inline)) _vec256_fabs_ps(__m256 &v)
{
  __m256i vec_sgn = _mm256_load_si256((__m256i*)sign_int_256);
  v = _mm256_andnot_ps((__m256)vec_sgn, v);
}




inline __m256 __attribute__((always_inline)) _vec256_atan_ps(__m256 x)
{
  __m256i vec_sgn;
  __m256 x_sign;
  __m256 gr1;
  __m256 gr2;
  __m256 z1;
  __m256 z2;
  __m256 x2;
  
  //extract the sign of x and make x positive
  vec_sgn = _mm256_load_si256((__m256i*)sign_int_256);
  x_sign = _mm256_and_ps((__m256)vec_sgn, x);
  x = _mm256_xor_ps(x_sign, x);
  
  
  //if x > sqrt(2)+1, then set x = -1/x
  //else if x > sqrt(2)-1, then set x = (x-1)/(x+1)
  gr1 = _mm256_cmpgt_ps(x, sqrt2_plus_1_256);
  gr2 = _mm256_cmpgt_ps(x, sqrt2_minus_1_256);
  z1 = _mm256_div_ps(negone_256, x);
  z1 = _mm256_and_ps(z1, gr1);
  z2 = _mm256_add_ps(x, negone_256);
  x2 = _mm256_sub_ps(x, negone_256);
  z2 = _mm256_div_ps(z2, x2);
  x = _mm256_andnot_ps(gr2, x);
  gr2 = _mm256_xor_ps(gr2, gr1);
  z2 = _mm256_and_ps(z2, gr2);
  x = _mm256_or_ps(x, z1);
  x = _mm256_or_ps(x, z2);
  
  
  //compute z1 = atan(x) from a Chebyshev polynomial (in monomial form) using Horner's scheme
  x2 = _mm256_mul_ps(x, x);
  z1 = _mm256_mul_ps(x2, atanC3_256);
  z1 = _mm256_add_ps(z1, atanC2_256);
  z2 = _mm256_mul_ps(x2, z1);
  z2 = _mm256_add_ps(z2, atanC1_256);
  z1 = _mm256_mul_ps(x2, z2);
  z1 = _mm256_add_ps(z1, atanC0_256);
  z1 = _mm256_mul_ps(z1, x);
  
  //add either pi/4 or pi/2 to z1 and set result to x2, depending on the initial value of x
  x = _mm256_and_ps(gr1, pi_over_two_256);
  x2 = _mm256_and_ps(gr2, pi_over_four_256);
  x2 = _mm256_or_ps(x2, x);
  x2 = _mm256_add_ps(x2, z1);
  
  //recover the original sign of x
  x2 = _mm256_xor_ps(x2, x_sign);
  
  return x2;
}


#endif
