#include <xmmintrin.h>
#include <emmintrin.h>

static const __m128 twoTo23 = (const __m128) { 0x1.0p23f, 0x1.0p23f, 0x1.0p23f, 0x1.0p23f };

inline __m128 __attribute__((always_inline)) _vec_floor_ps(__m128 v)
{
  // b = fabs(v)
  __m128 b = (__m128) _mm_srli_epi32( _mm_slli_epi32( (__m128i) v, 1 ), 1 );
  // The essence of the floor routine
  __m128 d = _mm_sub_ps( _mm_add_ps( _mm_add_ps( _mm_sub_ps( v, twoTo23 ), twoTo23 ), twoTo23 ), twoTo23 );
  // ?1 if v >= 2**23
  __m128 largeMaskE = _mm_cmpgt_ps( b, twoTo23 );
  // Check for possible off by one error
  __m128 g = _mm_cmplt_ps( v, d );
  // Convert positive check result to -1.0, negative to 0.0
  __m128 h = _mm_cvtepi32_ps( (__m128i) g );
  // Add in the error if there is one
  __m128 t = _mm_add_ps( d, h );
  //Select between output result and input value based on v >= 2**23
  v = _mm_and_ps( v, largeMaskE );
  t = _mm_andnot_ps( largeMaskE, t );
  return _mm_or_ps( t, v );
}


static const __m128 pi = {0x3.243f6a8885a308d4p0f, 0x3.243f6a8885a308d4p0f, 0x3.243f6a8885a308d4p0f, 0x3.243f6a8885a308d4p0f};
static const __m128 twopi = {0x6.487ed5110b4611a8p0f, 0x6.487ed5110b4611a8p0f, 0x6.487ed5110b4611a8p0f, 0x6.487ed5110b4611a8p0f};
static const __m128 pi_over_two = {0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f, 0x1.921fb54442d1846ap0f};
static const __m128 pi_over_four = {0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f, 0xc.90fdaa22168c2350p-4f};

static const __m128 sqrt2_minus_1 = {0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f, 0x6.a09e667f3bcc9080p-4f};
static const __m128 sqrt2_plus_1 = {0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f, 0x2.6a09e667f3bcc9080p0f};

static const __m128 zero = {0., 0., 0., 0.};
static const __m128 one = {1., 1., 1., 1.};
static const __m128 three = {3., 3., 3., 3.};
static const __m128 negone = {-1., -1., -1., -1.};
static const __m128 one_o_2 = {0.5,0.5,0.5,0.5};

static const unsigned int sign_int[4] __attribute__((aligned(16))) = {0x80000000,0x80000000,0x80000000,0x80000000};

static const __m128 atanC0 = {0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f, 0xf.fffb771eba87d370p-4f};
static const __m128 atanC1 = {-0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f, -0x5.542eef19db937268p-4f};
static const __m128 atanC2 = {0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f, 0x3.1dcf607e2808c0d4p-4f};
static const __m128 atanC3 = {-0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f, -0x1.ab85dd26f5264feep-4f};


inline __m128 __attribute__((always_inline)) _vec_rsqrt_ps(__m128 x)
{
  __m128 x0 = _mm_rsqrt_ps(x);
  return _mm_mul_ps( one_o_2, _mm_mul_ps( x0, _mm_sub_ps( three, _mm_mul_ps( x0, _mm_mul_ps(x0, x) ) ) ) );
}


inline __m128 __attribute__((always_inline)) _vec_sqrt_ps(__m128 x)
{
  __m128 x0 = _mm_rsqrt_ps(x);
  return _mm_mul_ps(_mm_mul_ps( one_o_2, _mm_mul_ps( x0, _mm_sub_ps( three, _mm_mul_ps( x0, _mm_mul_ps(x0, x) ) ) ) ), x);
}


inline __m128 __attribute__((always_inline)) _vec_rec_ps(__m128 x)
{
  __m128 a = _vec_rsqrt_ps(x);
  return _mm_mul_ps(a, a);
}


inline __m128 __attribute__((always_inline)) _vec_atan_ps(__m128 x)
{
  __m128i vec_sgn;
  __m128 x_sign;
  __m128 gr1;
  __m128 gr2;
  __m128 z1;
  __m128 z2;
  __m128 x2;
  
  //extract the sign of x and make x positive
  vec_sgn = _mm_load_si128((__m128i*)sign_int);
  x_sign = _mm_and_ps((__m128)vec_sgn, x);
  x = _mm_xor_ps(x_sign, x);
  
  
  //if x > sqrt(2)+1, then set x = -1/x
  //else if x > sqrt(2)-1, then set x = (x-1)/(x+1)
  gr1 = _mm_cmpgt_ps(x, sqrt2_plus_1);
  gr2 = _mm_cmpgt_ps(x, sqrt2_minus_1);
  z1 = _mm_div_ps(negone, x);
  z1 = _mm_and_ps(z1, gr1);
  z2 = _mm_add_ps(x, negone);
  x2 = _mm_sub_ps(x, negone);
  z2 = _mm_div_ps(z2, x2);
  x = _mm_andnot_ps(gr2, x);
  gr2 = _mm_xor_ps(gr2, gr1);
  z2 = _mm_and_ps(z2, gr2);
  x = _mm_or_ps(x, z1);
  x = _mm_or_ps(x, z2);
  
  
  //compute z1 = atan(x) from a Chebyshev polynomial (in monomial form) using Horner's scheme
  x2 = _mm_mul_ps(x, x);
  z1 = _mm_mul_ps(x2, atanC3);
  z1 = _mm_add_ps(z1, atanC2);
  z2 = _mm_mul_ps(x2, z1);
  z2 = _mm_add_ps(z2, atanC1);
  z1 = _mm_mul_ps(x2, z2);
  z1 = _mm_add_ps(z1, atanC0);
  z1 = _mm_mul_ps(z1, x);
  
  //add either pi/4 or pi/2 to z1 and set result to x2, depending on the initial value of x
  x = _mm_and_ps(gr1, pi_over_two);
  x2 = _mm_and_ps(gr2, pi_over_four);
  x2 = _mm_or_ps(x2, x);
  x2 = _mm_add_ps(x2, z1);
  
  //recover the original sign of x
  x2 = _mm_xor_ps(x2, x_sign);
  
  return x2;
}


inline __m128 __attribute__((always_inline)) _vec_atan2_ps(__m128 y, __m128 x)
{
  __m128 eq0 = _mm_cmpeq_ps(x, zero);
  __m128 ratio = _mm_div_ps(y, x);
  __m128 atanval = _vec_atan_ps(ratio);
  __m128 zero_pio2 = _mm_and_ps(eq0, pi_over_two);
  __m128i vec_sgn = _mm_load_si128((__m128i*)sign_int);
  __m128 y_sign = _mm_and_ps((__m128)vec_sgn, y);
  __m128 less0 = _mm_cmplt_ps(x, zero);
  zero_pio2 = _mm_xor_ps(y_sign, zero_pio2);
  __m128 zero_pi = _mm_and_ps(less0, pi);
  zero_pi = _mm_xor_ps(y_sign, zero_pi);
  atanval = _mm_andnot_ps(eq0, atanval);
  atanval = _mm_add_ps(zero_pio2, atanval);
  atanval = _mm_add_ps(zero_pi, atanval);
  
  return atanval;
}


inline void __attribute__((always_inline)) _vec_fabs_ps(__m128 &v)
{
  __m128i vec_sgn = _mm_load_si128((__m128i*)sign_int);
  v = _mm_andnot_ps((__m128)vec_sgn, v);
}

static const __m128 one_over_twopi = {0x2.8be60db9391054ap-4f, 0x2.8be60db9391054ap-4f, 0x2.8be60db9391054ap-4f, 0x2.8be60db9391054ap-4f};

static const __m128 fourth = {0x0.4p0f, 0x0.4p0f, 0x0.4p0f, 0x0.4p0f};
static const __m128 half = {0x0.8p0f, 0x0.8p0f, 0x0.8p0f, 0x0.8p0f};

static const __m128 sinC0 = {0x6.487c58e5205cd791p0f, 0x6.487c58e5205cd791p0f, 0x6.487c58e5205cd791p0f, 0x6.487c58e5205cd791p0f};
static const __m128 sinC1 = {-0x2.955c385f44a6765fp4f, -0x2.955c385f44a6765fp4f, -0x2.955c385f44a6765fp4f, -0x2.955c385f44a6765fp4f};
static const __m128 sinC2 = {0x5.145d3f35fa67830ep4f, 0x5.145d3f35fa67830ep4f, 0x5.145d3f35fa67830ep4f, 0x5.145d3f35fa67830ep4f};
static const __m128 sinC3 = {-0x4.65f4793b5cd9628fp4f, -0x4.65f4793b5cd9628fp4f, -0x4.65f4793b5cd9628fp4f, -0x4.65f4793b5cd9628fp4f};


inline void __attribute__((always_inline)) _vec_sin_cos_ps(__m128 x, __m128 &s, __m128 &c)
{
  //we will compute via a polynomial for sin(2*pi*x)
  //sin(x) = sin(2*pi*(x/2*pi)), so first we divide x by 2*pi
  x = _mm_mul_ps(x, one_over_twopi);
  //set x to the fractional part of x
  x = _mm_sub_ps(x, _vec_floor_ps(x));
  //subtract 1/2 from x if x>1/2
  __m128 mod_half = _mm_cmpgt_ps(x, half);
  __m128 t1 = _mm_andnot_ps(mod_half, zero);
  __m128 t2 = _mm_and_ps(mod_half, half);
  t1 = _mm_xor_ps(t1, t2);
  x = _mm_sub_ps(x, t1);
  //subtract 1/4 from x if x>1/4
  __m128 mod_fourth = _mm_cmpgt_ps(x, fourth);
  t1 = _mm_andnot_ps(mod_fourth, zero);
  t2 = _mm_and_ps(mod_fourth, fourth);
  t1 = _mm_xor_ps(t1, t2);
  x = _mm_sub_ps(x, t1);
  //z = 1/4 - x
  __m128 z = _mm_sub_ps(fourth, x);
  //compute t1 = sin(x) from a Chebyshev polynomial (in monomial form) using Horner's scheme
  //if we were using higher precision we would want to calculate sqrt(1 - [sin(x)]^2), but for this low precision we will just use the polynomial again
  //compute k1 = sin(z) from a Chebyshev polynomial (in monomial form) using Horner's scheme
  __m128 k1, k2, x2, z2;
  x2 = _mm_mul_ps(x, x);
  z2 = _mm_mul_ps(z, z);
  t1 = _mm_mul_ps(x2, sinC3);
  k1 = _mm_mul_ps(z2, sinC3);
  t1 = _mm_add_ps(t1, sinC2);
  k1 = _mm_add_ps(k1, sinC2);
  t2 = _mm_mul_ps(x2, t1);
  k2 = _mm_mul_ps(z2, k1);
  t2 = _mm_add_ps(t2, sinC1);
  k2 = _mm_add_ps(k2, sinC1);
  t1 = _mm_mul_ps(x2, t2);
  k1 = _mm_mul_ps(z2, k2);
  t1 = _mm_add_ps(t1, sinC0);
  k1 = _mm_add_ps(k1, sinC0);
  t1 = _mm_mul_ps(t1, x);
  k1 = _mm_mul_ps(k1, z);
  
  //set s and c to the appropriate magnitudes of sin and cosine
  s = _mm_andnot_ps(mod_fourth, t1);
  t2 = _mm_and_ps(mod_fourth, k1);
  s = _mm_xor_ps(s, t2);
  
  c = _mm_andnot_ps(mod_fourth, k1);
  k2 = _mm_and_ps(mod_fourth, t1);
  c = _mm_xor_ps(c, k2);
  
  //set the appropriate signs of s and c
  __m128i vec_sgn = _mm_load_si128((__m128i*)sign_int);
  __m128 x_sign = _mm_and_ps(mod_half, (__m128)vec_sgn);
  s = _mm_xor_ps(s, x_sign);
  
  __m128 mod_cos = _mm_xor_ps(mod_fourth, mod_half);
  x_sign = _mm_and_ps(mod_cos, (__m128)vec_sgn);
  c = _mm_xor_ps(c, x_sign);
}


static const __m128d pid = {0x3.243f6a8885a308d4p0, 0x3.243f6a8885a308d4p0};
static const __m128d twopid = {0x6.487ed5110b4611a8p0, 0x6.487ed5110b4611a8p0};
static const __m128d pi_over_twod = {0x1.921fb54442d1846ap0, 0x1.921fb54442d1846ap0};
static const __m128d pi_over_fourd = {0xc.90fdaa22168c2350p-4, 0xc.90fdaa22168c2350p-4};

static const unsigned int sign_int_d[4] __attribute__((aligned(16))) = {0x80000000,0x00000000,0x80000000,0x00000000};

static const __m128d atanD0 = {0xf.fffffffffffffffffab9f803d2af0fb0p-4, 0xf.fffffffffffffffffab9f803d2af0fb0p-4};
static const __m128d atanD1 = {-0x5.5555555555555540c86e6b5fd8e468b0p-4, -0x5.5555555555555540c86e6b5fd8e468b0p-4};
static const __m128d atanD2 = {0x3.3333333333331b3d286002f2c41b81e0p-4, 0x3.3333333333331b3d286002f2c41b81e0p-4};
static const __m128d atanD3 = {-0x2.4924924924851ef4ced41f651e628f40p-4, -0x2.4924924924851ef4ced41f651e628f40p-4};
static const __m128d atanD4 = {0x1.c71c71c7184d82416ebb498f58e3a070p-4, 0x1.c71c71c7184d82416ebb498f58e3a070p-4};
static const __m128d atanD5 = {-0x1.745d1744fcf2617d04fabeb866a4f6dep-4, -0x1.745d1744fcf2617d04fabeb866a4f6dep-4};
static const __m128d atanD6 = {0x1.3b13b11e1c3da8ad9c8d73ccfc6c3722p-4, 0x1.3b13b11e1c3da8ad9c8d73ccfc6c3722p-4};
static const __m128d atanD7 = {-0x1.11110e44194b2e57596c2fdb8bfd10a0p-4, -0x1.11110e44194b2e57596c2fdb8bfd10a0p-4};
static const __m128d atanD8 = {0xf.0f0be69f64251bcf8af8470c615a7c40p-8, 0xf.0f0be69f64251bcf8af8470c615a7c40p-8};
static const __m128d atanD9 = {-0xd.79192575eea6d23daa4764613eb39850p-8, -0xd.79192575eea6d23daa4764613eb39850p-8};
static const __m128d atanD10 = {0xc.2f1d52fbd8638e0fd3d27cdf0e6e13c0p-8, 0xc.2f1d52fbd8638e0fd3d27cdf0e6e13c0p-8};
static const __m128d atanD11 = {-0xb.1518a42d3671c2ee1bf46f36650357a0p-8, -0xb.1518a42d3671c2ee1bf46f36650357a0p-8};
static const __m128d atanD12 = {0x9.f9678bbe523adaab81aff27ebf0ec070p-8, 0x9.f9678bbe523adaab81aff27ebf0ec070p-8};
static const __m128d atanD13 = {-0x8.68d5692c1b536ea0b6afe3a59bdae0b0p-8, -0x8.68d5692c1b536ea0b6afe3a59bdae0b0p-8};
static const __m128d atanD14 = {0x5.c3234a24f6727d6cacaf7b5c9647c2b8p-8, 0x5.c3234a24f6727d6cacaf7b5c9647c2b8p-8};
static const __m128d atanD15 = {-0x2.46516323aab74d114a581091dd99ed84p-8, -0x2.46516323aab74d114a581091dd99ed84p-8};

static const __m128d sqrt2_minus_1d = {0x6.a09e667f3bcc9080p-4, 0x6.a09e667f3bcc9080p-4};
static const __m128d sqrt2_plus_1d = {0x2.6a09e667f3bcc9080p0, 0x2.6a09e667f3bcc9080p0};

static const __m128d zerod = {0., 0.};
static const __m128d oned = {1., 1.};
static const __m128d negoned = {-1., -1.};



inline __m128d __attribute__((always_inline)) _vec_atan_dfull_ps(__m128d x)
{
   __m128i vec_sgn;
  __m128d x_sign;
  __m128d gr1;
  __m128d gr2;
  __m128d z1;
  __m128d z2;
  __m128d x2;
  
  //extract the sign of x and make x positive
  vec_sgn = _mm_load_si128((__m128i*)sign_int_d);
  x_sign = _mm_and_pd((__m128d)vec_sgn, x);
  x = _mm_xor_pd(x_sign, x);
  
  //if x > sqrt(2)+1, then set x = -1/x
  //else if x > sqrt(2)-1, then set x = (x-1)/(x+1)
  gr1 = _mm_cmpgt_pd(x, sqrt2_plus_1d);
  gr2 = _mm_cmpgt_pd(x, sqrt2_minus_1d);
  z1 = _mm_div_pd(negoned, x);
  z1 = _mm_and_pd(z1, gr1);
  z2 = _mm_add_pd(x, negoned);
  x2 = _mm_sub_pd(x, negoned);
  z2 = _mm_div_pd(z2, x2);
  x = _mm_andnot_pd(gr2, x);
  gr2 = _mm_xor_pd(gr2, gr1);
  z2 = _mm_and_pd(z2, gr2);
  x = _mm_or_pd(x, z1);
  x = _mm_or_pd(x, z2);
  
  //compute z1 = atan(x) from a Chebyshev polynomial (in monomial form) using Horner's scheme
  x2 = _mm_mul_pd(x, x);
  z1 = _mm_mul_pd(x2, atanD15);
  z1 = _mm_add_pd(z1, atanD14);
  z2 = _mm_mul_pd(x2, z1);
  z2 = _mm_add_pd(z2, atanD13);
  z1 = _mm_mul_pd(x2, z2);
  z1 = _mm_add_pd(z1, atanD12);
  z2 = _mm_mul_pd(x2, z1);
  z2 = _mm_add_pd(z2, atanD11);
  z1 = _mm_mul_pd(x2, z2);
  z1 = _mm_add_pd(z1, atanD10);
  z2 = _mm_mul_pd(x2, z1);
  z2 = _mm_add_pd(z2, atanD9);
  z1 = _mm_mul_pd(x2, z2);
  z1 = _mm_add_pd(z1, atanD8);
  z2 = _mm_mul_pd(x2, z1);
  z2 = _mm_add_pd(z2, atanD7);
  z1 = _mm_mul_pd(x2, z2);
  z1 = _mm_add_pd(z1, atanD6);
  z2 = _mm_mul_pd(x2, z1);
  z2 = _mm_add_pd(z2, atanD5);
  z1 = _mm_mul_pd(x2, z2);
  z1 = _mm_add_pd(z1, atanD4);
  z2 = _mm_mul_pd(x2, z1);
  z2 = _mm_add_pd(z2, atanD3);
  z1 = _mm_mul_pd(x2, z2);
  z1 = _mm_add_pd(z1, atanD2);
  z2 = _mm_mul_pd(x2, z1);
  z2 = _mm_add_pd(z2, atanD1);
  z1 = _mm_mul_pd(x2, z2);
  z1 = _mm_add_pd(z1, atanD0);
  z1 = _mm_mul_pd(z1, x);
  
  //add either pi/4 or pi/2 to z1 and set result to x2, depending on the initial value of x
  x = _mm_and_pd(gr1, pi_over_twod);
  x2 = _mm_and_pd(gr2, pi_over_fourd);
  x2 = _mm_or_pd(x2, x);
  x2 = _mm_add_pd(x2, z1);
  
  //recover the original sign of x
  x2 = _mm_xor_pd(x2, x_sign);
  
  return x2;
}










