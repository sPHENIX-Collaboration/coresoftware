#include "sPHENIXTracker.h"
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace SeamStress;

static inline double sign(double x) {
  return ((double)(x > 0.)) - ((double)(x < 0.));
}

void sPHENIXTracker::projectToLayer(SimpleTrack3D& seed, unsigned int layer,
                                    float& x, float& y, float& z) {
  float phi = seed.phi;
  float d = seed.d;
  float k = seed.kappa;
  float z0 = seed.z0;
  float dzdl = seed.dzdl;

  float hitx = seed.hits.back().get_x();
  float hity = seed.hits.back().get_y();

  float rad_det = detector_radii[layer];

  float cosphi = cos(phi);
  float sinphi = sin(phi);

  k = fabs(k);

  float kd = (d * k + 1.);
  float kcx = kd * cosphi;
  float kcy = kd * sinphi;
  float kd_inv = 1. / kd;
  float R2 = rad_det * rad_det;
  float a = 0.5 * (k * R2 + (d * d * k + 2. * d)) * kd_inv;
  float tmp1 = a * kd_inv;
  float P2x = kcx * tmp1;
  float P2y = kcy * tmp1;

  float h = sqrt(R2 - a * a);

  float ux = -kcy * kd_inv;
  float uy = kcx * kd_inv;

  float x1 = P2x + ux * h;
  float y1 = P2y + uy * h;
  float x2 = P2x - ux * h;
  float y2 = P2y - uy * h;
  float diff1 = (x1 - hitx) * (x1 - hitx) + (y1 - hity) * (y1 - hity);
  float diff2 = (x2 - hitx) * (x2 - hitx) + (y2 - hity) * (y2 - hity);
  float signk = 0.;
  if (diff1 < diff2) {
    signk = 1.;
  } else {
    signk = -1.;
  }
  x = P2x + signk * ux * h;
  y = P2y + signk * uy * h;

  double sign_dzdl = sign(dzdl);
  double onedzdl2_inv = 1. / (1. - dzdl * dzdl);
  double startx = d * cosphi;
  double starty = d * sinphi;
  double D = sqrt((startx - x) * (startx - x) + (starty - y) * (starty - y));
  double D_inv = 1. / D;
  double v = 0.5 * k * D;
  z = 0.;
  if (v > 0.1) {
    if (v >= 0.999999) {
      v = 0.999999;
    }
    double s = 2. * asin(v) / k;
    double s_inv = 1. / s;
    double sqrtvv = sqrt(1 - v * v);
    double dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    z = z0 + sign_dzdl * dz;
  } else {
    double s = 0.;
    double temp1 = k * D * 0.5;
    temp1 *= temp1;
    double temp2 = D * 0.5;
    s += 2. * temp2;
    temp2 *= temp1;
    s += temp2 / 3.;
    temp2 *= temp1;
    s += (3. / 20.) * temp2;
    temp2 *= temp1;
    s += (5. / 56.) * temp2;
    double s_inv = 1. / s;
    double dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    z = z0 + sign_dzdl * dz;
  }
}
