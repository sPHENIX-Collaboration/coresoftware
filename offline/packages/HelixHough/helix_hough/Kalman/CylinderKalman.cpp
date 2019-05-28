#include "CylinderKalman.h"

#include "HelixKalmanState.h"
#include "SimpleHit3D.h"

#include <Eigen/Core>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>

using namespace std;
using namespace Eigen;

static inline float sign(float x) {
  return ((float)(x > 0.)) - ((float)(x < 0.));
}

CylinderKalman::CylinderKalman(vector<float>& detector_radii,
                               vector<float>& detector_material, float B)
    : HelixKalman(B), nlayers(detector_radii.size()), signk_store(0.) {
  det_rad = detector_radii;
  det_scatter_variance = detector_material;
  for (unsigned int i = 0; i < det_scatter_variance.size(); ++i) {
    det_scatter_variance[i] = det_scatter_variance[i] * 0.0136 * 0.0136;
  }
}

CylinderKalman::~CylinderKalman() {}

void CylinderKalman::calculate_dxda(SimpleHit3D& hit, HelixKalmanState& state,
                                    Matrix<float, 3, 5>& dxda, float& x,
                                    float& y, float& z) {
  float phi = state.phi;
  float d = state.d;
  float k = state.kappa;
  float z0 = state.z0;
  float dzdl = state.dzdl;

  float rad_det = det_rad[hit.get_layer()];

  float cosphi = cos(phi);
  float sinphi = sin(phi);

  //   float signk = (float)(sign(k));
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
  float diff1 = (x1 - hit.get_x()) * (x1 - hit.get_x()) +
                (y1 - hit.get_y()) * (y1 - hit.get_y());
  float diff2 = (x2 - hit.get_x()) * (x2 - hit.get_x()) +
                (y2 - hit.get_y()) * (y2 - hit.get_y());
  float signk = 0.;
  if (diff1 < diff2) {
    signk = 1.;
  } else {
    signk = -1.;
  }
  signk_store = signk;
  x = P2x + signk * ux * h;
  y = P2y + signk * uy * h;

  // dx/dphi

  // (d/dA)[ (d*k +1.)*cos(phi) ]
  float dkcxdA = -kd * sinphi;
  // (d/dA)[ (d*k +1.)*sin(phi) ]
  float dkcydA = kd * cosphi;
  // (d/dA)[ 1/(k*d + 1) ]
  float dkd_invdA = 0.;
  // 0.5*[ (k*R2 + ( d*d*k + 2.*d ))*dkd_inv/dA + kd_inv*(R2*dk/dA +
  // 2*dd/dA*(k*d + 1) + (d^2)*dk/dA) ]
  float dadA = 0.;
  float dtmp1dA = dadA * kd_inv + a * dkd_invdA;
  float dP2xdA = tmp1 * dkcxdA + kcx * dtmp1dA;
  float duxdA = -kd_inv * dkcydA - kcy * dkd_invdA;
  float dhdA = (0.5 / sqrt(R2 - a * a)) * (-2. * a * dadA);
  dxda(0, 0) = dP2xdA + signk * (duxdA * h + ux * dhdA);

  // dy/dphi

  float duydA = kd_inv * dkcxdA + kcx * dkd_invdA;
  float dP2ydA = tmp1 * dkcydA + kcy * dtmp1dA;
  dxda(1, 0) = dP2ydA + signk * (duydA * h + uy * dhdA);

  // dx/d

  //   (d/dA)[ (d*k +1.)*cos(phi) ]
  dkcxdA = cosphi * k;
  //   (d/dA)[ (d*k +1.)*sin(phi) ]
  dkcydA = sinphi * k;
  //   (d/dA)[ 1/(k*d + 1) ]
  dkd_invdA = -kd_inv * kd_inv * k;
  // 0.5*[ (k*R2 + ( d*d*k + 2.*d ))*dkd_inv/dA + kd_inv*(R2*dk/dA +
  // 2*dd/dA*(k*d + 1) + (d^2)*dk/dA) ]
  dadA =
      0.5 * ((k * R2 + (d * d * k + 2. * d)) * dkd_invdA + kd_inv * (2. * kd));
  dtmp1dA = dadA * kd_inv + a * dkd_invdA;
  dP2xdA = tmp1 * dkcxdA + kcx * dtmp1dA;
  duxdA = -kd_inv * dkcydA - kcy * dkd_invdA;
  dhdA = (0.5 / sqrt(R2 - a * a)) * (-2. * a * dadA);
  dxda(0, 1) = dP2xdA + signk * (duxdA * h + ux * dhdA);

  // dy/d

  duydA = kd_inv * dkcxdA + kcx * dkd_invdA;
  dP2ydA = tmp1 * dkcydA + kcy * dtmp1dA;
  dxda(1, 1) = dP2ydA + signk * (duydA * h + uy * dhdA);

  // dx/dk

  //   (d/dA)[ (d*k +1.)*cos(phi) ]
  dkcxdA = cosphi * d;
  //   (d/dA)[ (d*k +1.)*sin(phi) ]
  dkcydA = sinphi * d;
  //   (d/dA)[ 1/(k*d + 1) ]
  dkd_invdA = -kd_inv * kd_inv * d;
  // 0.5*[ (k*R2 + ( d*d*k + 2.*d ))*dkd_inv/dA + kd_inv*(R2*dk/dA +
  // 2*dd/dA*(k*d + 1) + (d^2)*dk/dA) ]
  dadA = 0.5 *
         ((k * R2 + (d * d * k + 2. * d)) * dkd_invdA + kd_inv * (R2 + d * d));
  dtmp1dA = dadA * kd_inv + a * dkd_invdA;
  dP2xdA = tmp1 * dkcxdA + kcx * dtmp1dA;
  duxdA = -kd_inv * dkcydA - kcy * dkd_invdA;
  dhdA = (0.5 / sqrt(R2 - a * a)) * (-2. * a * dadA);
  dxda(0, 2) = dP2xdA + signk * (duxdA * h + ux * dhdA);

  // dy/k

  duydA = kd_inv * dkcxdA + kcx * dkd_invdA;
  dP2ydA = tmp1 * dkcydA + kcy * dtmp1dA;
  dxda(1, 2) = dP2ydA + signk * (duydA * h + uy * dhdA);

  // now for the z direction

  float sign_dzdl = sign(dzdl);
  float onedzdl2_inv = 1. / (1. - dzdl * dzdl);
  float startx = d * cosphi;
  float starty = d * sinphi;
  float D = sqrt((startx - x) * (startx - x) + (starty - y) * (starty - y));
  float D_inv = 1. / D;
  float v = 0.5 * k * D;
  z = 0.;
  if (v > 0.1) {
    if (v >= 0.999999) {
      v = 0.999999;
    }
    float s = 2. * asin(v) / k;
    float s_inv = 1. / s;
    float sqrtvv = sqrt(1 - v * v);
    float dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    z = z0 + sign_dzdl * dz;
    float dz_inv = 1. / dz;
    float dz2 = dz * dz;

    // phi
    float dstartxdA = -d * sinphi;
    float dstartydA = d * cosphi;
    float dkdA = 0.;
    float ddzdldA = 0.;
    float dz0dA = 0.;

    float dDdA = 0.5 * D_inv * (2. * (startx - x) * dstartxdA +
                                2. * (starty - y) * dstartydA);
    float dvdA = 0.5 * (k * dDdA + D * dkdA);
    float dsdA = (2. / (k * sqrtvv)) * dvdA - (s / k) * dkdA;
    float ddzdA = 0.5 * dz_inv *
                  (2. * (dsdA)*dz2 * s_inv +
                   s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 0) = dz0dA + sign_dzdl * ddzdA;

    // d
    dstartxdA = cosphi;
    dstartydA = sinphi;
    dkdA = 0.;
    ddzdldA = 0.;
    dz0dA = 0.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dvdA = 0.5 * (k * dDdA + D * dkdA);
    dsdA = (2. / (k * sqrtvv)) * dvdA - (s / k) * dkdA;
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 1) = dz0dA + sign_dzdl * ddzdA;

    // k
    dstartxdA = 0.;
    dstartydA = 0.;
    dkdA = 1.;
    ddzdldA = 0.;
    dz0dA = 0.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dvdA = 0.5 * (k * dDdA + D * dkdA);
    dsdA = (2. / (k * sqrtvv)) * dvdA - (s / k) * dkdA;
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 2) = dz0dA + sign_dzdl * ddzdA;

    // z0
    dstartxdA = 0.;
    dstartydA = 0.;
    dkdA = 0.;
    ddzdldA = 0.;
    dz0dA = 1.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dvdA = 0.5 * (k * dDdA + D * dkdA);
    dsdA = (2. / (k * sqrtvv)) * dvdA - (s / k) * dkdA;
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 3) = dz0dA + sign_dzdl * ddzdA;

    // dzdl
    dstartxdA = 0.;
    dstartydA = 0.;
    dkdA = 0.;
    ddzdldA = 1.;
    dz0dA = 0.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dvdA = 0.5 * (k * dDdA + D * dkdA);
    dsdA = (2. / (k * sqrtvv)) * dvdA - (s / k) * dkdA;
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 4) = dz0dA + sign_dzdl * ddzdA;
  } else {
    float s = 0.;
    float temp1 = k * D * 0.5;
    temp1 *= temp1;
    float temp2 = D * 0.5;
    s += 2. * temp2;
    temp2 *= temp1;
    s += temp2 / 3.;
    temp2 *= temp1;
    s += (3. / 20.) * temp2;
    temp2 *= temp1;
    s += (5. / 56.) * temp2;
    float s_inv = 1. / s;
    float dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    z = z0 + sign_dzdl * dz;
    float dz_inv = 1. / dz;
    float dz2 = dz * dz;

    float k2 = k * k;
    float k3 = k2 * k;
    float k4 = k3 * k;
    float k5 = k4 * k;
    float k6 = k5 * k;
    float D2 = D * D;
    float D3 = D2 * D;
    float D4 = D3 * D;
    float D5 = D4 * D;
    float D6 = D5 * D;
    float D7 = D6 * D;

    // phi
    float dstartxdA = -d * sinphi;
    float dstartydA = d * cosphi;
    float dkdA = 0.;
    float ddzdldA = 0.;
    float dz0dA = 0.;

    float dDdA = 0.5 * D_inv * (2. * (startx - x) * dstartxdA +
                                2. * (starty - y) * dstartydA);
    float dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    float ddzdA = 0.5 * dz_inv *
                  (2. * (dsdA)*dz2 * s_inv +
                   s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 0) = dz0dA + sign_dzdl * ddzdA;

    // d
    dstartxdA = cosphi;
    dstartydA = sinphi;
    dkdA = 0.;
    ddzdldA = 0.;
    dz0dA = 0.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 1) = dz0dA + sign_dzdl * ddzdA;

    // k
    dstartxdA = 0.;
    dstartydA = 0.;
    dkdA = 1.;
    ddzdldA = 0.;
    dz0dA = 0.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 2) = dz0dA + sign_dzdl * ddzdA;

    // z0
    dstartxdA = 0.;
    dstartydA = 0.;
    dkdA = 0.;
    ddzdldA = 0.;
    dz0dA = 1.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 3) = dz0dA + sign_dzdl * ddzdA;

    // dzdl
    dstartxdA = 0.;
    dstartydA = 0.;
    dkdA = 0.;
    ddzdldA = 1.;
    dz0dA = 0.;

    dDdA = 0.5 * D_inv *
           (2. * (startx - x) * dstartxdA + 2. * (starty - y) * dstartydA);
    dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dxda(2, 4) = dz0dA + sign_dzdl * ddzdA;
  }

  // we want dx/dnu instead of dx/dk
  // nu = sqrt(k)
  // dx/dnu = dx/dk * dk/dnu
  // k = nu^2 , dk/dnu = 2.*nu
  for (int i = 0; i < 3; ++i) {
    dxda(i, 2) *= 2. * state.nu;
  }
}

void CylinderKalman::subtractProjections(Matrix<float, 2, 1>& m,
                                         Matrix<float, 2, 1>& ha,
                                         Matrix<float, 2, 1>& diff) {
  diff = m - ha;
  if (diff(0, 0) > M_PI) {
    diff(0, 0) -= (2. * M_PI);
  }
  if (diff(0, 0) < (-M_PI)) {
    diff(0, 0) += (2. * M_PI);
  }
}

void CylinderKalman::calculateProjections(SimpleHit3D& hit,
                                          HelixKalmanState& state,
                                          Matrix<float, 2, 5>& H,
                                          Matrix<float, 2, 1>& ha) {
  det_rad[hit.get_layer()] =
      sqrt(hit.get_x() * hit.get_x() + hit.get_y() * hit.get_y());

  // calculate dx/dA, A are the helix params, x is (x,y,z)
  // first calculate for x-y plane
  Matrix<float, 3, 5> dxda = Matrix<float, 3, 5>::Zero(3, 5);
  float x = 0.;
  float y = 0.;
  float z = 0.;
  calculate_dxda(hit, state, dxda, x, y, z);

  // calculate dm/dx , m is (phi, z)
  Matrix<float, 2, 3> dmdx = Matrix<float, 2, 3>::Zero(2, 3);
  // phi = atan2(y, x);
  float r2_inv = 1. / (x * x + y * y);
  dmdx(0, 0) = -y * r2_inv;
  dmdx(0, 1) = x * r2_inv;
  dmdx(1, 2) = 1.;

  H = dmdx * dxda;

  ha(0, 0) = atan2(y, x);
  if (ha(0, 0) < 0.) {
    ha(0, 0) += 2. * M_PI;
  }
  ha(1, 0) = z;
}

void CylinderKalman::calculateMeasurements(SimpleHit3D& hit,
                                           Matrix<float, 2, 1>& m,
                                           Matrix<float, 2, 2>& G) {
  Matrix<float, 2, 2> V = Matrix<float, 2, 2>::Zero(2, 2);

  V(0, 0) = 3.33333333333333426e-01 *
            (
	     (2.0*sqrt(hit.get_size(0,0))) *
	     (2.0*sqrt(hit.get_size(0,0))) +
	     (2.0*sqrt(hit.get_size(1,1))) *
	     (2.0*sqrt(hit.get_size(1,1)))
	     ) /
            ((hit.get_x()) * (hit.get_x()) + (hit.get_y()) * (hit.get_y()));
  
  V(1, 1) = 3.33333333333333426e-01 *
    (2.0*sqrt(hit.get_size(2,2))) *
    (2.0*sqrt(hit.get_size(2,2)));

  G = V.fullPivLu().inverse();

  m = Matrix<float, 2, 1>::Zero(2, 1);
  m(0) = atan2(hit.get_y(), hit.get_x());
  if (m(0) < 0.) {
    m(0) += 2. * M_PI;
  }
  m(1) = hit.get_z();
}

void CylinderKalman::updateIntersection(HelixKalmanState& state, int layer) {
  float phi = state.phi;
  float d = state.d;
  float k = state.kappa;
  float z0 = state.z0;
  float dzdl = state.dzdl;

  float rad_det = det_rad[layer];

  float cosphi = cos(phi);
  float sinphi = sin(phi);

  //   float signk = (float)(sign(k));
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

  float signk = signk_store;
  state.x_int = P2x + signk * ux * h;
  state.y_int = P2y + signk * uy * h;

  float sign_dzdl = sign(dzdl);
  float startx = d * cosphi;
  float starty = d * sinphi;
  float D = sqrt((startx - state.x_int) * (startx - state.x_int) +
                 (starty - state.y_int) * (starty - state.y_int));
  float v = 0.5 * k * D;
  if (v > 0.1) {
    if (v >= 0.999999) {
      v = 0.999999;
    }
    float s = 2. * asin(v) / k;
    float dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    state.z_int = z0 + sign_dzdl * dz;
  } else {
    float s = 0.;
    float temp1 = k * D * 0.5;
    temp1 *= temp1;
    float temp2 = D * 0.5;
    s += 2. * temp2;
    temp2 *= temp1;
    s += temp2 / 3.;
    temp2 *= temp1;
    s += (3. / 20.) * temp2;
    temp2 *= temp1;
    s += (5. / 56.) * temp2;
    float dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    state.z_int = z0 + sign_dzdl * dz;
  }

  state.position = (layer + 1);
}

bool CylinderKalman::calculateScatteringVariance(HelixKalmanState& state,
                                                 float& var) {
  if ((state.position == 0) || (state.position > nlayers)) {
    var = 0.;
    return false;
  } else {
    float k = state.kappa;
    float p_inv = 3.33333333333333314e+02 * k * Bfield_inv *
                  sqrt(1. - state.dzdl * state.dzdl);
    var = 2. * p_inv * p_inv * det_scatter_variance[state.position - 1];
    return true;
  }
}
