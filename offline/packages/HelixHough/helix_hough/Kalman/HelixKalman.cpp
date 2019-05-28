#include "HelixKalman.h"
#include "HelixKalmanState.h"
#include "SimpleHit3D.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <cmath>

using namespace std;
using namespace Eigen;

static inline float sign(float x) {
  return ((float)(x > 0.)) - ((float)(x < 0.));
}

HelixKalman::HelixKalman(float B) : Bfield(B), Bfield_inv(1. / B) {}

HelixKalman::~HelixKalman() {}

void HelixKalman::subtractProjections(Matrix<float, 2, 1>& m,
                                      Matrix<float, 2, 1>& ha,
                                      Matrix<float, 2, 1>& diff) {
  diff = m - ha;
}

void HelixKalman::addHit(SimpleHit3D& hit, HelixKalmanState& state) {
  Matrix<float, 2, 5> H = Matrix<float, 2, 5>::Zero(2, 5);
  Matrix<float, 2, 1> ha = Matrix<float, 2, 1>::Zero(2, 1);
  calculateProjections(hit, state, H, ha);

  Matrix<float, 2, 1> m = Matrix<float, 2, 1>::Zero(2, 1);
  Matrix<float, 2, 2> G = Matrix<float, 2, 2>::Zero(2, 2);
  calculateMeasurements(hit, m, G);

  Matrix<float, 5, 5> Q = Matrix<float, 5, 5>::Zero(5, 5);
  calculateMSCovariance(state, Q);

  Matrix<float, 5, 2> Ht = H.transpose();

  Matrix<float, 5, 5> C_proj = state.C + Q;

  Matrix<float, 5, 5> C_proj_inv = C_proj.fullPivLu().inverse();

  Matrix<float, 5, 5> Cadd = Ht * G * H;

  Matrix<float, 5, 5> Cinv = (C_proj_inv + Cadd);

  state.C = Cinv.fullPivLu().inverse();

  Matrix<float, 5, 2> K = state.C * Ht * G;

  Matrix<float, 2, 1> proj_diff = Matrix<float, 2, 1>::Zero(2, 1);
  subtractProjections(m, ha, proj_diff);

  Matrix<float, 5, 1> state_add = K * proj_diff;

  Matrix<float, 5, 1> old_state = Matrix<float, 5, 1>::Zero(5, 1);
  old_state(0, 0) = state.phi;
  old_state(1, 0) = state.d;
  //   old_state(2,0) = state.kappa;
  old_state(2, 0) = state.nu;
  old_state(3, 0) = state.z0;
  old_state(4, 0) = state.dzdl;

  Matrix<float, 5, 1> new_state = old_state + state_add;

  state.phi = new_state(0, 0);
  state.d = new_state(1, 0);
  state.nu = new_state(2, 0);
  state.kappa = state.nu * state.nu;
  state.z0 = new_state(3, 0);
  state.dzdl = new_state(4, 0);

  int layer = hit.get_layer();
  updateIntersection(state, layer);
  Matrix<float, 5, 1> state_diff = (new_state - old_state);
  Matrix<float, 1, 5> state_diff_T = state_diff.transpose();

  Matrix<float, 1, 1> chi2_add = (state_diff_T * C_proj_inv * state_diff);
  state.chi2 += chi2_add(0, 0);
}

// dA/dA' , where A are the global helix parameters with pivot (0,0,0), and A'
// are the 3 helix parameters phi,k,dzdl with the pivot being the intersection
// at the current layer
void HelixKalman::calculate_dAdAp(HelixKalmanState& state,
                                  Matrix<float, 5, 3>& dAdAp, float& phi_p,
                                  float& cosphi_p, float& sinphi_p) {
  float phi = state.phi;
  float d = state.d;
  float k = state.kappa;
  //   float z0 = state.z0;
  float dzdl = state.dzdl;

  float cosphi = cos(phi);
  float sinphi = sin(phi);

  //   float signk = (float)(sign(k));
  k = fabs(k);

  // calculate primed variables after translating (x0,y0,z0) -> (0,0,0)

  float x0 = state.x_int;
  float y0 = state.y_int;

  phi_p = atan2((1. + k * d) * sinphi - k * y0, (1. + k * d) * cosphi - k * x0);
  cosphi_p = cos(phi_p);
  sinphi_p = sin(phi_p);

  // dA/dA' :

  //   tx = cosphi' + k*x0;
  //   ty = sinphi' + k*y0;
  //   phi = atan2(ty, tx);
  //
  //   dphi/dtx = -ty/(tx^2+ty^2)
  //   dphi/dty = tx/(tx^2+ty^2)

  //   d = (1/k)*(sqrt( tx*tx + ty*ty ) - 1.)
  //   dd/dtx = (1/k)*( (tx*tx + ty*ty)^(-1/2) )*( tx )
  //   dd/dty = (1/k)*( (tx*tx + ty*ty)^(-1/2) )*( ty )

  // The A params are phi,d,k,z0,dzdl
  // The A' params are phi',k',dzdl'

  // first the x-y plane

  float tx = cosphi_p + k * x0;
  float ty = sinphi_p + k * y0;
  float tx2ty2_inv = 1. / (tx * tx + ty * ty);
  float dphidtx = -ty * tx2ty2_inv;
  float dphidty = tx * tx2ty2_inv;
  float dtxdA_p = -sinphi_p;
  float dtydA_p = cosphi_p;
  dAdAp(0, 0) = dphidtx * dtxdA_p + dphidty * dtydA_p;

  dtxdA_p = x0;
  dtydA_p = y0;
  dAdAp(0, 1) = dphidtx * dtxdA_p + dphidty * dtydA_p;

  if (k != 0.) {
    float k_inv = 1. / k;
    float tx2ty2sqrtinv = sqrt(tx2ty2_inv);
    float dddtx = tx2ty2sqrtinv * tx * k_inv;
    float dddty = tx2ty2sqrtinv * ty * k_inv;

    dtxdA_p = -sinphi_p;
    dtydA_p = cosphi_p;
    dAdAp(1, 0) = dddtx * dtxdA_p + dddty * dtydA_p;

    dtxdA_p = x0;
    dtydA_p = y0;
    dAdAp(1, 1) = dphidtx * dtxdA_p + dphidty * dtydA_p -
                  (sqrt(tx * tx + ty * ty) - 1.) * k_inv * k_inv;
  } else {
    // d = x0*cosphi_p + y0*sinphi_p
    dAdAp(1, 0) = y0 * cosphi_p - x0 * sinphi_p;
    dAdAp(1, 1) = (dAdAp(1, 0)) * (dAdAp(1, 0)) * 0.5;
  }

  dAdAp(2, 1) = 1.;

  // now the z direction

  float sign_dzdl = sign(dzdl);

  float dx = d * cosphi;
  float dy = d * sinphi;
  float D = sqrt((x0 - dx) * (x0 - dx) + (y0 - dy) * (y0 - dy));
  float D_inv = 1. / D;
  float v = 0.5 * k * D;
  if (v > 0.1) {
    float k_inv = 1. / k;
    float s = 2. * asin(v) * k_inv;
    float s_inv = 1. / s;
    float dsdv = 2. / (k * sqrt(1 - v * v));
    float dsdk = -s * k_inv;
    float tmp1 = s * s * dzdl * dzdl;
    float tmp2 = dzdl * tmp1;
    float tmp3 = 1. / tmp2;
    float dz = sqrt(tmp1 / (1. - dzdl * dzdl));
    float dz3 = dz * dz * dz;

    // phi'
    float ddxdA = -d * sinphi * dAdAp(0, 0) + cosphi * dAdAp(1, 0);
    float ddydA = d * cosphi * dAdAp(0, 0) + sinphi * dAdAp(1, 0);
    float dkdA = 0.;
    float ddzdldA = 0.;

    float dDdA =
        0.5 * D_inv * (2. * (dx - x0) * ddxdA + 2. * (dy - y0) * ddydA);
    float dvdA = 0.5 * (k * dDdA + D * dkdA);
    float dsdA = dsdv * dvdA + dsdk * dkdA;
    float ddzdA = (dz * s_inv) * dsdA + (dz3 * tmp3) * ddzdldA;
    dAdAp(3, 0) = -sign_dzdl * ddzdA;

    // k'
    ddxdA = 0.;
    ddydA = 0.;
    dkdA = 1.;
    ddzdldA = 0.;

    dDdA = 0.5 * D_inv * (2. * (dx - x0) * ddxdA + 2. * (dy - y0) * ddydA);
    dvdA = 0.5 * (k * dDdA + D * dkdA);
    dsdA = dsdv * dvdA + dsdk * dkdA;
    ddzdA = (dz * s_inv) * dsdA + (dz3 * tmp3) * ddzdldA;
    dAdAp(3, 1) = -sign_dzdl * ddzdA;

    // dzdl'
    ddxdA = 0.;
    ddydA = 0.;
    dkdA = 0.;
    ddzdldA = 1.;

    dDdA = 0.5 * D_inv * (2. * (dx - x0) * ddxdA + 2. * (dy - y0) * ddydA);
    dvdA = 0.5 * (k * dDdA + D * dkdA);
    dsdA = dsdv * dvdA + dsdk * dkdA;
    ddzdA = (dz * s_inv) * dsdA + (dz3 * tmp3) * ddzdldA;
    dAdAp(3, 2) = -sign_dzdl * ddzdA;
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
    float onedzdl2_inv = 1. / (1. - dzdl * dzdl);
    float dz = sqrt(s * s * dzdl * dzdl * onedzdl2_inv);
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

    // phi'
    float ddxdA = -d * sinphi * dAdAp(0, 0) + cosphi * dAdAp(1, 0);
    float ddydA = d * cosphi * dAdAp(0, 0) + sinphi * dAdAp(1, 0);
    float dkdA = 0.;
    float ddzdldA = 0.;

    float dDdA =
        0.5 * D_inv * (2. * (dx - x0) * ddxdA + 2. * (dy - y0) * ddydA);
    float dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    float ddzdA = 0.5 * dz_inv *
                  (2. * (dsdA)*dz2 * s_inv +
                   s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dAdAp(3, 0) = -sign_dzdl * ddzdA;

    // k'
    ddxdA = 0.;
    ddydA = 0.;
    dkdA = 1.;
    ddzdldA = 0.;

    dDdA = 0.5 * D_inv * (2. * (dx - x0) * ddxdA + 2. * (dy - y0) * ddydA);
    dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dAdAp(3, 1) = -sign_dzdl * ddzdA;

    // dzdl'
    ddxdA = 0.;
    ddydA = 0.;
    dkdA = 0.;
    ddzdldA = 1.;

    dDdA = 0.5 * D_inv * (2. * (dx - x0) * ddxdA + 2. * (dy - y0) * ddydA);
    dsdA = dDdA;
    dsdA += (1. / 24.) * (3. * k2 * D2 * dDdA + 2. * k * D3 * dkdA);
    dsdA += (3. / 640.) * (5. * D4 * k4 * dDdA + 4. * k3 * D5 * dkdA);
    dsdA += (5. / 7168.) * (7. * D6 * k6 * dDdA + 6. * k5 * D7 * dkdA);
    ddzdA = 0.5 * dz_inv *
            (2. * (dsdA)*dz2 * s_inv +
             s * s * (2. * dzdl * ddzdldA * onedzdl2_inv * onedzdl2_inv));
    dAdAp(3, 2) = -sign_dzdl * ddzdA;
  }

  dAdAp(4, 2) = 1.;

  // we want to substitute nu for k
  // dnu/dnu' = 1 , just as for kappa
  // for the rest of the variables, dx/dnu = dx/dk * dk/dnu
  dAdAp(0, 1) *= 2. * state.nu;
  dAdAp(1, 1) *= 2. * state.nu;
  dAdAp(3, 1) *= 2. * state.nu;
  dAdAp(4, 1) *= 2. * state.nu;
}

// dA'/dp , where p is the momentum vector
void HelixKalman::calculate_dApdp(HelixKalmanState& state,
                                  Matrix<float, 3, 3>& dApdp, Vector3f& p,
                                  float phi, float cosphi, float sinphi) {
  float k = state.kappa;
  float dzdl = state.dzdl;

  // kappa = 0.003*B/sqrt(px*px + py*py)
  // phi = atan2(-px, py)
  // dzdl' = pz/sqrt(px^2 + py^2 + pz^2)

  float pT = 0.003 * Bfield / k;
  float pT_inv = 1. / pT;
  float pT_inv2 = pT_inv * pT_inv;

  float px = -sinphi * pT;
  float py = cosphi * pT;
  float pz = dzdl * pT / sqrt(1. - dzdl * dzdl);

  // dphi/dpx
  dApdp(0, 0) = -py * pT_inv2;
  // dphi/dpy
  dApdp(0, 1) = px * pT_inv2;

  // dk/dpx = -0.003*B*(px*px+py*py)^(-3/2)*px
  float pTneg3 = (pT_inv * pT_inv2);
  dApdp(1, 0) = -0.003 * Bfield * px * pTneg3;
  // dk/dpy
  dApdp(1, 1) = -0.003 * Bfield * py * pTneg3;

  float pmag = sqrt(pT * pT + pz * pz);
  float pmag_inv = 1. / pmag;
  float pmag_inv_3 = pmag_inv * pmag_inv * pmag_inv;
  // ddzdl/dpx = -pz*(px^2 + py^2 + pz^2)^(-3/2)*px
  dApdp(2, 0) = -pz * pmag_inv_3 * px;
  dApdp(2, 1) = -pz * pmag_inv_3 * py;
  dApdp(2, 2) = -pz * pmag_inv_3 * pz + pmag_inv;

  p(0) = px;
  p(1) = py;
  p(2) = pz;

  // dnu/dx = dk/dx * dnu/dk
  // TODO put this directly into the calculation above
  dApdp(1, 0) *= 0.5 / state.nu;
  dApdp(1, 1) *= 0.5 / state.nu;
}

// dp/db : b is defined by p = R*s*b , with R being a rotation matrix rotating
// the unit vector in the z-direction to the direction of p, and S is a scalar
// with a value of the magnitude of the current value of p
void HelixKalman::calculate_dpdb(Vector3f& p, Matrix<float, 3, 3>& dpdb) {
  Vector3f b = Vector3f::Zero(3);
  b(2) = 1.;

  float p_mag = sqrt(p.dot(p));

  Vector3f p_unit = p / p_mag;

  Vector3f axis = b.cross(p_unit);
  float angle = acos(p_unit.dot(b));

  axis /= sqrt(axis.dot(axis));

  Matrix<float, 3, 3> rot_p = Matrix<float, 3, 3>::Zero(3, 3);
  float cos_p = cos(angle);
  float sin_p = sin(angle);
  rot_p(0, 0) = cos_p + axis(0) * axis(0) * (1. - cos_p);
  rot_p(0, 1) = axis(0) * axis(1) * (1. - cos_p) - axis(2) * sin_p;
  rot_p(0, 2) = axis(0) * axis(2) * (1. - cos_p) + axis(1) * sin_p;
  rot_p(1, 0) = axis(1) * axis(0) * (1. - cos_p) + axis(2) * sin_p;
  rot_p(1, 1) = cos_p + axis(1) * axis(1) * (1. - cos_p);
  rot_p(1, 2) = axis(1) * axis(2) * (1. - cos_p) - axis(0) * sin_p;
  rot_p(2, 0) = axis(2) * axis(0) * (1. - cos_p) - axis(1) * sin_p;
  rot_p(2, 1) = axis(2) * axis(1) * (1. - cos_p) + axis(0) * sin_p;
  rot_p(2, 2) = cos_p + axis(2) * axis(2) * (1. - cos_p);

  dpdb = rot_p * p_mag;
}

// db/dt , with b a vector, and t a being the two rotation angles, one around
// the x axis and one around the y axis.  The derivative db/dt is calculated at
// the value of b being the unit vector in the z direction
void HelixKalman::calculate_dbdt(Matrix<float, 3, 2>& dbdt_out) {
  // bz = (1 + [tan(t1)]^2 + [tan(t2)]^2)^(-1/2)
  // bx = bz*tan(t1)
  // by = bz*tan(t2)

  // dbz/dt1 = -[(1 + [tan(t1)]^2 + [tan(t2)]^2)^(-3/2)]*( tan(t1)*( 1 +
  // [tan(t1)]^2 ) )

  // dbx/dt1 = tan(t1)*(dbz/dt1) + bz*( 1 + [tan(t1)]^2 )
  // dbx/dt2 = tan(t1)*(dbz/dt2)

  // dby/dt1 = tan(t2)*(dbz/dt1)
  // dby/dt2 = tan(t2)*(dbz/dt2) + bz*( 1 + [tan(t2)]^2 )

  // t1 = t2 = 0

  //   float tant1 = tan(t1);
  //   float tant1_2 = tant1*tant1;
  //   float tant2 = tan(t2);
  //   float tant2_2 = tant2*tant2;
  //   float bz = 1./sqrt(1. + tant1_2 + tant2_2);
  //   dbdt_out(2,0) = -bz*bz*bz*tant1*(1. + tant1_2);
  //   dbdt_out(2,1) = -bz*bz*bz*tant2*(1. + tant2_2);
  //
  //   dbdt_out(0,0) = tant1*dbdt_out(2,0) + bz*(1. + tant1_2);
  //   dbdt_out(0,1) = tant1*dbdt_out(2,1);
  //
  //   dbdt_out(1,0) = tant2*dbdt_out(2,0);
  //   dbdt_out(1,1) = tant2*dbdt_out(2,1) + bz*(1. + tant2_2);

  dbdt_out(2, 0) = 0.;
  dbdt_out(2, 1) = 0.;

  dbdt_out(0, 0) = 1.;
  dbdt_out(0, 1) = 0.;

  dbdt_out(1, 0) = 0.;
  dbdt_out(1, 1) = 1.;
}

void HelixKalman::calculateMSCovariance(HelixKalmanState& state,
                                        Matrix<float, 5, 5>& Q) {
  Matrix<float, 5, 3> dAdAp = Matrix<float, 5, 3>::Zero(5, 3);
  float phi_p = 0.;
  float cosphi_p = 0.;
  float sinphi_p = 0.;

  float var = 0.;
  if (calculateScatteringVariance(state, var) == false) {
    return;
  }

  calculate_dAdAp(state, dAdAp, phi_p, cosphi_p, sinphi_p);

  Matrix<float, 3, 3> dApdp = Matrix<float, 3, 3>::Zero(3, 3);
  Vector3f p = Vector3f::Zero(3);
  calculate_dApdp(state, dApdp, p, phi_p, cosphi_p, sinphi_p);

  Matrix<float, 3, 3> dpdb = Matrix<float, 3, 3>::Zero(3, 3);
  calculate_dpdb(p, dpdb);

  Matrix<float, 3, 2> dbdt = Matrix<float, 3, 2>::Zero(3, 2);
  calculate_dbdt(dbdt);

  Matrix<float, 5, 2> dAdt = dAdAp * dApdp * dpdb * dbdt;

  for (unsigned int j = 0; j < 5; ++j) {
    for (unsigned int i = 0; i < 5; ++i) {
      Q(i, j) = (dAdt(i, 0) * dAdt(j, 0) + dAdt(i, 1) * dAdt(j, 1)) * var;
    }
  }
}
