#include "Track3D.h"

#include <Eigen/Core>                       // for MatrixXf, Product, Transpose
#include <Eigen/Dense>

#include <cfloat>
#include <climits>
#include <cmath>                           // for sqrt
#include <memory>                           // for allocator_traits<>::value...

using namespace std;
using namespace Eigen;

Track3D::Track3D()
  : hits(std::vector<Cluster3D> ()),
    phi(FLT_MAX), 
    d(FLT_MAX), 
    kappa(FLT_MAX),
    dzdl(FLT_MAX), 
    z0(FLT_MAX),
    index(UINT_MAX),
    vertex_id(UINT_MAX)
{
    
}	

void Track3D::reset() {
	hits.clear();
	cluster_ids.clear();
}

float Track3D::fit_track(float scale) {
/*
  cout << "Track3D: "
      << "id: " << index 
      << "(phi,d,kappa,dzdl,z0) = (" << phi << "," << d << "," << kappa <<","<<dzdl<<","<<z0 << ") "
      << endl;
*/
  std::vector<float> chi2_hit;
  chi2_hit.clear();
  std::vector<float> xy_res;
  std::vector<float> xy_res_inv;
  std::vector<float> z_res;
  std::vector<float> z_res_inv;

  for (unsigned int i = 0; i < hits.size(); i++) {

    float ex = (2.0*sqrt(hits[i].get_size(0,0))) * scale;
    float ey = (2.0*sqrt(hits[i].get_size(1,1))) * scale;
    float ez = (2.0*sqrt(hits[i].get_size(2,2))) * scale;

    if (hits[i].get_layer() < 0) {
      ex = 0.0001 * scale;
      ey = 0.0001 * scale;
      ez = 0.0001 * scale;
    }

    xy_res.push_back(sqrt( ex * ex + ey * ey ));
    xy_res_inv.push_back(1. / xy_res.back());
    z_res.push_back( ez );
    z_res_inv.push_back(1. / z_res.back());
  }

  chi2_hit.resize(hits.size(), 0.);

  MatrixXf y = MatrixXf::Zero(hits.size(), 1);
  for (unsigned int i = 0; i <hits.size(); i++) {
    y(i, 0) = (pow(hits[i].get_x(), 2) + pow(hits[i].get_y(), 2));
    y(i, 0) *= xy_res_inv[i];
  }

  MatrixXf X = MatrixXf::Zero(hits.size(), 3);
  for (unsigned int i = 0; i < hits.size(); i++) {
    X(i, 0) = hits[i].get_x();
    X(i, 1) = hits[i].get_y();
    X(i, 2) = -1.;
    X(i, 0) *= xy_res_inv[i];
    X(i, 1) *= xy_res_inv[i];
    X(i, 2) *= xy_res_inv[i];
  }

  MatrixXf Xt = X.transpose();

  MatrixXf prod = Xt * X;

  MatrixXf Xty = Xt * y;
  MatrixXf beta = prod.ldlt().solve(Xty);

  float cx = beta(0, 0) * 0.5;
  float cy = beta(1, 0) * 0.5;
  float r = sqrt(cx * cx + cy * cy - beta(2, 0));

  phi = atan2(cy, cx);
  d = sqrt(cx * cx + cy * cy) - r;
  kappa = 1. / r;

  MatrixXf diff = y - (X * beta);
  MatrixXf chi2 = (diff.transpose()) * diff;

  float dx = d * cos(phi);
  float dy = d * sin(phi);

  MatrixXf y2 = MatrixXf::Zero(hits.size(), 1);
  for (unsigned int i = 0; i < hits.size(); i++) {
    y2(i, 0) = hits[i].get_z();
    y2(i, 0) *= z_res_inv[i];
  }

  MatrixXf X2 = MatrixXf::Zero(hits.size(), 2);
  for (unsigned int i = 0; i < hits.size(); i++) {
    float D = sqrt(pow(dx - hits[i].get_x(), 2) + pow(dy - hits[i].get_y(), 2));
    float s = 0.0;

    if (0.5 * kappa * D > 0.1) {
      float v = 0.5 * kappa * D;
      if (v >= 0.999999) {
        v = 0.999999;
      }
      s = 2. * asin(v) / kappa;
    } else {
      float temp1 = kappa * D * 0.5;
      temp1 *= temp1;
      float temp2 = D * 0.5;
      s += 2. * temp2;
      temp2 *= temp1;
      s += temp2 / 3.;
      temp2 *= temp1;
      s += (3. / 20.) * temp2;
      temp2 *= temp1;
      s += (5. / 56.) * temp2;
    }

    X2(i, 0) = s;
    X2(i, 1) = 1.0;

    X2(i, 0) *= z_res_inv[i];
    X2(i, 1) *= z_res_inv[i];
  }

  MatrixXf Xt2 = X2.transpose();
  MatrixXf prod2 = Xt2 * X2;

  MatrixXf Xty2 = Xt2 * y2;
  MatrixXf beta2 = prod2.ldlt().solve(Xty2);

  MatrixXf diff2 = y2 - (X2 * beta2);
  MatrixXf chi2_z = (diff2.transpose()) * diff2;

  z0 = beta2(1, 0);
  dzdl = beta2(0, 0) / sqrt(1. + beta2(0, 0) * beta2(0, 0));
 

  if (kappa != 0.) {
    r = 1. / kappa;
  } else {
    r = 1.0e10;
  }

  cx = (d + r) * cos(phi);
  cy = (d + r) * sin(phi);

  float chi2_tot = 0.;
  for (unsigned int h = 0; h < hits.size(); h++) {
    float dx1 = hits[h].get_x() - cx;
    float dy1 = hits[h].get_y() - cy;

    float dx2 = hits[h].get_x() + cx;
    float dy2 = hits[h].get_y() + cy;

    float xydiff1 = sqrt(dx1 * dx1 + dy1 * dy1) - r;
    float xydiff2 = sqrt(dx2 * dx2 + dy2 * dy2) - r;
    float xydiff = xydiff2;
    if (fabs(xydiff1) < fabs(xydiff2)) {
      xydiff = xydiff1;
    }

    float ls_xy = xy_res[h];

    chi2_hit[h] = 0.;
    chi2_hit[h] += xydiff * xydiff / (ls_xy * ls_xy);
    chi2_hit[h] += diff2(h, 0) * diff2(h, 0);

    chi2_tot += chi2_hit[h];
  }

  unsigned int deg_of_freedom = 2 * hits.size() - 5;

  return (chi2_tot) / ((float)(deg_of_freedom));

}

