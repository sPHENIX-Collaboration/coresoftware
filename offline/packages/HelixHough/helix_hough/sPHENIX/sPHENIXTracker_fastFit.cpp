#include "sPHENIXTracker.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

float sPHENIXTracker::fitTrack(SimpleTrack3D& track) {
  vector<float> chi2_hit;
  return sPHENIXTracker::fitTrack(track, chi2_hit);
}

float sPHENIXTracker::fitTrack(SimpleTrack3D& track, vector<float>& chi2_hit) {
  chi2_hit.clear();
  vector<float> xyres;
  vector<float> xyres_inv;
  vector<float> zres;
  vector<float> zres_inv;
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    xyres.push_back(sqrt(track.hits[i].get_ex() * track.hits[i].get_ex() +
                         track.hits[i].get_ey() * track.hits[i].get_ey()));
    xyres_inv.push_back(1. / xyres.back());
    zres.push_back(track.hits[i].get_ez());
    zres_inv.push_back(1. / zres.back());
  }

  chi2_hit.resize(track.hits.size(), 0.);

  MatrixXf y = MatrixXf::Zero(track.hits.size(), 1);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    y(i, 0) = (pow(track.hits[i].get_x(), 2) + pow(track.hits[i].get_y(), 2));
    y(i, 0) *= xyres_inv[i];
  }

  MatrixXf X = MatrixXf::Zero(track.hits.size(), 3);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    X(i, 0) = track.hits[i].get_x();
    X(i, 1) = track.hits[i].get_y();
    X(i, 2) = -1.;
    X(i, 0) *= xyres_inv[i];
    X(i, 1) *= xyres_inv[i];
    X(i, 2) *= xyres_inv[i];
  }

  MatrixXf Xt = X.transpose();

  MatrixXf prod = Xt * X;

  MatrixXf Xty = Xt * y;
  MatrixXf beta = prod.ldlt().solve(Xty);

  float cx = beta(0, 0) * 0.5;
  float cy = beta(1, 0) * 0.5;
  float r = sqrt(cx * cx + cy * cy - beta(2, 0));

  float phi = atan2(cy, cx);
  float d = sqrt(cx * cx + cy * cy) - r;
  float k = 1. / r;

  MatrixXf diff = y - (X * beta);
  MatrixXf chi2 = (diff.transpose()) * diff;

  float dx = d * cos(phi);
  float dy = d * sin(phi);

  MatrixXf y2 = MatrixXf::Zero(track.hits.size(), 1);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    y2(i, 0) = track.hits[i].get_z();
    y2(i, 0) *= zres_inv[i];
  }

  MatrixXf X2 = MatrixXf::Zero(track.hits.size(), 2);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    float D = sqrt(pow(dx - track.hits[i].get_x(), 2) + pow(dy - track.hits[i].get_y(), 2));
    float s = 0.0;

    if (0.5 * k * D > 0.1) {
      float v = 0.5 * k * D;
      if (v >= 0.999999) {
        v = 0.999999;
      }
      s = 2. * asin(v) / k;
    } else {
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
    }

    X2(i, 0) = s;
    X2(i, 1) = 1.0;

    X2(i, 0) *= zres_inv[i];
    X2(i, 1) *= zres_inv[i];
  }

  MatrixXf Xt2 = X2.transpose();
  MatrixXf prod2 = Xt2 * X2;

  MatrixXf Xty2 = Xt2 * y2;
  MatrixXf beta2 = prod2.ldlt().solve(Xty2);

  MatrixXf diff2 = y2 - (X2 * beta2);
  MatrixXf chi2_z = (diff2.transpose()) * diff2;

  float z0 = beta2(1, 0);
  float dzdl = beta2(0, 0) / sqrt(1. + beta2(0, 0) * beta2(0, 0));

  track.phi = phi;
  track.d = d;
  track.kappa = k;
  track.dzdl = dzdl;
  track.z0 = z0;

  if (track.kappa != 0.) {
    r = 1. / track.kappa;
  } else {
    r = 1.0e10;
  }

  cx = (track.d + r) * cos(track.phi);
  cy = (track.d + r) * sin(track.phi);

  float chi2_tot = 0.;
  for (unsigned int h = 0; h < track.hits.size(); h++) {
    float dx1 = track.hits[h].get_x() - cx;
    float dy1 = track.hits[h].get_y() - cy;

    float dx2 = track.hits[h].get_x() + cx;
    float dy2 = track.hits[h].get_y() + cy;

    float xydiff1 = sqrt(dx1 * dx1 + dy1 * dy1) - r;
    float xydiff2 = sqrt(dx2 * dx2 + dy2 * dy2) - r;
    float xydiff = xydiff2;
    if (fabs(xydiff1) < fabs(xydiff2)) {
      xydiff = xydiff1;
    }

    float ls_xy = xyres[h];

    chi2_hit[h] = 0.;
    chi2_hit[h] += xydiff * xydiff / (ls_xy * ls_xy);
    chi2_hit[h] += diff2(h, 0) * diff2(h, 0);

    chi2_tot += chi2_hit[h];
  }

  unsigned int deg_of_freedom = 2 * track.hits.size() - 5;

  return (chi2_tot) / ((float)(deg_of_freedom));
}

float sPHENIXTracker::fitTrack_3(SimpleTrack3D& track) {
  vector<double> chi2_hit;

  vector<double> xyres;
  vector<double> xyres_inv;
  vector<double> zres;
  vector<double> zres_inv;
  for (unsigned int i = 0; i < 3; i++) {
    xyres.push_back(sqrt(track.hits[i].get_ex() * track.hits[i].get_ex() +
                         track.hits[i].get_ey() * track.hits[i].get_ey()));
    xyres_inv.push_back(1. / xyres.back());
    zres.push_back(track.hits[i].get_ez());
    zres_inv.push_back(1. / zres.back());
  }

  chi2_hit.resize(3, 0.);

  Matrix<double, 3, 1> y;
  for (unsigned int i = 0; i < 3; i++) {
    y(i, 0) = (pow(track.hits[i].get_x(), 2) + pow(track.hits[i].get_y(), 2));
    y(i, 0) *= xyres_inv[i];
  }

  Matrix<double, 3, 3> X;
  for (unsigned int i = 0; i < 3; i++) {
    X(i, 0) = track.hits[i].get_x();
    X(i, 1) = track.hits[i].get_y();
    X(i, 2) = -1.;
    X(i, 0) *= xyres_inv[i];
    X(i, 1) *= xyres_inv[i];
    X(i, 2) *= xyres_inv[i];
  }

  Matrix<double, 3, 3> Xt = X.transpose();

  Matrix<double, 3, 3> prod = Xt * X;

  Matrix<double, 3, 1> Xty = Xt * y;
  Matrix<double, 3, 1> beta = prod.ldlt().solve(Xty);

  double cx = beta(0, 0) * 0.5;
  double cy = beta(1, 0) * 0.5;
  double r = sqrt(cx * cx + cy * cy - beta(2, 0));

  double phi = atan2(cy, cx);
  double d = sqrt(cx * cx + cy * cy) - r;
  double k = 1. / r;

  double dx = d * cos(phi);
  double dy = d * sin(phi);

  Matrix<double, 3, 1> y2;
  for (unsigned int i = 0; i < 3; i++) {
    y2(i, 0) = track.hits[i].get_z();
    y2(i, 0) *= zres_inv[i];
  }

  Matrix<double, 3, 2> X2;
  for (unsigned int i = 0; i < 3; i++) {
    double D =
        sqrt(pow(dx - track.hits[i].get_x(), 2) + pow(dy - track.hits[i].get_y(), 2));
    double s = 0.0;

    if (0.5 * k * D > 0.1) {
      double v = 0.5 * k * D;
      if (v >= 0.999999) {
        v = 0.999999;
      }
      s = 2. * asin(v) / k;
    } else {
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
    }

    X2(i, 0) = s;
    X2(i, 1) = 1.0;

    X2(i, 0) *= zres_inv[i];
    X2(i, 1) *= zres_inv[i];
  }

  Matrix<double, 2, 3> Xt2 = X2.transpose();
  Matrix<double, 2, 2> prod2 = Xt2 * X2;

  Matrix<double, 2, 1> Xty2 = Xt2 * y2;
  Matrix<double, 2, 1> beta2 = prod2.ldlt().solve(Xty2);

  Matrix<double, 3, 1> diff2 = y2 - (X2 * beta2);

  double z0 = beta2(1, 0);
  double dzdl = beta2(0, 0) / sqrt(1. + beta2(0, 0) * beta2(0, 0));

  track.phi = phi;
  track.d = d;
  track.kappa = k;
  track.dzdl = dzdl;
  track.z0 = z0;

  if (track.kappa != 0.) {
    r = 1. / track.kappa;
  } else {
    r = 1.0e10;
  }

  cx = (track.d + r) * cos(track.phi);
  cy = (track.d + r) * sin(track.phi);

  double chi2_tot = 0.;
  for (unsigned int h = 0; h < 3; h++) {
    double dx1 = track.hits[h].get_x() - cx;
    double dy1 = track.hits[h].get_y() - cy;

    double dx2 = track.hits[h].get_x() + cx;
    double dy2 = track.hits[h].get_y() + cy;

    double xydiff1 = sqrt(dx1 * dx1 + dy1 * dy1) - r;
    double xydiff2 = sqrt(dx2 * dx2 + dy2 * dy2) - r;
    double xydiff = xydiff2;
    if (fabs(xydiff1) < fabs(xydiff2)) {
      xydiff = xydiff1;
    }

    double ls_xy = xyres[h];

    chi2_hit[h] = 0.;
    chi2_hit[h] += xydiff * xydiff / (ls_xy * ls_xy);
    chi2_hit[h] += diff2(h, 0) * diff2(h, 0);

    chi2_tot += chi2_hit[h];
  }

  return (chi2_tot);
}
