#pragma once

#include <vector>

class Tpc_FittingTools
{
 public:
  enum FitMode
  {
    FIT_LINEAR  = 0,
    FIT_SAGITTA = 1
  };

  struct FitPoint
  {
    FitPoint();
    FitPoint(double x_, double y_, double w_ = 1.0);

    double x;
    double y;
    double w;
  };

  struct LineFit
  {
    LineFit();

    bool ok;
    double slope;
    double intercept;
    double chi2;
    int ndof;
  };

  struct SagittaFit
  {
    SagittaFit();

    bool ok;
    double S;
    double x0;
    double invR;
    double theta;
    double b;
    double chi2;
    int ndof;
  };

  struct Point
  {
    double x {0.0};
    double y {0.0};
    double z {0.0};
  };

  struct FitResult
  {
    bool ok {false};
    bool is_line {false};
    double d0 {0.0};
    double z0 {0.0};
    double phi0 {0.0};
    double theta {0.0};
    double curvature {0.0};
    double line_x {0.0};
    double line_y {0.0};
    double line_z {0.0};
    double line_dx {0.0};
    double line_dy {0.0};
    double line_dz {0.0};
    double chi2_xy {0.0};
    double chi2_z {0.0};
    int ndof_xy {0};
    int ndof_z {0};
  };

  static double adcWeight(double adc, double maxadc, double power, double floor_frac);

  static bool weightedLineFit(const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& w,
                              double& m,
                              double& b,
                              double& chi2,
                              int& ndof);

  static LineFit fitLine(const std::vector<FitPoint>& points);

  static double sagittaModel(double xrot, double S, double x0, double invR);

  static bool weightedSagittaFit(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 const std::vector<double>& w,
                                 double& S,
                                 double& x0,
                                 double& invR,
                                 double& theta,
                                 double& bline,
                                 double& chi2,
                                 int& ndof);

  static SagittaFit fitSagitta(const std::vector<FitPoint>& points);

  static bool fit(const std::vector<Point>& points, FitResult& fit);
  static bool fitLine3D(const std::vector<Point>& points, FitResult& fit);
};
