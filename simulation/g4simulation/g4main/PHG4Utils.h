// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4UTILS_H
#define G4MAIN_PHG4UTILS_H

#include <string>
#include <utility>  // for pair

class G4VisAttributes;

class PHG4Utils
{
 public:
  static double GetLengthForRapidityCoverage(const double radius, const double eta);
  static double GetLengthForRapidityCoverage(const double radius);
  static void SetPseudoRapidityCoverage(const double eta);
  static void SetColour(G4VisAttributes *att, const std::string &mat);
  static double get_theta(const double eta);
  static double get_eta(const double theta);
  static std::pair<double, double> get_etaphi(const double x, const double y, const double z);
  static double get_eta(const double radius, const double z);
  static std::pair<bool, double> line_and_rectangle_intersect(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy);
  static std::pair<bool, std::pair<double, double>> lines_intersect(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy);
  static double circle_rectangle_intersection(double x1, double y1, double x2, double y2, double mx, double my, double r);
  static double sA(double r, double x, double y);
  static std::string md5sum(const std::string &filename);

 private:
  static double _eta_coverage;
};

#endif  // G4MAIN_PHG4UTILS_H
