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

 private:
  static double _eta_coverage;
};

#endif  // G4MAIN_PHG4UTILS_H
