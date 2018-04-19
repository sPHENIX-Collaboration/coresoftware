#ifndef PHG4Utils__H
#define PHG4Utils__H

#ifndef __CINT__
#include <Geant4/G4RotationMatrix.hh>
#endif

#include <string>

class G4LogicalVolume;
class G4VisAttributes;
class G4VSolid;

class PHG4Utils
{
 public:
  static double GetLengthForRapidityCoverage( const double radius, const double eta );
  static double GetLengthForRapidityCoverage( const double radius);
  static void SetPseudoRapidityCoverage( const double eta);
  static void SetColour(G4VisAttributes* att, const std::string &mat);
  static double get_theta(const double eta);
  static double get_eta(const double theta);
  static std::pair<double, double> get_etaphi(const double x, const double y, const double z);
  static double get_eta(const double radius, const double z);
#ifndef __CINT__
  static void DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm = nullptr);
  static void DisplayVolume(G4LogicalVolume *checksolid, G4LogicalVolume *logvol, G4RotationMatrix *rotm = nullptr);
#endif

 private:
  static double _eta_coverage;

};

#endif
