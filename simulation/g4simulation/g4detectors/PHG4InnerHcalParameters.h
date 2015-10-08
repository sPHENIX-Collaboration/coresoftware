#ifndef PHG4INNERHCALPARAMETERS_H
#define PHG4INNERHCALPARAMETERS_H

#include <Geant4/globals.hh>

// contains parameters in G4 internal units
class PHG4InnerHcalParameters
{
 public:
  PHG4InnerHcalParameters();
  virtual ~PHG4InnerHcalParameters() {}
  void print() const;

  G4double inner_radius;
  G4double outer_radius;
  G4double size_z;
  G4double scinti_gap;
  G4double tilt_angle;
  G4int n_scinti_plates;
  G4int n_scinti_tiles;
  G4double scinti_tile_thickness;
  G4double scinti_gap_neighbor;
  G4double scinti_eta_coverage;
  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double x_rot;
  G4double y_rot;
  G4double z_rot;
  G4int active;
  G4int absorberactive;
  G4int ncross;
  G4int blackhole;
  G4String material;
  G4double steplimits;
  G4int absorbertruth;
};

#endif
