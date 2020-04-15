// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $                                                                                             

/*!
 * \file PHG4CylinderCellGeom_Sectored.h
 * \brief Cylindrical detector geom for detectors that tile phi in sectors, such that there are dead areas periodically in phi
 * \author Ross Corliss <ross.corliss@stonybrook.edu>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4DETECTORS_PHG4CYLINDERCELLGEOMSECTORED_H
#define G4DETECTORS_PHG4CYLINDERCELLGEOMSECTORED_H

#include "PHG4CylinderCellGeom.h"

#include <phool/phool.h>
#include <cmath>
#include <cstdlib>
#include <iostream>                // for cout, ostream
#include <map>
#include <utility>                 // for pair

/*!
 * \brief PHG4CylinderCellGeom_Sectored
 */
class PHG4CylinderCellGeom_Sectored : public PHG4CylinderCellGeom
{
public:
  PHG4CylinderCellGeom_Sectored();
  virtual
  ~PHG4CylinderCellGeom_Sectored();

  //from base class:
  void identify(std::ostream& os = std::cout) const;
  virtual std::pair<double, double>    get_phibounds(const int ibin) const;
  virtual double  get_phicenter(const int ibin) const;
  virtual int  get_phibin(const double phi) const;

  
  //sectored-cylinder specific:
  void set_sectors(const int i); //set the number of sectors in azimuth
  void set_phibins_per_sector(const int i); //set the number of phibins per sector
  void set_sectormargin(const double phi); //the dead area on each edge of 1/nsectors azimuth region
  //void set_sector_inset(const double phi); //the distance from the edge of the dead area to the center of the nearest bin (allows the edge bins to have a different effective geometry)

  //sector-wise getters:
  std::pair<int,int> get_sector_phibin_bounds(const int s); //get the first phibin in range, and the one just past the end.
  std::pair<double, double> get_sector_phibounds(const int s); //get the double edges of the active area of the sector
  int get_sector(const double phi);
  int get_sector(const int i);//gets the sector number for phibin i;
  int get_sectors();//returns number of sectors, per 'get_zbins' etc.
  double get_sectormargin();
  double get_sectorstep();
  double get_phistep();//for thoroughness.  technically the base class is fine.
  int get_phibins_per_sector();
  

  //bin getters:
  bool check_phi_active(const double phi); //returns true if this is an active area

  //don't need?
  typedef std::pair<double, double> bound_t; //ordered pair: lower and upper bounds of a cell
  //typedef std::map<int, bound_t> bound_map_t; //pair of pad index and bounds pair
 protected:
  int nsectors;
  int nphipersector;
  double sectorstep;
  double sectormargin;
  //double sectorinset;

ClassDef(PHG4CylinderCellGeom_Sectored,1)

};

#endif /* G4DETECTORS_PHG4CYLINDERCELLGEOMSECTORED_H */
