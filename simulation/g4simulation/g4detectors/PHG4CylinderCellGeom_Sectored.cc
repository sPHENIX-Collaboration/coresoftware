// $Id: $                                                                                             

/*!
 * \file PHG4CylinderCellGeom_Sectored.cc
 * \brief 
 * \author Ross Corliss <ross.corliss@stonybrook.edu>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4CylinderCellGeom_Sectored.h"
#include "PHG4CylinderCellDefs.h"

#include <phool/phool.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;

PHG4CylinderCellGeom_Sectored::PHG4CylinderCellGeom_Sectored():
  nsectors(-1),
  nphipersector(-1),
  sectorstep(NAN)
  sectormargin(NAN)
{
  return;
}

PHG4CylinderCellGeom_Sectored::~PHG4CylinderCellGeom_Sectored()
{
  
}

void
PHG4CylinderCellGeom_Sectored::identify(std::ostream& os) const{
  os << "PHG4CylinderCellGeom::identify - layer: " << layer
     << ", radius: " << radius
     << ", thickness: " << thickness;

  os << ", total phibins: " <<nphibins;

  os << ", phisectors: " << nsectors
     << ", sectorstep: " <<sectorstep
     << ", sectorbins: " << nphipersector
     << ", sectormargin: " << sectormargin
     << ", phistep: " << phistep;

  os << endl;

  return;
}

//from base class that need to be overrided to handle sector boundaries:
std::pair<double, double>
PHG4CylinderCellGeom_Sectored:: get_phibounds(const int ibin) const{
  if (ibin < 0 || ibin > nphibins)    {
    cout << PHWHERE << " Asking for invalid bin in phi: " << ibin << endl;
    exit(1);
  }
  int s=ibin/nphipersector;//sector
  int p=ibin%nphipersector;//phi relative to the sector
  double philow=s*sectorstep+sectormargin+p*phistep;
  double phihigh=philow+phistep;
  
  return make_pair(philow,phihigh);
}
  
double  PHG4CylinderCellGeom_Sectored::get_phicenter(const int ibin) const{
  if (ibin < 0 || ibin > nphibins)    {
    cout << PHWHERE << " Asking for invalid bin in phi: " << ibin << endl;
    exit(1);
  }
  int s=ibin/nphipersector;//sector
  int p=ibin%nphipersector;//phi relative to the sector
  return s*sectorstep+sectormargin+(p+0.5)*phistep;
}

int  PHG4CylinderCellGeom_Sectored::get_phibin(const double phi) const{
  double norm_phi = phi;
  if(phi < 0 || phi > (2*M_PI)) {
    int nwraparound = -floor((phi) * 0.5 / M_PI);
    norm_phi += 2*M_PI*nwraparound;
  }
  int s=floor(phi/sectorstep);
  double phi_rel=norm_phi-s*sectorstep-sectormargin;
  if (phi_rel<0) return -1; //not a valid phibin: in the low side margin
  int p=floor(phi_rel/phistep);
  if (!(p<nphipersector)) return -1;//not a valid phibi: in the high side margin
  
  return p;
}

  
//sectored-cylinder specific:
void PHG4CylinderCellGeom_Sectored::set_sectors(const int i){
  //set the number of sectors in azimuth
  nsectors=i;
  //and update the dependent variables.
  sectorstep=2*M_PI/nsectors;
  nphibins=nsectors*nphipersector;
  phistep=(sectorstep-2*sectormargin)/nphipersector;
  return;
}
  
void PHG4CylinderCellGeom_Sectored::set_phibins_per_sector(const int i){
  //set the number of phibins per sector
  nphipersector=i;
  //and update the dependent variables.
  nphibins=nsectors*nphipersector;
  phistep=(sectorstep-2*sectormargin)/nphipersector;
  return;
}
  
void PHG4CylinderCellGeom_Sectored::set_sectormargin(const double phi){
  //the dead phi area on each edge of each sector
  //note this does not check whether these values are sane!  Expect weird behavior if margin exceeds the total step.
  sectormargin=phi;
  //and update the dependent variables.
  phistep=(sectorstep-2*sectormargin)/nphipersector;
  return;
}

//void PHG4CylinderCellGeom_Sectored::set_sector_inset(const double phi){
  //the distance from the edge of the dead area to the center of the nearest bin (allows the edge bins to have a different effective geometry)
//  sectorinset=phi;
//  return;
//}


//sector-wise getters:
//

std::pair<int,int> PHG4CylinderCellGeom_Sectored::get_sector_phibin_bounds(const int is){
  //get the phibins that mark the included lower bound and excluded upper bound of the sector.
   if (is < 0 || is > nsectors)    {
    cout << PHWHERE << " Asking for invalid sector in phi: " << is << endl;
    exit(1);
  }
   return make_pair(is*nphipersector,(is+1)*nphipersector);
}
   
std::pair<double, double> PHG4CylinderCellGeom_Sectored::get_sector_phibounds(const int is){
  //get the double edges of the active area of the sector
   if (is < 0 || is > nsectors)    {
    cout << PHWHERE << " Asking for invalid sector in phi: " << is << endl;
    exit(1);
  }
   return make_pair(is*sectorstep+sectormargin,(is+1)*sectorstep-sectormargin);
}

  
int PHG4CylinderCellGeom_Sectored::get_sector(const double phi){
  double norm_phi = phi;
  if(phi < 0 || phi > (2*M_PI)) {
    int nwraparound = -floor((phi) * 0.5 / M_PI);
    norm_phi += 2*M_PI*nwraparound;
  }
  return floor(phi/sectorstep);
}
  
int PHG4CylinderCellGeom_Sectored::get_sector(const int ibin){
  //gets the sector number for phibin ibin;
   if (ibin < 0 || ibin > nphibins)    {
    cout << PHWHERE << " Asking for invalid bin in phi: " << ibin << endl;
    exit(1);
  }
  return ibin/nphipersector;//sector
}

double PHG4CylinderCellGeom_Sectored::get_sectormargin(){
  return sectormargin;
}

int PHG4CylinderCellGeom_Sectored::get_sectors(){
  return nsectors;
}

int PHG4CylinderCellGeom_Sectored::get_phibins_per_sector(){
  return nphipersector;
}

double PHG4CylinderCellGeom_Sectored::get_sectorstep(){
  return sectorstep;
}
double PHG4CylinderCellGeom_Sectored::get_phistep(){
  return phistep;
}

bool PHG4CylinderCellGeom_Sectored::check_phi_active(const double phi){
  //returns true if this is an active area
  double norm_phi = phi;
  if(phi < 0 || phi > (2*M_PI)) {
    int nwraparound = -floor((phi) * 0.5 / M_PI);
    norm_phi += 2*M_PI*nwraparound;
  }
  int s=floor(phi/sectorstep);
  double phi_rel=norm_phi-s*sectorstep-sectormargin;
  if (phi_rel<0) return false; //not a valid phibin: in the low side margin
  int p=floor(phi_rel/phistep);
  if (!(p<nphipersector)) return false;//not a valid phibi: in the high side margin
  
  return true;
}







