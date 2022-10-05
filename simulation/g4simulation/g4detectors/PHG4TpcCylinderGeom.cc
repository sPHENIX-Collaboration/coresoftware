#include "PHG4TpcCylinderGeom.h"
#include "PHG4CylinderCellDefs.h"

#include <phool/phool.h>
#include <cstdlib>

void PHG4TpcCylinderGeom::set_zbins(const int i)
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  nzbins = i;
}

void PHG4TpcCylinderGeom::set_zmin(const double z)
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  zmin = z;
}

int PHG4TpcCylinderGeom::get_zbins() const
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return nzbins;
}

double
PHG4TpcCylinderGeom::get_zmin() const
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return zmin;
}

double
PHG4TpcCylinderGeom::get_zstep() const
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return zstep;
}

void PHG4TpcCylinderGeom::set_zstep(const double z)
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  zstep = z;
}

int PHG4TpcCylinderGeom::get_phibins() const
{
  check_binning_method_phi("PHG4TpcCylinderGeom::get_phibins");
  return nphibins;
}

double
PHG4TpcCylinderGeom::get_phistep() const
{
  check_binning_method_phi("PHG4TpcCylinderGeom::get_phistep");
  return phistep;
}

double
PHG4TpcCylinderGeom::get_phimin() const
{
  check_binning_method_phi("PHG4TpcCylinderGeom::get_phimin");
  return phimin;
}

void PHG4TpcCylinderGeom::set_phibins(const int i)
{
  check_binning_method_phi("PHG4TpcCylinderGeom::set_phibins");
  nphibins = i;
}

void PHG4TpcCylinderGeom::set_phistep(const double r)
{
  check_binning_method_phi("PHG4TpcCylinderGeom::set_phistep");
  phistep = r;
}

void PHG4TpcCylinderGeom::set_phimin(const double r)
{
  check_binning_method_phi("PHG4TpcCylinderGeom::set_phimin");
  phimin = r;
}

int PHG4TpcCylinderGeom::get_etabins() const
{
  check_binning_method_eta("PHG4TpcCylinderGeom::get_etabins");
  return nzbins;
}

double
PHG4TpcCylinderGeom::get_etastep() const
{
  check_binning_method_eta("PHG4TpcCylinderGeom::get_etastep");
  return zstep;
}
double
PHG4TpcCylinderGeom::get_etamin() const
{
  check_binning_method_eta("PHG4TpcCylinderGeom::get_etamin");
  return zmin;
}

void PHG4TpcCylinderGeom::set_etamin(const double z)
{
  check_binning_method_eta("PHG4TpcCylinderGeom::set_etamin");
  zmin = z;
}

void PHG4TpcCylinderGeom::set_etastep(const double z)
{
  check_binning_method_eta("PHG4TpcCylinderGeom::set_etastep");
  zstep = z;
}

void PHG4TpcCylinderGeom::set_etabins(const int i)
{
  check_binning_method_eta("PHG4TpcCylinderGeom::set_etabins");
  nzbins = i;
}

void PHG4TpcCylinderGeom::identify(std::ostream& os) const
{
  os << "PHG4TpcCylinderGeom::identify - layer: " << layer
     << ", radius: " << radius
     << ", thickness: " << thickness;
  switch (binning)
  {
  case PHG4CylinderCellDefs::sizebinning:
    os << ", zbins: " << nzbins
       << ", zmin: " << zmin
       << ", zstepsize: " << zstep;
    break;
  case PHG4CylinderCellDefs::etaphibinning:
    os << ", etabins: " << nzbins
       << ", etamin: " << zmin
       << ", etastepsize: " << zstep;
    break;
  case PHG4CylinderCellDefs::etaslatbinning:
    os << ", etabins: " << nzbins
       << ", etamin: " << zmin
       << ", etastepsize: " << zstep;
    break;
  case PHG4CylinderCellDefs::spacalbinning:
    os << ", etabins: " << nzbins << " for Spacal";
    break;
  default:
    os << "no valid binning method: " << binning << std::endl;
    return;
    break;
  }
  os << ", phimin: " << phimin
     << ", phibins: " << nphibins
     << ", phistep: " << phistep
     << std::endl;
  return;
}

std::pair<double, double>
PHG4TpcCylinderGeom::get_zbounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
  {
    std::cout << PHWHERE << " Asking for invalid bin in z: " << ibin << std::endl;
    exit(1);
  }
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  double zlow = zmin + ibin * zstep;
  double zhigh = zlow + zstep;
  return std::make_pair(zlow, zhigh);
}

std::pair<double, double>
PHG4TpcCylinderGeom::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
  {
    std::cout << PHWHERE << " Asking for invalid bin in z: " << ibin << std::endl;
    exit(1);
  }
  check_binning_method_eta("PHG4TpcCylinderGeom::get_etabounds");
  //  check_binning_method(PHG4CylinderCellDefs::etaphibinning);
  double zlow = zmin + ibin * zstep;
  double zhigh = zlow + zstep;
  return std::make_pair(zlow, zhigh);
}

std::pair<double, double>
PHG4TpcCylinderGeom::get_phibounds(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
  {
    std::cout << PHWHERE << "Asking for invalid bin in phi: " << ibin << std::endl;
    exit(1);
  }

  double philow = phimin + ibin * phistep;
  double phihigh = philow + phistep;
  return std::make_pair(philow, phihigh);
}

int PHG4TpcCylinderGeom::get_zbin(const double z) const
{
  if (z < zmin || z > (zmin + nzbins * zstep))
  {
    //    cout << PHWHERE << "Asking for bin for z outside of z range: " << z << endl;
    return -1;
  }

  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return floor((z - zmin) / zstep);
}

int PHG4TpcCylinderGeom::get_etabin(const double eta) const
{
  if (eta < zmin || eta > (zmin + nzbins * zstep))
  {
    //    cout << "Asking for bin for eta outside of eta range: " << eta << endl;
    return -1;
  }
  check_binning_method_eta();
  return floor((eta - zmin) / zstep);
}

int PHG4TpcCylinderGeom::get_phibin(const double phi) const
{
  double norm_phi = phi;
  if (phi < phimin || phi > (phimin + nphibins * phistep))
  {
    int nwraparound = -floor((phi - phimin) * 0.5 / M_PI);
    norm_phi += 2 * M_PI * nwraparound;
  }
  check_binning_method_phi();
  return floor((norm_phi - phimin) / phistep);
  
}

int PHG4TpcCylinderGeom::get_phibin_new(const double phi) const
{
  const int side = 0 ;
  //double norm_phi = phi;
  int phi_bin = -1;
  if (phi >= phimin && phi <= -1*phimin){
    for(std::size_t s=0;s<sector_max_Phi[side].size();s++){
      if(phi < sector_max_Phi[side][s] && phi > sector_min_Phi[side][s]){
        phi_bin = (int) ( std::abs(sector_max_Phi[side][s] - phi)/phistep  + nphibins/12 * s);
        break;
      }
      if (s==11){
        if(phi < sector_max_Phi[side][s] && phi >= -M_PI){
          phi_bin = int( std::abs(sector_max_Phi[side][s] - phi)/phistep ) + nphibins/12 * s;
          break;
        }
        if(phi > sector_min_Phi[side][s]+2*M_PI && phi <= M_PI){
          phi_bin = int( std::abs(sector_max_Phi[side][s]+2*M_PI - phi)/phistep ) + nphibins/12 * s;
          break;
        }
      }
    }
  }
  else{
    std::cout << PHWHERE << "Asking for bin for phi outside of phi range: " << phi << std::endl;
    exit(1);
  }
  if (phi_bin<0){
    for(std::size_t s=0;s<sector_max_Phi[side].size();s++){    
      std::cout << sector_max_Phi[side][s]<< " > phi > " << sector_min_Phi[side][s] << std::endl;
    }
    std::cout << PHWHERE << "PHG4TpcCylinderGeom::get_phibin - Phi-bin is negative: " << phi << std::endl;
    exit(1);
  }
  return phi_bin;
}

double
PHG4TpcCylinderGeom::get_zcenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
  {
    std::cout << PHWHERE << "Asking for invalid bin in z: " << ibin << std::endl;
    exit(1);
  }
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return zmin + (ibin + 0.5) * zstep;
}

double
PHG4TpcCylinderGeom::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
  {
    std::cout << PHWHERE << "Asking for invalid bin in eta: " << ibin << std::endl;
    std::cout << "minbin: 0, maxbin " << nzbins << std::endl;
    exit(1);
  }
  check_binning_method_eta();
  return zmin + (ibin + 0.5) * zstep;
}

double
PHG4TpcCylinderGeom::get_phicenter(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
  {
    std::cout << PHWHERE << "Asking for invalid bin in phi: " << ibin << std::endl;
    exit(1);
  }

  check_binning_method_phi();

  return (phimin + (ibin + 0.5) * phistep);
}


double
PHG4TpcCylinderGeom::get_phicenter_new(const int ibin) const
{
  //double phi_center = -999;
  if (ibin < 0 || ibin > nphibins)
  {
    std::cout << PHWHERE << "Asking for invalid bin in phi: " << ibin << std::endl;
    exit(1);
  }

  check_binning_method_phi();

  const int side = 0 ;
  unsigned int pads_per_sector = nphibins / 12;
  unsigned int sector = ibin / pads_per_sector;  
  //std::cout << "PHG4TpcCylinderGeom::get_phicenter : \n ibin = " << ibin 
  //          << " sector:" << sector << " nphibins = " << nphibins << " phistep = " << phistep << std::endl;
  //std::cout << sector_max_Phi[side][sector]  << " - " 
  //          <<sector_min_Phi[side][sector] << std::endl; 

  return (sector_max_Phi[side][sector] - (ibin + 0.5 - sector * pads_per_sector) * phistep);
}

std::string
PHG4TpcCylinderGeom::methodname(const int i) const
{
  switch (i)
  {
  case PHG4CylinderCellDefs::sizebinning:
    return "Bins in cm";
    break;
  case PHG4CylinderCellDefs::etaphibinning:
    return "Eta/Phi bins";
    break;
  case PHG4CylinderCellDefs::etaslatbinning:
    return "Eta/numslat bins";
    break;
  case PHG4CylinderCellDefs::spacalbinning:
    return "SPACAL Tower bins";
    break;
  default:
    break;
  }
  return "Unknown";
}

void PHG4TpcCylinderGeom::check_binning_method(const int i) const
{
  if (binning != i)
  {
    std::cout << "different binning method used " << methodname(binning)
              << ", not : " << methodname(i)
              << std::endl;
    exit(1);
  }
  return;
}

void PHG4TpcCylinderGeom::check_binning_method_eta(const std::string& src) const
{
  if (binning != PHG4CylinderCellDefs::etaphibinning &&
      binning != PHG4CylinderCellDefs::etaslatbinning &&
      binning != PHG4CylinderCellDefs::spacalbinning)
  {
    if (src.size())
      std::cout << src << " : ";

    std::cout << "different binning method used " << methodname(binning)
              << ", not : " << methodname(PHG4CylinderCellDefs::etaphibinning)
              << " or " << methodname(PHG4CylinderCellDefs::etaslatbinning)
              << " or " << methodname(PHG4CylinderCellDefs::spacalbinning)
              << std::endl;
    exit(1);
  }
  return;
}

void PHG4TpcCylinderGeom::check_binning_method_phi(const std::string& src) const
{
  if (binning != PHG4CylinderCellDefs::etaphibinning &&
      binning != PHG4CylinderCellDefs::sizebinning &&
      binning != PHG4CylinderCellDefs::etaslatbinning &&
      binning != PHG4CylinderCellDefs::spacalbinning)
  {
    if (src.size())
      std::cout << src << " : ";

    std::cout << "different binning method used " << methodname(binning)
              << ", not : " << methodname(PHG4CylinderCellDefs::etaphibinning)
              << " or " << methodname(PHG4CylinderCellDefs::sizebinning)
              << " or " << methodname(PHG4CylinderCellDefs::etaslatbinning)
              << " or " << methodname(PHG4CylinderCellDefs::spacalbinning)
              << std::endl;
    exit(1);
  }
  return;
}
