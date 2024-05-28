#include "PHG4TpcCylinderGeom.h"
#include "PHG4CylinderCellDefs.h"

#include <phool/phool.h>

#include <cstdlib>

//#include <cmath> // for std::abs
#include <algorithm> // for std::min_element


float roundTo1e5(float num)
{ 
  // Scale the number up to avoid floating-point imprecision
  float scale = 1e5;
  // Round the number to the nearest integer
  float rounded = std::round(num * scale);
  // Scale it back down to the original range
  return rounded / scale;
}

int PHG4TpcCylinderGeom::nearest_element(const double phi, int s, int side) const
{

  //=======
  // Initialize variables to track the minimum difference and the index
  double min_diff = std::abs((sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 + pow(-1,side)*layer_pad_phi[0] - phi);
  int closest_index = 0;
  
  // Iterate through the vector to find the closest element
  for (int i = 1; i < int(layer_pad_phi.size()); ++i) {
      double current_diff = std::abs((sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 - pow(-1,side)*layer_pad_phi[i] - phi);
      if (current_diff < min_diff) {
          min_diff = current_diff;
          closest_index = i;
      }
  }

  return closest_index;

}


namespace
{
  // streamer for internal 2dimensional arrays
  using array_t = std::array<std::vector<double>, PHG4TpcCylinderGeom::NSides>;
  std::ostream& operator<<(std::ostream& out, const array_t& array)
  {
    out << "{ ";
    for (const auto& iside : array)
    {
      out << "{";
      bool first = true;
      for (const auto& value : iside)
      {
        if (!first)
        {
          out << ", ";
        }
        first = false;
        out << value;
      }
      out << "} ";
    }
    out << " }";
    return out;
  }
}  // namespace

std::ostream& operator<<(std::ostream& out, const PHG4TpcCylinderGeom& geom)
{
  out << "PHG4TpcCylinderGeom - layer: " << geom.layer << std::endl;
  out
      << "  binnig: " << geom.binning
      << ", radius: " << geom.radius
      << ", nzbins: " << geom.nzbins
      << ", zmin: " << geom.zmin
      << ", zstep: " << geom.zstep
      << ", nphibins: " << geom.nphibins
      << ", phimin: " << geom.phimin
      << ", phistep: " << geom.phistep
      << ", thickness: " << geom.thickness
      << std::endl;

  out << "  sector_R_bias: " << geom.sector_R_bias << std::endl;
  out << "  sector_Phi_bias: " << geom.sector_Phi_bias << std::endl;
  out << "  sector_min_Phi: " << geom.sector_min_Phi << std::endl;
  out << "  sector_max_Phi: " << geom.sector_max_Phi << std::endl;

  return out;
}

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


int PHG4TpcCylinderGeom::find_phibin(const double phi, int side) const
{
  //std::cout<< "PHG4TpcCylinderGeom::find_phibin 0 side = " << side <<std::endl;
  //std::cout<<"PHG4TpcCylinderGeom::find_phibin: phi = "<<phi<<std::endl;
  //double dphi = layer_pad_phi[2]-layer_pad_phi[1];

  double norm_phi = roundTo1e5(phi);
  if (phi < phimin || phi > (phimin + nphibins * phistep))
  {
    int nwraparound = -floor((phi - phimin) * 0.5 / M_PI);
    norm_phi += roundTo1e5(2 * M_PI * nwraparound);
  }
   if (phi >  M_PI){
     norm_phi = roundTo1e5(phi - 2* M_PI);
   }
   if (phi < phimin){
     norm_phi = roundTo1e5(phi + 2* M_PI);
   }

  int phi_bin = -1;
  //int sector_number = -1;

  //std::cout<< "PHG4TpcCylinderGeom::find_phibin 1" <<std::endl;
  for (std::size_t s = 0; s < sector_max_Phi[side].size(); s++)
  {
    //std::cout<< "s = "<< s << std::endl;
    if (s == 11 && phi_bin<0)
    {
      //std::cout << "norm_phi = " << norm_phi << " sector_min_Phi[side][s] = " << sector_min_Phi[side][s] << " sector_max_Phi[side][s] = " << sector_max_Phi[side][s] << std::endl;
      if(norm_phi >= sector_min_Phi[side][s] + 2 * M_PI){  
        phi_bin =  nearest_element(-1*norm_phi, s, side);
        //for(long unsigned int i=0;i<layer_pad_phi.size();i++){
        //  double set_pad_phi = (sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 - pow(-1,side)*layer_pad_phi[i];
        //  if(i<layer_pad_phi.size()-1) dphi=layer_pad_phi[i+1]-layer_pad_phi[i];
        //  double phi_min = roundTo1e5(set_pad_phi-dphi/2);
        //  double phi_max = roundTo1e5(set_pad_phi+dphi/2);
        //  //std::cout<< " i= "<< i <<": " << phi_min << "<" << -1.0*norm_phi << "<" << phi_max << " phi_min - phi =" << phi_min-norm_phi << " phi_max - phi =" << phi_max-norm_phi  << std::endl;
        //  if(-1.0*norm_phi<=phi_max && -1.0*norm_phi>=phi_min)
        //  {
        //    phi_bin = abs(int(i) - side * nphibins / 12) + nphibins / 12 * s;
        //    //std::cout<< i << " vs. nearest_element = " << nearest_element(-1*norm_phi, s, side)<< std::endl;
        //    int new_phi_bin =  nearest_element(-1*norm_phi, s, side);
        //    std::cout<< i << " vs. nearest_element = " << new_phi_bin << " dN="<<new_phi_bin-i << " nphibins / 12 = " << nphibins / 12 << std::endl;
        //    break;
        //  }
        //  if(i==layer_pad_phi.size()-1){
        //    if(-1.0*norm_phi<= roundTo1e5(sector_min_Phi[side][s]) && -1.0*norm_phi>=phi_min ){
        //      phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
        //      break;
        //    }          
        //  }
        //  if(i==0){
        //    //std::cout<< "i= "<< i <<": " << sector_max_Phi[side][s] << "<" << -1.0*norm_phi << "<" << phi_max << " phi_min - phi =" << phi_min-norm_phi << " phi_max - phi =" << phi_max-norm_phi  << std::endl;
        //    if(-1.0*norm_phi<=phi_max  &&  -1.0*norm_phi>=roundTo1e5(sector_max_Phi[side][s]) ){
        //      phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
        //      break;
        //    }          
        //  }
        //}
      }
      if (-1.0*norm_phi <= sector_max_Phi[side][s] && -1.0*norm_phi >= -M_PI){
        phi_bin = std::abs(nphibins / 12 - nearest_element(-1*norm_phi, s, side));
        //for(long unsigned int i=0;i<layer_pad_phi.size();i++){
        //  double set_pad_phi = (sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 + pow(-1,side)*layer_pad_phi[i];
        //  if(i<layer_pad_phi.size()-1) dphi=layer_pad_phi[i+1]-layer_pad_phi[i];
        //  double phi_min = roundTo1e5(set_pad_phi-dphi/2);
        //  double phi_max = roundTo1e5(set_pad_phi+dphi/2);
        //  //std::cout<< "i= "<< i <<": " << phi_min << "<" << -1.0*norm_phi << "<" << phi_max << " phi_min - phi =" << phi_min-norm_phi << " phi_max - phi =" << phi_max-norm_phi  << std::endl;
        //  if(-1.0*norm_phi<=phi_max && -1.0*norm_phi>=phi_min){
        //    phi_bin = abs(int(i) - side * nphibins / 12) + nphibins / 12 * s;
        //    //int new_phi_bin = std::abs(nphibins / 12 - nearest_element(-1*norm_phi, s, side));
        //    //std::cout<< i << " vs. nearest_element = " << new_phi_bin) << " dN="<<new_phi_bin-i << " nphibins / 12 = " << nphibins / 12 << std::endl;
        //    break;
        //  }  
        //  if(i==layer_pad_phi.size()-1){
        //    //std::cout<< "i= "<< i <<": " << phi_min << "<" << -1.0*norm_phi << "<" << sector_min_Phi[side][s] << " phi_min - phi =" << phi_min-norm_phi << " phi_max - phi =" << phi_max-norm_phi  << std::endl;
        //    if(-1.0*norm_phi<= roundTo1e5(sector_min_Phi[side][s]) && -1.0*norm_phi>=phi_min ){
        //      phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
        //      break;
        //    }          
        //  }
        //  if(i==0){
        //    //std::cout<< "i= "<< i <<": " << sector_max_Phi[side][s] << "<" << -1.0*norm_phi << "<" << phi_max << " phi_min - phi =" << phi_min-norm_phi << " phi_max - phi =" << phi_max-norm_phi  << std::endl;
        //    if(-1.0*norm_phi<=phi_max  && -1.0*norm_phi>=sector_max_Phi[side][s]  ){
        //      phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
        //      break;
        //    }          
        //  }
        //}   
      }
    }
    //std::cout<< sector_min_Phi[side][s] << " < " << norm_phi << " < " << sector_max_Phi[side][s] << std::endl;
    if (norm_phi < sector_max_Phi[side][s] && norm_phi >= sector_min_Phi[side][s])
    {
      //if(sector_number>-1 || phi_bin>-1) break;
      // NOLINTNEXTLINE(bugprone-integer-division)
      //sector_number = s;
      phi_bin = abs(nearest_element(norm_phi, s, side) - side*nphibins / 12) + nphibins / 12 * s;
      //for(long unsigned int i=0;i<layer_pad_phi.size();i++){
      //  double set_pad_phi = (sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 - pow(-1,side)*layer_pad_phi[i];
      //  if(i<layer_pad_phi.size()-1) dphi=layer_pad_phi[i+1]-layer_pad_phi[i];
      //  double phi_min = roundTo1e5(set_pad_phi-dphi/2);
      //  double phi_max = roundTo1e5(set_pad_phi+dphi/2);
      //  if(i==layer_pad_phi.size()-1){
      //    phi_min = set_pad_phi-10*dphi;
      //  }
      //  if(i==0){
      //    phi_max = set_pad_phi+10*dphi;
      //  }
      //  //std::cout<< "i= "<< i <<": " << phi_min << "<" << norm_phi << "<" << phi_max << " phi_min - phi =" << phi_min-norm_phi << " phi_max - phi =" << phi_max-norm_phi  << std::endl;
      //  if(norm_phi<=phi_max && norm_phi>=phi_min){
      //    phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
      //    //std::cout<< "i = " << i << " vs. nearest_element = " << nearest_element(norm_phi, s, side) << std::endl;
      //    break;
      //  }
      //
      //  //if(i==layer_pad_phi.size()-1){
      //  //  std::cout<< "i==layer_pad_phi.size()-1 " << set_pad_phi+dphi/2 << "<" << norm_phi << "<" << sector_min_Phi[side][s] << " phi_min - phi =" << phi_min-norm_phi << " phi_max - phi =" << phi_max-norm_phi  << std::endl;
      //  //  if(norm_phi<0){
      //  //    if(norm_phi<= roundTo1e5(sector_min_Phi[side][s]) && norm_phi>=phi_max){
      //  //      phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
      //  //      break;
      //  //    }      
      //  //  } else{
      //  //    if(norm_phi<= roundTo1e5(sector_max_Phi[side][s]) && norm_phi>=phi_max){
      //  //      phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
      //  //      break;
      //  //    }      
      //  //  }   
      //  //}
      // //if(i==0){
      // //   std::cout<< "i==0 " << sector_max_Phi[side][s] << "<" << norm_phi << "<" << set_pad_phi+dphi/2 << std::endl;
      // //   if(norm_phi>0){
      // //     if(norm_phi<=phi_max  && norm_phi>=roundTo1e5(sector_max_Phi[side][s])){
      // //       phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
      // //       break;
      // //     }
      // //   }else{
      // //     if(norm_phi<=phi_max  && norm_phi>=roundTo1e5(sector_min_Phi[side][s])){
      // //       phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s;
      // //       break;
      // //     }            
      // //   }          
      // // }
      //}

    }
    
    //if (phi_bin>-1) break;

  }
  //std::cout<< "PHG4TpcCylinderGeom::find_phibin 2" <<std::endl;
    //  if (norm_phi < sector_max_Phi[side][s] && norm_phi >= -M_PI)
    //  {
    //    // NOLINTNEXTLINE(bugprone-integer-division)
    //    phi_bin = floor(std::abs(sector_max_Phi[side][s] - norm_phi) / phistep) + nphibins / 12 * s;
    //    break;
    //  }
    //  if (norm_phi > sector_min_Phi[side][s] + 2 * M_PI)
    //  {
    //    // NOLINTNEXTLINE(bugprone-integer-division)
    //    phi_bin = floor(std::abs(sector_max_Phi[side][s] - (norm_phi - 2 * M_PI)) / phistep) + nphibins / 12 * s;
    //    break;
    //  }
    //}
  //std::cout<< "PHG4TpcCylinderGeom::find_phibin PhiBin = " << phi_bin <<std::endl;

  
  return phi_bin;
}

int PHG4TpcCylinderGeom::find_vecbin(const double phi, int side) const
{

  //std::cout<<"PHG4TpcCylinderGeom::find_phibin: phi = "<<phi<<std::endl;
  double dphi = layer_pad_phi[2]-layer_pad_phi[1];

  double norm_phi = phi;
  //if (phi < phimin || phi > (phimin + nphibins * phistep))
  //{
  //  int nwraparound = -floor((phi - phimin) * 0.5 / M_PI);
  //  norm_phi += 2 * M_PI * nwraparound;
  //}

  int phi_bin = -1;
  int sector_number = -1;

  //std::cout<<"PHG4TpcCylinderGeom::find_phibin: for loop for sectors "<<std::endl;
  for (std::size_t s = 0; s < sector_max_Phi[side].size(); s++)
  {
    if (norm_phi < sector_max_Phi[side][s] && norm_phi > sector_min_Phi[side][s])
    {
      if(sector_number>-1 || phi_bin>-1) break;
      // NOLINTNEXTLINE(bugprone-integer-division)
      sector_number = s;
      //phi_bin_old = (floor(std::abs(sector_max_Phi[side][s] - norm_phi) / phistep) + nphibins / 12 * s);
      for(long unsigned int i=0;i<layer_pad_phi.size();i++){
        double set_pad_phi = (sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 - pow(-1,side)*layer_pad_phi[i];
        if(i<layer_pad_phi.size()-1) dphi=layer_pad_phi[i+1]-layer_pad_phi[i];

        //std::cout<<"1 s = "<<s<<" :"<<" i = "<< i<< " layer_pad_phi[i] = " << layer_pad_phi[i] << ": " << set_pad_phi-dphi/2 << "<" << norm_phi<< "<" << set_pad_phi+dphi/2 <<std::endl;

        if(norm_phi<=set_pad_phi+dphi/2 && norm_phi>=set_pad_phi-dphi/2){
          //std::cout<<"1.1 s = "<<s<<" :"<<" i = "<< i<< " " << set_pad_phi-dphi/2 << "<" << norm_phi<< "<" << set_pad_phi+dphi/2 <<std::endl;
          //std::cout<<"i = "<< i <<std::endl;
          phi_bin = i;
          //std::cout<< "phi_bin = abs(int(i) - side*nphibins / 12) + nphibins / 12 * s : "<<
          //             phi_bin << " = abs("  << int(i)  << " - " << side << "*" << nphibins 
          //             << " / " << 12 << ") + " << nphibins << " / " << 12  << "*" <<  s << std::endl;
          //std::cout<<"0 phi_bin set to: "<<phi_bin<<std::endl;
          break;
        }
        if(i==layer_pad_phi.size()-1){
          if(norm_phi<= sector_min_Phi[side][s] && norm_phi>=set_pad_phi-dphi/2){
            phi_bin = i;
            break;
          }          
        }
       if(i==0){
          if(norm_phi<=set_pad_phi+dphi/2  && norm_phi>=sector_max_Phi[side][s]){
            //std::cout<<"7.1 s = "<< s <<" side = "<< side <<" :"<< sector_max_Phi[side][s] << "<" << norm_phi<< ">" << set_pad_phi+dphi/2 <<std::endl;
            //std::cout<<"7.2 s = "<< s <<" side = "<< side <<" :"<< sector_min_Phi[side][s] << "<" << norm_phi<< ">" << set_pad_phi+dphi/2 <<std::endl;
            phi_bin = i;
            //std::cout<<"2 phi_bin set to: "<<phi_bin<<std::endl;
            break;
          }          
        }
      }
      //std::cout<<"3 phi_bin set to: "<<phi_bin<<std::endl;

    }
    
    if (phi_bin>-1) break;

    if (s == 11 && phi_bin<0)
    {

      if(norm_phi > sector_min_Phi[side][s] + 2 * M_PI){  
        //std::cout<<"norm_phi = " << norm_phi <<"> sector_min_Phi[side][s] + 2 * M_PI = "   << sector_min_Phi[side][s] + 2 * M_PI<< " sector_min_Phi[side][s]="<< sector_min_Phi[side][s]<< " sector_max_Phi[side][s]="<<sector_max_Phi[side][s] <<std::endl;
        //phi_bin_old = floor(std::abs(sector_max_Phi[side][s] - (norm_phi - 2 * M_PI)) / phistep) + nphibins / 12 * s;

        for(long unsigned int i=0;i<layer_pad_phi.size();i++){
          double set_pad_phi = (sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 - pow(-1,side)*layer_pad_phi[i];
          //std::cout<< "nphibins*phistep*s = "<<nphibins<<"*"<<phistep<<"*"<<s<<std::endl;
          if(i<layer_pad_phi.size()-1) dphi=layer_pad_phi[i+1]-layer_pad_phi[i];

          //std::cout << set_pad_phi+dphi/2 << "<" << -1.0*norm_phi << "<" << set_pad_phi-dphi/2 << "dphi/2=" << dphi/2 << " set_pad_phi = "<<set_pad_phi<< std::endl;
      
          //if(-1.0*norm_phi<set_pad_phi+phistep/2 && -1.0*norm_phi>set_pad_phi-phistep/2){
          if(-1.0*norm_phi<set_pad_phi+dphi/2 && -1.0*norm_phi>=set_pad_phi-dphi/2){
            phi_bin = i;
            //std::cout <<"phi_bin = "<< phi_bin << ": " << set_pad_phi+dphi/2 << "<" << -1.0*norm_phi << "<" << set_pad_phi-dphi/2 << "dphi/2=" << dphi/2 << " set_pad_phi = "<<set_pad_phi<< std::endl;
            break;
          }
          if(i==layer_pad_phi.size()-1){
            //std::cout<<"4 s = "<<s<<" :"<< set_pad_phi-dphi/2 << "<" << -1.0*norm_phi<< "<" << sector_max_Phi[side][s] <<std::endl;
            if(-1.0*norm_phi<= sector_max_Phi[side][s] && -1.0*norm_phi>=set_pad_phi-dphi/2){
              //std::cout<<"5 s = "<<s<<" :"<< set_pad_phi-dphi/2 << "<" << norm_phi<< ">" << sector_max_Phi[side][s] <<std::endl;
              phi_bin = i;
              //std::cout<<"1 phi_bin set to: "<<phi_bin<<std::endl;
              break;
            }          
          }
          if(i==0){
            //std::cout<<"6 s = "<<s<<" :"<< sector_min_Phi[side][s] << "<" << -1.0*norm_phi<< "<" << set_pad_phi+dphi/2 <<std::endl;
            if(-1.0*norm_phi<=set_pad_phi+dphi/2  && -1.0*norm_phi>=sector_min_Phi[side][s]){
              //std::cout<<"7 s = "<<s<<" :"<< sector_min_Phi[side][s] << "<" << norm_phi<< ">" << set_pad_phi+dphi/2 <<std::endl;
              phi_bin = i;
              //std::cout<<"2 phi_bin set to: "<<phi_bin<<std::endl;
              break;
            }          
          }
        }
      }
      if (-1.0*norm_phi < sector_max_Phi[side][s] && -1.0*norm_phi >= -M_PI){
        //std::cout<<"norm_phi = " << norm_phi <<"< sector_max_Phi[side][s] = "   << sector_max_Phi[side][s] <<std::endl;
        //phi_bin_old = floor(std::abs(sector_max_Phi[side][s] - norm_phi) / phistep) + nphibins / 12 * s;
        for(long unsigned int i=0;i<layer_pad_phi.size();i++){
          double set_pad_phi = (sector_max_Phi[side][s]+sector_min_Phi[side][s])/2 + pow(-1,side)*layer_pad_phi[i];
          if(i<layer_pad_phi.size()-1) dphi=layer_pad_phi[i+1]-layer_pad_phi[i];
          //if(-1.0*norm_phi<set_pad_phi+phistep/2 && -1.0*norm_phi>set_pad_phi-phistep/2){
          if(-1.0*norm_phi<set_pad_phi+dphi/2 && -1.0*norm_phi>=set_pad_phi-dphi/2){
            phi_bin = i;
            break;
          }  
          if(i==layer_pad_phi.size()-1){
            //std::cout<<"14 s = "<<s<<" :"<< set_pad_phi-dphi/2 << "<" << norm_phi<< ">" << sector_max_Phi[side][s] <<std::endl;
            if(-1.0*norm_phi<= sector_max_Phi[side][s] && -1.0*norm_phi>=set_pad_phi-dphi/2){
              //std::cout<<"15 s = "<<s<<" :"<< set_pad_phi-dphi/2 << "<" << norm_phi<< ">" << sector_max_Phi[side][s] <<std::endl;
              phi_bin = i;
              //std::cout<<"11 phi_bin set to: "<<phi_bin<<std::endl;
              break;
            }          
          }
          if(i==0){
            //std::cout<<"16 s = "<<s<<" :"<< sector_min_Phi[side][s] << "<" << norm_phi<< ">" << set_pad_phi+dphi/2 <<std::endl;
            if(-1.0*norm_phi<=set_pad_phi+dphi/2  && -1.0*norm_phi>=sector_min_Phi[side][s]){
              //std::cout<<"17 s = "<<s<<" :"<< sector_min_Phi[side][s] << "<" << norm_phi<< ">" << set_pad_phi+dphi/2 <<std::endl;
              phi_bin = i;
              //std::cout<<"12 phi_bin set to: "<<phi_bin<<std::endl;
              break;
            }          
          }
        }   
      }
    }
    //  if (norm_phi < sector_max_Phi[side][s] && norm_phi >= -M_PI)
    //  {
    //    // NOLINTNEXTLINE(bugprone-integer-division)
    //    phi_bin = floor(std::abs(sector_max_Phi[side][s] - norm_phi) / phistep) + nphibins / 12 * s;
    //    break;
    //  }
    //  if (norm_phi > sector_min_Phi[side][s] + 2 * M_PI)
    //  {
    //    // NOLINTNEXTLINE(bugprone-integer-division)
    //    phi_bin = floor(std::abs(sector_max_Phi[side][s] - (norm_phi - 2 * M_PI)) / phistep) + nphibins / 12 * s;
    //    break;
    //  }
    

  }
  //std::cout<< "PHG4TpcCylinderGeom::get_phibin PhiBin = " << phi_bin <<std::endl;
  return phi_bin;
}

float PHG4TpcCylinderGeom::get_pad_float(const double phi, int side) const
{
  double norm_phi = phi;
  if (phi < phimin || phi > (phimin + nphibins * phistep))
  {
    int nwraparound = -floor((phi - phimin) * 0.5 / M_PI);
    norm_phi += 2 * M_PI * nwraparound;
  }
  // if (phi >  M_PI){
  //   norm_phi = phi - 2* M_PI;
  // }
  // if (phi < phimin){
  //   norm_phi = phi + 2* M_PI;
  // }
  side = 0;

  float phi_bin = -1;

  for (std::size_t s = 0; s < sector_max_Phi[side].size(); s++)
  {
    if (norm_phi < sector_max_Phi[side][s] && norm_phi > sector_min_Phi[side][s])
    {
      // NOLINTNEXTLINE(bugprone-integer-division)
      phi_bin = (std::abs(sector_max_Phi[side][s] - norm_phi) / phistep) + nphibins / 12 * s;
      break;
    }
    if (s == 11)
    {
      if (norm_phi < sector_max_Phi[side][s] && norm_phi >= -M_PI)
      {
        // NOLINTNEXTLINE(bugprone-integer-division)
        phi_bin = (std::abs(sector_max_Phi[side][s] - norm_phi) / phistep) + nphibins / 12 * s;
        break;
      }
      if (norm_phi > sector_min_Phi[side][s] + 2 * M_PI)
      {
        // NOLINTNEXTLINE(bugprone-integer-division)
        phi_bin = (std::abs(sector_max_Phi[side][s] - (norm_phi - 2 * M_PI)) / phistep) + nphibins / 12 * s;
        break;
      }
    }
  }
  return phi_bin - 0.5;
}

float PHG4TpcCylinderGeom::get_tbin_float(const double z) const
{
  if (z < zmin || z > (zmin + nzbins * zstep))
  {
    //    cout << PHWHERE << "Asking for bin for z outside of z range: " << z << endl;
    return -1;
  }

  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return ((z - zmin) / zstep) - 0.5;
}

int PHG4TpcCylinderGeom::get_phibin(const double phi, int side) const
{
  //std::cout<< "PHG4TpcCylinderGeom::get_phibin 0" <<std::endl;

  double new_phi = phi;

  if (new_phi > 2 * M_PI) {
    // Use fmod to get the remainder of angle divided by 2*pi
    new_phi = fmod(new_phi, 2 * M_PI);
    if(phi>2*M_PI){
      std::cout << "PHG4TpcCylinderGeom::get_phibin: Asking for bin for phi outside of phi range: " << phi << " cahnging to: " << new_phi << std::endl;
    }
  }
  // Adjust if the angle is outside the range [-pi, pi]
  if (new_phi > M_PI) {
    new_phi -= 2 * M_PI;
  } else if (new_phi < -M_PI) {
    new_phi += 2 * M_PI;
  }



  // Get phi-bin number
  int phi_bin = find_phibin(new_phi,side);
  //std::cout<< "PHG4TpcCylinderGeom::get_phibin 1" << phi_bin <<std::endl;

  //side = 0;
  // If phi-bin is not defined, check that it is in the dead area and put it to the edge of sector
  //std::cout<< "PHG4TpcCylinderGeom::get_phibin 2" <<std::endl;
  //if (phi_bin < 0)
  //{
  //  for (std::size_t s = 0; s < sector_max_Phi[side].size(); s++)
  //  {
  //    double daPhi = 0;
  //    if (s == 0)
  //    {
  //      daPhi = fabs(sector_min_Phi[side][11] + 2 * M_PI - sector_max_Phi[side][s]);
  //    }
  //    else
  //    {
  //      daPhi = fabs(sector_min_Phi[side][s - 1] - sector_max_Phi[side][s]);
  //    }

  //    double min_phi = sector_max_Phi[side][s];
  //    double max_phi = sector_max_Phi[side][s] + daPhi;
  //    if (new_phi <= max_phi && new_phi >= min_phi)
  //    {
  //      if (fabs(max_phi - new_phi) > fabs(new_phi - min_phi))
  //      {
  //        new_phi = min_phi - phistep / 5;
  //      }
  //      else
  //      {
  //        new_phi = max_phi + phistep / 5;
  //      }
  //    }
  //  }
  //  // exit(1);

  //  //phi_bin = find_phibin(new_phi,side);
  //  //std::cout<< "PHG4TpcCylinderGeom::get_phibin 3 "<< phi_bin <<std::endl;
  //  if (phi_bin < 0)
  //  {
  //    std::cout << PHWHERE << "PHG4TpcCylinderGeom::get_phibin: Asking for bin for phi outside of phi range: " << phi << " new_phi = " << new_phi << std::endl;
  //    //for (std::size_t s = 0; s < sector_max_Phi[side].size(); s++)
  //    //{
  //    //  std::cout<<sector_min_Phi[side][s]<< "<"<< phi <<"<"<< sector_max_Phi[side][s] <<" dphi -> max-phi = "<< sector_max_Phi[side][s] - phi << " min-phi = " <<sector_min_Phi[side][s] - phi <<std::endl;
  //    //}
  //    exit(1);
  //    // phi_bin=0;
  //  }
  //}

  if (phi_bin < 0)
  {
    std::cout << PHWHERE << "PHG4TpcCylinderGeom::get_phibin: Asking for bin for phi outside of phi range: " << phi << " new_phi = " << new_phi << std::endl;
    //for (std::size_t s = 0; s < sector_max_Phi[side].size(); s++)
    //{
    //  std::cout<<sector_min_Phi[side][s]<< "<"<< phi <<"<"<< sector_max_Phi[side][s] <<" dphi -> max-phi = "<< sector_max_Phi[side][s] - phi << " min-phi = " <<sector_min_Phi[side][s] - phi <<std::endl;
    //}
    exit(1);
    // phi_bin=0;
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
PHG4TpcCylinderGeom::get_phicenter(const int ibin, int side) const
{
  // double phi_center = -999;
  if (ibin < 0 || ibin > nphibins)
  {
    std::cout << PHWHERE << "Asking for invalid bin in phi: " << ibin << std::endl;
    exit(1);
  }

  check_binning_method_phi();

  unsigned int pads_per_sector = nphibins / 12;
  unsigned int sector = ibin / pads_per_sector;
  //double phi_center = (sector_max_Phi[side][sector] - (ibin + 0.5 - sector * pads_per_sector) * phistep);
  int vbin = ibin -  pads_per_sector * sector;
  if (vbin < 0 || vbin > int(pads_per_sector))
  {
    std::cout << PHWHERE << "Asking for invalid bin in layer_pad_phi: " << vbin << " for bin in phi: "<< ibin << std::endl;
    exit(1);
  }
  //double phi_center = (sector_max_Phi[side][sector]+sector_min_Phi[side][sector])/2 - pow(-1,side)*layer_pad_phi[vbin];
  double phi_center = (sector_max_Phi[side][sector]+sector_min_Phi[side][sector])/2 - layer_pad_phi[vbin];
  if (phi_center <= -M_PI)
  {
    phi_center += 2 * M_PI;
  }
  return phi_center;
}

double
PHG4TpcCylinderGeom::get_phi(const float ibin, int side) const
{

  // double phi_center = -999;
  if (ibin < 0 || ibin > nphibins)
  {
    std::cout << PHWHERE << "Asking for invalid bin in phi: " << ibin << std::endl;
    exit(1);
  }

  check_binning_method_phi();

  //const int side = 0;
  unsigned int pads_per_sector = nphibins / 12;
  unsigned int sector = ibin / pads_per_sector;
  //double phi_old = (sector_max_Phi[side][sector] - (ibin + 0.5 - sector * pads_per_sector) * phistep);
  int vbin = ibin -  pads_per_sector * sector;
  //int vbin1 = layer_pad_phi.size()-vbin;
  //double phi = (sector_max_Phi[side][sector]+sector_min_Phi[side][sector])/2 - pow(-1,side)*layer_pad_phi[vbin];
  double phi = (sector_max_Phi[side][sector]+sector_min_Phi[side][sector])/2 - layer_pad_phi[vbin];
  //for(size_t k=0; k<layer_pad_phi.size();k++){
  //  std::cout << k << "; "<< layer_pad_phi[k] << "| ";
  //}
  //std::cout <<std::endl;
  //std::cout << "PHG4TpcCylinderGeom::get_phi sector = " << sector << 
  //            " side = " << side <<
  //            "sector_max_Phi[side][sector] = " << sector_max_Phi[side][sector] <<
  //            "sector_min_Phi[side][sector]" << sector_min_Phi[side][sector] <<
  //            "(sector_max_Phi[side][sector]+sector_min_Phi[side][sector])/2" << (sector_max_Phi[side][sector]+sector_min_Phi[side][sector])/2 << std::endl;

  if (phi <= -M_PI)
  {
    phi += 2 * M_PI;
  }
  return phi;
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
    {
      std::cout << src << " : ";
    }

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
    {
      std::cout << src << " : ";
    }

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
