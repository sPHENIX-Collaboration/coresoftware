#include <phool/phool.h>
#include <TFile.h>
#include <TNtuple.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include <PHBFieldMapObj.h>

using namespace std;

const double PHBFieldMapObj::CELL_WIDTH_R = 2.0;    // cm
const double PHBFieldMapObj::CELL_WIDTH_PHI = 180.0;  // degrees
const double PHBFieldMapObj::CELL_WIDTH_Z = 2.0;    // cm

long PHBFieldMapObj::MIN_R = LONG_MAX;
long PHBFieldMapObj::MAX_R = LONG_MIN;
long PHBFieldMapObj::MIN_PHI = LONG_MAX;
long PHBFieldMapObj::MAX_PHI = LONG_MIN;
long PHBFieldMapObj::MIN_Z = LONG_MAX;
long PHBFieldMapObj::MAX_Z = LONG_MIN;

const double PHBFieldMapObj::RAD_TO_DEG = 180.0/M_PI;
const double PHBFieldMapObj::DEG_TO_RAD = M_PI/180.0;

template <typename T>
static inline T
sqr (const T& x)
{
  return x * x;
}

PHBFieldMapObj::~PHBFieldMapObj()
{
  // clear b-field lookup table
  for (unsigned int i=0; i<_lookup.size(); i++)
    {
      if ( _lookup[i]!=null_fieldpoint() ) delete _lookup[i];
    }
}

void 
PHBFieldMapObj::initialize()
{
  cout << PHWHERE << " Magnetic field lookup initialization" << endl;
  
  TFile *f = 0; 
  
  // check if the map file exists
  ifstream fp(_filename.c_str(),ios::in);

  if(!fp || fp.bad()){

    string fpath("/phenix/upgrades/decadal/fieldmaps/");

    cerr << PHWHERE
	 << ": Unable to find local file " << _filename << endl;
    cerr << PHWHERE
	 << ": Trying " << fpath << " area...";

    ifstream fp((fpath+_filename).c_str(),ios::in);
    if(!fp || fp.bad()){
      cerr << endl;
      cerr << PHWHERE
	   << ": ERROR - Unable to find field map file " << _filename << endl;
      return;
    }
    else{
      cerr << "OK." << endl;
    }

    fp.close();

    f = new TFile((fpath+_filename).c_str());

  }  
  else{
    cerr << PHWHERE
	 << ": Using local copy of " << _filename << endl;
  }
  
  if(!f) {
    fp.close();
    f = new TFile(_filename.c_str());
  } 
  
  TNtuple* nt = (TNtuple*)gDirectory->Get("fieldmap");
  if(!nt){
    throw std::runtime_error(DESCRIPTION("can't find input file ntuple: fieldmap"));
  }

  // Setup locals for nt access
  //
  Float_t z=0;
  Float_t r=0;
  Float_t phi=0;
  Float_t bz=0;
  Float_t br=0;
  Float_t bphi=0;
  
  nt->SetBranchAddress("z",&z);
  nt->SetBranchAddress("r",&r);
  nt->SetBranchAddress("bz",&bz);
  nt->SetBranchAddress("br",&br);

  bool is_3d = nt->FindBranch("phi") && nt->FindBranch("bphi");
  if (is_3d){
  nt->SetBranchAddress("phi",&phi);
  nt->SetBranchAddress("bphi",&bphi);
  }  

  // determine the bounds of the lookup
  Stat_t nentries = nt->GetEntries();
  for (int j=0; j<nentries;j++) {
    nt->GetEntry(j);   
    if(is_3d){ update_bounds(r,phi,z);
    } else {
     update_bounds(r,180.,z);
     update_bounds(r,0.,z);
    }
  }

  // initialize lookup table
  cout << PHWHERE << " initializing lookup" << endl;
  initialize_lookup();
  
  // set the root field point read from field map
  for (int j=0; j<nentries;j++) {
    nt->GetEntry(j);   
    if(is_3d) add_point(r,phi,z,br,bphi,bz);    
    else {
     add_point(r,180.,z,br,0.,bz);
     add_point(r,0.,z,br,0.,bz);
    }
  }
  _initialized = true;

  // fill out the rest of the voxel if voxel mode
  if(_voxel_mode){
    cout << PHWHERE << " initializing voxel lookup" << endl;
    for (int jentry=0; jentry<nentries;jentry++) {
      nt->GetEntry(jentry);   
      add_neighbours(r,phi,z);    
    }
    _voxel_initialized = true;
  }
  f->Close();
  cout << PHWHERE << " magnetic field initialization done" << endl;
}

void 
PHBFieldMapObj::print_voxel(const double& x, const double& y, const double& z) 
{
  // convert input to cylindrical coords
  double r = std::sqrt(sqr(x) + sqr(y));
  double phi = std::atan2(y,x);
  
  // translate [-pi,pi] to [0,2pi] and convert to degrees
  phi = phi < 0 ? PHBFieldMapObj::RAD_TO_DEG*(phi+2*M_PI) : PHBFieldMapObj::RAD_TO_DEG*phi; 
  
  // print all the field points in a given voxel
  size_t key = hash_function(r,phi,z);  
  _lookup[key]->print();
}

const PHBFieldMapObj::BField& 
PHBFieldMapObj::get_field(const double& x, const double& y, const double& z) 
{

  // convert input to cylindrical coords
  double r = std::sqrt(sqr(x) + sqr(y));
  double phi = std::atan2(y,x);
  
  // translate [-pi,pi] to [0,2pi] and convert to degrees
  phi = phi < 0 ? PHBFieldMapObj::RAD_TO_DEG*(phi+2*M_PI) : PHBFieldMapObj::RAD_TO_DEG*phi; 
  
  // cache FieldPoints if previous lookup is in the same cell
  static FieldPoint *r_lo_phi_lo_z_lo = 0;
  static FieldPoint *r_lo_phi_lo_z_hi = 0;
  static FieldPoint *r_lo_phi_hi_z_lo = 0;
  static FieldPoint *r_lo_phi_hi_z_hi = 0;
  static FieldPoint *r_hi_phi_lo_z_lo = 0;
  static FieldPoint *r_hi_phi_lo_z_hi = 0;
  static FieldPoint *r_hi_phi_hi_z_lo = 0;
  static FieldPoint *r_hi_phi_hi_z_hi = 0;

  // prev_key is how we determine whether we are in the same cell
  static size_t prev_key = _lookup.size();

  size_t key = hash_function(r, phi, z);

  // only do lookup if we are in a different b-field cell
  if ( key!=prev_key )
    {
      prev_key = key;

      if (key>=_lookup.size()) {
        prev_key = _lookup.size();
        return null_field();
      }
      r_lo_phi_lo_z_lo = _lookup[key];

      key = hash_function(r, phi, z, 0, 0, 1);
      if (key>=_lookup.size()) return null_field();
      r_lo_phi_lo_z_hi = _lookup[key];
 
      key = hash_function(r, phi, z, 0, 1, 0);
      if (key>=_lookup.size()) return null_field();
      r_lo_phi_hi_z_lo = _lookup[key];
 
      key = hash_function(r, phi, z, 0, 0, 1);
      if (key>=_lookup.size()) return null_field();
      r_lo_phi_hi_z_hi = _lookup[key];
  
      key = hash_function(r, phi, z, 1, 0, 0);
      if (key>=_lookup.size()) return null_field();  
      r_hi_phi_lo_z_lo = _lookup[key];
  
      key = hash_function(r, phi, z, 1, 0, 1);
      if (key>=_lookup.size()) return null_field();
      r_hi_phi_lo_z_hi = _lookup[key];
  
      key = hash_function(r, phi, z, 1, 1, 0);
      if (key>=_lookup.size()) return null_field();
      r_hi_phi_hi_z_lo = _lookup[key];
  
      key = hash_function(r, phi, z, 1, 1, 1);
      if (key>=_lookup.size()) return null_field();
      r_hi_phi_hi_z_hi = _lookup[key];

    }

  double t = (r - r_lo_phi_lo_z_lo->r)/CELL_WIDTH_R;
  double u = (phi - r_lo_phi_lo_z_lo->phi)/CELL_WIDTH_PHI;
  double v = (z - r_lo_phi_lo_z_lo->z)/CELL_WIDTH_Z;
  
  // bx interpolation
  {
    // bilinear interpolation in z phi plane at r_lo
    double z1 = (1-v) * r_lo_phi_lo_z_lo->bx + v * r_lo_phi_lo_z_hi->bx;
    double z2 = (1-v) * r_lo_phi_hi_z_lo->bx + v * r_lo_phi_hi_z_hi->bx;              
    double z_phi_1 = (1-u) * z1 + u * z2;
    
    // bilinear interpolation in z phi plane at r_hi
    double z3 = (1-v) * r_hi_phi_lo_z_lo->bx + v * r_hi_phi_lo_z_hi->bx;
    double z4 = (1-v) * r_hi_phi_hi_z_lo->bx + v * r_hi_phi_hi_z_hi->bx;              
    double z_phi_2 = (1-u) * z3 + u * z4;
    
    // now interpolate between the two r planes 
    _field.bx = (1-t) * z_phi_1 + t * z_phi_2;
  }
  
  // by interpolation
  {
    double z1 = (1-v) * r_lo_phi_lo_z_lo->by + v * r_lo_phi_lo_z_hi->by;
    double z2 = (1-v) * r_lo_phi_hi_z_lo->by + v * r_lo_phi_hi_z_hi->by;              
    double z_phi_1 = (1-u) * z1 + u * z2;
    
    double z3 = (1-v) * r_hi_phi_lo_z_lo->by + v * r_hi_phi_lo_z_hi->by;
    double z4 = (1-v) * r_hi_phi_hi_z_lo->by + v * r_hi_phi_hi_z_hi->by;              
    double z_phi_2 = (1-u) * z3 + u * z4;
    
    _field.by = (1-t) * z_phi_1 + t * z_phi_2;
  }
  
  // bz interpolation
  {
    double z1 = (1-v) * r_lo_phi_lo_z_lo->bz + v * r_lo_phi_lo_z_hi->bz;
    double z2 = (1-v) * r_lo_phi_hi_z_lo->bz + v * r_lo_phi_hi_z_hi->bz;              
    double z_phi_1 = (1-u) * z1 + u * z2;
    
    double z3 = (1-v) * r_hi_phi_lo_z_lo->bz + v * r_hi_phi_lo_z_hi->bz;
    double z4 = (1-v) * r_hi_phi_hi_z_lo->bz + v * r_hi_phi_hi_z_hi->bz;              
    double z_phi_2 = (1-u) * z3 + u * z4;
    
    _field.bz = (1-t) * z_phi_1 + t * z_phi_2;
  }

  _field.error = false; 

  return _field;
}

const PHBFieldMapObj::BField& 
PHBFieldMapObj::get_field_voxel(const double& x, const double& y, const double& z) 
{
#ifdef VOXEL_MODE

  double r = std::sqrt(sqr(x) + sqr(y));
  double phi = std::atan2(y,x);
  
  phi = phi < 0 ? PHBFieldMapObj::RAD_TO_DEG*(phi+2*M_PI) : PHBFieldMapObj::RAD_TO_DEG*phi; 
  
  long key = hash_function(r, phi, z);

  if(key>=_lookup.size()) return null_field();

  const FieldPoint& fp = _lookup[key];
    
  double t = (r - fp.r)/CELL_WIDTH_R;
  double u = (phi - fp.phi)/CELL_WIDTH_PHI;
  double v = (z - fp.z)/CELL_WIDTH_Z;
  
  // bx interpolation
  {
    // Bilinear interpolation in z phi plane at r_lo
    double z1 = (1-v) * fp.bx + v * fp.bx_rlo_phihi_zhi;
    double z2 = (1-v) * fp.bx_rlo_phihi_zlo + v * fp.bx_rlo_phihi_zhi;              
    double z_phi_1 = (1-u) * z1 + u * z2;
    
    // Bilinear interpolation in z phi plane at r_hi
    double z3 = (1-v) * fp.bx_rhi_philo_zlo + v * fp.bx_rhi_philo_zhi;
    double z4 = (1-v) * fp.bx_rhi_phihi_zlo + v * fp.bx_rhi_phihi_zhi;              
    double z_phi_2 = (1-u) * z3 + u * z4;
    
    // Now interpolate between the two r planes 
    _field.bx = (1-t) * z_phi_1 + t * z_phi_2;
  }
  
  // by interpolation
  {
    double z1 = (1-v) * fp.by + v * fp.by_rlo_philo_zhi;
    double z2 = (1-v) * fp.by_rlo_phihi_zlo + v * fp.by_rlo_phihi_zhi;              
    double z_phi_1 = (1-u) * z1 + u * z2;
    
    double z3 = (1-v) * fp.by_rhi_philo_zlo + v * fp.by_rhi_philo_zhi;
    double z4 = (1-v) * fp.by_rhi_phihi_zlo + v * fp.by_rhi_phihi_zhi;              
    double z_phi_2 = (1-u) * z3 + u * z4;
    
    _field.by = (1-t) * z_phi_1 + t * z_phi_2;
  }
  
  // bz interpolation
  {
    double z1 = (1-v) * fp.bz + v * fp.bz_rlo_philo_zhi;
    double z2 = (1-v) * fp.bz_rlo_phihi_zlo + v * fp.bz_rlo_phihi_zhi;              
    double z_phi_1 = (1-u) * z1 + u * z2;
    
    double z3 = (1-v) * fp.bz_rhi_philo_zlo + v * fp.bz_rhi_philo_zhi;
    double z4 = (1-v) * fp.bz_rhi_phihi_zlo + v * fp.bz_rhi_phihi_zhi;              
    double z_phi_2 = (1-u) * z3 + u * z4;
    
    _field.bz = (1-t) * z_phi_1 + t * z_phi_2;
  }

#endif
#ifndef VOXEL_MODE
  throw std::logic_error(DESCRIPTION("PHBFieldMapObj not compiled for voxel lookup"));
#endif 
  return _field;
}

void 
PHBFieldMapObj::add_point(const double r, const double phi, const double z,
		       const double br, const double bphi, const double bz)
{
  // Convert to cartesian
  double sin_phi = std::sin(PHBFieldMapObj::DEG_TO_RAD*phi);
  double cos_phi = std::cos(PHBFieldMapObj::DEG_TO_RAD*phi);
  
  double bx = br*cos_phi - bphi*sin_phi;
  double by = br*sin_phi + bphi*cos_phi;

  FieldPoint *point = new FieldPoint();
  point->r = r;
  point->phi = phi;
  point->z = z;
  point->bx = bx;
  point->by = by;
  point->bz = bz;

  // Use at() during initialization 
  size_t key = hash_function(r,phi,z);
  if ( key<0 || key>=_lookup.size() )
    {
      cout << PHWHERE << "Invalid key to Field Lookup" << endl;
      return;
    }

  _lookup[hash_function(r,phi,z)] = point;
}

void 
PHBFieldMapObj::add_neighbours(const double& r, const double& phi, const double& z)
{
#ifdef VOXEL_MODE
  // fill in the 7 other field points associated with a given voxel.  The hash_function is
  // overloaded to handle the requested offsets.
  _lookup[hash_function(r,phi,z)]->bx_rlo_phihi_zlo = _lookup[hash_function(r,phi,z,0,1,0)]->bx;
  _lookup[hash_function(r,phi,z)]->by_rlo_phihi_zlo = _lookup[hash_function(r,phi,z,0,1,0)]->by;
  _lookup[hash_function(r,phi,z)]->bz_rlo_phihi_zlo = _lookup[hash_function(r,phi,z,0,1,0)]->bz;

  _lookup[hash_function(r,phi,z)]->bx_rlo_philo_zhi = _lookup[hash_function(r,phi,z,0,0,1)]->bx;
  _lookup[hash_function(r,phi,z)]->by_rlo_philo_zhi = _lookup[hash_function(r,phi,z,0,0,1)]->by;
  _lookup[hash_function(r,phi,z)]->bz_rlo_philo_zhi = _lookup[hash_function(r,phi,z,0,0,1)]->bz;

  _lookup[hash_function(r,phi,z)]->bx_rlo_phihi_zhi = _lookup[hash_function(r,phi,z,0,1,1)]->bx;
  _lookup[hash_function(r,phi,z)]->by_rlo_phihi_zhi = _lookup[hash_function(r,phi,z,0,1,1)]->by;
  _lookup[hash_function(r,phi,z)]->bz_rlo_phihi_zhi = _lookup[hash_function(r,phi,z,0,1,1)]->bz;

  _lookup[hash_function(r,phi,z)]->bx_rhi_philo_zlo = _lookup[hash_function(r,phi,z,1,0,0)]->bx;
  _lookup[hash_function(r,phi,z)]->by_rhi_philo_zlo = _lookup[hash_function(r,phi,z,1,0,0)]->by;
  _lookup[hash_function(r,phi,z)]->bz_rhi_philo_zlo = _lookup[hash_function(r,phi,z,1,0,0)]->bz;

  _lookup[hash_function(r,phi,z)]->bx_rhi_phihi_zlo = _lookup[hash_function(r,phi,z,1,1,0)]->bx;
  _lookup[hash_function(r,phi,z)]->by_rhi_phihi_zlo = _lookup[hash_function(r,phi,z,1,1,0)]->by;
  _lookup[hash_function(r,phi,z)]->bz_rhi_phihi_zlo = _lookup[hash_function(r,phi,z,1,1,0)]->bz;

  _lookup[hash_function(r,phi,z)]->bx_rhi_philo_zhi = _lookup[hash_function(r,phi,z,1,0,1)]->bx;
  _lookup[hash_function(r,phi,z)]->by_rhi_philo_zhi = _lookup[hash_function(r,phi,z,1,0,1)]->by;
  _lookup[hash_function(r,phi,z)]->bz_rhi_philo_zhi = _lookup[hash_function(r,phi,z,1,0,1)]->bz;

  _lookup[hash_function(r,phi,z)]->bx_rhi_phihi_zhi = _lookup[hash_function(r,phi,z,1,1,1)]->bx;
  _lookup[hash_function(r,phi,z)]->by_rhi_phihi_zhi = _lookup[hash_function(r,phi,z,1,1,1)]->by;
  _lookup[hash_function(r,phi,z)]->bz_rhi_phihi_zhi = _lookup[hash_function(r,phi,z,1,1,1)]->bz;
#endif
}

void 
PHBFieldMapObj::update_bounds(const double& r, const double& phi, const double& z)
{
  long i_r = static_cast<int>(std::floor(r/CELL_WIDTH_R));
  long i_phi = static_cast<int>(std::floor(phi/CELL_WIDTH_PHI));
  long i_z = static_cast<int>(std::floor(z/CELL_WIDTH_Z));

  MIN_R = std::min(i_r,MIN_R);
  MIN_PHI = std::min(i_phi,MIN_PHI);
  MIN_Z = std::min(i_z,MIN_Z);
  MAX_R = std::max(i_r-MIN_R,MAX_R);
  MAX_PHI = std::max(i_phi-MIN_PHI,MAX_PHI);
  MAX_Z = std::max(i_z,MAX_Z);

}

void 
PHBFieldMapObj::initialize_lookup()
{
  // Offset 2 is (1 for distance calc + 1 for empty border)
  _size_r = (MAX_R-MIN_R + 2); 
  _size_phi = (MAX_PHI-MIN_PHI + 2); 
  _size_z = (MAX_Z-MIN_Z + 2); 

  // allocate memory for pointers (set to null_field initially so that 
  // field outside mapped region is automatically set to 0)
  size_t n_cells = _size_r*_size_phi*_size_z;
  _lookup.reserve( n_cells );
  for (size_t i=0; i<n_cells; i++)
    {
      _lookup.push_back( null_fieldpoint() );
    }
}

inline size_t PHBFieldMapObj::hash_function(const double& r, const double& phi, const double& z,
					 size_t r_offset, size_t phi_offset, size_t z_offset)  
{  
  int i_r = static_cast<int>(std::floor(r/CELL_WIDTH_R)) - MIN_R;
  
  int i_phi = static_cast<int>(std::floor(phi/CELL_WIDTH_PHI)) - MIN_PHI;
  
  int i_z = static_cast<int>(std::floor(z/CELL_WIDTH_Z)) - MIN_Z;

  if( (i_r<=MAX_R) && (i_phi<=MAX_PHI) && (i_z<=MAX_Z) && (i_r>=0) && (i_phi>=0) && (i_z>=0))
    return (i_z+z_offset)*_size_r*_size_phi + (i_phi+phi_offset)*_size_r + i_r+r_offset;
  else
    return (_lookup.size() + 1);

}

inline size_t PHBFieldMapObj::hash_function(const double& r, const double& phi, const double& z)
{  
  int i_r = static_cast<int>(std::floor(r/CELL_WIDTH_R)) - MIN_R;
  
  int i_phi = static_cast<int>(std::floor(phi/CELL_WIDTH_PHI)) - MIN_PHI;
  
  int i_z = static_cast<int>(std::floor(z/CELL_WIDTH_Z)) - MIN_Z;
  
  if( (i_r<=(MAX_R-MIN_R)) && (i_phi<=(MAX_PHI-MIN_PHI)) && (i_z<=(MAX_Z-MIN_Z)) && (i_r>=0) && (i_phi>=0) && (i_z>=0))
    return i_z*_size_r*_size_phi + i_phi*_size_r + i_r;
  else
    return (_lookup.size() + 1);
  
}

