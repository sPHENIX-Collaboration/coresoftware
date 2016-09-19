
#ifndef __PHBFIELDMAPOBJ_H__
#define __PHBFIELDMAPOBJ_H__

#include <TDataType.h>
#include <PHException.h>
#include <vector>
#include <iomanip>
#include <climits>

#include <boost/smart_ptr.hpp>
#include <boost/array.hpp>

//! 3D B Field lookup with linear interpolation.

/*!
  Upon construction PHBFieldMapObj reads the sPHENIX 2D magnetic field and populates
  an lookup table. The input B field map is a ROOT format file that is specified
  at construction time.  If voxel_mode (construction parameter) is true (the default)
  then lookup table looks up and caches the 8 field points that make up a voxel.  
  This mode takes longer to initialize but reduces subsequent lookup times by about
  a factor of five.
  <br>
  <br>
  The input grid pitch is set to delta r = 2cm, delta z = 2cm and delta phi = 180 degrees.
  The phi pitch is defined and given an artificial value to accomodate full 3D field map.
  If the input files is of a different grid the code will not work -- ie the input
  grid spacing is hardwired not determined at runtime from the input file.
  <br>
*/

class PHBFieldMapObj
{
 public: 

  PHBFieldMapObj(const std::string& filename, bool voxel_mode=false) : 
    _filename(filename),
    _initialized(false),
    _voxel_mode(voxel_mode),
    _voxel_initialized(false)
    {
      initialize();
    }
  
  virtual ~PHBFieldMapObj();

  struct BField {        
    BField() : bx(0),by(0),bz(0),error(true){;}
    double bx;
    double by;
    double bz;
    bool error;    
  };
  

  const BField& get_field(const double& x, const double& y, const double& z);
  const BField& get_field_voxel(const double& x, const double& y, const double& z);
  
  void print_voxel(const double& x, const double& y, const double& z);
  
 private:
  
  size_t hash_function(const double& r, const double& phi, const double& z);
  size_t hash_function(const double& r, const double& phi, const double& z,
		       size_t r_offset, size_t phi_offset, size_t z_offset);

  void initialize();  
  void update_bounds(const double& r, const double& phi, const double& z);
  void initialize_lookup();
  void add_point(const double r, const double phi, const double z,
		 const double br, const double bphi, const double bz);
  void add_neighbours(const double& r, const double& phi, const double& z);


  struct FieldPoint 
  {    
    float r;
    float phi;
    float z;

    float bx;
    float by;
    float bz;

#ifdef VOXEL_MODE

    float bx_rlo_phihi_zlo;
    float by_rlo_phihi_zlo;
    float bz_rlo_phihi_zlo;

    float bx_rlo_philo_zhi;
    float by_rlo_philo_zhi;
    float bz_rlo_philo_zhi;

    float bx_rlo_phihi_zhi;
    float by_rlo_phihi_zhi;
    float bz_rlo_phihi_zhi;

    float bx_rhi_philo_zlo;
    float by_rhi_philo_zlo;
    float bz_rhi_philo_zlo;

    float bx_rhi_phihi_zlo;
    float by_rhi_phihi_zlo;
    float bz_rhi_phihi_zlo;

    float bx_rhi_philo_zhi;
    float by_rhi_philo_zhi;
    float bz_rhi_philo_zhi;

    float bx_rhi_phihi_zhi;
    float by_rhi_phihi_zhi;
    float bz_rhi_phihi_zhi;

#endif

    void print(std::ostream& os=std::cout) const {
      // for conversion to cylindrical 
      os << std::setw(10) << std::setprecision(5) << std::setiosflags(std::ios::showpoint);

      // print point cylindrical coordinates
      os << "x_cyln = {" << r << ", " << phi << ", " << z << "}" << std::endl;

      // print the voxel
      os << "B 0-0-0 = {" << bx << ", " << by << ", " << bz << "}" << std::endl;

#ifdef VOXEL_MODE

      os << "B 0-1-0 = {" << bx_rlo_phihi_zlo << ", " << by_rlo_phihi_zlo <<  	 
	", " << bz_rlo_phihi_zlo << "}" << std::endl;

      os << "B 0-0-1 = {" << bx_rlo_philo_zhi << ", " << by_rlo_philo_zhi <<  	 
	", " << bz_rlo_philo_zhi << "}" << std::endl;

      os << "B 0-1-1 = {" << bx_rlo_phihi_zhi << ", " << by_rlo_phihi_zhi <<  	 
	", " << bz_rlo_phihi_zhi << "}" << std::endl;
       
      os << "B 1-0-0 = {" << bx_rhi_philo_zlo << ", " << by_rhi_philo_zlo <<  	 
	", " << bz_rhi_philo_zlo << "}" << std::endl;

      os << "B 1-1-0 = {" << bx_rhi_phihi_zlo << ", " << by_rhi_phihi_zlo <<  	 
	", " << bz_rhi_phihi_zlo << "}" << std::endl;

      os << "B 1-0-1 = {" << bx_rhi_philo_zhi << ", " << by_rhi_philo_zhi <<  	 
	", " << bz_rhi_philo_zhi << "}" << std::endl;

      os << "B 1-1-1 = {" << bx_rhi_phihi_zhi << ", " << by_rhi_phihi_zhi <<  	 
	", " << bz_rhi_phihi_zhi << "}" << std::endl;       

#endif

    }
  };
  
  static const double CELL_WIDTH_R;
  static const double CELL_WIDTH_PHI;
  static const double CELL_WIDTH_Z;

  static long MIN_R;
  static long MIN_PHI;
  static long MIN_Z;
  static long MAX_R;
  static long MAX_PHI;
  static long MAX_Z;

  static const double RAD_TO_DEG;
  static const double DEG_TO_RAD;

  size_t _size_r;
  size_t _size_phi;
  size_t _size_z;

  std::vector<FieldPoint*> _lookup;
  
  BField _field;    
  
  const BField& null_field(){
    static BField null;
    return null;
  }

  FieldPoint *null_fieldpoint() {
    static FieldPoint null;
    return &null;
  }

  std::string _filename;
  
  bool _initialized;
  bool _voxel_mode;
  bool _voxel_initialized;
};

#endif

