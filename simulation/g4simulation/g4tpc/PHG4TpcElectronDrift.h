// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H

#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>

#include <cmath>
#include <memory>
#include <phparameter/PHParameterInterface.h>

#include <gsl/gsl_rng.h>
#include <string>                              // for string

class PHG4TpcPadPlane;
class PHCompositeNode;
class TH1;
class TH3;
class TNtuple;
class TFile;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcElectronDrift : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcElectronDrift(const std::string &name = "PHG4TpcElectronDrift");
  virtual ~PHG4TpcElectronDrift() = default;
  virtual int Init(PHCompositeNode*);
  virtual int InitRun(PHCompositeNode*);
  virtual int process_event(PHCompositeNode*);
  virtual int End(PHCompositeNode*);

  class DistortionStruct
  {
    
    public:
    using List = std::vector<DistortionStruct>;
    
    /// constructor
    DistortionStruct() = default;
    
    float _r = 0;
    float _phi = 0;
    float _z = 0;
    
    float _dr = 0;
    float _dphi = 0;
    float _dz = 0;
  };
 
  /// track container
  class Container: public PHObject
  {

    public:

    /// constructor
    explicit Container() = default;

    /// copy constructor
    explicit Container(const Container &) = delete;

    /// assignment operator
    Container& operator = ( const Container& ) = delete;

    /// reset
    virtual void Reset();

    /// distrotions
    const DistortionStruct::List& distortions() const
    { return _distortions; }
    
    /// add distortion
    void addDistortion( const DistortionStruct& distortion )
    { _distortions.push_back( distortion ); }

    private:

    /// event struct
    DistortionStruct::List _distortions;

    ClassDef(Container,1)

  };
  
  void SetDefaultParameters();

  //! detector name
  void Detector(const std::string &d) 
  { detector = d; }

  //! detector name
  std::string Detector() const 
  { return detector; }
  
  //! random seed
  void set_seed(const unsigned int iseed);
  
  //! space charge distortions
  void set_enable_distortions( bool value )
  { m_enable_distortions = value; }
  
  //! distortion filename
  void set_distortion_filename( const std::string& value )
  { m_distortion_filename = value; }
  
  //! setup readout plane
  void registerPadPlane(PHG4TpcPadPlane *padplane);

 private:
  
  //! map a given x,y,z coordinates to plane hits
  void MapToPadPlane(const double x, const double y, const double z, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);

  TrkrHitSetContainer *hitsetcontainer = nullptr;
  TrkrHitTruthAssoc *hittruthassoc = nullptr;
  std::unique_ptr<TrkrHitSetContainer> temp_hitsetcontainer;
  std::unique_ptr<PHG4TpcPadPlane> padplane;

  //! evaluation node
  Container* m_container = nullptr;

  //! space charge distortion file name
  bool m_enable_distortions = false;
  std::string m_distortion_filename;
  TFile *m_distortion_tfile = nullptr;

  //!@name space charge distortion histograms
  //@{
  TH3 *hDRint = nullptr;
  TH3 *hDPint = nullptr;
  TH3 *hDZint = nullptr;
  //@}
  
  TH1 *dlong = nullptr;
  TH1 *dtrans = nullptr;
  TFile *m_outf = nullptr;
  TNtuple *nt = nullptr;
  TNtuple *nthit = nullptr;
  TNtuple *ntfinalhit = nullptr;
  TNtuple *ntpad = nullptr;
  std::string detector;
  std::string hitnodename;
  std::string seggeonodename;
  double diffusion_trans = NAN;
  double added_smear_sigma_trans = NAN;
  double diffusion_long = NAN;
  double added_smear_sigma_long = NAN;
  double drift_velocity = NAN;
  double tpc_length = NAN;
  double electrons_per_gev = NAN;
  double min_active_radius = NAN;
  double max_active_radius = NAN;
  double min_time = NAN;
  double max_time = NAN;

  //! rng de-allocator
  class Deleter
  {
    public:
    //! deletion operator
    void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };
  std::unique_ptr<gsl_rng, Deleter> RandomGenerator;
  
};

#endif  // G4TPC_PHG4TPCELECTRONDRIFT_H
