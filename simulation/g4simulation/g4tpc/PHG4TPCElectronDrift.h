#ifndef PHG4TPCElectronDrift_h
#define PHG4TPCElectronDrift_h

#include <fun4all/SubsysReco.h>

#include <g4detectors/PHG4ParameterInterface.h>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

class PHCompositeNode;

class  PHG4TPCElectronDrift: public SubsysReco, public PHG4ParameterInterface
{
public:
  PHG4TPCElectronDrift(const std::string& name = "PHG4TPCElectronDrift");
  virtual ~PHG4TPCElectronDrift() {}
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d) {detector = d;}
  std::string Detector() const {return detector;}
  void set_seed(const unsigned int iseed);
private:
  std::string detector;
  std::string hitnodename;
  unsigned int seed;
  double diffusion_trans;
  double diffusion_long;
  double drift_velocity;
  double electrons_per_gev;
#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif // PHG4TPCElectronDrift_h
