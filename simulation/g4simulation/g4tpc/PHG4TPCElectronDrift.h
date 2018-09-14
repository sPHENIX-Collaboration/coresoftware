#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameterInterface.h>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <vector>

class PHG4CellContainer;
class PHG4TPCPadPlane;
class PHCompositeNode;
class TH1;
class TNtuple;

class  PHG4TPCElectronDrift: public SubsysReco, public PHParameterInterface
{
public:
  PHG4TPCElectronDrift(const std::string& name = "PHG4TPCElectronDrift");
  virtual ~PHG4TPCElectronDrift();
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d) {detector = d;}
  std::string Detector() const {return detector;}
  void set_seed(const unsigned int iseed);
//  void Amplify(const double x, const double y, const double z);
  void MapToPadPlane(const double x, const double y, const double t);
  void registerPadPlane(PHG4TPCPadPlane *padplane);

private:
  PHG4CellContainer *g4cells;
  TH1 *dlong;
  TH1 *dtrans;
  TNtuple *nt;
  TNtuple *nthit;
  TNtuple *ntpad;
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  unsigned int seed;
  double diffusion_trans;
  double diffusion_long;
  double drift_velocity;
  double electrons_per_gev;
  double min_active_radius;
  double max_active_radius;
  double min_time;
  double max_time;
  std::vector<PHG4TPCPadPlane *> tpcpadplane;
#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif // G4TPC_PHG4TPCELECTRONDRIFT_H
