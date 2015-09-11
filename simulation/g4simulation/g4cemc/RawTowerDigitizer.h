#ifndef RawTowerDigitizer_H__
#define RawTowerDigitizer_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeom;

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

//! simple tower digitizer which sum all cell to produce photon yield and pedstal noises
class RawTowerDigitizer : public SubsysReco
{

public:
  RawTowerDigitizer(const std::string& name = "RawTowerDigitizer");
  virtual
  ~RawTowerDigitizer();

  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);
  void
  Detector(const std::string &d)
  {
    detector = d;
  }
  void
  EminCut(const double e)
  {
    emin = e;
  }

  void
  set_seed(const unsigned int iseed);
  unsigned int
  get_seed() const
  {
    return seed;
  }

protected:
  void
  CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _sim_towers;
  RawTowerContainer* _raw_towers;
  RawTowerGeom *rawtowergeom;

  std::string detector;
  std::string SimTowerNodeName;
  std::string RawTowerNodeName;
  std::string TowerGeomNodeName;

  int _cell_binning;
  double emin;

  PHTimeServer::timer _timer;

  unsigned int seed;
#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif /* RawTowerDigitizer_H__ */
