#ifndef PHG4TPCPadPlane_h
#define PHG4TPCPadPlane_h

#include <phparameter/PHParameterInterface.h>
#include <fun4all/SubsysReco.h>

class PHG4CellContainer;
class PHCompositeNode;

class PHG4TPCPadPlane: public SubsysReco, public PHParameterInterface
{
public:
  PHG4TPCPadPlane(const std::string &name="PHG4TPCPadPlane"); 

virtual ~PHG4TPCPadPlane(){}

#ifndef __CINT__
  int process_event(PHCompositeNode *) final {return 0;}
#else
  int process_event(PHCompositeNode *) {return 0;}
#endif
  int InitRun(PHCompositeNode *topNode);
  virtual void UpdateInternalParameters() {return;}
  virtual void MapToPadPlane(PHG4CellContainer *g4cells, const double x_gem, const double y_gem, const double t_gem) {}
  void Detector(const std::string &name) {detector = name;}
protected:

  std::string detector;
};

#endif
