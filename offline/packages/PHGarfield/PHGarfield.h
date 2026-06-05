#ifndef PHGARFIELD__H
#define PHGARFIELD__H

#include <cmath>
#include <iostream>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class CdbUrlSave;
class CDBInterface;
class CDBTTree;
class PHCompositeNode;
class TpcMap;
class PHField3DCartesian;
class ComponentUser;
class TPolyLine3D;

namespace Garfield
{
  class ComponentUser;
  class MediumMagboltz;
}


class PHGarfield : public SubsysReco
{
 public:
  PHGarfield(const std::string &name = "PHGarfield");
  ~PHGarfield() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  bool StopHere(const double x, const double y, const double z, const double zPrevious);

  void PrintMaps();
  void PrintGarfield(double x, double y, double z);

  //  These are left in public namespace for easy plotting macros...
  //  The user is encouraged to add more routine to fit their analysis goals...
  TPolyLine3D *ReverseDrift (double x_cm, double y_cm, double z_cm, double step_ns=50.0); // Drifts electrons from some initial point until they hit a detector boundary...
  double radii[48];  // Radius on each layer just for test purposes...need to be cm!
  
 private:
  CDBInterface             *m_cdb            {nullptr}; // Access to all thiungs CDB...
  CDBTTree                 *m_cdbTPCMAPttree {nullptr}; // Locations of the pads from CDB...
  PHField3DCartesian       *m_field          {nullptr}; // The stanards sPHENIX field holding container.
  Garfield::ComponentUser  *m_component      {nullptr}; // This handles the interface of the electric and magnetic fields as handed to Garfield
  Garfield::MediumMagboltz *m_gas            {nullptr}; // This is the pre-tabulated gas properties required by Garfield...

  void GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double& bx_t,   double& by_t,   double& bz_t) ;   // Feeds magnetic field to Garfield
  void GetElectricFieldVcm  (double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm) ; // Feeds electric field to Garfield
  void InitializeGas        (std::string dir);
  void FillRadii();

  //  These are utilities for a spot check of the overall routine:
  //std::string calibdir;
  //std::string m_DiodeContainerName;
  double bounder(double phi, double phi_min);
  double PHI_MIN;

};

#endif
