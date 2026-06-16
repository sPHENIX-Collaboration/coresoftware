#ifndef PHGARFIELD__H
#define PHGARFIELD__H

#include <fun4all/SubsysReco.h>

#include <array>
#include <numbers>
#include <string>

class CDBTTree;
class PHField3DCartesian;
class TPolyLine3D;

namespace Garfield
{
  class ComponentUser;
  class MediumMagboltz;
}  // namespace Garfield

class PHGarfield : public SubsysReco
{
 public:
  PHGarfield(const std::string &name = "PHGarfield");
  ~PHGarfield() override = default;

  int InitRun(PHCompositeNode *) override;
  int process_event(PHCompositeNode * /*topNode*/) override;

  bool StopHere(const double x, const double y, const double z, const double zPrevious);

  void PrintMaps() const;
  void PrintGarfield(double x, double y, double z) const;

  //  These are left in public namespace for easy plotting macros...
  //  The user is encouraged to add more routine to fit their analysis goals...
  TPolyLine3D *ReverseDrift(double x_cm, double y_cm, double z_cm, double step_ns = 50.0);  // Drifts electrons from some initial point until they hit a detector boundary...

  double GetRadius(size_t index) const {return radii.at(index);}

 private:
  void GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double &bx_t, double &by_t, double &bz_t) const;      // Feeds magnetic field to Garfield
  void GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double &ex_vcm, double &ey_vcm, double &ez_vcm) const;  // Feeds electric field to Garfield
  void InitializeGas(const std::string &dir);
  void FillRadii();
  static double bounder(double phi, double phi_min);

  CDBTTree *m_cdbTPCMAPttree{nullptr};            // Locations of the pads from CDB...
  PHField3DCartesian *m_field{nullptr};           // The stanards sPHENIX field holding container.
  Garfield::ComponentUser *m_component{nullptr};  // This handles the interface of the electric and magnetic fields as handed to Garfield
  Garfield::MediumMagboltz *m_gas{nullptr};       // This is the pre-tabulated gas properties required by Garfield...

  //  These are utilities for a spot check of the overall routine:
  // std::string calibdir;
  // std::string m_DiodeContainerName;
  double PHI_MIN{-std::numbers::pi};
  std::array<double, 48> radii{};  // Radius on each layer just for test purposes...need to be cm!
};

#endif
