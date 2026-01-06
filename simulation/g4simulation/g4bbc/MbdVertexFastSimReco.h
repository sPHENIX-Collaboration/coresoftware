#ifndef G4MBD_MBDVERTEXFASTSIMRECO_H
#define G4MBD_MBDVERTEXFASTSIMRECO_H

//===========================================================
/// \file MbdVertexFastSimReco.h
/// \brief simple truth vertex smearing algorithm
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <limits>
#include <string>  // for string

class PHCompositeNode;

/// \class MbdVertexFastSimReco
///
/// \brief simple truth vertex smearing algorithm
///
class MbdVertexFastSimReco : public SubsysReco
{
 public:
  MbdVertexFastSimReco(const std::string &name = "MbdVertexFastSimReco");
  ~MbdVertexFastSimReco() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_t_smearing(const float t_smear) { m_T_Smear = t_smear; }
  void set_z_smearing(const float z_smear) { m_Z_Smear = z_smear; }

 private:
  int CreateNodes(PHCompositeNode *topNode);

  float m_T_Smear{std::numeric_limits<float>::quiet_NaN()};
  float m_Z_Smear{std::numeric_limits<float>::quiet_NaN()};

  gsl_rng *RandomGenerator;
};

#endif  // G4MBD_MBDVERTEXFASTSIMRECO_H
