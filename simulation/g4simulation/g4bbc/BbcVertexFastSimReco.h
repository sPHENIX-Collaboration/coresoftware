#ifndef G4BBC_BBCVERTEXFASTSIMRECO_H
#define G4BBC_BBCVERTEXFASTSIMRECO_H

//===========================================================
/// \file BbcVertexFastSimReco.h
/// \brief simple truth vertex smearing algorithm
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>                // for string

class PHCompositeNode;

/// \class BbcVertexFastSimReco
///
/// \brief simple truth vertex smearing algorithm
///
class BbcVertexFastSimReco : public SubsysReco
{
 public:
  BbcVertexFastSimReco(const std::string &name = "BbcVertexFastSimReco");
  ~BbcVertexFastSimReco() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_t_smearing(const float t_smear) { m_T_Smear = t_smear; }
  void set_z_smearing(const float z_smear) { m_Z_Smear = z_smear; }

 private:
  int CreateNodes(PHCompositeNode *topNode);

  float m_T_Smear;
  float m_Z_Smear;

  gsl_rng *RandomGenerator;
};

#endif  // G4BBC_BBCVERTEXFASTSIMRECO_H
