#ifndef G4BBC_BBCVERTEXFASTSIMRECO_H
#define G4BBC_BBCVERTEXFASTSIMRECO_H

//===========================================================
/// \file BbcVertexFastSimReco.h
/// \brief simple truth vertex smearing algorithm
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

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
  virtual ~BbcVertexFastSimReco();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_t_smearing(const float t_smear) { m_T_Smear = t_smear; }
  void set_z_smearing(const float z_smear) { m_Z_Smear = z_smear; }

 private:
  int CreateNodes(PHCompositeNode *topNode);

  float m_T_Smear;
  float m_Z_Smear;

#if !defined(__CINT__) || defined(__CLING__)
  gsl_rng *RandomGenerator;
#endif
};

#endif  // G4BBC_BBCVERTEXFASTSIMRECO_H
