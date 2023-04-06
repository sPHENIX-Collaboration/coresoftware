// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4VERTEX_GLOBALVERTEXFASTSIMRECO_H
#define G4VERTEX_GLOBALVERTEXFASTSIMRECO_H

//===========================================================
/// \file GlobalVertexFastSimReco.h
/// \brief simple truth vertex smearing algorithm
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <limits>
#include <string>  // for string

class PHCompositeNode;

/// \class GlobalVertexFastSimReco
///
/// \brief simple truth vertex smearing algorithm
///
class GlobalVertexFastSimReco : public SubsysReco
{
 public:
  GlobalVertexFastSimReco(const std::string &name = "GlobalVertexFastSimReco");
  ~GlobalVertexFastSimReco() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_x_smearing(const float x_smear) { _x_smear = x_smear; }
  void set_y_smearing(const float y_smear) { _y_smear = y_smear; }
  void set_z_smearing(const float z_smear) { _z_smear = z_smear; }
  void set_t_smearing(const float t_smear) { _t_smear = t_smear; }

 private:
  int CreateNodes(PHCompositeNode *topNode);

  float _x_smear = std::numeric_limits<float>::quiet_NaN();
  float _y_smear = std::numeric_limits<float>::quiet_NaN();
  float _z_smear = std::numeric_limits<float>::quiet_NaN();
  float _t_smear = std::numeric_limits<float>::quiet_NaN();
  gsl_rng *RandomGenerator = nullptr;
};

#endif  // G4VERTEX_GLOBALVERTEXFASTSIMRECO_H
