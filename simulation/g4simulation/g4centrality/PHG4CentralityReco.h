#ifndef G4CENTRALITY_PHG4CENTRALITYRECO_H
#define G4CENTRALITY_PHG4CENTRALITYRECO_H

//===========================================================
/// \file PHG4CentralityReco.h
/// \brief Centrality quantity construction & calibration
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <string>

class PHCompositeNode;

class PHG4CentralityReco : public SubsysReco
{
 public:
  PHG4CentralityReco(const std::string &name = "PHG4CentralityReco");
  ~PHG4CentralityReco() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  void CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  float _epd_N = NAN;
  float _epd_S = NAN;
  float _epd_NS = NAN;

  float _mbd_N = NAN;
  float _mbd_S = NAN;
  float _mbd_NS = NAN;
};

#endif
