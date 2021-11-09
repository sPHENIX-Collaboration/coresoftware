#ifndef G4CENTRALITY_PHG4CENTRALITYRECO_H
#define G4CENTRALITY_PHG4CENTRALITYRECO_H

//===========================================================
/// \file PHG4CentralityReco.h
/// \brief Centrality quantity construction & calibration
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <map>
#include <utility>               // for pair

#include <centrality/CentralityInfov1.h>

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

  float _epd_NS;

  float _mbd_NS;

};

#endif
