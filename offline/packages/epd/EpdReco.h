// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef __EPDRECO_H__
#define __EPDRECO_H__

//===========================================================
/// \author Ejiro Umaka
//===========================================================

#include <cdbobjects/CDBTTree.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_const_cgsm.h>

#include <string> // for string
#include <vector> // for vector

class PHCompositeNode;

class EpdReco : public SubsysReco {
public:
  EpdReco(const std::string &name = "EpdReco");
  ~EpdReco() override = default;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode * /*topNode*/) override;

private:
  CDBTTree *cdbttree{nullptr};
  std::string m_fieldname;
  std::string m_calibName;
  bool m_overrideCalibName{false};
  bool m_overrideFieldName{false};
  std::string m_Detector = "SEPD";
  std::string m_TowerInfoNodeName_calib = "TOWERINFO_CALIB_SEPD";
  int Getrmap(int rindex);
  int Getphimap(int phiindex);
  float GetTilePhi(int thisphi);
  float GetTilePhi0(int thisphi0);
  float GetTileR(int thisr);
  float GetTileZ(int thisz);
  void CreateNodes(PHCompositeNode *topNode);
};

#endif // EPDRECO_H