// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANE_EPFINDERRECO_H
#define EVENTPLANE_EPFINDERRECO_H

#include <fun4all/SubsysReco.h>

#include <string>

//Forward declarations
class PHCompositeNode;
class EpFinder;
class EpInfo;
class RawTowerContainer;
class RawTowerGeomContainer;
class PHG4HitContainer;
class TH1D;

class EpFinderReco : public SubsysReco
{
 public:
  EpFinderReco(const std::string &name = "EpFinderReco");
  ~EpFinderReco() override;

  int Init(PHCompositeNode *) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  int ResetEvent(PHCompositeNode * /*topNode*/) override;

  void set_algo_node(const std::string &algonode) { _algonode = algonode; }

  void set_ep_mode(int do_ep)
  {
    _do_ep = do_ep;
  }

  void set_sEPD_calib(bool do_sepd_calib)
  {
    _do_sepd_calib = do_sepd_calib;
  }

  void set_sEPD_MPV_in_GeV(const double sepdmpv)
  {
    _sepdmpv = sepdmpv;
  }

  void Detector(const std::string &d)
  {
    detector = d;
  }

 private:
  TH1D *heta = nullptr;
  TH1D *hcent = nullptr;

  void GetEventPlanes(PHCompositeNode *);
  int GetNodes(PHCompositeNode *);
  int CreateNodes(PHCompositeNode *);

  int GetPhiBin(float tphi, int numPhiDivisions);
  float GetMeanPhi(int iphi, int numPhiDivisions);

  std::string _algonode = "EVENT_PLANE";
  int _do_ep = 0;
  bool _do_sepd_calib = false;
  double _sepdmpv = 1.;

  RawTowerContainer *_calib_towers = nullptr;
  RawTowerGeomContainer *rawtowergeom = nullptr;
  RawTowerContainer *cemctowers = nullptr;
  RawTowerGeomContainer *cemctowergeom = nullptr;
  RawTowerContainer *hcalotowers = nullptr;
  RawTowerGeomContainer *hcalotowergeom = nullptr;
  RawTowerContainer *hcalitowers = nullptr;
  RawTowerGeomContainer *hcalitowergeom = nullptr;
  PHG4HitContainer *e_hit_container = nullptr;
  PHG4HitContainer *b_hit_container = nullptr;

  std::string detector = "CEMC";

  EpFinder *EpFinder_1 = nullptr;
  EpFinder *EpFinder_2 = nullptr;
  EpFinder *EpFinder_3 = nullptr;
  EpFinder *EpFinder_4 = nullptr;

  EpInfo *_CALO_EpInfo = nullptr;
  EpInfo *_BBC_EpInfoS = nullptr;
  EpInfo *_BBC_EpInfoN = nullptr;
  EpInfo *_EPD_EpInfoS = nullptr;
  EpInfo *_EPD_EpInfoN = nullptr;
  EpInfo *_CEMCHCAL_EpInfo = nullptr;
  EpInfo *_EPD_EpInfoS_calib = nullptr;
  EpInfo *_EPD_EpInfoN_calib = nullptr;

  std::string CaliTowerNodeName;
  std::string TowerGeomNodeName;
  std::string EPNodeName;
};

#endif  //* EVENTPLANE_EPFINDERRECO_H *//
