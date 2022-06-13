#ifndef __EpFinderReco_H__
#define __EpFinderReco_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>
#include <iostream>

//Forward declerations
class PHCompositeNode;
class EpFinder;
class EpInfo;
class PHG4TruthInfoContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class PHG4HitContainer;
class TH1D;

class EpFinderReco: public SubsysReco
{
 public:
  
  EpFinderReco(const std::string &name="EpFinderReco");

  int Init(PHCompositeNode *);

  int process_event(PHCompositeNode *);

  int End(PHCompositeNode *);

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
  

  void GetEventPlanes(PHCompositeNode*);
  int GetNodes(PHCompositeNode *);
  int CreateNodes(PHCompositeNode *);
 
  int GetPhiBin(float tphi, int numPhiDivisions);
  float GetMeanPhi(int iphi, int numPhiDivisions);
    
  
  RawTowerContainer *_calib_towers;
  RawTowerGeomContainer *rawtowergeom;

  RawTowerContainer *cemctowers;
  RawTowerContainer *hcalotowers;
  RawTowerContainer *hcalitowers;
  RawTowerGeomContainer *cemctowergeom;
  RawTowerGeomContainer *hcalotowergeom;
  RawTowerGeomContainer *hcalitowergeom;
 
  PHG4HitContainer* e_hit_container;
  PHG4HitContainer* b_hit_container;

  EpInfo *_calo_EpInfo;
  EpInfo *_CEMCHCAL_EpInfo;
  EpInfo *_BBC_EpInfoS;
  EpInfo *_BBC_EpInfoN;
  EpInfo *_EPD_EpInfoS;
  EpInfo *_EPD_EpInfoN;
  EpInfo *_EPD_EpInfoS_calib;
  EpInfo *_EPD_EpInfoN_calib;
    
  EpFinder *EpFinder_1;
  EpFinder *EpFinder_2;
  EpFinder *EpFinder_3;
  EpFinder *EpFinder_4;
 
  std::string _algonode;
    
  int _do_ep;
  bool _do_sepd_calib;
  double _sepdmpv;
  std::string detector;
  std::string CaliTowerNodeName;
  std::string TowerGeomNodeName;
  std::string EPNodeName;
    
};

#endif //* __EpFinderReco_H__ *//

