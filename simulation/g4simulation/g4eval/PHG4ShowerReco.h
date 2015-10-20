#ifndef __PHG4SHOWERRECO_H__
#define __PHG4SHOWERRECO_H__

//===============================================
/// \file PHG4ShowerReco.h
/// \brief Replaces disk space used by g4hit & g4cells with g4showers
/// \author Michael P. McCumber
//===============================================

#include "CaloTruthEval.h"
#include "CaloRawTowerEval.h"

#include <g4main/PHG4ShowerMap.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4cemc/RawTowerContainer.h>

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>

#include <string>
#include <set>
#include <map>

/// \class PHG4ShowerReco
///
/// \brief Replaces disk space used by g4hit & g4cells with g4showers
///
/// Plan: This module will used the evaluator objects to trace the
/// connections between RawTowers and the PHG4Particles and build a
/// new space-saving ancestry through the PHG4Shower objects and then
/// delete the calorimeter-only g4hit, g4cells, and secondaries
///
class PHG4ShowerReco : public SubsysReco {
  
public:
 
  PHG4ShowerReco(const std::string &name = "PHG4SHOWERRECO");
  virtual ~PHG4ShowerReco() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

  int CreateNodes(PHCompositeNode* topNode);
  int GetNodes(PHCompositeNode* topNode);

  unsigned int _ievent;
  
  PHG4TruthInfoContainer* _truth_info;
  
  std::map<PHG4Shower::VOLUME,std::string>        _volume_names;
  std::map<PHG4Shower::VOLUME,CaloTruthEval*>     _volume_truthevals;
  std::map<PHG4Shower::VOLUME,PHG4HitContainer*>  _volume_g4hits;
  std::map<PHG4Shower::VOLUME,RawTowerContainer*> _volume_simtowers;
  std::map<PHG4Shower::VOLUME,RawTowerContainer*> _volume_rawtowers;
  std::map<PHG4Shower::VOLUME,RawTowerContainer*> _volume_calibtowers;
  std::map<PHG4Shower::VOLUME,CaloRawTowerEval*>  _volume_towerevals;
  
  PHG4ShowerMap* _shower_map;
};

#endif // __PHG4SHOWERRECO_H__
