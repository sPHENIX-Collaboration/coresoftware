#ifndef __CALOEVALUATOR_H__
#define __CALOEVALUATOR_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <string>
#include <stdint.h>

class PHCompositeNode;
class PHG4HitContainer;
class Fun4AllHistoManager;
class TH1F;
class TTree;
class SvtxEvalStack;
class PHG4Particle;
class RawTowerGeom;
class RawTowerContainer;
class SvtxTrack;

/// \class QAG4SimulationCalorimeter
class QAG4SimulationCalorimeter : public SubsysReco
{

public:

  enum enu_flags
  {
    kProcessSF = 1 << 1, kProcessTower = 1 << 2, kProcessMCPhoton = 1 << 3,

    kDefaultFlag = kProcessSF | kProcessTower
  };

  QAG4SimulationCalorimeter(enu_flags flags = kDefaultFlag);
  virtual
  ~QAG4SimulationCalorimeter();

  int
  Init(PHCompositeNode *topNode);
  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);

  uint32_t
  get_flags() const
  {
    return _flags;
  }

  void
  set_flags(enu_flags flags)
  {
    _flags = (uint32_t) flags;
  }

  void
  set_flag(enu_flags flag)
  {
    _flags |= (uint32_t) flag;
  }

  bool
  flag(enu_flags flag)
  {
    return _flags & flag;
  }

  void
  reset_flag(enu_flags flag)
  {
    _flags &= ~(uint32_t) flag;
  }

  float
  get_mag_field() const
  {
    return _magfield;
  }
  void
  set_mag_field(float magfield)
  {
    _magfield = magfield;
  }

private:

  int
  Init_SF(PHCompositeNode *topNode);
  int
  process_event_SF(PHCompositeNode *topNode);

  int
  Init_Tower(PHCompositeNode *topNode);
  int
  process_event_Tower(PHCompositeNode *topNode);

//  int
//  Init_MCPhoton(PHCompositeNode *topNode);
//  int
//  process_event_MCPhoton(PHCompositeNode *topNode);

  enum enu_calo
  {
    kCEMC, kHCALIN, kHCALOUT
  };

  SvtxEvalStack * _eval_stack;

  double _magfield;

  uint32_t _flags;
  unsigned long _ievent;

  PHG4HitContainer* _hcalout_hit_container;
  PHG4HitContainer* _hcalin_hit_container;
  PHG4HitContainer* _cemc_hit_container;
  PHG4HitContainer* _hcalout_abs_hit_container;
  PHG4HitContainer* _hcalin_abs_hit_container;
  PHG4HitContainer* _cemc_abs_hit_container;
  PHG4HitContainer* _magnet_hit_container;
  PHG4HitContainer* _bh_hit_container;
  PHG4HitContainer* _bh_plus_hit_container;
  PHG4HitContainer* _bh_minus_hit_container;
  PHG4HitContainer* _cemc_electronics_hit_container;
  PHG4HitContainer* _hcalin_spt_hit_container;
  PHG4HitContainer* _svtx_hit_container;
  PHG4HitContainer* _svtx_support_hit_container;
};

#endif // __CALOEVALUATOR_H__
