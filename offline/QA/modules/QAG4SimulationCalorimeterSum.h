#ifndef QA_QAG4SIMULATIONCALORIMETERSUM_H
#define QA_QAG4SIMULATIONCALORIMETERSUM_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <memory>
#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4Particle;
class CaloEvalStack;
class SvtxEvalStack;
class SvtxTrack;

/// \class QAG4SimulationCalorimeterSum
class QAG4SimulationCalorimeterSum : public SubsysReco
{
 public:
  enum enu_flags
  {
    //    kProcessTower = 1 << 1,
    kProcessCluster = 1 << 2,    // histograms of best cluster matched to truth particle
    kProcessTrackProj = 1 << 3,  // histograms of tower/tower sums VS track projections

    kDefaultFlag = kProcessCluster | kProcessTrackProj
  };

  QAG4SimulationCalorimeterSum(enu_flags flags = kDefaultFlag);

  virtual ~QAG4SimulationCalorimeterSum() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

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

  //! common prefix for QA histograms
  std::string
  get_histo_prefix();

  std::string
  get_calo_name_cemc() const
  {
    return _calo_name_cemc;
  }

  void
  set_calo_name_cemc(const std::string &caloNameCemc)
  {
    _calo_name_cemc = caloNameCemc;
  }

  std::string
  get_calo_name_hcalin() const
  {
    return _calo_name_hcalin;
  }

  void
  set_calo_name_hcalin(const std::string &caloNameHcalin)
  {
    _calo_name_hcalin = caloNameHcalin;
  }

  std::string
  get_calo_name_hcalout() const
  {
    return _calo_name_hcalout;
  }

  void
  set_calo_name_hcalout(const std::string &caloNameHcalout)
  {
    _calo_name_hcalout = caloNameHcalout;
  }

  float get_mag_field() const { return _magField; }
  void set_mag_field(float magField) { _magField = magField; }

 private:
  //  int
  //  Init_Tower(PHCompositeNode *topNode);
  //  int
  //  process_event_Tower(PHCompositeNode *topNode);

  int Init_Cluster(PHCompositeNode *topNode);
  int process_event_Cluster(PHCompositeNode *topNode);

  int Init_TrackProj(PHCompositeNode *topNode);
  int process_event_TrackProj(PHCompositeNode *topNode);

  std::shared_ptr<CaloEvalStack> _caloevalstack_cemc;
  std::shared_ptr<CaloEvalStack> _caloevalstack_hcalin;
  std::shared_ptr<CaloEvalStack> _caloevalstack_hcalout;
  std::shared_ptr<SvtxEvalStack> _svtxevalstack;

  uint32_t _flags;

  std::string _calo_name_cemc;
  std::string _calo_name_hcalin;
  std::string _calo_name_hcalout;

  PHG4TruthInfoContainer *_truth_container;

  //! fetch the truth particle to be analyzed. By default it is the last primary particle in truth container (therefore works in single particle embedding)
  PHG4Particle *
  get_truth_particle();

  //! fetch tower around track and histogram energy distributions
  bool
  eval_trk_proj(const std::string &detector, SvtxTrack *track,
                PHCompositeNode *topNode);

  //! central magnetic field strength in T
  float _magField;

  enum
  {
    //! max number of tower row/column to process around a track projection
    Max_N_Tower = 11
  };
};

#endif  // QA_QAG4SIMULATIONCALORIMETERSUM_H
