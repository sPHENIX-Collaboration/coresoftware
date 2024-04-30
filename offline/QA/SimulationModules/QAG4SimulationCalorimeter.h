#ifndef QA_QAG4SIMULATIONCALORIMETER_H
#define QA_QAG4SIMULATIONCALORIMETER_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <memory>
#include <string>

class CaloEvalStack;
class PHCompositeNode;
class PHG4HitContainer;
class PHG4TruthInfoContainer;

/// \class QAG4SimulationCalorimeter
class QAG4SimulationCalorimeter : public SubsysReco
{
 public:
  enum enu_flags
  {
    kProcessG4Hit = 1 << 1,
    kProcessTower = 1 << 2,
    kProcessCluster = 1 << 3,

    kDefaultFlag = kProcessG4Hit | kProcessTower | kProcessCluster
  };

  QAG4SimulationCalorimeter(const std::string &calo_name, enu_flags flags =
                                                              kDefaultFlag);
  virtual ~QAG4SimulationCalorimeter() {}

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

 private:
  int Init_G4Hit(PHCompositeNode *topNode);
  int process_event_G4Hit(PHCompositeNode *topNode);

  int Init_Tower(PHCompositeNode *topNode);
  int process_event_Tower(PHCompositeNode *topNode);

  int Init_Cluster(PHCompositeNode *topNode);
  int process_event_Cluster(PHCompositeNode *topNode);

  std::shared_ptr<CaloEvalStack> _caloevalstack;

  std::string _calo_name;
  uint32_t _flags;

  PHG4HitContainer *_calo_hit_container;
  PHG4HitContainer *_calo_abs_hit_container;
  PHG4TruthInfoContainer *_truth_container;
};

#endif  // QA_QAG4SIMULATIONCALORIMETER_H
