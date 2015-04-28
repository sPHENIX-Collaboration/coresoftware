#ifndef TOWERDIGITIZATION_H__
#define TOWERDIGITIZATION_H__

#include "RawTowerv2.h"
#include "RawTowerContainer.h"

#include <fun4all/SubsysReco.h>

#include <phool/PHTimeServer.h>

#include <string>

class PHCompositeNode;
class RawTowerv2;
class RawTowerContainer;

/**
 * \brief SubsysReco module adding post-Geant4 detector effects to calorimeter tower objects.
 *
 * \author Nils Feege <nils.feege@stonybrook.edu>
 *
 */
class TowerDigitization : public SubsysReco {

public:
  TowerDigitization(const std::string& name="TowerDigitization", const std::string& nameRaw="DEFAULT_RAW", const std::string& nameDigi="DEFAULT_DIGI", int randSeed=1234);

  virtual ~TowerDigitization(){}

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  void SetApplyPoissonSmearing( bool switchOn=true , int poissonMean=200 )
  { _applyPoissonSmearing = switchOn;
    _poissonMean = poissonMean; }

protected:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "DIGI_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  void ApplyPoissonSmearing( RawTowerv2& tower );

  RawTowerContainer* _towersDigi;

  std::string _nodeNameTowerRaw;
  std::string _nodeNameTowerDigi;

  int _randSeed;

  bool _applyPoissonSmearing;

  double _poissonMean;

  PHTimeServer::timer _timer;

};

#endif /* RAWTOWERBUILDER_H__ */
