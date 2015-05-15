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
class CrystalCalorimeterDigitization : public SubsysReco {

public:
  CrystalCalorimeterDigitization(const std::string& name="CrystalCalorimeterDigitization", const std::string& nameRaw="DEFAULT_RAW", const std::string& nameDigi="DEFAULT_DIGI", int randSeed=1234);

  virtual ~CrystalCalorimeterDigitization(){}

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  void SetMeanLightYield( double meanLY )
  {
    _meanLY = meanLY;
  }

  void SetApplyPhotonStatistic( bool switchOn=true )
  {
    _applyPhotonStatistic = switchOn;
  }


protected:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "DIGI_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  void ApplyPhotonStatistic( RawTowerv2& tower );

  RawTowerContainer* _towersDigi;

  std::string _nodeNameTowerRaw;
  std::string _nodeNameTowerDigi;

  double _meanLY;

  bool _applyPhotonStatistic;
  int _randSeed;

  PHTimeServer::timer _timer;

};

#endif /* RAWTOWERBUILDER_H__ */
