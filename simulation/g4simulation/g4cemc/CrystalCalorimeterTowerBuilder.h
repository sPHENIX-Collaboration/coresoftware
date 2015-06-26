#ifndef _CRYSTAL_CALORIMETER_TOWER_BUILDER_H__
#define _CRYSTAL_CALORIMETER_TOWER_BUILDER_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class PHG4HitContainer;

/**
 * \brief SubsysReco module creating calorimeter tower objects (RawTowerv2) from hits
 * (PHG4Hit) using j,k indeces of these hits
 *
 * \author Nils Feege <nils.feege@stonybrook.edu>
 *
 */
class CrystalCalorimeterTowerBuilder : public SubsysReco {

public:
  CrystalCalorimeterTowerBuilder(const std::string& name="CrystalCalorimeterTowerBuilder");
  virtual ~CrystalCalorimeterTowerBuilder(){}

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  /** Name of the detector node the G4Hits should be taken from.
   */
  void Detector(const std::string &d) {detector = d;}

  /** Define minimum tower energy. After processing an event, towers with lower energy
   * are will be deleted.
   */
  void EminCut(const double e) {emin = e;}

protected:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "TOWER_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _towers;

  std::string detector;
  std::string hitnodename;
  std::string TowerNodeName;

  double emin;

  PHTimeServer::timer _timer;

};

#endif
