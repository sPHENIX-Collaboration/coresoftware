#ifndef RAWTOWERBUILDERCONE_H__
#define RAWTOWERBUILDERCONE_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeom;

/**
 * \brief SubsysReco module creating calorimeter tower objects (RawTowerv2) from hits
 * (PHG4Hit) in cone-shaped sensitive volumes like forward-calorimeter.
 *
 * \author Nils Feege <nils.feege@stonybrook.edu>
 *
 * Based on RawTowerBuilder and PHG4CylinderCellReco.
 *
 * Tower binning: All towers have the same size in theta (\b not eta) and phi.
 *
 * Example of usage in Fun4All macro:
 *
 * \code{.cc}
 * RawTowerBuilderCone *TowerBuilderCone = new RawTowerBuilderCone("RawTowerBuilderCone");
 * TowerBuilderCone->Detector("FHCAL");
 * TowerBuilderCone->SetZMin(350);
 * TowerBuilderCone->SetZMax(450);
 * TowerBuilderCone->SetEtaMin(1.15);
 * TowerBuilderCone->SetEtaMax(5.0);
 * TowerBuilderCone->SetEtaNBins(30);
 * TowerBuilderCone->SetPhiNBins(36);
 * se->registerSubsystem( TowerBuilderCone );
 * \endcode
 *
 * \see RawTowerBuilder (module doing the same for cylinder detectors)
 * \see PHG4CylinderCellReco
 * \see PHG4ConeDetector
 */
class RawTowerBuilderCone : public SubsysReco {

public:
  RawTowerBuilderCone(const std::string& name="RawTowerBuilderCone");
  virtual ~RawTowerBuilderCone(){}

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

  /** @name Set Parameters
   *  Group of functions to set geometry and tower binning parameters.
   */
  ///@{
  /** Unit: cm. */
  void SetZMin( double set ) { _zmin = set; }
  /** Unit: cm. */
  void SetZMax( double set ) { _zmax = set; }

  void SetEtaMin( double set ) { _etamin = set; }
  void SetEtaMax( double set ) { _etamax = set; }
  void SetEtaNBins( int set ) { _netabins = set; }

  /** Default: 0. */
  void SetPhiMin( double set ) { _phimin = set; }
  /** Default: 2*Pi. */
  void SetPhiMax( double set ) { _phimax = set; }
  void SetPhiNBins( int set ) { _nphibins = set; }


  void SetGroupID( std::string & id ) { GroupID = id; }
  ///@}

protected:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "TOWER_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  /** Calculate eta and phi for a hit based on hit position x, y, z. */
  std::pair<double, double> get_etaphi(const double x, const double y, const double z);

  /** Convert eta to theta. */
  double eta2theta(const double eta);

  /** Convert theta to eta. */
  double theta2eta(const double theta);

  /** Get eta / theta bin (tower) for given angle theta.
   * \param[in] theta Theta angle (of a hit).
   * \return Corresponding tower bin in theta direction.
   */
  double get_thetabin(const double theta);

  /** Get phi bin (tower) for given angle phi.
   * \param[in] phi Phi angle (of a hit).
   * \return Corresponding tower bin in phi direction.
   */
  double get_phibin(const double);

  RawTowerContainer* _towers;
  RawTowerGeom *rawtowergeom;

  std::string detector;
  std::string hitnodename;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;
  std::string GroupID;

  double emin;
  int _nlayers;
  int _nphibins;
  int _netabins;
  double _zmin;
  double _zmax;
  double _etamin;
  double _etamax;
  double _phimin;
  double _phimax;
  double _thetamin;
  double _thetamax;

  double _thetastep;
  double _phistep;

  PHTimeServer::timer _timer;

};

#endif /* RAWTOWERBUILDER_H__ */
