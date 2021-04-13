#ifndef G4CALO__RAWTOWERBUILDERDRCALO_H
#define G4CALO__RAWTOWERBUILDERDRCALO_H

#include <calobase/RawTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

/**
 * \brief SubsysReco module creating calorimeter tower objects (RawTowerv1) from hits
 * (PHG4Hit) using j,k indeces of these hits
 *
 */
class RawTowerBuilderDRCALO : public SubsysReco
{
 public:
  RawTowerBuilderDRCALO(const std::string &name = "RawTowerBuilderDRCALO");
  virtual ~RawTowerBuilderDRCALO() {}

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  /** Name of the detector node the G4Hits should be taken from.
   */
  void Detector(const std::string &d);
  /** Specifiy text-file with table for tower mapping
   */
  void GeometryTableFile(const std::string &d)
  {
    m_MappingTowerFile = d;
  }
  /** Define minimum tower energy. After processing an event, towers with lower energy
   * are will be deleted.
   */
  void EminCut(const double e) { m_Emin = e; }

  /** Get prefix for tower collection to identify simulated towers
   * before digitization.
   */
  std::string
  get_sim_tower_node_prefix() const
  {
    return m_SimTowerNodePrefix;
  }

  /** Set prefix for tower collection to identify simulated towers
   * before digitization.
   */
  void
  set_sim_tower_node_prefix(const std::string &simTowerNodePrefix)
  {
    m_SimTowerNodePrefix = simTowerNodePrefix;
  }

 private:
  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "TOWER_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  /** Read geometry information from table stored in text-file
   */
  bool ReadGeometryFromTable();

  RawTowerContainer *m_Towers;
  RawTowerGeomContainer *m_Geoms;

  std::string m_Detector;
  std::string m_SimTowerNodePrefix;

  std::string m_MappingTowerFile;


  RawTowerDefs::CalorimeterId m_CaloId;

  double m_GlobalPlaceInX;
  double m_GlobalPlaceInY;
  double m_GlobalPlaceInZ;

  double m_RotInX;
  double m_RotInY;
  double m_RotInZ;

  double m_Emin;

  std::map<std::string, double> m_GlobalParameterMap;
};

#endif
