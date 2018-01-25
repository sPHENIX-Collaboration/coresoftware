#ifndef _RAW_TOWER_BUILDER_BY_HIT_INDEX__
#define _RAW_TOWER_BUILDER_BY_HIT_INDEX__

#include <calobase/RawTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <phool/PHTimeServer.h>

#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;
class PHG4HitContainer;

/**
 * \brief SubsysReco module creating calorimeter tower objects (RawTowerv1) from hits
 * (PHG4Hit) using j,k indeces of these hits
 *
 */
class RawTowerBuilderByHitIndex : public SubsysReco {

public:

  RawTowerBuilderByHitIndex( const std::string& name="RawTowerBuilderByHitIndex" );
  virtual ~RawTowerBuilderByHitIndex(){}

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  /** Name of the detector node the G4Hits should be taken from.
   */
  void Detector( const std::string &d );

  /** Specifiy text-file with table for tower mapping
   */
  void GeometryTableFile( const std::string &d )
  { mapping_tower_file_ = d; }

  /** Define minimum tower energy. After processing an event, towers with lower energy
   * are will be deleted.
   */
  void EminCut(const double e) {emin_ = e;}

  /** Get prefix for tower collection to identify simulated towers
   * before digitization.
   */
  std::string
  get_sim_tower_node_prefix() const
  {
    return sim_tower_node_prefix_;
  }

  /** Set prefix for tower collection to identify simulated towers
   * before digitization.
   */
  void
  set_sim_tower_node_prefix(std::string simTowerNodePrefix)
  {
    sim_tower_node_prefix_ = simTowerNodePrefix;
  }



protected:

  /** Create nodes for output.
   *
   * Name of output node for RawTowerContainer: "TOWER_" + detector;
   */
  void CreateNodes(PHCompositeNode *topNode);

  /** Read geometry information from table stored in text-file
   */
  bool ReadGeometryFromTable();

  RawTowerContainer* towers_;
  RawTowerGeomContainer* geoms_;

  std::string detector_;
  std::string node_name_hits_;
  std::string node_name_towers_;
  std::string node_name_tower_geometries_;
  std::string sim_tower_node_prefix_;

  std::string mapping_tower_file_;

  RawTowerDefs::CalorimeterId calo_id_;

  double global_place_in_x_;
  double global_place_in_y_;
  double global_place_in_z_;

  double rot_in_x_;
  double rot_in_y_;
  double rot_in_z_;

  double emin_;

  std::map< std::string, double > map_global_parameter_;

  PHTimeServer::timer timer_;

};

#endif
