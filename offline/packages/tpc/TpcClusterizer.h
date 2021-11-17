#ifndef TPC_TPCCLUSTERIZER_H
#define TPC_TPCCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>

#include <map> 
#include <vector>
#include <string>

class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4CylinderCellGeom;
class PHG4CylinderCellGeomContainer;

//typedef std::pair<int, int> iphiz;
//typedef std::pair<double, iphiz> ihit;
typedef std::pair<unsigned short, unsigned short> iphiz;
typedef std::pair<unsigned short, iphiz> ihit;

class TpcClusterizer : public SubsysReco
{
 public:
  TpcClusterizer(const std::string &name = "TpcClusterizer");
  ~TpcClusterizer() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_sector_fiducial_cut(const double cut){SectorFiducialCut = cut; }
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}
  void set_max_cluster_half_size_phi(unsigned short size) { MaxClusterHalfSizePhi = size ;}
  void set_max_cluster_half_size_z(unsigned short size) { MaxClusterHalfSizeZ = size ;}

  void set_drift_velocity_scale(double value) { m_drift_velocity_scale = value; }
  
 private:
  bool is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom) const;

  TrkrHitSetContainer *m_hits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;
  ActsSurfaceMaps *m_surfMaps = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;

  bool do_hit_assoc = true;
  double pedestal = 74.4;
  double SectorFiducialCut = 0.5;
  unsigned short MaxClusterHalfSizePhi = 3;
  unsigned short MaxClusterHalfSizeZ = 5;

  /// drift velocity scale factor
  /** 
   * represents the ratio vdrift_measured/vdrift_true
   * it is used to get cluster z from its timebin
   **/
  double m_drift_velocity_scale = 1.0;
  
  // TPC shaping offset correction parameters
  // From Tony Frawley May 13, 2021
  double par0_neg = 0.0503;
  double par0_pos = -0.0503;
  
};

#endif
