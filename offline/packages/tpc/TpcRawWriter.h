#ifndef TPC_TPCRAWWRITER_H
#define TPC_TPCRAWWRITER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsGeometry.h>

#include <map> 
#include <vector>
#include <string>

class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class RawHitSetContainerv1;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4CylinderCellGeom;
class PHG4TpcCylinderGeomContainer;

//typedef std::pair<int, int> iphiz;
//typedef std::pair<double, iphiz> ihit;
typedef std::pair<unsigned short, unsigned short> iphiz;
typedef std::pair<unsigned short, iphiz> ihit;

class TpcRawWriter : public SubsysReco
{
 public:
  TpcRawWriter(const std::string &name = "TpcRawWriter");
  ~TpcRawWriter() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  
  void set_sector_fiducial_cut(const double cut){SectorFiducialCut = cut; }
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}
  void set_do_wedge_emulation(bool do_wedge){ do_wedge_emulation = do_wedge;}
  void set_do_sequential(bool do_seq){ do_sequential = do_seq;}
  void set_max_cluster_half_size_phi(unsigned short size) { MaxClusterHalfSizePhi = size ;}
  void set_max_cluster_half_size_z(unsigned short size) { MaxClusterHalfSizeZ = size ;}
  void set_drift_velocity_scale(double value) { m_drift_velocity_scale = value; }
  void set_cluster_version(int value) { cluster_version = value; }
  
 private:
  
  TrkrHitSetContainer *m_hits = nullptr;
  RawHitSetContainerv1 *m_rawhits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;
  //  ActsSurfaceMaps *m_surfMaps = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  bool do_hit_assoc = true;
  bool do_wedge_emulation = true;
  bool do_sequential = false;
  double pedestal = 74.4;
  double SectorFiducialCut = 0.5;
  unsigned short MaxClusterHalfSizePhi = 3;
  unsigned short MaxClusterHalfSizeZ = 5;
  int cluster_version = 3;
  /// drift velocity scale factor

  double m_drift_velocity_scale = 1.0;
  
  // TPC shaping offset correction parameters
  // From Tony Frawley May 13, 2021
  //  double par0_neg = 0.0503;
  // double par0_pos = -0.0503;
  
};

#endif
