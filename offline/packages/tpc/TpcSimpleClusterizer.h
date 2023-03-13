#ifndef TPC_TPCSIMPLECLUSTERIZER_H
#define TPC_TPCSIMPLECLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsGeometry.h>

#include <map> 
#include <vector>
#include <string>

class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;

//typedef std::pair<int, int> iphiz;
//typedef std::pair<double, iphiz> ihit;
typedef std::pair<unsigned short, unsigned short> iphiz;
typedef std::pair<unsigned short, iphiz> ihit;

class TpcSimpleClusterizer : public SubsysReco
{
 public:
  TpcSimpleClusterizer(const std::string &name = "TpcSimpleClusterizer");
  ~TpcSimpleClusterizer() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_sector_fiducial_cut(const double cut){SectorFiducialCut = cut; }
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}

 private:
  bool is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom) const;

  TrkrHitSetContainer *m_hits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;
  ActsGeometry *m_tGeometry = nullptr;

  bool do_hit_assoc = true;
  double pedestal = 74.4;
  double SectorFiducialCut = 0.5;
  
  // TPC shaping offset correction parameters
  // From Tony Frawley May 13, 2021
  double par0_neg = 0.0503;
  double par0_pos = -0.0503;

};

#endif
