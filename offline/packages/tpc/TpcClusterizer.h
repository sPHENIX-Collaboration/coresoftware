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
  virtual ~TpcClusterizer(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_sector_fiducial_cut(const double cut){SectorFiducialCut = cut; }
  void set_search_bins(const int bins){NSearch = bins;}
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}

 private:
  bool is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom);

  TrkrHitSetContainer *m_hits;
  TrkrClusterContainer *m_clusterlist;
  TrkrClusterHitAssoc *m_clusterhitassoc;
  ActsSurfaceMaps *m_surfMaps;
  ActsTrackingGeometry *m_tGeometry;

  bool do_hit_assoc;
  double zz_shaping_correction;
  double pedestal;
  double SectorFiducialCut;

  int NSearch;
  int NZBinsMax;

};

#endif
