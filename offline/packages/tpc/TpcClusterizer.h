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
  ~TpcClusterizer() override{}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

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
  double pedestal;
  double SectorFiducialCut;

  // TPC shaping offset correction parameters
  // From Tony Frawley May 13, 2021
  std::pair<double,double> par0_neg = std::make_pair(0.0538913, 0.000252096);
  std::pair<double,double> par0_pos = std::make_pair(-0.0647731, 0.000296734);
  std::pair<double,double> par1_neg = std::make_pair(-0.000208279, 1.9205e-06);
  std::pair<double,double> par1_pos = std::make_pair(-0.000195514, 2.26467e-06);
  
  int NSearch;
  int NZBinsMax;

};

#endif
