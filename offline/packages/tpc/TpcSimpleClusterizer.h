#ifndef TPC_TPCSIMPLECLUSTERIZER_H
#define TPC_TPCSIMPLECLUSTERIZER_H

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
  bool is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom) const;

  TrkrHitSetContainer *m_hits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;
  ActsSurfaceMaps *m_surfMaps = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;

  bool do_hit_assoc = true;
  double pedestal = 74.4;
  double SectorFiducialCut = 0.5;

  // TPC shaping offset correction parameters
  // From Tony Frawley May 13, 2021
//   std::pair<double,double> par0_neg = std::make_pair(0.0538913, 0.000252096);
//   std::pair<double,double> par0_pos = std::make_pair(-0.0647731, 0.000296734);
//   std::pair<double,double> par1_neg = std::make_pair(-0.000208279, 1.9205e-06);
//   std::pair<double,double> par1_pos = std::make_pair(-0.000195514, 2.26467e-06);
  
  // revisited by Hugo May 28, 2021
  // should check with Tony, in particular the inversion of the par0_pos layer slope
  std::pair<double,double> par0_neg = std::make_pair(0.05465077, 0.000252096);
  std::pair<double,double> par0_pos = std::make_pair(-0.05392444, -0.000296734);
  std::pair<double,double> par1_neg = std::make_pair(-0.000208279, 1.9205e-06);
  std::pair<double,double> par1_pos = std::make_pair(-0.000195514, 2.26467e-06);

};

#endif
