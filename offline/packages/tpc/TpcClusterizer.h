#ifndef TPC_TPCCLUSTERIZER_H
#define TPC_TPCCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsGeometry.h>

#include <map> 
#include <vector>
#include <string>

class ClusHitsVerbosev1;
class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class RawHitSet;
class RawHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrainingHitsContainer;
class PHG4TpcCylinderGeom;
class PHG4TpcCylinderGeomContainer;

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
  void set_store_hits(bool store_hits) { _store_hits = store_hits; }
  void set_use_nn(bool use_nn) { _use_nn = use_nn; }
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}
  void set_do_wedge_emulation(bool do_wedge){ do_wedge_emulation = do_wedge;}
  void set_do_sequential(bool do_seq){ do_sequential = do_seq;}
  void set_do_split(bool split){ do_split = split;}
  void set_pedestal(float val) { pedestal = val;}
  void set_seed_threshold(float val) { seed_threshold = val;}
  void set_edge_threshold(float val) { edge_threshold = val;}
  void set_min_err_squared(float val) { min_err_squared = val;}
  void set_min_clus_size(float val) { min_clus_size = val;}
  void set_min_adc_sum(float val) { min_adc_sum = val;}
  void set_remove_singles(bool do_sing){ do_singles = do_sing;}
  void set_read_raw(bool read_raw){ do_read_raw = read_raw;}
  void set_max_cluster_half_size_phi(unsigned short size) { MaxClusterHalfSizePhi = size ;}
  void set_max_cluster_half_size_z(unsigned short size) { MaxClusterHalfSizeT = size ;}
  
  void set_ClusHitsVerbose(bool set=true) { record_ClusHitsVerbose = set; };
  void set_rawdata_reco() {
    set_do_hit_association(false);
    set_do_split(false);
    set_pedestal(0);
    set_seed_threshold(5);
    set_edge_threshold(3);
    set_min_err_squared(0);
    set_min_clus_size(0);
    set_min_adc_sum(5);
    set_remove_singles(true);
    set_max_cluster_half_size_phi(5);
    set_max_cluster_half_size_z(8);
  };
  ClusHitsVerbosev1* mClusHitsVerbose { nullptr };
  
 private:
  bool is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom) const;
  bool record_ClusHitsVerbose { false };

  TrkrHitSetContainer *m_hits = nullptr;
  RawHitSetContainer *m_rawhits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  bool _store_hits = false;
  bool _use_nn = false;
  bool do_hit_assoc = true;
  bool do_wedge_emulation = false;
  bool do_sequential = false;
  bool do_read_raw = false;
  bool do_singles = false;
  bool do_split = true;
  double pedestal = 74.4;
  double seed_threshold = 5;
  double edge_threshold = 0;
  double min_err_squared = 0.01;
  double min_clus_size = 1;
  double min_adc_sum = 10;
  double SectorFiducialCut = 0.5;
  unsigned short MaxClusterHalfSizePhi = 3;
  unsigned short MaxClusterHalfSizeT = 5;
 
  double m_tdriftmax = 0;
  double AdcClockPeriod = 53.0;   // ns 

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  double m_sampa_tbias = 39.6;  // ns  

  TrainingHitsContainer *m_training;
};

#endif
