// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TpcPadBaselineShift_H
#define G4TPC_PHG4TpcPadBaselineShift_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>

#include <map> 
#include <vector>
#include <string>

class PHCompositeNode;

class TTree;
class TFile;

class TrkrHitSet;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4CylinderCellGeom;
class PHG4CylinderCellGeomContainer;

typedef std::pair<unsigned short, unsigned short> iphiz;
typedef std::pair<unsigned short, iphiz> ihit;

class PHG4TpcPadBaselineShift : public SubsysReco
{
 public:

  PHG4TpcPadBaselineShift(const std::string &name = "PHG4TpcPadBaselineShift");

  virtual ~PHG4TpcPadBaselineShift();
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override; 
  //int ResetEvent(PHCompositeNode *topNode) override;

  //int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  //int Reset(PHCompositeNode * /*topNode*/) override;

  //void Print(const std::string &what = "ALL") const override;

   void setScale(float CScale);
   void setFileName(std::string filename);
   void writeTree(int f_writeTree);
 

 private:
   bool is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom);
   float _hit_z  ;
   float _hit_r  ;
   float _hit_phi;
   float _hit_e;
   int _hit_adc;
   int _hit_adc_bls;
   int _hit_layer;
   int _hit_sector;

   TrkrHitSetContainer *m_hits;
   TrkrClusterContainer *m_clusterlist;
   TrkrClusterHitAssoc *m_clusterhitassoc;
   ActsSurfaceMaps *m_surfMaps;
   ActsTrackingGeometry *m_tGeometry;

   bool do_hit_assoc;
   double pedestal;
   int _writeTree;
   double SectorFiducialCut;

   int NSearch;
   int NZBinsMax;
   float _CScale;

   TFile *outfile;
   std::string _filename;

   TTree *_rawHits;

};

#endif // PHG4TpcPadBaselineShift_H
