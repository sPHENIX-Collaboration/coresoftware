// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file               TrackToCalo.h
 *  \brief              Track To Calo matching, for TPC drift velocity calibration
 *  \author Xudong Yu <xyu3@bnl.gov>
 */

#ifndef TRACKTOCALO_H
#define TRACKTOCALO_H

#include <fun4all/SubsysReco.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <calobase/RawClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;
class TFile;
class TTree;

class TrackToCalo : public SubsysReco
{
 public:

  TrackToCalo(const std::string &name = "TrackToCalo", const std::string &file = "output.root");

  ~TrackToCalo() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void ResetTreeVectors();

  void EMcalRadiusUser(bool use) {m_use_emcal_radius = use;}
  void setEMcalRadius(float r) {m_emcal_radius_user = r;}

  void setRawClusContEMName(const std::string& name) {m_RawClusCont_EM_name = name;}

  void setTrackPtLowCut(float pt) {m_track_pt_low_cut = pt;}
  void setEmcalELowCut(float e) {m_emcal_e_low_cut = e;}
  void setnTpcClusters(int n) {m_ntpc_low_cut = n;}

 private:
   std::string _outfilename;
   TFile *_outfile = nullptr;
   TTree *_tree = nullptr;

   std::vector<float> _track_phi_emc;
   std::vector<float> _track_z_emc;
   std::vector<float> _emcal_phi;
   std::vector<float> _emcal_z;

   std::string m_RawClusCont_EM_name = "TOPOCLUSTER_EMCAL";

   SvtxTrackMap *svtxTrackMap = nullptr;
   RawClusterContainer *rawClusterContainer = nullptr;
   TrkrClusterContainer *trkrClusterContainer = nullptr;
   RawTowerGeomContainer *rawTowerGeomContainer = nullptr;
   SvtxTrackState *thisState = nullptr;
   SvtxTrack *track = nullptr;
   TrackSeed *tpc_seed = nullptr;
   TrkrCluster *trkrCluster = nullptr;
   RawCluster *cluster = nullptr;

   int cnt = 0;
   bool m_use_emcal_radius = false;
   float m_emcal_radius_user = 93.5;
   float radius_scale = 1;
   float m_track_pt_low_cut = 0.5;
   float m_emcal_e_low_cut = 0.2;
   int m_ntpc_low_cut = 22;
};

#endif // TRACKTOCALO_H
