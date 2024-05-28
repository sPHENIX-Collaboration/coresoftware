// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GEANTTESTER_H
#define GEANTTESTER_H

#include <fun4all/SubsysReco.h>

#include <decayfinder/DecayFinderContainer_v1.h>  // for DecayFinderContainer_v1

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4Particle;
class PHG4VtxPoint;

class truthDecayTester : public SubsysReco
{
 public:
  truthDecayTester(const std::string &name = "truthDecayTester");

  virtual ~truthDecayTester();

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

  std::string get_histo_prefix();

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void setMotherPDG(int PDGID) { m_decay_pdg_id = PDGID; }
  void setMinPT(float value) { m_min_pt = value; }
  void setEtaRange(float min, float max)
  {
    m_min_eta = min;
    m_max_eta = max;
  }
  void setDFNodeName(const std::string &name) { m_df_module_name = name; }
  void setOutputName(const std::string &name) { m_outfile_name = name; }
  void writeTuple(bool write) { m_write_nTuple = write; }

 private:
  unsigned int m_nTracks;
  float m_min_pt;
  float m_min_eta;
  ;
  float m_max_eta;
  bool m_write_nTuple;
  int m_decay_pdg_id;

  PHG4TruthInfoContainer *m_truth_info;
  PHG4Particle *m_g4particle;
  std::string m_df_module_name;
  std::string m_outfile_name;
  TFile *m_outfile;
  TTree *m_tree;
  bool m_write_QAHists;
  DecayFinderContainer_v1 *m_decayMap = nullptr;

  void initializeBranches();
  void getMotherPDG(PHCompositeNode *topNode);
  std::vector<int> getDecayFinderMothers(PHCompositeNode *topNode);
  bool isInRange(float min, float value, float max);
  void resetValues();

  static const int max_tracks = 20;
  unsigned int m_event_number = 0;
  float m_mother_mass = -99;
  float m_daughter_sum_mass = 0;
  float m_mother_decayLength = -99;
  float m_mother_decayTime = -99;
  int m_mother_pdg_id = -99;
  float m_mother_px = 0;
  float m_mother_py = 0;
  float m_mother_pz = 0;
  float m_mother_pE = 0;
  float m_mother_pT = 0;
  float m_mother_eta = 0;
  int m_mother_barcode = -99;
  int m_track_pdg_id[max_tracks] = {0};
  float m_track_px[max_tracks] = {0};
  float m_track_py[max_tracks] = {0};
  float m_track_pz[max_tracks] = {0};
  float m_track_pE[max_tracks] = {0};
  float m_track_pT[max_tracks] = {0};
  float m_track_eta[max_tracks] = {0};
  float m_track_mass[max_tracks] = {0};
  int m_track_mother_barcode[max_tracks] = {0};
  float m_delta_px = 0;
  float m_delta_py = 0;
  float m_delta_pz = 0;
  float m_delta_pE = 0;
  bool m_accept_px_1percent = false;
  bool m_accept_py_1percent = false;
  bool m_accept_pz_1percent = false;
  bool m_accept_pE_1percent = false;
  bool m_accept_px_5percent = false;
  bool m_accept_py_5percent = false;
  bool m_accept_pz_5percent = false;
  bool m_accept_pE_5percent = false;
  bool m_accept_px_15percent = false;
  bool m_accept_py_15percent = false;
  bool m_accept_pz_15percent = false;
  bool m_accept_pE_15percent = false;
  bool m_accept_eta = true;
  bool m_accept_pT = true;
};

#endif  // GEANTTESTER_H
