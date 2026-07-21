#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

class IdealPadMap;

#include <mutex>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class TH1D;
class TH2D;

class Tpc_ModuleTrackContainer;
class Tpc_AssembledTrackContainer;
class TrkrHitSetContainer;

class Tpc_AssembledTrackReco : public SubsysReco
{
 public:
  explicit Tpc_AssembledTrackReco(const std::string& name = "Tpc_AssembledTrackReco",
                              const std::string& filename = "Tpc_AssembledTracks.root");
  ~Tpc_AssembledTrackReco() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n)  { m_inputNodeName  = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setDebugOutputFileName(const std::string& n) { m_debugOutputFileName = n; }

  void setConnectMaxLayerGap(unsigned int n)   { m_connectMaxLayerGap = n; }

  // Connection window in global phi radians and tbin units at the match radius.
  void setConnectWindow(double dphi, double dtbin)
  {
    m_connect_dphi  = dphi;
    m_connect_dtbin = dtbin;
  }

  // Slope windows: d(phi)/d(radius) and d(tbin)/d(radius).
  void setConnectSlopeWindow(double dphi_slope, double dtbin_slope)
  {
    m_connect_dphi_slope  = dphi_slope;
    m_connect_dtbin_slope = dtbin_slope;
  }

  void setUseSagittaPhiFit(bool v) { m_useSagittaPhiFit = v; }

  void setSeedCovarianceDiagonal(double sigmaX, double sigmaY, double sigmaZ,
                                 double sigmaPx, double sigmaPy, double sigmaPz)
  {
    m_seedSigmaX = sigmaX;
    m_seedSigmaY = sigmaY;
    m_seedSigmaZ = sigmaZ;
    m_seedSigmaPx = sigmaPx;
    m_seedSigmaPy = sigmaPy;
    m_seedSigmaPz = sigmaPz;
  }

 public:
  // ---------------------------------------------------------------
  // Piece: one Tpc_ModuleTrack refit from stored hit indices for temporary matching.
  // Fits are in ideal (radius, global phi) and (radius, timebin) coordinates.
  // ---------------------------------------------------------------
  struct Piece
  {
    Piece();

    unsigned int source_index;
    unsigned int source_track_id;
    unsigned int event;
    unsigned int region;
    unsigned int sector;
    int          side;

    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int nblobs;
    unsigned int nrawhits;

    double phi_slope;
    double phi_intercept;
    double phi_S;
    double phi_x0;
    double phi_invR;
    double phi_theta;
    double phi_bline;
    bool phi_sagitta_ok;

    double tbin_slope;
    double tbin_intercept;

    std::vector<double> radius_values;
    std::vector<double> phi_values;
    std::vector<double> tbin_values;
    std::vector<double> weights;

    std::vector<TrkrDefs::hitsetkey> hitsetkeys;
    std::vector<TrkrDefs::hitkey>    hitkeys;
  };

  // ---------------------------------------------------------------
  // Candidate: a growing assembled track built from one or more Pieces.
  // Fits are temporary and are not copied into the output Tpc_AssembledTrack object.
  // ---------------------------------------------------------------
  struct SeedParameters
  {
    bool ok {false};
    double x {0.0};
    double y {0.0};
    double z {0.0};
    double px {0.0};
    double py {0.0};
    double pz {0.0};
    double cov[6][6] {};
  };

  struct Candidate
  {
    Candidate();

    unsigned int event;
    int          side;
    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int first_sector;
    unsigned int last_sector;
    unsigned int first_region;
    unsigned int last_region;
    unsigned int nsegments;
    unsigned int nblobs;
    unsigned int nrawhits;

    double phi_slope;
    double phi_intercept;
    double phi_S;
    double phi_x0;
    double phi_invR;
    double phi_theta;
    double phi_bline;
    bool phi_sagitta_ok;

    double tbin_slope_r;
    double tbin_intercept_r;
    double chi2_phi;
    double chi2_tbin;
    int    ndof_phi;
    int    ndof_tbin;

    std::vector<unsigned int>        piece_indices;
    std::vector<TrkrDefs::hitsetkey> hitsetkeys;
    std::vector<TrkrDefs::hitkey>    hitkeys;
  };

 private:
  int  getNodes(PHCompositeNode*);
  int  createNodes(PHCompositeNode*);
  void reset_tree_vars();
  void create_debug_histograms();
  void write_debug_histograms();

  bool make_piece(unsigned int source_index, Piece& p) const;

  double predict_phi(const Piece& p, double radius) const;
  double predict_phi(const Candidate& c, double radius) const;
  double predict_phi_slope(const Piece& p, double radius) const;
  double predict_phi_slope(const Candidate& c, double radius) const;

  // Refit a growing candidate for temporary matching only.
  bool refit_candidate(const std::vector<Piece>& pieces,
                       const std::vector<unsigned int>& piece_indices,
                       Candidate& c) const;

  // Check whether piece b can be appended to candidate a.
  bool candidates_can_connect(const Candidate& a,
                               const Piece&     b,
                               double&          score,
                               double&          b_phi_intercept_shifted) const;
  bool candidates_can_connect(const Candidate& a,
                               const Candidate& b,
                               double&          score) const;

  void connect_sector_pieces(const std::vector<Piece>& pieces,
                              int side,
                              unsigned int sector,
                              std::vector<Candidate>& output) const;
  void connect_side_candidates(const std::vector<Piece>& pieces,
                                const std::vector<Candidate>& seeds,
                                int side,
                                std::vector<Candidate>& output) const;
  SeedParameters make_seed_parameters(const Candidate& c) const;

  // ---------------------------------------------------------------
  std::string m_outputFileName;
  std::string m_debugOutputFileName;
  std::string m_inputNodeName;
  std::string m_outputNodeName;

  TFile* m_outputFile;
  TFile* m_debugOutputFile;
  TTree* m_tree;

  Tpc_ModuleTrackContainer* m_tpcModuleTrackContainer;
  Tpc_AssembledTrackContainer*     m_assembledTrackContainer;
  TrkrHitSetContainer*    m_hits;

  int          m_event;
  IdealPadMap* m_idealPadMap;

  unsigned int m_connectMaxLayerGap;
  double       m_connect_dphi;
  double       m_connect_dtbin;
  double       m_connect_dphi_slope;
  double       m_connect_dtbin_slope;
  bool         m_useSagittaPhiFit;
  double       m_seedSigmaX;
  double       m_seedSigmaY;
  double       m_seedSigmaZ;
  double       m_seedSigmaPx;
  double       m_seedSigmaPy;
  double       m_seedSigmaPz;

  // Debug matching histograms
  TH1D* m_h_dphi;
  TH1D* m_h_dtbin;
  TH1D* m_h_dmphi;
  TH1D* m_h_dmtbin;
  TH1D* m_h_score;
  TH2D* m_h_dphi_vs_dtbin;
  TH2D* m_h_dmphi_vs_dmtbin;
  TH2D* m_h_dphi_vs_dmphi;
  TH2D* m_h_tbin_slope_vs_first_tbin;
  TH2D* m_h_tbin_slope_vs_last_tbin;
  TH2D* m_h_track_tbin_slope_vs_tbin_span_3modules;
  TH2D* m_h_track_tbin_slope_vs_first_tbin_3modules;
  TH2D* m_h_track_tbin_slope_vs_last_tbin_3modules;
  TH1D* m_h_layer_gap;
  TH1D* m_h_nsegments;
  TH1D* m_h_matched_sector_delta;

  mutable std::mutex m_debugMutex;

  // ---------------------------------------------------------------
  // TTree branches
  // ---------------------------------------------------------------
  int                        m_tree_event;
  std::vector<unsigned int>  m_tree_track_id;
  std::vector<int>           m_tree_side;
  std::vector<unsigned int>  m_tree_nsegments;
  std::vector<unsigned int>  m_tree_nblobs;
  std::vector<unsigned int>  m_tree_nrawhits;
  std::vector<unsigned int>  m_tree_first_layer;
  std::vector<unsigned int>  m_tree_last_layer;
  std::vector<unsigned int>  m_tree_first_sector;
  std::vector<unsigned int>  m_tree_last_sector;
  std::vector<unsigned int>  m_tree_first_region;
  std::vector<unsigned int>  m_tree_last_region;

  std::vector<unsigned int>  m_tree_source_assembled_track_id;
  std::vector<unsigned int>  m_tree_source_inmodule_track_id;
  std::vector<unsigned int>  m_tree_source_region;
  std::vector<unsigned int>  m_tree_source_sector;
  std::vector<int>           m_tree_source_side;

  std::vector<unsigned int>       m_tree_hit_assembled_track_id;
  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;
};