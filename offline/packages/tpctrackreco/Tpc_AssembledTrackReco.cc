#include "Tpc_AssembledTrackReco.h"

#include "Tpc_FittingTools.h"
#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackv1.h"
#include "Tpc_AssembledTrackContainer.h"
#include "Tpc_AssembledTrackContainerv1.h"
#include "IdealPadMap.h"

#include "Tpc_ModuleTrack.h"
#include "Tpc_ModuleTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

namespace
{
  struct PieceStartSort
  {
    const std::vector<Tpc_AssembledTrackReco::Piece>* pieces;
    explicit PieceStartSort(const std::vector<Tpc_AssembledTrackReco::Piece>* p) : pieces(p) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const Tpc_AssembledTrackReco::Piece& pa = (*pieces)[a];
      const Tpc_AssembledTrackReco::Piece& pb = (*pieces)[b];
      if (pa.first_layer != pb.first_layer) return pa.first_layer < pb.first_layer;
      if (pa.last_layer != pb.last_layer) return pa.last_layer < pb.last_layer;
      if (pa.sector != pb.sector) return pa.sector < pb.sector;
      return pa.source_track_id < pb.source_track_id;
    }
  };

  struct CandidateStartSort
  {
    const std::vector<Tpc_AssembledTrackReco::Candidate>* candidates;
    explicit CandidateStartSort(const std::vector<Tpc_AssembledTrackReco::Candidate>* c) : candidates(c) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const Tpc_AssembledTrackReco::Candidate& ca = (*candidates)[a];
      const Tpc_AssembledTrackReco::Candidate& cb = (*candidates)[b];
      if (ca.first_layer != cb.first_layer) return ca.first_layer < cb.first_layer;
      if (ca.last_layer != cb.last_layer) return ca.last_layer < cb.last_layer;
      if (ca.first_sector != cb.first_sector) return ca.first_sector < cb.first_sector;
      return ca.nsegments < cb.nsegments;
    }
  };

  struct RadiusSort
  {
    const std::vector<double>* radius;
    explicit RadiusSort(const std::vector<double>* r) : radius(r) {}
    bool operator()(unsigned int a, unsigned int b) const { return (*radius)[a] < (*radius)[b]; }
  };

  double unwrap_phi_to_reference(double phi, const double ref)
  {
    while (phi - ref > M_PI) phi -= 2.0 * M_PI;
    while (phi - ref < -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  double wrap_to_pi(double phi)
  {
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  int wrapped_sector_delta(const unsigned int sector_a, const unsigned int sector_b)
  {
    int d = static_cast<int>(sector_b) - static_cast<int>(sector_a);
    while (d > 6) d -= 12;
    while (d < -6) d += 12;
    return d;
  }

  double sagitta_model_derivative(double xrot, double x0, double invR)
  {
    const double dx = xrot - x0;
    const double dx2 = dx * dx;
    const double invR2 = invR * invR;
    const double invR3 = invR2 * invR;
    const double invR5 = invR3 * invR2;
    return -invR * dx - 0.5 * invR3 * dx2 * dx - 0.375 * invR5 * dx2 * dx2 * dx;
  }

  double predict_sagitta_phi(double radius, double S, double x0, double invR, double theta, double bline)
  {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    double yy = std::tan(theta) * radius;

    for (unsigned int iter = 0; iter < 25; ++iter)
    {
      const double xrot = c * radius + s * yy;
      const double yrot = -s * radius + c * yy;
      const double f = Tpc_FittingTools::sagittaModel(xrot, S, x0, invR);
      const double g = yrot - f;
      const double df = sagitta_model_derivative(xrot, x0, invR);
      const double dg = c - df * s;
      if (std::fabs(dg) < 1.0e-12) break;
      const double step = g / dg;
      yy -= step;
      if (std::fabs(step) < 1.0e-10) break;
    }

    return bline + yy;
  }

  bool fit_points(const std::vector<double>& radius,
                  const std::vector<double>& phi,
                  const std::vector<double>& tbin,
                  const std::vector<double>& weight,
                  bool use_sagitta,
                  double& phi_slope,
                  double& phi_intercept,
                  double& phi_S,
                  double& phi_x0,
                  double& phi_invR,
                  double& phi_theta,
                  double& phi_bline,
                  bool& phi_sagitta_ok,
                  double& tbin_slope,
                  double& tbin_intercept,
                  double& chi2_phi,
                  double& chi2_tbin,
                  int& ndof_phi,
                  int& ndof_tbin)
  {
    if (radius.size() < 2 || radius.size() != phi.size() || radius.size() != tbin.size() || radius.size() != weight.size()) return false;

    std::vector<Tpc_FittingTools::FitPoint> phi_points;
    std::vector<Tpc_FittingTools::FitPoint> tbin_points;
    phi_points.reserve(radius.size());
    tbin_points.reserve(radius.size());
    for (unsigned int i = 0; i < radius.size(); ++i)
    {
      phi_points.push_back(Tpc_FittingTools::FitPoint(radius[i], phi[i], weight[i]));
      tbin_points.push_back(Tpc_FittingTools::FitPoint(radius[i], tbin[i], weight[i]));
    }

    const Tpc_FittingTools::LineFit phi_line = Tpc_FittingTools::fitLine(phi_points);
    if (!phi_line.ok) return false;

    phi_slope = phi_line.slope;
    phi_intercept = phi_line.intercept;
    chi2_phi = phi_line.chi2;
    ndof_phi = phi_line.ndof;
    phi_S = 0.0;
    phi_x0 = 0.0;
    phi_invR = 0.0;
    phi_theta = std::atan(phi_slope);
    phi_bline = phi_intercept;
    phi_sagitta_ok = false;

    if (use_sagitta && radius.size() >= 3)
    {
      const Tpc_FittingTools::SagittaFit phi_sagitta = Tpc_FittingTools::fitSagitta(phi_points);
      if (phi_sagitta.ok)
      {
        phi_S = phi_sagitta.S;
        phi_x0 = phi_sagitta.x0;
        phi_invR = phi_sagitta.invR;
        phi_theta = phi_sagitta.theta;
        phi_bline = phi_sagitta.b;
        chi2_phi = phi_sagitta.chi2;
        ndof_phi = phi_sagitta.ndof;
        phi_sagitta_ok = true;
      }
    }

    const Tpc_FittingTools::LineFit tbin_line = Tpc_FittingTools::fitLine(tbin_points);
    if (!tbin_line.ok) return false;

    tbin_slope = tbin_line.slope;
    tbin_intercept = tbin_line.intercept;
    chi2_tbin = tbin_line.chi2;
    ndof_tbin = tbin_line.ndof;
    return true;
  }
}

Tpc_AssembledTrackReco::Piece::Piece()
  : source_index(0)
  , source_track_id(0)
  , event(0)
  , region(0)
  , sector(0)
  , side(0)
  , first_layer(0)
  , last_layer(0)
  , nblobs(0)
  , nrawhits(0)
  , phi_slope(0.0)
  , phi_intercept(0.0)
  , phi_S(0.0)
  , phi_x0(0.0)
  , phi_invR(0.0)
  , phi_theta(0.0)
  , phi_bline(0.0)
  , phi_sagitta_ok(false)
  , tbin_slope(0.0)
  , tbin_intercept(0.0)
{
}

Tpc_AssembledTrackReco::Candidate::Candidate()
  : event(0)
  , side(0)
  , first_layer(0)
  , last_layer(0)
  , first_sector(0)
  , last_sector(0)
  , first_region(0)
  , last_region(0)
  , nsegments(0)
  , nblobs(0)
  , nrawhits(0)
  , phi_slope(0.0)
  , phi_intercept(0.0)
  , phi_S(0.0)
  , phi_x0(0.0)
  , phi_invR(0.0)
  , phi_theta(0.0)
  , phi_bline(0.0)
  , phi_sagitta_ok(false)
  , tbin_slope_r(0.0)
  , tbin_intercept_r(0.0)
  , chi2_phi(0.0)
  , chi2_tbin(0.0)
  , ndof_phi(0)
  , ndof_tbin(0)
{
}

Tpc_AssembledTrackReco::Tpc_AssembledTrackReco(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_debugOutputFileName("Tpc_AssembledTrackRecoDebug.root")
  , m_inputNodeName("TPC_MODULETRACKS")
  , m_outputNodeName("TPC_ASSEMBLEDTRACKS")
  , m_outputFile(nullptr)
  , m_debugOutputFile(nullptr)
  , m_tree(nullptr)
  , m_tpcModuleTrackContainer(nullptr)
  , m_assembledTrackContainer(nullptr)
  , m_hits(nullptr)
  , m_event(0)
  , m_idealPadMap(nullptr)
  , m_connectMaxLayerGap(16)
  , m_connect_dphi(0.03)
  , m_connect_dtbin(8.0)
  , m_connect_dphi_slope(0.01)
  , m_connect_dtbin_slope(2.0)
  , m_useSagittaPhiFit(true)
  , m_seedSigmaX(5.0)
  , m_seedSigmaY(5.0)
  , m_seedSigmaZ(10.0)
  , m_seedSigmaPx(1.0)
  , m_seedSigmaPy(1.0)
  , m_seedSigmaPz(1.0)
  , m_h_dphi(nullptr)
  , m_h_dtbin(nullptr)
  , m_h_dmphi(nullptr)
  , m_h_dmtbin(nullptr)
  , m_h_score(nullptr)
  , m_h_dphi_vs_dtbin(nullptr)
  , m_h_dmphi_vs_dmtbin(nullptr)
  , m_h_dphi_vs_dmphi(nullptr)
  , m_h_tbin_slope_vs_first_tbin(nullptr)
  , m_h_tbin_slope_vs_last_tbin(nullptr)
  , m_h_track_tbin_slope_vs_tbin_span_3modules(nullptr)
  , m_h_track_tbin_slope_vs_first_tbin_3modules(nullptr)
  , m_h_track_tbin_slope_vs_last_tbin_3modules(nullptr)
  , m_h_layer_gap(nullptr)
  , m_h_nsegments(nullptr)
  , m_h_matched_sector_delta(nullptr)
{
}

Tpc_AssembledTrackReco::~Tpc_AssembledTrackReco()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;

  if (m_debugOutputFile)
  {
    delete m_debugOutputFile;
    m_debugOutputFile = nullptr;
  }

  if (m_outputFile)
  {
    delete m_outputFile;
    m_outputFile = nullptr;
  }
}

int Tpc_AssembledTrackReco::Init(PHCompositeNode*)
{
  if (Verbosity() <= 0) return Fun4AllReturnCodes::EVENT_OK;

  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_debugOutputFile = new TFile(m_debugOutputFileName.c_str(), "RECREATE");
  if (!m_debugOutputFile || m_debugOutputFile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot create debug file " << m_debugOutputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  create_debug_histograms();

  m_tree = new TTree("Tpc_AssembledTracks", "Assembled tracks connected from Tpc_ModuleTrackReco");
  m_tree->Branch("event", &m_tree_event, "event/I");
  m_tree->Branch("track_id", &m_tree_track_id);
  m_tree->Branch("side", &m_tree_side);
  m_tree->Branch("nsegments", &m_tree_nsegments);
  m_tree->Branch("nblobs", &m_tree_nblobs);
  m_tree->Branch("nrawhits", &m_tree_nrawhits);
  m_tree->Branch("first_layer", &m_tree_first_layer);
  m_tree->Branch("last_layer", &m_tree_last_layer);
  m_tree->Branch("first_sector", &m_tree_first_sector);
  m_tree->Branch("last_sector", &m_tree_last_sector);
  m_tree->Branch("first_region", &m_tree_first_region);
  m_tree->Branch("last_region", &m_tree_last_region);

  m_tree->Branch("source_assembled_track_id", &m_tree_source_assembled_track_id);
  m_tree->Branch("source_inmodule_track_id", &m_tree_source_inmodule_track_id);
  m_tree->Branch("source_region", &m_tree_source_region);
  m_tree->Branch("source_sector", &m_tree_source_sector);
  m_tree->Branch("source_side", &m_tree_source_side);

  m_tree->Branch("hit_assembled_track_id", &m_tree_hit_assembled_track_id);
  m_tree->Branch("hit_hitsetkey", &m_tree_hit_hitsetkey);
  m_tree->Branch("hit_hitkey", &m_tree_hit_hitkey);

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - cannot load IdealPadMap from CDB" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::End(PHCompositeNode*)
{
  if (m_outputFile)
  {
    m_outputFile->cd();
    if (m_tree) m_tree->Write();
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = nullptr;
  }

  write_debug_histograms();
  if (m_debugOutputFile)
  {
    m_debugOutputFile->Close();
    delete m_debugOutputFile;
    m_debugOutputFile = nullptr;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::getNodes(PHCompositeNode* topNode)
{
  m_tpcModuleTrackContainer = findNode::getClass<Tpc_ModuleTrackContainer>(topNode, m_inputNodeName);
  if (!m_tpcModuleTrackContainer)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_assembledTrackContainer = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_outputNodeName);
  if (!m_assembledTrackContainer)
  {
    m_assembledTrackContainer = new Tpc_AssembledTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_assembledTrackContainer, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void Tpc_AssembledTrackReco::create_debug_histograms()
{
  if (!m_debugOutputFile) return;
  m_debugOutputFile->cd();

  m_h_dphi = new TH1D("h_dphi", Form("#Delta#phi at match point, cut=%.4g;|#Delta#phi| [rad];tested pairs", m_connect_dphi), 200, 0.0, std::max(0.2, 5.0 * m_connect_dphi));
  m_h_dtbin = new TH1D("h_dtbin", Form("#Deltatbin at match point, cut=%.4g;|#Deltatbin|;tested pairs", m_connect_dtbin), 200, 0.0, std::max(50.0, 5.0 * m_connect_dtbin));
  m_h_dmphi = new TH1D("h_dmphi", Form("#Delta(d#phi/dr), cut=%.4g;|#Delta(d#phi/dr)| [rad/cm];tested pairs", m_connect_dphi_slope), 200, 0.0, std::max(0.08, 5.0 * m_connect_dphi_slope));
  m_h_dmtbin = new TH1D("h_dmtbin", Form("#Delta(dtbin/dr), cut=%.4g;|#Delta(dtbin/dr)| [tbin/cm];tested pairs", m_connect_dtbin_slope), 200, 0.0, std::max(20.0, 5.0 * m_connect_dtbin_slope));
  m_h_score = new TH1D("h_score", "accepted connection score;score;accepted connections", 200, 0.0, 20.0);

  m_h_dphi_vs_dtbin = new TH2D("h_dphi_vs_dtbin", Form("#Delta#phi vs #Deltatbin, cuts #Delta#phi<%.4g #Deltatbin<%.4g;|#Delta#phi| [rad];|#Deltatbin|", m_connect_dphi, m_connect_dtbin), 160, 0.0, std::max(0.2, 5.0 * m_connect_dphi), 160, 0.0, std::max(50.0, 5.0 * m_connect_dtbin));
  m_h_dmphi_vs_dmtbin = new TH2D("h_dmphi_vs_dmtbin", Form("slope residuals, cuts #Delta(d#phi/dr)<%.4g #Delta(dtbin/dr)<%.4g;|#Delta(d#phi/dr)| [rad/cm];|#Delta(dtbin/dr)| [tbin/cm]", m_connect_dphi_slope, m_connect_dtbin_slope), 160, 0.0, std::max(0.08, 5.0 * m_connect_dphi_slope), 160, 0.0, std::max(20.0, 5.0 * m_connect_dtbin_slope));
  m_h_dphi_vs_dmphi = new TH2D("h_dphi_vs_dmphi", Form("#phi position vs slope residual, cuts #Delta#phi<%.4g #Delta(d#phi/dr)<%.4g;|#Delta#phi| [rad];|#Delta(d#phi/dr)| [rad/cm]", m_connect_dphi, m_connect_dphi_slope), 160, 0.0, std::max(0.2, 5.0 * m_connect_dphi), 160, 0.0, std::max(0.08, 5.0 * m_connect_dphi_slope));

  m_h_tbin_slope_vs_first_tbin = new TH2D("h_tbin_slope_vs_first_tbin", "tested connection dtbin/dr vs first timebin;first timebin;dtbin/dr [tbin/cm]", 200, 0.0, 600.0, 200, -20.0, 20.0);
  m_h_tbin_slope_vs_last_tbin = new TH2D("h_tbin_slope_vs_last_tbin", "tested connection dtbin/dr vs last timebin;last timebin;dtbin/dr [tbin/cm]", 200, 0.0, 600.0, 200, -20.0, 20.0);
  m_h_track_tbin_slope_vs_tbin_span_3modules = new TH2D("h_track_tbin_slope_vs_tbin_span_3modules", "3-module tracks dtbin/dr vs last-first timebin;last timebin - first timebin;dtbin/dr [tbin/cm]", 200, -600.0, 600.0, 200, -20.0, 20.0);
  m_h_track_tbin_slope_vs_first_tbin_3modules = new TH2D("h_track_tbin_slope_vs_first_tbin_3modules", "3-module tracks dtbin/dr vs first timebin;first timebin;dtbin/dr [tbin/cm]", 200, 0.0, 600.0, 200, -20.0, 20.0);
  m_h_track_tbin_slope_vs_last_tbin_3modules = new TH2D("h_track_tbin_slope_vs_last_tbin_3modules", "3-module tracks dtbin/dr vs last timebin;last timebin;dtbin/dr [tbin/cm]", 200, 0.0, 600.0, 200, -20.0, 20.0);

  m_h_layer_gap = new TH1D("h_layer_gap", "accepted connection layer gap;b.first_layer - a.last_layer - 1;accepted connections", 16, -0.5, 15.5);
  m_h_nsegments = new TH1D("h_nsegments", "pieces per assembled track;nsegments;assembled tracks", 16, -0.5, 15.5);
  m_h_matched_sector_delta = new TH1D("h_matched_sector_delta", "accepted matched sector difference;wrapped #Delta sector;accepted connections", 25, -12.5, 12.5);
}

void Tpc_AssembledTrackReco::write_debug_histograms()
{
  if (!m_debugOutputFile) return;
  m_debugOutputFile->cd();

  if (m_h_dphi) m_h_dphi->Write();
  if (m_h_dtbin) m_h_dtbin->Write();
  if (m_h_dmphi) m_h_dmphi->Write();
  if (m_h_dmtbin) m_h_dmtbin->Write();
  if (m_h_score) m_h_score->Write();
  if (m_h_dphi_vs_dtbin) m_h_dphi_vs_dtbin->Write();
  if (m_h_dmphi_vs_dmtbin) m_h_dmphi_vs_dmtbin->Write();
  if (m_h_dphi_vs_dmphi) m_h_dphi_vs_dmphi->Write();
  if (m_h_tbin_slope_vs_first_tbin) m_h_tbin_slope_vs_first_tbin->Write();
  if (m_h_tbin_slope_vs_last_tbin) m_h_tbin_slope_vs_last_tbin->Write();
  if (m_h_track_tbin_slope_vs_tbin_span_3modules) m_h_track_tbin_slope_vs_tbin_span_3modules->Write();
  if (m_h_track_tbin_slope_vs_first_tbin_3modules) m_h_track_tbin_slope_vs_first_tbin_3modules->Write();
  if (m_h_track_tbin_slope_vs_last_tbin_3modules) m_h_track_tbin_slope_vs_last_tbin_3modules->Write();
  if (m_h_layer_gap) m_h_layer_gap->Write();
  if (m_h_nsegments) m_h_nsegments->Write();
  if (m_h_matched_sector_delta) m_h_matched_sector_delta->Write();

}

void Tpc_AssembledTrackReco::reset_tree_vars()
{
  m_tree_event = m_event;
  m_tree_track_id.clear();
  m_tree_side.clear();
  m_tree_nsegments.clear();
  m_tree_nblobs.clear();
  m_tree_nrawhits.clear();
  m_tree_first_layer.clear();
  m_tree_last_layer.clear();
  m_tree_first_sector.clear();
  m_tree_last_sector.clear();
  m_tree_first_region.clear();
  m_tree_last_region.clear();
  m_tree_source_assembled_track_id.clear();
  m_tree_source_inmodule_track_id.clear();
  m_tree_source_region.clear();
  m_tree_source_sector.clear();
  m_tree_source_side.clear();
  m_tree_hit_assembled_track_id.clear();
  m_tree_hit_hitsetkey.clear();
  m_tree_hit_hitkey.clear();
}

bool Tpc_AssembledTrackReco::make_piece(unsigned int source_index, Piece& p) const
{
  const Tpc_ModuleTrack* trk = m_tpcModuleTrackContainer->get_track(source_index);
  if (!trk || !trk->isValid()) return false;
  if (trk->get_last_layer() < trk->get_first_layer()) return false;
  if (!m_idealPadMap || !m_idealPadMap->is_loaded() || !m_hits) return false;

  p.source_index = source_index;
  p.source_track_id = trk->get_track_id();
  p.event = trk->get_event();
  p.region = trk->get_region();
  p.sector = trk->get_sector();
  p.side = trk->get_side();
  p.first_layer = trk->get_first_layer();
  p.last_layer = trk->get_last_layer();
  p.nblobs = trk->get_nblobs();
  p.nrawhits = trk->get_nrawhits();

  p.radius_values.clear();
  p.phi_values.clear();
  p.tbin_values.clear();
  p.weights.clear();
  p.hitsetkeys.clear();
  p.hitkeys.clear();

  double maxadc = 0.0;
  for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
  {
    const Tpc_ModuleTrack::HitIndex hi = trk->get_hit_index(ih);
    TrkrHitSet* hitset = m_hits->findHitSet(hi.first);
    if (!hitset) continue;
    TrkrHit* hit = hitset->getHit(hi.second);
    if (!hit) continue;
    const double adc = static_cast<double>(hit->getAdc());
    if (adc > maxadc) maxadc = adc;
  }

  for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
  {
    const Tpc_ModuleTrack::HitIndex hi = trk->get_hit_index(ih);
    TrkrHitSet* hitset = m_hits->findHitSet(hi.first);
    if (!hitset) continue;
    TrkrHit* hit = hitset->getHit(hi.second);
    if (!hit) continue;

    const unsigned int layer = TrkrDefs::getLayer(hi.first);
    const unsigned int pad = TpcDefs::getPad(hi.second);
    const unsigned int tbin = TpcDefs::getTBin(hi.second);

    const double radius = m_idealPadMap->get_radius(layer);
    const double phi = wrap_to_pi(m_idealPadMap->get_phi(static_cast<unsigned int>(p.side), layer, pad));
    if (!std::isfinite(radius) || !std::isfinite(phi)) continue;

    p.radius_values.push_back(radius);
    p.phi_values.push_back(phi);
    p.tbin_values.push_back(static_cast<double>(tbin));
    p.weights.push_back(Tpc_FittingTools::adcWeight(static_cast<double>(hit->getAdc()), maxadc, 0.5, 0.15));
    p.hitsetkeys.push_back(hi.first);
    p.hitkeys.push_back(hi.second);
  }

  if (p.radius_values.size() < 2) return false;

  std::vector<unsigned int> order;
  order.reserve(p.radius_values.size());
  for (unsigned int i = 0; i < p.radius_values.size(); ++i) order.push_back(i);
  std::sort(order.begin(), order.end(), RadiusSort(&p.radius_values));

  std::vector<double> r_sorted, phi_sorted, tbin_sorted, w_sorted;
  std::vector<TrkrDefs::hitsetkey> hsk_sorted;
  std::vector<TrkrDefs::hitkey> hk_sorted;
  r_sorted.reserve(order.size());
  phi_sorted.reserve(order.size());
  tbin_sorted.reserve(order.size());
  w_sorted.reserve(order.size());
  hsk_sorted.reserve(order.size());
  hk_sorted.reserve(order.size());

  for (unsigned int io = 0; io < order.size(); ++io)
  {
    const unsigned int i = order[io];
    double phi = p.phi_values[i];
    if (!phi_sorted.empty()) phi = unwrap_phi_to_reference(phi, phi_sorted.back());
    r_sorted.push_back(p.radius_values[i]);
    phi_sorted.push_back(phi);
    tbin_sorted.push_back(p.tbin_values[i]);
    w_sorted.push_back(p.weights[i]);
    hsk_sorted.push_back(p.hitsetkeys[i]);
    hk_sorted.push_back(p.hitkeys[i]);
  }

  p.radius_values.swap(r_sorted);
  p.phi_values.swap(phi_sorted);
  p.tbin_values.swap(tbin_sorted);
  p.weights.swap(w_sorted);
  p.hitsetkeys.swap(hsk_sorted);
  p.hitkeys.swap(hk_sorted);

  double chi2_phi = 0.0;
  double chi2_tbin = 0.0;
  int ndof_phi = 0;
  int ndof_tbin = 0;
  return fit_points(p.radius_values, p.phi_values, p.tbin_values, p.weights, m_useSagittaPhiFit,
                    p.phi_slope, p.phi_intercept, p.phi_S, p.phi_x0, p.phi_invR,
                    p.phi_theta, p.phi_bline, p.phi_sagitta_ok,
                    p.tbin_slope, p.tbin_intercept, chi2_phi, chi2_tbin, ndof_phi, ndof_tbin);
}

double Tpc_AssembledTrackReco::predict_phi(const Piece& p, double radius) const
{
  if (m_useSagittaPhiFit && p.phi_sagitta_ok) return predict_sagitta_phi(radius, p.phi_S, p.phi_x0, p.phi_invR, p.phi_theta, p.phi_bline);
  return p.phi_slope * radius + p.phi_intercept;
}

double Tpc_AssembledTrackReco::predict_phi(const Candidate& c, double radius) const
{
  if (m_useSagittaPhiFit && c.phi_sagitta_ok) return predict_sagitta_phi(radius, c.phi_S, c.phi_x0, c.phi_invR, c.phi_theta, c.phi_bline);
  return c.phi_slope * radius + c.phi_intercept;
}

double Tpc_AssembledTrackReco::predict_phi_slope(const Piece& p, double radius) const
{
  if (!(m_useSagittaPhiFit && p.phi_sagitta_ok)) return p.phi_slope;
  const double eps = 1.0e-3;
  return (predict_phi(p, radius + eps) - predict_phi(p, radius - eps)) / (2.0 * eps);
}

double Tpc_AssembledTrackReco::predict_phi_slope(const Candidate& c, double radius) const
{
  if (!(m_useSagittaPhiFit && c.phi_sagitta_ok)) return c.phi_slope;
  const double eps = 1.0e-3;
  return (predict_phi(c, radius + eps) - predict_phi(c, radius - eps)) / (2.0 * eps);
}

bool Tpc_AssembledTrackReco::refit_candidate(const std::vector<Piece>& pieces, const std::vector<unsigned int>& piece_indices, Candidate& c) const
{
  if (piece_indices.empty()) return false;

  std::vector<double> radius, phi, tbin, weight;
  c = Candidate();
  c.piece_indices = piece_indices;

  for (unsigned int ii = 0; ii < piece_indices.size(); ++ii)
  {
    const Piece& p = pieces[piece_indices[ii]];

    for (unsigned int ih = 0; ih < p.radius_values.size(); ++ih)
    {
      double phiv = p.phi_values[ih];
      if (!phi.empty()) phiv = unwrap_phi_to_reference(phiv, phi.back());
      radius.push_back(p.radius_values[ih]);
      phi.push_back(phiv);
      tbin.push_back(p.tbin_values[ih]);
      weight.push_back(p.weights[ih]);
    }

    if (ii == 0)
    {
      c.event = p.event;
      c.side = p.side;
      c.first_layer = p.first_layer;
      c.last_layer = p.last_layer;
      c.first_sector = p.sector;
      c.last_sector = p.sector;
      c.first_region = p.region;
      c.last_region = p.region;
    }
    else
    {
      if (p.first_layer < c.first_layer)
      {
        c.first_layer = p.first_layer;
        c.first_sector = p.sector;
        c.first_region = p.region;
      }
      if (p.last_layer > c.last_layer)
      {
        c.last_layer = p.last_layer;
        c.last_sector = p.sector;
        c.last_region = p.region;
      }
    }

    c.nblobs += p.nblobs;
    c.nrawhits += p.nrawhits;
    for (unsigned int ih = 0; ih < p.hitsetkeys.size(); ++ih)
    {
      c.hitsetkeys.push_back(p.hitsetkeys[ih]);
      c.hitkeys.push_back(p.hitkeys[ih]);
    }
  }

  c.nsegments = static_cast<unsigned int>(piece_indices.size());
  return fit_points(radius, phi, tbin, weight, m_useSagittaPhiFit,
                    c.phi_slope, c.phi_intercept, c.phi_S, c.phi_x0, c.phi_invR,
                    c.phi_theta, c.phi_bline, c.phi_sagitta_ok,
                    c.tbin_slope_r, c.tbin_intercept_r,
                    c.chi2_phi, c.chi2_tbin, c.ndof_phi, c.ndof_tbin);
}

bool Tpc_AssembledTrackReco::candidates_can_connect(const Candidate& a, const Piece& b, double& score, double& b_phi_intercept_shifted) const
{
  score = std::numeric_limits<double>::max();
  b_phi_intercept_shifted = b.phi_intercept;

  if (a.side != b.side) return false;
  if (a.last_layer >= b.first_layer) return false;

  const unsigned int gap = b.first_layer - a.last_layer - 1;
  if (gap > m_connectMaxLayerGap) return false;

  const double ra = m_idealPadMap->get_radius(a.last_layer);
  const double rb = m_idealPadMap->get_radius(b.first_layer);
  if (!std::isfinite(ra) || !std::isfinite(rb) || ra <= 0.0 || rb <= 0.0) return false;

  const double rmatch = 0.5 * (ra + rb);
  const double phi_a = predict_phi(a, rmatch);
  const double phi_b_raw = predict_phi(b, rmatch);
  const double phi_b = unwrap_phi_to_reference(phi_b_raw, phi_a);
  b_phi_intercept_shifted = b.phi_intercept + (phi_b - phi_b_raw);

  const double tbin_a = a.tbin_slope_r * rmatch + a.tbin_intercept_r;
  const double tbin_b = b.tbin_slope * rmatch + b.tbin_intercept;
  const double dphi = std::fabs(phi_a - phi_b);
  const double dtbin = std::fabs(tbin_a - tbin_b);
  const double dmphi = std::fabs(predict_phi_slope(a, rmatch) - predict_phi_slope(b, rmatch));
  const double dmtbin = std::fabs(a.tbin_slope_r - b.tbin_slope);

  {
    std::lock_guard<std::mutex> lock(m_debugMutex);
    if (m_h_dphi) m_h_dphi->Fill(dphi);
    if (m_h_dtbin) m_h_dtbin->Fill(dtbin);
    if (m_h_dmphi) m_h_dmphi->Fill(dmphi);
    if (m_h_dmtbin) m_h_dmtbin->Fill(dmtbin);
    if (m_h_dphi_vs_dtbin) m_h_dphi_vs_dtbin->Fill(dphi, dtbin);
    if (m_h_dmphi_vs_dmtbin) m_h_dmphi_vs_dmtbin->Fill(dmphi, dmtbin);
    if (m_h_dphi_vs_dmphi) m_h_dphi_vs_dmphi->Fill(dphi, dmphi);
    if (m_h_tbin_slope_vs_last_tbin) m_h_tbin_slope_vs_last_tbin->Fill(a.tbin_slope_r * ra + a.tbin_intercept_r, a.tbin_slope_r);
    if (m_h_tbin_slope_vs_first_tbin && !b.tbin_values.empty()) m_h_tbin_slope_vs_first_tbin->Fill(b.tbin_values.front(), b.tbin_slope);
  }

  if (dphi > m_connect_dphi) return false;
  if (dtbin > m_connect_dtbin) return false;
  if (dmphi > m_connect_dphi_slope) return false;
  if (dmtbin > m_connect_dtbin_slope) return false;

  constexpr double w_phi   = 1.0;   // φ position — highest weight
  constexpr double w_mphi  = 1.0;   // φ slope
  constexpr double w_tbin  = 1.0;   // tbin position — softer penalty
  constexpr double w_mtbin = 2.0;   // tbin slope

  score = w_phi   * (dphi / m_connect_dphi) * (dphi / m_connect_dphi)
        + w_tbin  * (dtbin / m_connect_dtbin) * (dtbin / m_connect_dtbin)
        + w_mphi  * (dmphi / m_connect_dphi_slope) * (dmphi / m_connect_dphi_slope)
        + w_mtbin * (dmtbin / m_connect_dtbin_slope) * (dmtbin / m_connect_dtbin_slope)
        + 0.05 * static_cast<double>(gap);

  return true;
}

void Tpc_AssembledTrackReco::connect_sector_pieces(const std::vector<Piece>& pieces, int side, unsigned int sector, std::vector<Candidate>& output) const
{
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < pieces.size(); ++i)
  {
    if (pieces[i].side == side && pieces[i].sector == sector) order.push_back(i);
  }
  if (order.empty()) return;

  std::sort(order.begin(), order.end(), PieceStartSort(&pieces));
  std::vector<int> used(pieces.size(), 0);

  for (unsigned int io = 0; io < order.size(); ++io)
  {
    const unsigned int iseed = order[io];
    if (used[iseed]) continue;

    std::vector<unsigned int> current_indices;
    current_indices.push_back(iseed);
    used[iseed] = 1;

    Candidate current;
    if (!refit_candidate(pieces, current_indices, current)) continue;

    bool merged_any = true;
    while (merged_any)
    {
      merged_any = false;
      int best_j = -1;
      double best_score = std::numeric_limits<double>::max();

      for (unsigned int jo = 0; jo < order.size(); ++jo)
      {
        const unsigned int j = order[jo];
        if (used[j]) continue;

        double score = 0.0;
        double shifted_intercept = 0.0;
        if (!candidates_can_connect(current, pieces[j], score, shifted_intercept)) continue;

        if (score < best_score)
        {
          best_score = score;
          best_j = static_cast<int>(j);
        }
      }

      if (best_j >= 0)
      {
        std::vector<unsigned int> trial_indices = current_indices;
        trial_indices.push_back(static_cast<unsigned int>(best_j));

        Candidate refit;
        if (refit_candidate(pieces, trial_indices, refit))
        {
          const Piece& accepted_piece = pieces[static_cast<unsigned int>(best_j)];
          const unsigned int accepted_gap = accepted_piece.first_layer - current.last_layer - 1;
          {
            std::lock_guard<std::mutex> lock(m_debugMutex);
            if (m_h_score) m_h_score->Fill(best_score);
            if (m_h_layer_gap) m_h_layer_gap->Fill(static_cast<double>(accepted_gap));
            if (m_h_matched_sector_delta)
            {
              m_h_matched_sector_delta->Fill(static_cast<double>(wrapped_sector_delta(current.last_sector, accepted_piece.sector)));
            }
          }

          current = refit;
          current_indices.swap(trial_indices);
          used[best_j] = 1;
          merged_any = true;
        }
      }
    }

    output.push_back(current);
  }
}

bool Tpc_AssembledTrackReco::candidates_can_connect(const Candidate& a, const Candidate& b, double& score) const
{
  score = std::numeric_limits<double>::max();

  if (a.side != b.side) return false;
  if (a.last_layer >= b.first_layer) return false;
  if (a.last_sector == b.first_sector) return false;

  const unsigned int gap = b.first_layer - a.last_layer - 1;
  if (gap > m_connectMaxLayerGap) return false;

  const double ra = m_idealPadMap->get_radius(a.last_layer);
  const double rb = m_idealPadMap->get_radius(b.first_layer);
  if (!std::isfinite(ra) || !std::isfinite(rb) || ra <= 0.0 || rb <= 0.0) return false;

  const double rmatch = 0.5 * (ra + rb);
  const double phi_a = predict_phi(a, rmatch);
  const double phi_b = unwrap_phi_to_reference(predict_phi(b, rmatch), phi_a);
  const double tbin_a = a.tbin_slope_r * rmatch + a.tbin_intercept_r;
  const double tbin_b = b.tbin_slope_r * rmatch + b.tbin_intercept_r;
  const double dphi = std::fabs(phi_a - phi_b);
  const double dtbin = std::fabs(tbin_a - tbin_b);
  const double dmphi = std::fabs(predict_phi_slope(a, rmatch) - predict_phi_slope(b, rmatch));
  const double dmtbin = std::fabs(a.tbin_slope_r - b.tbin_slope_r);

  {
    std::lock_guard<std::mutex> lock(m_debugMutex);
    if (m_h_dphi) m_h_dphi->Fill(dphi);
    if (m_h_dtbin) m_h_dtbin->Fill(dtbin);
    if (m_h_dmphi) m_h_dmphi->Fill(dmphi);
    if (m_h_dmtbin) m_h_dmtbin->Fill(dmtbin);
    if (m_h_dphi_vs_dtbin) m_h_dphi_vs_dtbin->Fill(dphi, dtbin);
    if (m_h_dmphi_vs_dmtbin) m_h_dmphi_vs_dmtbin->Fill(dmphi, dmtbin);
    if (m_h_dphi_vs_dmphi) m_h_dphi_vs_dmphi->Fill(dphi, dmphi);
    if (m_h_tbin_slope_vs_last_tbin) m_h_tbin_slope_vs_last_tbin->Fill(a.tbin_slope_r * ra + a.tbin_intercept_r, a.tbin_slope_r);
    if (m_h_tbin_slope_vs_first_tbin) m_h_tbin_slope_vs_first_tbin->Fill(b.tbin_slope_r * rb + b.tbin_intercept_r, b.tbin_slope_r);
  }

  if (dphi > m_connect_dphi) return false;
  if (dtbin > m_connect_dtbin) return false;
  if (dmphi > m_connect_dphi_slope) return false;
  if (dmtbin > m_connect_dtbin_slope) return false;

  constexpr double w_phi   = 1.0;
  constexpr double w_mphi  = 1.0;
  constexpr double w_tbin  = 1.0;
  constexpr double w_mtbin = 2.0;

  score = w_phi   * (dphi / m_connect_dphi) * (dphi / m_connect_dphi)
        + w_tbin  * (dtbin / m_connect_dtbin) * (dtbin / m_connect_dtbin)
        + w_mphi  * (dmphi / m_connect_dphi_slope) * (dmphi / m_connect_dphi_slope)
        + w_mtbin * (dmtbin / m_connect_dtbin_slope) * (dmtbin / m_connect_dtbin_slope)
        + 0.05 * static_cast<double>(gap);

  return true;
}

void Tpc_AssembledTrackReco::connect_side_candidates(const std::vector<Piece>& pieces,
                                                 const std::vector<Candidate>& seeds,
                                                 int side,
                                                 std::vector<Candidate>& output) const
{
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < seeds.size(); ++i)
  {
    if (seeds[i].side == side) order.push_back(i);
  }
  if (order.empty()) return;

  std::sort(order.begin(), order.end(), CandidateStartSort(&seeds));
  std::vector<int> used(seeds.size(), 0);

  for (unsigned int io = 0; io < order.size(); ++io)
  {
    const unsigned int iseed = order[io];
    if (used[iseed]) continue;

    std::vector<unsigned int> current_indices = seeds[iseed].piece_indices;
    used[iseed] = 1;

    Candidate current;
    if (!refit_candidate(pieces, current_indices, current)) continue;

    bool merged_any = true;
    while (merged_any)
    {
      merged_any = false;
      int best_j = -1;
      double best_score = std::numeric_limits<double>::max();

      for (unsigned int jo = 0; jo < order.size(); ++jo)
      {
        const unsigned int j = order[jo];
        if (used[j]) continue;

        double score = 0.0;
        if (!candidates_can_connect(current, seeds[j], score)) continue;

        if (score < best_score)
        {
          best_score = score;
          best_j = static_cast<int>(j);
        }
      }

      if (best_j >= 0)
      {
        std::vector<unsigned int> trial_indices = current_indices;
        const Candidate& accepted_seed = seeds[static_cast<unsigned int>(best_j)];
        trial_indices.insert(trial_indices.end(), accepted_seed.piece_indices.begin(), accepted_seed.piece_indices.end());

        Candidate refit;
        if (refit_candidate(pieces, trial_indices, refit))
        {
          const unsigned int accepted_gap = accepted_seed.first_layer - current.last_layer - 1;
          {
            std::lock_guard<std::mutex> lock(m_debugMutex);
            if (m_h_score) m_h_score->Fill(best_score);
            if (m_h_layer_gap) m_h_layer_gap->Fill(static_cast<double>(accepted_gap));
            if (m_h_matched_sector_delta)
            {
              m_h_matched_sector_delta->Fill(static_cast<double>(wrapped_sector_delta(current.last_sector, accepted_seed.first_sector)));
            }
          }

          current = refit;
          current_indices.swap(trial_indices);
          used[best_j] = 1;
          merged_any = true;
        }
      }
    }

    output.push_back(current);
  }
}


Tpc_AssembledTrackReco::SeedParameters
Tpc_AssembledTrackReco::make_seed_parameters(const Candidate& c) const
{
  SeedParameters seed;
  if (!m_idealPadMap) return seed;

  const double r_first = m_idealPadMap->get_radius(c.first_layer);
  const double r_last = m_idealPadMap->get_radius(c.last_layer);
  if (!std::isfinite(r_first) || !std::isfinite(r_last)) return seed;

  const double radius = 0.5 * (r_first + r_last);
  const double phi = wrap_to_pi(predict_phi(c, radius));
  const double dphi_dr = predict_phi_slope(c, radius);
  const double dtbin_dr = c.tbin_slope_r;
  const double tbin = c.tbin_slope_r * radius + c.tbin_intercept_r;
  if (!std::isfinite(radius) || !std::isfinite(phi) ||
      !std::isfinite(dphi_dr) || !std::isfinite(dtbin_dr) || !std::isfinite(tbin))
  {
    return seed;
  }

  seed.x = radius * std::cos(phi);
  seed.y = radius * std::sin(phi);
  seed.z = tbin;

  const double dx_dr = std::cos(phi) - radius * std::sin(phi) * dphi_dr;
  const double dy_dr = std::sin(phi) + radius * std::cos(phi) * dphi_dr;
  const double dz_dr = dtbin_dr;
  const double norm = std::sqrt(dx_dr * dx_dr + dy_dr * dy_dr + dz_dr * dz_dr);
  if (!std::isfinite(norm) || norm <= 0.0) return seed;

  seed.px = dx_dr / norm;
  seed.py = dy_dr / norm;
  seed.pz = dz_dr / norm;

  for (unsigned int i = 0; i < 6; ++i)
  {
    for (unsigned int j = 0; j < 6; ++j) seed.cov[i][j] = 0.0;
  }
  seed.cov[0][0] = m_seedSigmaX * m_seedSigmaX;
  seed.cov[1][1] = m_seedSigmaY * m_seedSigmaY;
  seed.cov[2][2] = m_seedSigmaZ * m_seedSigmaZ;
  seed.cov[3][3] = m_seedSigmaPx * m_seedSigmaPx;
  seed.cov[4][4] = m_seedSigmaPy * m_seedSigmaPy;
  seed.cov[5][5] = m_seedSigmaPz * m_seedSigmaPz;
  seed.ok = std::isfinite(seed.x) && std::isfinite(seed.y) && std::isfinite(seed.z) &&
            std::isfinite(seed.px) && std::isfinite(seed.py) && std::isfinite(seed.pz);
  return seed;
}

int Tpc_AssembledTrackReco::process_event(PHCompositeNode*)
{
  reset_tree_vars();
  if (m_assembledTrackContainer) m_assembledTrackContainer->Reset();

  std::vector<Piece> pieces;
  if (m_tpcModuleTrackContainer)
  {
    const unsigned int n = m_tpcModuleTrackContainer->size();
    pieces.reserve(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      Piece p;
      if (make_piece(i, p)) pieces.push_back(p);
    }
  }

  std::vector<std::vector<Candidate> > sector_outputs(24);
  std::vector<std::thread> workers;
  workers.reserve(24);
  for (int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      const unsigned int index = static_cast<unsigned int>(side) * 12 + sector;
      workers.push_back(std::thread(&Tpc_AssembledTrackReco::connect_sector_pieces, this, std::cref(pieces), side, sector, std::ref(sector_outputs[index])));
    }
  }
  for (std::thread& worker : workers) worker.join();

  std::vector<Candidate> sector_tracks;
  for (unsigned int i = 0; i < sector_outputs.size(); ++i)
  {
    sector_tracks.insert(sector_tracks.end(), sector_outputs[i].begin(), sector_outputs[i].end());
  }

  std::vector<Candidate> assembled_tracks;
  connect_side_candidates(pieces, sector_tracks, 0, assembled_tracks);
  connect_side_candidates(pieces, sector_tracks, 1, assembled_tracks);

  for (unsigned int it = 0; it < assembled_tracks.size(); ++it)
  {
    const Candidate& c = assembled_tracks[it];
    const unsigned int assembled_id = m_assembledTrackContainer ? m_assembledTrackContainer->size() : it;
    if (m_h_nsegments) m_h_nsegments->Fill(static_cast<double>(c.nsegments));
    if (c.nsegments == 3 && m_idealPadMap)
    {
      const double rfirst = m_idealPadMap->get_radius(c.first_layer);
      const double rlast = m_idealPadMap->get_radius(c.last_layer);
      if (std::isfinite(rfirst) && std::isfinite(rlast))
      {
        const double first_tbin = c.tbin_slope_r * rfirst + c.tbin_intercept_r;
        const double last_tbin = c.tbin_slope_r * rlast + c.tbin_intercept_r;
        const double dtbin_track = last_tbin - first_tbin;
        if (m_h_track_tbin_slope_vs_tbin_span_3modules) m_h_track_tbin_slope_vs_tbin_span_3modules->Fill(dtbin_track, c.tbin_slope_r);
        if (m_h_track_tbin_slope_vs_first_tbin_3modules) m_h_track_tbin_slope_vs_first_tbin_3modules->Fill(first_tbin, c.tbin_slope_r);
        if (m_h_track_tbin_slope_vs_last_tbin_3modules) m_h_track_tbin_slope_vs_last_tbin_3modules->Fill(last_tbin, c.tbin_slope_r);
      }
    }

    Tpc_AssembledTrackv1* out = new Tpc_AssembledTrackv1();
    out->set_event(static_cast<unsigned int>(m_event));
    out->set_track_id(assembled_id);
    out->set_side(c.side);
    out->set_nsegments(c.nsegments);
    out->set_nblobs(c.nblobs);
    out->set_nrawhits(c.nrawhits);
    out->set_first_layer(c.first_layer);
    out->set_last_layer(c.last_layer);
    out->set_first_sector(c.first_sector);
    out->set_last_sector(c.last_sector);
    out->set_first_region(c.first_region);
    out->set_last_region(c.last_region);

    const SeedParameters seed = make_seed_parameters(c);
    if (seed.ok)
    {
      out->set_seed_valid(1);
      out->set_seed_x(seed.x);
      out->set_seed_y(seed.y);
      out->set_seed_z(seed.z);
      out->set_seed_px(seed.px);
      out->set_seed_py(seed.py);
      out->set_seed_pz(seed.pz);
      for (unsigned int iseed = 0; iseed < 6; ++iseed)
      {
        for (unsigned int jseed = 0; jseed < 6; ++jseed)
        {
          out->set_seed_cov(iseed, jseed, seed.cov[iseed][jseed]);
        }
      }
    }

    m_tree_track_id.push_back(assembled_id);
    m_tree_side.push_back(c.side);
    m_tree_nsegments.push_back(c.nsegments);
    m_tree_nblobs.push_back(c.nblobs);
    m_tree_nrawhits.push_back(c.nrawhits);
    m_tree_first_layer.push_back(c.first_layer);
    m_tree_last_layer.push_back(c.last_layer);
    m_tree_first_sector.push_back(c.first_sector);
    m_tree_last_sector.push_back(c.last_sector);
    m_tree_first_region.push_back(c.first_region);
    m_tree_last_region.push_back(c.last_region);

    for (unsigned int ip = 0; ip < c.piece_indices.size(); ++ip)
    {
      const Piece& p = pieces[c.piece_indices[ip]];
      out->add_source_track(p.source_track_id, p.region, p.sector);
      m_tree_source_assembled_track_id.push_back(assembled_id);
      m_tree_source_inmodule_track_id.push_back(p.source_track_id);
      m_tree_source_region.push_back(p.region);
      m_tree_source_sector.push_back(p.sector);
      m_tree_source_side.push_back(p.side);
    }

    for (unsigned int ih = 0; ih < c.hitsetkeys.size(); ++ih)
    {
      out->add_hit_index(c.hitsetkeys[ih], c.hitkeys[ih]);
      m_tree_hit_assembled_track_id.push_back(assembled_id);
      m_tree_hit_hitsetkey.push_back(static_cast<unsigned long long>(c.hitsetkeys[ih]));
      m_tree_hit_hitkey.push_back(static_cast<unsigned long long>(c.hitkeys[ih]));
    }

    if (m_assembledTrackContainer) m_assembledTrackContainer->add_track(out);
    else delete out;
  }

  if (m_tree) m_tree->Fill();

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " input pieces=" << pieces.size()
              << " assembled_tracks=" << assembled_tracks.size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
