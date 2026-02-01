#include "EventPlaneRecov2.h"

#include "EventplaneinfoMapv1.h"
#include "Eventplaneinfov2.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <ffamodules/CDBInterface.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// -- event
#include <ffaobjects/EventHeader.h>

// -- Centrality
#include <centrality/CentralityInfo.h>

// -- sEPD
#include <epd/EpdGeom.h>

// -- root includes --
#include <TFile.h>
#include <TTree.h>

// c++ includes --
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <format>
#include <algorithm>

//____________________________________________________________________________..
EventPlaneRecov2::EventPlaneRecov2(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
EventPlaneRecov2::~EventPlaneRecov2()
{
  std::cout << "EventPlaneRecov2::~EventPlaneRecov2() Calling dtor" << std::endl;
}

bool EventPlaneRecov2::hasValidTree(const std::string &filePath)
{
  // 1. Attempt to open the file
  // "READ" is the default, but being explicit is good practice
  std::unique_ptr<TFile> file(TFile::Open(filePath.c_str(), "READ"));

  // 2. Validate the file pointer and check if the file is "Zombie" (corrupt/unreadable)
  if (!file || file->IsZombie())
  {
    std::cout << "Error: Could not open file: " << filePath << std::endl;
    return false;
  }

  // 3. Attempt to get the object by name
  TObject *obj = file->Get("Multiple");

  // 4. Validate existence and check if it actually inherits from TTree
  if (obj && obj->InheritsFrom(TTree::Class()))
  {
    return true;
  }

  std::cout << "Error: Object 'Multiple' not found or is not a TTree." << std::endl;
  return false;
}

//____________________________________________________________________________..
int EventPlaneRecov2::Init(PHCompositeNode *topNode)
{
  std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);

  if (!m_directURL_EventPlaneCalib.empty() && hasValidTree(m_directURL_EventPlaneCalib))
  {
    m_cdbttree = std::make_unique<CDBTTree>(m_directURL_EventPlaneCalib);
    std::cout << PHWHERE << " Custom Event Plane Calib Found: " << m_directURL_EventPlaneCalib << std::endl;
  }
  else if (!calibdir.empty())
  {
    m_cdbttree = std::make_unique<CDBTTree>(calibdir);
    std::cout << PHWHERE << " Event Plane Calib Found: " << calibdir << std::endl;
  }
  else if (m_doAbortNoEventPlaneCalib)
  {
    std::cout << PHWHERE << " Error: No Event Plane Calib Found and m_doAbortNoEventPlaneCalib is true. Aborting." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  else
  {
    std::cout << PHWHERE << " Error: No Event Plane Calib Found. Skipping Event Plane Calibrations." << std::endl;
    m_doNotCalib = true;
  }

  if (!m_doNotCalib)
  {
    LoadCalib();
  }

  if (Verbosity() > 0)
  {
    print_correction_data();
  }

  CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

std::array<std::array<double, 2>, 2> EventPlaneRecov2::calculate_flattening_matrix(double xx, double yy, double xy, int n, int cent_bin, const std::string& det_label)
{
  std::array<std::array<double, 2>, 2> mat{};

  double D_arg = (xx * yy) - (xy * xy);
  if (D_arg <= 0)
  {
    std::cout << PHWHERE << "Invalid D-term " << D_arg << " for n=" << n << ", cent bin=" << cent_bin << ", det=" << det_label << std::endl;
    // Return Identity Matrix to preserve Recentered vector
    mat[0][0] = 1.0;
    mat[1][1] = 1.0;
    return mat;
  }
  double D = std::sqrt(D_arg);

  double N_term = D * (xx + yy + (2 * D));
  if (N_term <= 0)
  {
    std::cout << PHWHERE << "Invalid N-term " << N_term << " for n=" << n << ", cent bin=" << cent_bin << ", det=" << det_label << std::endl;
    // Return Identity Matrix to preserve Recentered vector
    mat[0][0] = 1.0;
    mat[1][1] = 1.0;
    return mat;
  }
  double inv_sqrt_N = 1.0 / std::sqrt(N_term);

  mat[0][0] = inv_sqrt_N * (yy + D);
  mat[0][1] = -inv_sqrt_N * xy;
  mat[1][0] = mat[0][1];
  mat[1][1] = inv_sqrt_N * (xx + D);
  return mat;
}

//____________________________________________________________________________..
void EventPlaneRecov2::LoadCalib()
{
  size_t south_idx = static_cast<size_t>(Subdetector::S);
  size_t north_idx = static_cast<size_t>(Subdetector::N);
  size_t ns_idx = static_cast<size_t>(Subdetector::NS);

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];

    std::string S_x_avg_name = std::format("Q_S_x_{}_avg", n);
    std::string S_y_avg_name = std::format("Q_S_y_{}_avg", n);
    std::string N_x_avg_name = std::format("Q_N_x_{}_avg", n);
    std::string N_y_avg_name = std::format("Q_N_y_{}_avg", n);

    std::string S_xx_avg_name = std::format("Q_S_xx_{}_avg", n);
    std::string S_yy_avg_name = std::format("Q_S_yy_{}_avg", n);
    std::string S_xy_avg_name = std::format("Q_S_xy_{}_avg", n);
    std::string N_xx_avg_name = std::format("Q_N_xx_{}_avg", n);
    std::string N_yy_avg_name = std::format("Q_N_yy_{}_avg", n);
    std::string N_xy_avg_name = std::format("Q_N_xy_{}_avg", n);

    std::string NS_xx_avg_name = std::format("Q_NS_xx_{}_avg", n);
    std::string NS_yy_avg_name = std::format("Q_NS_yy_{}_avg", n);
    std::string NS_xy_avg_name = std::format("Q_NS_xy_{}_avg", n);

    for (size_t cent_bin = 0; cent_bin < m_bins_cent; ++cent_bin)
    {
      int key = static_cast<int>(cent_bin);

      // South
      auto& dataS = m_correction_data[h_idx][cent_bin][south_idx];
      dataS.avg_Q.x = m_cdbttree->GetDoubleValue(key, S_x_avg_name);
      dataS.avg_Q.y = m_cdbttree->GetDoubleValue(key, S_y_avg_name);

      dataS.avg_Q_xx = m_cdbttree->GetDoubleValue(key, S_xx_avg_name);
      dataS.avg_Q_yy = m_cdbttree->GetDoubleValue(key, S_yy_avg_name);
      dataS.avg_Q_xy = m_cdbttree->GetDoubleValue(key, S_xy_avg_name);

      dataS.X_matrix = calculate_flattening_matrix(dataS.avg_Q_xx, dataS.avg_Q_yy, dataS.avg_Q_xy, n, cent_bin, "South");

      // North
      auto& dataN = m_correction_data[h_idx][cent_bin][north_idx];
      dataN.avg_Q.x = m_cdbttree->GetDoubleValue(key, N_x_avg_name);
      dataN.avg_Q.y = m_cdbttree->GetDoubleValue(key, N_y_avg_name);

      dataN.avg_Q_xx = m_cdbttree->GetDoubleValue(key, N_xx_avg_name);
      dataN.avg_Q_yy = m_cdbttree->GetDoubleValue(key, N_yy_avg_name);
      dataN.avg_Q_xy = m_cdbttree->GetDoubleValue(key, N_xy_avg_name);

      dataN.X_matrix = calculate_flattening_matrix(dataN.avg_Q_xx, dataN.avg_Q_yy, dataN.avg_Q_xy, n, cent_bin, "North");

      // North South
      // Note: We do NOT load avg_Q (x,y) for NS because NS is recentered by summing the recentered S and N vectors.
      auto& dataNS = m_correction_data[h_idx][cent_bin][ns_idx];

      dataNS.avg_Q_xx = m_cdbttree->GetDoubleValue(key, NS_xx_avg_name);
      dataNS.avg_Q_yy = m_cdbttree->GetDoubleValue(key, NS_yy_avg_name);
      dataNS.avg_Q_xy = m_cdbttree->GetDoubleValue(key, NS_xy_avg_name);

      dataNS.X_matrix = calculate_flattening_matrix(dataNS.avg_Q_xx, dataNS.avg_Q_yy, dataNS.avg_Q_xy, n, cent_bin, "NorthSouth");
    }
  }
}

//____________________________________________________________________________..
void EventPlaneRecov2::print_correction_data()
{
  std::cout << std::format("\n{:=>60}\n", "");
  std::cout << std::format("{:^60}\n", "EVENT PLANE CORRECTION DATA SUMMARY");
  std::cout << std::format("{:=>60}\n", "");

  // Iterate through harmonics {2, 3, 4}
  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];
    std::cout << std::format("\n>>> HARMONIC n = {} <<<\n", n);

    // Iterate through Centrality Bins (0-7)
    for (size_t cent = 0; cent < m_bins_cent; ++cent)
    {
      std::cout << std::format("\n  Centrality Bin: {}\n", cent);
      std::cout << std::format("  {:->30}\n", "");

      // Header with fixed column widths
      std::cout << std::format("    {:<12} {:>10} {:>10} {:>10} {:>10} {:>10}\n",
                               "Detector", "Avg Qx", "Avg Qy", "Avg Qxx", "Avg Qyy", "Avg Qxy");

      // Iterate through Subdetectors {S, N}
      for (size_t det_idx = 0; det_idx < 3; ++det_idx)
      {
        std::string det_name;
        if (det_idx == 0)
        {
          det_name = "South";
        }
        else if (det_idx == 1)
        {
          det_name = "North";
        }
        else
        {
          det_name = "NorthSouth";
        }

        const auto& data = m_correction_data[h_idx][cent][det_idx];

        // For NS, Avg Qx/Qy will be 0.0 because they are not loaded from CDB.
        // This is expected behavior.
        std::cout << std::format("    {:<12} {:>10.6f} {:>10.6f} {:>10.6f} {:>10.6f} {:>10.6f}\n",
                                 det_name,
                                 data.avg_Q.x, data.avg_Q.y,
                                 data.avg_Q_xx, data.avg_Q_yy, data.avg_Q_xy);

        // Print X-Matrix in a bracketed layout
        std::cout << std::format("      X-Matrix: [ {:>8.6f}, {:>8.6f} ]\n",
                                 data.X_matrix[0][0], data.X_matrix[0][1]);
        std::cout << std::format("                [ {:>8.6f}, {:>8.6f} ]\n",
                                 data.X_matrix[1][0], data.X_matrix[1][1]);
      }
    }
  }
  std::cout << std::format("\n{:=>60}\n", "");
}

int EventPlaneRecov2::CreateNodes(PHCompositeNode *topNode) {
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    auto global_ptr = std::make_unique<PHCompositeNode>("GLOBAL");
    globalNode = global_ptr.get();
    dstNode->addNode(global_ptr.release());
  }

  EventplaneinfoMap *eps = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!eps)
  {
    auto eps_ptr = std::make_unique<EventplaneinfoMapv1>();
    auto epMapNode_ptr = std::make_unique<PHIODataNode<PHObject>>(eps_ptr.release(), "EventplaneinfoMap", "PHObject");
    globalNode->addNode(epMapNode_ptr.release());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventPlaneRecov2::process_centrality(PHCompositeNode *topNode)
{
  CentralityInfo* centInfo = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!centInfo)
  {
    std::cout << PHWHERE << " CentralityInfo is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_cent = centInfo->get_centile(CentralityInfo::PROP::mbd_NS) * 100;

  if (!std::isfinite(m_cent) || m_cent < 0)
  {
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << " Warning Centrality is out of range. Cent: " << m_cent << ". Cannot calibrate Q vector for this event." << std::endl;
    }
    m_doNotCalibEvent = true;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventPlaneRecov2::process_sEPD(PHCompositeNode* topNode)
{
  TowerInfoContainer* towerinfosEPD = findNode::getClass<TowerInfoContainer>(topNode, m_inputNode);
  if (!towerinfosEPD)
  {
    std::cout << PHWHERE << " TOWERINFO_CALIB_SEPD is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  EpdGeom* epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  if (!epdgeom)
  {
    std::cout << PHWHERE << " TOWERGEOM_EPD is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // sepd
  unsigned int nchannels_epd = towerinfosEPD->size();

  double sepd_total_charge_south = 0;
  double sepd_total_charge_north = 0;

  for (unsigned int channel = 0; channel < nchannels_epd; ++channel)
  {
    TowerInfo* tower = towerinfosEPD->get_tower_at_channel(channel);

    unsigned int key = TowerInfoDefs::encode_epd(channel);
    double charge = tower->get_energy();
    double phi = epdgeom->get_phi(key);

    // skip bad channels
    // skip channels with very low charge
    if (!tower->get_isGood() || charge < m_sepd_min_channel_charge)
    {
      continue;
    }

    // arm = 0: South
    // arm = 1: North
    unsigned int arm = TowerInfoDefs::get_epd_arm(key);

    // sepd charge sums
    double& sepd_total_charge = (arm == 0) ? sepd_total_charge_south : sepd_total_charge_north;

    // Compute total charge for the respective sEPD arm
    sepd_total_charge += charge;

    for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
    {
      int n = m_harmonics[h_idx];
      QVec q_n = {charge * std::cos(n * phi), charge * std::sin(n * phi)};
      m_Q_raw[h_idx][arm].x += q_n.x;
      m_Q_raw[h_idx][arm].y += q_n.y;
    }
  }

  // ensure both total charges are nonzero
  if (sepd_total_charge_south == 0 || sepd_total_charge_north == 0)
  {
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << " Error: Total sEPD Charge is Zero: "
                << "South = " << sepd_total_charge_south
                << ", North = " << sepd_total_charge_north << std::endl;
    }

    // ensure raw Q vec is reset
    m_Q_raw = {};
    m_doNotCalibEvent = true;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    m_Q_raw[h_idx][0].x /= sepd_total_charge_south;
    m_Q_raw[h_idx][0].y /= sepd_total_charge_south;

    m_Q_raw[h_idx][1].x /= sepd_total_charge_north;
    m_Q_raw[h_idx][1].y /= sepd_total_charge_north;

    // NEW: Calculate Raw NS (Sum of Raw S + Raw N)
    m_Q_raw[h_idx][2].x = m_Q_raw[h_idx][0].x + m_Q_raw[h_idx][1].x;
    m_Q_raw[h_idx][2].y = m_Q_raw[h_idx][0].y + m_Q_raw[h_idx][1].y;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void EventPlaneRecov2::correct_QVecs()
{
  size_t cent_bin = static_cast<size_t>(m_cent / 10.0);
  if (cent_bin >= m_bins_cent)
  {
    cent_bin = m_bins_cent - 1;  // Clamp max
  }

  size_t south_idx = static_cast<size_t>(Subdetector::S);
  size_t north_idx = static_cast<size_t>(Subdetector::N);
  size_t ns_idx = static_cast<size_t>(Subdetector::NS);

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    auto& dataS = m_correction_data[h_idx][cent_bin][south_idx];
    auto& dataN = m_correction_data[h_idx][cent_bin][north_idx];
    auto& dataNS = m_correction_data[h_idx][cent_bin][ns_idx];

    double Q_S_x_avg = dataS.avg_Q.x;
    double Q_S_y_avg = dataS.avg_Q.y;
    double Q_N_x_avg = dataN.avg_Q.x;
    double Q_N_y_avg = dataN.avg_Q.y;

    QVec q_S = m_Q_raw[h_idx][south_idx];
    QVec q_N = m_Q_raw[h_idx][north_idx];

    // Apply Recentering
    QVec q_S_recenter = {q_S.x - Q_S_x_avg, q_S.y - Q_S_y_avg};
    QVec q_N_recenter = {q_N.x - Q_N_x_avg, q_N.y - Q_N_y_avg};
    QVec q_NS_recenter = {q_S_recenter.x + q_N_recenter.x, q_S_recenter.y + q_N_recenter.y};

    m_Q_recentered[h_idx][south_idx] = q_S_recenter;
    m_Q_recentered[h_idx][north_idx] = q_N_recenter;
    m_Q_recentered[h_idx][ns_idx] = q_NS_recenter;

    // Flattening Matrix
    const auto &X_S = dataS.X_matrix;
    const auto &X_N = dataN.X_matrix;
    const auto &X_NS = dataNS.X_matrix;

    // Apply Flattening
    double Q_S_x_flat = X_S[0][0] * q_S_recenter.x + X_S[0][1] * q_S_recenter.y;
    double Q_S_y_flat = X_S[1][0] * q_S_recenter.x + X_S[1][1] * q_S_recenter.y;
    double Q_N_x_flat = X_N[0][0] * q_N_recenter.x + X_N[0][1] * q_N_recenter.y;
    double Q_N_y_flat = X_N[1][0] * q_N_recenter.x + X_N[1][1] * q_N_recenter.y;

    double Q_NS_x_flat = X_NS[0][0] * q_NS_recenter.x + X_NS[0][1] * q_NS_recenter.y;
    double Q_NS_y_flat = X_NS[1][0] * q_NS_recenter.x + X_NS[1][1] * q_NS_recenter.y;

    QVec q_S_flat = {Q_S_x_flat, Q_S_y_flat};
    QVec q_N_flat = {Q_N_x_flat, Q_N_y_flat};
    QVec q_NS_flat = {Q_NS_x_flat, Q_NS_y_flat};

    m_Q_flat[h_idx][south_idx] = q_S_flat;
    m_Q_flat[h_idx][north_idx] = q_N_flat;
    m_Q_flat[h_idx][ns_idx] = q_NS_flat;
  }
}

void EventPlaneRecov2::print_QVectors()
{
  std::string header_text = std::format("EVENT Q-VECTOR SUMMARY (Event: {}, CENTRALITY: {:.0f}%)", m_globalEvent, m_cent);

  std::cout << std::format("\n{:*>100}\n", "");
  std::cout << std::format("{:^100}\n", header_text);
  std::cout << std::format("{:*>100}\n", "");

  // Table Header
  std::cout << std::format("  {:<10} {:<10} | {:>21} | {:>21} | {:>21}\n",
                           "Harmonic", "Detector", "Raw (x, y)", "Recentered (x, y)", "Flattened (x, y)");
  std::cout << std::format("  {:-<100}\n", "");

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];

    for (size_t det_idx = 0; det_idx < 3; ++det_idx)
    {
      std::string det_name;
      if (det_idx == 0)
      {
        det_name = "South";
      }
      else if (det_idx == 1)
      {
        det_name = "North";
      }
      else
      {
        det_name = "NorthSouth";
      }

      const auto& raw = m_Q_raw[h_idx][det_idx];
      const auto& rec = m_Q_recentered[h_idx][det_idx];
      const auto& flat = m_Q_flat[h_idx][det_idx];

      std::string h_label = (det_idx == 0) ? std::format("n={}", n) : "";

      // Groups x and y into (val, val) pairs for better scannability
      std::string raw_str  = std::format("({:>8.5f}, {:>8.5f})", raw.x, raw.y);
      std::string rec_str  = std::format("({:>8.5f}, {:>8.5f})", rec.x, rec.y);
      std::string flat_str = std::format("({:>8.5f}, {:>8.5f})", flat.x, flat.y);

      std::cout << std::format("  {:<10} {:<10} | {:<21} | {:<21} | {:10}\n",
                               h_label, det_name, raw_str, rec_str, flat_str);
    }
    if (h_idx < m_harmonics.size() - 1)
    {
      std::cout << std::format("  {:.>100}\n", "");
    }
  }
  std::cout << std::format("{:*>100}\n\n", "");
}

int EventPlaneRecov2::FillNode(PHCompositeNode *topNode)
{
  EventplaneinfoMap *epmap = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!epmap)
  {
    std::cout << PHWHERE << " EventplaneinfoMap is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  size_t vec_size = static_cast<size_t>(*std::ranges::max_element(m_harmonics));

  std::vector<std::pair<double, double>> south_Qvec_raw(vec_size, {NAN, NAN});
  std::vector<std::pair<double, double>> south_Qvec_recentered(vec_size, {NAN, NAN});
  std::vector<std::pair<double, double>> south_Qvec(vec_size, {NAN, NAN});

  std::vector<std::pair<double, double>> north_Qvec_raw(vec_size, {NAN, NAN});
  std::vector<std::pair<double, double>> north_Qvec_recentered(vec_size, {NAN, NAN});
  std::vector<std::pair<double, double>> north_Qvec(vec_size, {NAN, NAN});

  std::vector<std::pair<double, double>> northsouth_Qvec_raw(vec_size, {NAN, NAN});
  std::vector<std::pair<double, double>> northsouth_Qvec_recentered(vec_size, {NAN, NAN});
  std::vector<std::pair<double, double>> northsouth_Qvec(vec_size, {NAN, NAN});

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];
    int idx = n - 1;

    // Fallback logic: Use raw if calibration failed or centrality is out of range
    const auto& Q_S = (m_doNotCalib || m_doNotCalibEvent) ? m_Q_raw[h_idx][0] : m_Q_flat[h_idx][0];
    const auto& Q_N = (m_doNotCalib || m_doNotCalibEvent) ? m_Q_raw[h_idx][1] : m_Q_flat[h_idx][1];
    const auto& Q_NS = (m_doNotCalib || m_doNotCalibEvent) ? m_Q_raw[h_idx][2] : m_Q_flat[h_idx][2];

    const auto& Q_S_raw = m_Q_raw[h_idx][0];
    const auto& Q_S_recentered = m_Q_recentered[h_idx][0];

    const auto& Q_N_raw = m_Q_raw[h_idx][1];
    const auto& Q_N_recentered = m_Q_recentered[h_idx][1];

    const auto& Q_NS_raw = m_Q_raw[h_idx][2];
    const auto& Q_NS_recentered = m_Q_recentered[h_idx][2];

    // South
    south_Qvec_raw[idx] = {Q_S_raw.x, Q_S_raw.y};
    south_Qvec_recentered[idx] = {Q_S_recentered.x, Q_S_recentered.y};
    south_Qvec[idx] = {Q_S.x, Q_S.y};

    // North
    north_Qvec_raw[idx] = {Q_N_raw.x, Q_N_raw.y};
    north_Qvec_recentered[idx] = {Q_N_recentered.x, Q_N_recentered.y};
    north_Qvec[idx] = {Q_N.x, Q_N.y};

    // Combined (North + South)
    northsouth_Qvec_raw[idx] = {Q_NS_raw.x, Q_NS_raw.y};
    northsouth_Qvec_recentered[idx] = {Q_NS_recentered.x, Q_NS_recentered.y};
    northsouth_Qvec[idx] = {Q_NS.x, Q_NS.y};
  }

  // Helper lambda to fill nodes using the class's GetPsi method
  auto create_and_fill = [&](const std::vector<std::pair<double, double>>& qvecs_raw, const std::vector<std::pair<double, double>>& qvecs_recentered, const std::vector<std::pair<double, double>>& qvecs) {
      auto node = std::make_unique<Eventplaneinfov2>();
      node->set_qvector_raw(qvecs_raw);
      node->set_qvector_recentered(qvecs_recentered);
      node->set_qvector(qvecs);

      std::vector<double> psi_vec(vec_size, NAN);
      for (int n : m_harmonics) {
          psi_vec[n-1] = node->GetPsi(qvecs[n-1].first, qvecs[n-1].second, n);
      }
      node->set_shifted_psi(psi_vec);
      return node;
  };

  epmap->insert(create_and_fill(south_Qvec_raw, south_Qvec_recentered, south_Qvec).release(), EventplaneinfoMap::sEPDS);
  epmap->insert(create_and_fill(north_Qvec_raw, north_Qvec_recentered, north_Qvec).release(), EventplaneinfoMap::sEPDN);
  epmap->insert(create_and_fill(northsouth_Qvec_raw, northsouth_Qvec_recentered, northsouth_Qvec).release(), EventplaneinfoMap::sEPDNS);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventPlaneRecov2::process_event(PHCompositeNode *topNode)
{
  EventHeader *eventInfo = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!eventInfo)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_globalEvent = eventInfo->get_EvtSequence();

  int ret = process_centrality(topNode);
  if (ret)
  {
    return ret;
  }

  ret = process_sEPD(topNode);
  if (ret)
  {
    return ret;
  }

  // Calibrate Q Vectors
  if (!m_doNotCalib && !m_doNotCalibEvent)
  {
    correct_QVecs();
  }

  ret = FillNode(topNode);
  if (ret)
  {
    return ret;
  }

  if (Verbosity() > 1)
  {
    print_QVectors();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventPlaneRecov2::ResetEvent([[maybe_unused]] PHCompositeNode *topNode)
{
  m_doNotCalibEvent = false;

  m_Q_raw = {};
  m_Q_recentered = {};
  m_Q_flat = {};

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventPlaneRecov2::End([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "EventPlaneRecov2::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
