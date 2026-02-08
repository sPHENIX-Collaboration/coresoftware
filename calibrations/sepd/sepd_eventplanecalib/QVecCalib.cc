#include "QVecCalib.h"
#include "EventPlaneData.h"

// ====================================================================
// sPHENIX Includes
// ====================================================================
#include <calobase/TowerInfoDefs.h>

// -- Fun4All
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

// -- Nodes
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// -- sEPD
#include <epd/EpdGeom.h>

// -- Run
#include <ffaobjects/RunHeader.h>

// -- CDBTTree
#include <cdbobjects/CDBTTree.h>

// ====================================================================
// Standard C++ Includes
// ====================================================================
#include <format>
#include <numbers>
#include <filesystem>

//____________________________________________________________________________..
QVecCalib::QVecCalib(const std::string &name):
 SubsysReco(name)
{
  std::cout << "QVecCalib::QVecCalib(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
int QVecCalib::Init([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "QVecCalib::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Print("NODETREE");

  int ret = process_QA_hist();
  if (ret)
  {
    return ret;
  }

  init_hists();

  if (m_pass == Pass::ApplyRecentering || m_pass == Pass::ApplyFlattening)
  {
    ret = load_correction_data();
    if (ret)
    {
      return ret;
    }
  }

  prepare_hists();

  return Fun4AllReturnCodes::EVENT_OK;
}

void QVecCalib::prepare_hists()
{
  if (m_pass == Pass::ComputeRecentering)
  {
    prepare_average_hists();
  }
  else if (m_pass == Pass::ApplyRecentering)
  {
    prepare_recenter_hists();
  }
  else if (m_pass == Pass::ApplyFlattening)
  {
    prepare_flattening_hists();
  }
}

int QVecCalib::process_QA_hist()
{
  TH1::AddDirectory(kFALSE);
  auto* file = TFile::Open(m_input_hist.c_str());

  // Check if the file was opened successfully.
  if (!file || file->IsZombie())
  {
    std::cout << PHWHERE << "Error! Cannot not open file: " << m_input_hist << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Get sEPD Total Charge Bounds as function of centrality
  int ret = process_sEPD_event_thresholds(file);
  if (ret)
  {
    return ret;
  }

  // Get List of Bad Channels
  ret = process_bad_channels(file);
  if (ret)
  {
    return ret;
  }

  // cleanup
  file->Close();
  delete file;

  return Fun4AllReturnCodes::EVENT_OK;
}

int QVecCalib::process_sEPD_event_thresholds(TFile* file)
{
  Fun4AllServer *se = Fun4AllServer::instance();

  std::string sepd_totalcharge_centrality = "h2SEPD_totalcharge_centrality";

  auto* hist = file->Get<TH2F>(sepd_totalcharge_centrality.c_str());

  // Check if the hist is stored in the file
  if (hist == nullptr)
  {
    std::cout << PHWHERE << "Error! Cannot find hist: " << sepd_totalcharge_centrality << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  h2SEPD_Charge = static_cast<TH2*>(hist->Clone("h2SEPD_Charge"));
  h2SEPD_Chargev2 = static_cast<TH2*>(hist->Clone("h2SEPD_Chargev2"));

  se->registerHisto(h2SEPD_Charge);
  se->registerHisto(h2SEPD_Chargev2);

  auto* h2SEPD_Charge_py = h2SEPD_Charge->ProfileY("h2SEPD_Charge_py", 1, -1, "s");

  int binsx = h2SEPD_Charge->GetNbinsX();
  int binsy = h2SEPD_Charge->GetNbinsY();
  double ymin = h2SEPD_Charge->GetYaxis()->GetXmin();
  double ymax = h2SEPD_Charge->GetYaxis()->GetXmax();

  hSEPD_Charge_Min = new TProfile("hSEPD_Charge_Min", "; Centrality [%]; sEPD Total Charge", binsy, ymin, ymax);
  hSEPD_Charge_Max = new TProfile("hSEPD_Charge_Max", "; Centrality [%]; sEPD Total Charge", binsy, ymin, ymax);

  se->registerHisto(hSEPD_Charge_Min);
  se->registerHisto(hSEPD_Charge_Max);

  for (int y = 1; y <= binsy; ++y)
  {
    double cent = h2SEPD_Charge_py->GetBinCenter(y);
    double mean = h2SEPD_Charge_py->GetBinContent(y);
    double sigma = h2SEPD_Charge_py->GetBinError(y);

    if (sigma == 0)
    {
      continue;
    }

    double charge_low  = mean - m_sEPD_sigma_threshold * sigma;
    double charge_high = mean + m_sEPD_sigma_threshold * sigma;

    hSEPD_Charge_Min->Fill(cent, charge_low);
    hSEPD_Charge_Max->Fill(cent, charge_high);

    for (int x = 1; x <= binsx; ++x)
    {
      double charge = h2SEPD_Charge->GetXaxis()->GetBinCenter(x);
      double zscore = (charge - mean) / sigma;

      if (std::abs(zscore) > m_sEPD_sigma_threshold)
      {
        h2SEPD_Chargev2->SetBinContent(x, y, 0);
        h2SEPD_Chargev2->SetBinError(x, y, 0);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QVecCalib::process_bad_channels(TFile* file)
{
  Fun4AllServer *se = Fun4AllServer::instance();

  std::string sepd_charge_hist = "hSEPD_Charge";

  auto* hSEPD_Charge = file->Get<TProfile>(sepd_charge_hist.c_str());

  // Check if the hist is stored in the file
  if (hSEPD_Charge == nullptr)
  {
    std::cout << PHWHERE << "Error! Cannot find hist: " << sepd_charge_hist << ", in file: " << file->GetName() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int rbins = 16;
  int bins_charge = 40;

  h2SEPD_South_Charge_rbin = new TH2F("h2SEPD_South_Charge_rbin",
                                      "sEPD South; r_{bin}; Avg Charge",
                                      rbins, -0.5, rbins - 0.5,
                                      bins_charge, 0, bins_charge);

  h2SEPD_North_Charge_rbin = new TH2F("h2SEPD_North_Charge_rbin",
                                      "sEPD North; r_{bin}; Avg Charge",
                                      rbins, -0.5, rbins - 0.5,
                                      bins_charge, 0, bins_charge);

  h2SEPD_South_Charge_rbinv2 = new TH2F("h2SEPD_South_Charge_rbinv2",
                                        "sEPD South; r_{bin}; Avg Charge",
                                        rbins, -0.5, rbins - 0.5,
                                        bins_charge, 0, bins_charge);

  h2SEPD_North_Charge_rbinv2 = new TH2F("h2SEPD_North_Charge_rbinv2",
                                        "sEPD North; r_{bin}; Avg Charge",
                                        rbins, -0.5, rbins - 0.5,
                                        bins_charge, 0, bins_charge);

  hSEPD_Bad_Channels = new TProfile("h_sEPD_Bad_Channels", "sEPD Bad Channels; Channel; Status", QVecShared::SEPD_CHANNELS, -0.5, QVecShared::SEPD_CHANNELS-0.5);

  se->registerHisto(h2SEPD_South_Charge_rbin);
  se->registerHisto(h2SEPD_North_Charge_rbin);
  se->registerHisto(h2SEPD_South_Charge_rbinv2);
  se->registerHisto(h2SEPD_North_Charge_rbinv2);
  se->registerHisto(hSEPD_Bad_Channels);

  auto* h2S = h2SEPD_South_Charge_rbin;
  auto* h2N = h2SEPD_North_Charge_rbin;

  auto* h2Sv2 = h2SEPD_South_Charge_rbinv2;
  auto* h2Nv2 = h2SEPD_North_Charge_rbinv2;

  auto* hBad = hSEPD_Bad_Channels;

  for (int channel = 0; channel < QVecShared::SEPD_CHANNELS; ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(channel);
    int rbin = TowerInfoDefs::get_epd_rbin(key);
    unsigned int arm = TowerInfoDefs::get_epd_arm(key);

    double avg_charge = hSEPD_Charge->GetBinContent(channel + 1);

    auto* h2 = (arm == 0) ? h2S : h2N;

    h2->Fill(rbin, avg_charge);
  }

  auto* hSpx = h2S->ProfileX("hSpx", 2, -1, "s");
  auto* hNpx = h2N->ProfileX("hNpx", 2, -1, "s");

  int ctr_dead = 0;
  int ctr_hot = 0;
  int ctr_cold = 0;

  for (int channel = 0; channel < QVecShared::SEPD_CHANNELS; ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(channel);
    int rbin = TowerInfoDefs::get_epd_rbin(key);
    unsigned int arm = TowerInfoDefs::get_epd_arm(key);

    auto* h2 = (arm == 0) ? h2Sv2 : h2Nv2;
    auto* hprof = (arm == 0) ? hSpx : hNpx;

    double charge = hSEPD_Charge->GetBinContent(channel + 1);
    double mean_charge = hprof->GetBinContent(rbin + 1);
    double sigma = hprof->GetBinError(rbin + 1);
    double zscore = 0.0;

    if (sigma > 0)
    {
      zscore = (charge - mean_charge) / sigma;
    }

    if (charge < m_sEPD_min_avg_charge_threshold || std::abs(zscore) > m_sEPD_sigma_threshold)
    {
      m_bad_channels.insert(channel);

      std::string type;
      int status_fill;

      // dead channel
      if (charge == 0)
      {
        type = "Dead";
        status_fill = static_cast<int>(QVecShared::ChannelStatus::Dead);
        ++ctr_dead;
      }
      // hot channel
      else if (zscore > m_sEPD_sigma_threshold)
      {
        type = "Hot";
        status_fill = static_cast<int>(QVecShared::ChannelStatus::Hot);
        ++ctr_hot;
      }
      // cold channel
      else
      {
        type = "Cold";
        status_fill = static_cast<int>(QVecShared::ChannelStatus::Cold);
        ++ctr_cold;
      }

      hBad->Fill(channel, status_fill);
      std::cout << std::format("{:4} Channel: {:3d}, arm: {}, rbin: {:2d}, Mean: {:5.2f}, Charge: {:5.2f}, Z-Score: {:5.2f}",
                                type, channel, arm, rbin, mean_charge, charge, zscore) << std::endl;
    }
    else
    {
      h2->Fill(rbin, charge);
    }
  }

  std::cout << std::format("Total Bad Channels: {}, Dead: {}, Hot: {}, Cold: {}", m_bad_channels.size(), ctr_dead, ctr_hot, ctr_cold) << std::endl;

  std::cout << "Finished processing Hot sEPD channels" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void QVecCalib::init_hists()
{
  unsigned int bins_psi = 126;
  double psi_low = -std::numbers::pi;
  double psi_high = std::numbers::pi;

  Fun4AllServer *se = Fun4AllServer::instance();

  hCentrality = new TH1F("hCentrality", "|z| < 10 cm and MB; Centrality [%]; Events", m_cent_bins, m_cent_low, m_cent_high);
  se->registerHisto(hCentrality);

  std::string pass_suffix;
  if (m_pass == Pass::ApplyRecentering)
  {
    pass_suffix = "_corr";
  }
  else if (m_pass == Pass::ApplyFlattening)
  {
    pass_suffix = "_corr2";
  }

  // n = 2, 3, 4, etc.
  for (int n : m_harmonics)
  {
    std::string name_S = std::format("h2_sEPD_Psi_S_{}{}", n, pass_suffix);
    std::string name_N = std::format("h2_sEPD_Psi_N_{}{}", n, pass_suffix);
    std::string name_NS = std::format("h2_sEPD_Psi_NS_{}{}", n, pass_suffix);

    std::string title_S = std::format("sEPD South #Psi (Order {0}); Centrality [%]; {0}#Psi^{{S}}_{{{0}}}", n);
    std::string title_N = std::format("sEPD North #Psi (Order {0}); Centrality [%]; {0}#Psi^{{N}}_{{{0}}}", n);
    std::string title_NS = std::format("sEPD North South #Psi (Order {0}); Centrality [%]; {0}#Psi^{{NS}}_{{{0}}}", n);

    m_hists2D[name_S] = new TH2F(name_S.c_str(), title_S.c_str(), m_cent_bins, m_cent_low, m_cent_high, bins_psi, psi_low, psi_high);
    m_hists2D[name_N] = new TH2F(name_N.c_str(), title_N.c_str(), m_cent_bins, m_cent_low, m_cent_high, bins_psi, psi_low, psi_high);
    m_hists2D[name_NS] = new TH2F(name_NS.c_str(), title_NS.c_str(), m_cent_bins, m_cent_low, m_cent_high, bins_psi, psi_low, psi_high);

    // South, North
    for (auto det : m_subdetectors)
    {
      std::string det_str = (det == QVecShared::Subdetector::S) ? "S" : "N";
      std::string det_name = (det == QVecShared::Subdetector::S) ? "South" : "North";

      std::string q_avg_sq_cross_name;
      std::string q_avg_sq_cross_title = std::format("sEPD {0}; Centrality [%]; <Q_{{{1},x}} Q_{{{1},y}}>", det_name, n);

      if (m_pass == Pass::ApplyRecentering)
      {
        q_avg_sq_cross_name = QVecShared::get_hist_name(det_str, "xy", n);
      }

      if (m_pass == Pass::ApplyFlattening)
      {
        q_avg_sq_cross_name = QVecShared::get_hist_name(det_str, "xy", n, "_corr");
      }

      if (!q_avg_sq_cross_name.empty())
      {
        m_profiles[q_avg_sq_cross_name] = new TProfile(q_avg_sq_cross_name.c_str(), q_avg_sq_cross_title.c_str(),
                                                                     m_cent_bins, m_cent_low, m_cent_high);
      }

      for (auto comp : m_components)
      {
        std::string comp_str = (comp == QVecShared::QComponent::X) ? "x" : "y";
        std::string name = QVecShared::get_hist_name(det_str, comp_str, n, pass_suffix);

        auto add_profile = [&](const std::string& prof_name, std::string_view label_suffix = "")
        {
          std::string title = std::format("sEPD {}; Centrality [%]; <Q_{{{},{}}}{}>", det_name, n, comp_str, label_suffix);
          m_profiles[prof_name] = new TProfile(prof_name.c_str(), title.c_str(), m_cent_bins, m_cent_low, m_cent_high);
        };

        add_profile(name);

        // 2. Only generate strings and profiles for the current pass
        switch (m_pass)
        {
          case Pass::ComputeRecentering:
          {
            break;
          }

          case Pass::ApplyRecentering:
          {
            std::string name_sq = QVecShared::get_hist_name(det_str, comp_str+comp_str, n);
            add_profile(name_sq, "^{2}");
            break;
          }

          case Pass::ApplyFlattening:
          {
            std::string name_sq_corr = QVecShared::get_hist_name(det_str, comp_str+comp_str, n, "_corr");
            add_profile(name_sq_corr, "^{2}");
            break;
          }
        }
      }
    }

    // Init for Combined NS Histograms
    if (m_pass == Pass::ApplyRecentering || m_pass == Pass::ApplyFlattening)
    {
      std::string det_str = "NS";

      // Initialize 2nd Moment Profiles for NS (needed to compute flattening)
      for (const auto* comp : {"xx", "yy", "xy"})
      {
        std::string name = QVecShared::get_hist_name(det_str, comp, n);
        std::string title = std::format("sEPD NS; Centrality [%]; <Q_{{{},{}}}>", n, comp);
        m_profiles[name] = new TProfile(name.c_str(), title.c_str(), m_cent_bins, m_cent_low, m_cent_high);
      }

      // Initialize Validation Profiles (Flattened NS)
      if (m_pass == Pass::ApplyFlattening)
      {
        for (const auto* comp : {"xx", "yy", "xy"})
        {
            std::string name = QVecShared::get_hist_name(det_str, comp, n, "_corr");
            std::string title = std::format("sEPD NS Corrected; Centrality [%]; <Q_{{{},{}}}^{{2}}>", n, comp);
            m_profiles[name] = new TProfile(name.c_str(), title.c_str(), m_cent_bins, m_cent_low, m_cent_high);
        }
      }
    }
  }
}

std::array<std::array<double, 2>, 2> QVecCalib::calculate_flattening_matrix(double xx, double yy, double xy, int n, int cent_bin, const std::string& det_label)
{
  double D_arg = (xx * yy) - (xy * xy);
  if (D_arg < 1e-12)
  {
    std::cout << "Warning: Near-zero determinant in bin " << cent_bin << ". Skipping matrix calc." << std::endl;
    return std::array<std::array<double, 2>, 2>{{{1, 0}, {0, 1}}};  // Return Identity
  }
  double D = std::sqrt(D_arg);

  double N_term = D * (xx + yy + (2 * D));
  if (N_term <= 0)
  {
    throw std::runtime_error(std::format(
        "Invalid N-term ({}) for n={}, cent={}, det={}", N_term, n, cent_bin, det_label));
  }
  double inv_sqrt_N = 1.0 / std::sqrt(N_term);

  std::array<std::array<double, 2>, 2> mat{};
  mat[0][0] = inv_sqrt_N * (yy + D);
  mat[0][1] = -inv_sqrt_N * xy;
  mat[1][0] = mat[0][1];
  mat[1][1] = inv_sqrt_N * (xx + D);
  return mat;
}

template <typename T>
T* QVecCalib::load_and_clone(TFile* file, const std::string& name) {
  auto* obj = file->Get<T>(name.c_str());
  if (!obj)
  {
    throw std::runtime_error(std::format("Could not find histogram '{}' in file '{}'", name, file->GetName()));
  }
  return static_cast<T*>(obj->Clone());
}

int QVecCalib::load_correction_data()
{
  Fun4AllServer *se = Fun4AllServer::instance();
  auto* file = TFile::Open(m_input_Q_calib.c_str());

  if (!file || file->IsZombie())
  {
    std::cout << PHWHERE << "Error! Cannot open: " << m_input_Q_calib << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  using SD = QVecShared::Subdetector;

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];

    // Helper to load and register histograms automatically
    auto load_reg = [&](const std::string& det, const std::string& var, const std::string& suffix = "")
    {
      std::string name = QVecShared::get_hist_name(det, var, n, suffix);
      m_profiles[name] = load_and_clone<TProfile>(file, name);
      se->registerHisto(m_profiles[name]);
      return name;
    };

    // Load standard Recentering averages for S and N
    std::string s_names[2][2];  // [det][comp]
    for (int d = 0; d < 2; ++d)
    {
      std::string det_str = (d == 0) ? "S" : "N";
      s_names[d][0] = load_reg(det_str, "x");
      s_names[d][1] = load_reg(det_str, "y");
    }

    // Load Flattening (2nd moment) data if needed
    if (m_pass == Pass::ApplyFlattening)
    {
      for (const auto& det_str : {"S", "N", "NS"})
      {
        for (const auto& var : {"xx", "yy", "xy"})
        {
          load_reg(det_str, var);
        }
      }
    }

    // Populate the CorrectionData matrix
    for (size_t cent_bin = 0; cent_bin < m_cent_bins; ++cent_bin)
    {
      int bin = static_cast<int>(cent_bin) + 1;

      // Populate Recentering (S, N)
      for (int d = 0; d < 2; ++d)
      {
        m_correction_data[cent_bin][h_idx][d].avg_Q = {m_profiles[s_names[d][0]]->GetBinContent(bin), m_profiles[s_names[d][1]]->GetBinContent(bin)};
      }

      if (m_pass == Pass::ApplyFlattening)
      {
        // Populate Flattening for S, N, and NS
        for (int d = 0; d < (int) SD::Count; ++d)
        {
          std::string det_str = (d == 0) ? "S" : (d == 1) ? "N" : "NS";
          double xx = m_profiles[QVecShared::get_hist_name(det_str, "xx", n)]->GetBinContent(bin);
          double yy = m_profiles[QVecShared::get_hist_name(det_str, "yy", n)]->GetBinContent(bin);
          double xy = m_profiles[QVecShared::get_hist_name(det_str, "xy", n)]->GetBinContent(bin);
          
          auto& data = m_correction_data[cent_bin][h_idx][d];
          data.X_matrix = calculate_flattening_matrix(xx, yy, xy, n, cent_bin, det_str);
          data.avg_Q_xx = xx;
          data.avg_Q_yy = yy;
          data.avg_Q_xy = xy;
        }
      }
    }
  }

  file->Close();
  delete file;
  return Fun4AllReturnCodes::EVENT_OK;
}

void QVecCalib::prepare_average_hists()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  for (int n : m_harmonics)
  {
    std::string S_x_avg_name = QVecShared::get_hist_name("S", "x", n);
    std::string S_y_avg_name = QVecShared::get_hist_name("S", "y", n);
    std::string N_x_avg_name = QVecShared::get_hist_name("N", "x", n);
    std::string N_y_avg_name = QVecShared::get_hist_name("N", "y", n);

    std::string psi_S_name = std::format("h2_sEPD_Psi_S_{}", n);
    std::string psi_N_name = std::format("h2_sEPD_Psi_N_{}", n);
    std::string psi_NS_name = std::format("h2_sEPD_Psi_NS_{}", n);

    AverageHists h;

    h.S_x_avg = m_profiles.at(S_x_avg_name);
    h.S_y_avg = m_profiles.at(S_y_avg_name);
    h.N_x_avg = m_profiles.at(N_x_avg_name);
    h.N_y_avg = m_profiles.at(N_y_avg_name);

    h.Psi_S = m_hists2D.at(psi_S_name);
    h.Psi_N = m_hists2D.at(psi_N_name);
    h.Psi_NS = m_hists2D.at(psi_NS_name);

    se->registerHisto(h.S_x_avg);
    se->registerHisto(h.S_y_avg);
    se->registerHisto(h.N_x_avg);
    se->registerHisto(h.N_y_avg);

    se->registerHisto(h.Psi_S);
    se->registerHisto(h.Psi_N);
    se->registerHisto(h.Psi_NS);

    m_average_hists.push_back(h);
  }
}

void QVecCalib::prepare_recenter_hists()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  for (int n : m_harmonics)
  {
    std::string S_x_corr_avg_name = QVecShared::get_hist_name("S", "x", n, "_corr");
    std::string S_y_corr_avg_name = QVecShared::get_hist_name("S", "y", n, "_corr");
    std::string N_x_corr_avg_name = QVecShared::get_hist_name("N", "x", n, "_corr");
    std::string N_y_corr_avg_name = QVecShared::get_hist_name("N", "y", n, "_corr");

    std::string S_xx_avg_name = QVecShared::get_hist_name("S", "xx", n);
    std::string S_yy_avg_name = QVecShared::get_hist_name("S", "yy", n);
    std::string S_xy_avg_name = QVecShared::get_hist_name("S", "xy", n);
    std::string N_xx_avg_name = QVecShared::get_hist_name("N", "xx", n);
    std::string N_yy_avg_name = QVecShared::get_hist_name("N", "yy", n);
    std::string N_xy_avg_name = QVecShared::get_hist_name("N", "xy", n);

    std::string NS_xx_avg_name = QVecShared::get_hist_name("NS", "xx", n);
    std::string NS_yy_avg_name = QVecShared::get_hist_name("NS", "yy", n);
    std::string NS_xy_avg_name = QVecShared::get_hist_name("NS", "xy", n);

    std::string psi_S_name = std::format("h2_sEPD_Psi_S_{}_corr", n);
    std::string psi_N_name = std::format("h2_sEPD_Psi_N_{}_corr", n);
    std::string psi_NS_name = std::format("h2_sEPD_Psi_NS_{}_corr", n);

    RecenterHists h;

    h.S_x_corr_avg = m_profiles.at(S_x_corr_avg_name);
    h.S_y_corr_avg = m_profiles.at(S_y_corr_avg_name);
    h.N_x_corr_avg = m_profiles.at(N_x_corr_avg_name);
    h.N_y_corr_avg = m_profiles.at(N_y_corr_avg_name);

    h.S_xx_avg = m_profiles.at(S_xx_avg_name);
    h.S_yy_avg = m_profiles.at(S_yy_avg_name);
    h.S_xy_avg = m_profiles.at(S_xy_avg_name);

    h.N_xx_avg = m_profiles.at(N_xx_avg_name);
    h.N_yy_avg = m_profiles.at(N_yy_avg_name);
    h.N_xy_avg = m_profiles.at(N_xy_avg_name);

    h.NS_xx_avg = m_profiles.at(NS_xx_avg_name);
    h.NS_yy_avg = m_profiles.at(NS_yy_avg_name);
    h.NS_xy_avg = m_profiles.at(NS_xy_avg_name);

    h.Psi_S_corr = m_hists2D.at(psi_S_name);
    h.Psi_N_corr = m_hists2D.at(psi_N_name);
    h.Psi_NS_corr = m_hists2D.at(psi_NS_name);

    se->registerHisto(h.S_x_corr_avg);
    se->registerHisto(h.S_y_corr_avg);
    se->registerHisto(h.N_x_corr_avg);
    se->registerHisto(h.N_y_corr_avg);

    se->registerHisto(h.S_xx_avg);
    se->registerHisto(h.S_yy_avg);
    se->registerHisto(h.S_xy_avg);

    se->registerHisto(h.N_xx_avg);
    se->registerHisto(h.N_yy_avg);
    se->registerHisto(h.N_xy_avg);

    se->registerHisto(h.NS_xx_avg);
    se->registerHisto(h.NS_yy_avg);
    se->registerHisto(h.NS_xy_avg);

    se->registerHisto(h.Psi_S_corr);
    se->registerHisto(h.Psi_N_corr);
    se->registerHisto(h.Psi_NS_corr);

    m_recenter_hists.push_back(h);
  }
}

void QVecCalib::prepare_flattening_hists()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  for (int n : m_harmonics)
  {
    std::string S_x_corr2_avg_name = QVecShared::get_hist_name("S", "x", n, "_corr2");
    std::string S_y_corr2_avg_name = QVecShared::get_hist_name("S", "y", n, "_corr2");
    std::string N_x_corr2_avg_name = QVecShared::get_hist_name("N", "x", n, "_corr2");
    std::string N_y_corr2_avg_name = QVecShared::get_hist_name("N", "y", n, "_corr2");

    std::string S_xx_corr_avg_name = QVecShared::get_hist_name("S", "xx", n, "_corr");
    std::string S_yy_corr_avg_name = QVecShared::get_hist_name("S", "yy", n, "_corr");
    std::string S_xy_corr_avg_name = QVecShared::get_hist_name("S", "xy", n, "_corr");
    std::string N_xx_corr_avg_name = QVecShared::get_hist_name("N", "xx", n, "_corr");
    std::string N_yy_corr_avg_name = QVecShared::get_hist_name("N", "yy", n, "_corr");
    std::string N_xy_corr_avg_name = QVecShared::get_hist_name("N", "xy", n, "_corr");

    std::string NS_xx_corr_avg_name = QVecShared::get_hist_name("NS", "xx", n, "_corr");
    std::string NS_yy_corr_avg_name = QVecShared::get_hist_name("NS", "yy", n, "_corr");
    std::string NS_xy_corr_avg_name = QVecShared::get_hist_name("NS", "xy", n, "_corr");

    std::string psi_S_name = std::format("h2_sEPD_Psi_S_{}_corr2", n);
    std::string psi_N_name = std::format("h2_sEPD_Psi_N_{}_corr2", n);
    std::string psi_NS_name = std::format("h2_sEPD_Psi_NS_{}_corr2", n);

    FlatteningHists h;

    h.S_x_corr2_avg = m_profiles.at(S_x_corr2_avg_name);
    h.S_y_corr2_avg = m_profiles.at(S_y_corr2_avg_name);
    h.N_x_corr2_avg = m_profiles.at(N_x_corr2_avg_name);
    h.N_y_corr2_avg = m_profiles.at(N_y_corr2_avg_name);

    h.S_xx_corr_avg = m_profiles.at(S_xx_corr_avg_name);
    h.S_yy_corr_avg = m_profiles.at(S_yy_corr_avg_name);
    h.S_xy_corr_avg = m_profiles.at(S_xy_corr_avg_name);

    h.N_xx_corr_avg = m_profiles.at(N_xx_corr_avg_name);
    h.N_yy_corr_avg = m_profiles.at(N_yy_corr_avg_name);
    h.N_xy_corr_avg = m_profiles.at(N_xy_corr_avg_name);

    h.NS_xx_corr_avg = m_profiles.at(NS_xx_corr_avg_name);
    h.NS_yy_corr_avg = m_profiles.at(NS_yy_corr_avg_name);
    h.NS_xy_corr_avg = m_profiles.at(NS_xy_corr_avg_name);

    h.Psi_S_corr2 = m_hists2D.at(psi_S_name);
    h.Psi_N_corr2 = m_hists2D.at(psi_N_name);
    h.Psi_NS_corr2 = m_hists2D.at(psi_NS_name);

    se->registerHisto(h.S_x_corr2_avg);
    se->registerHisto(h.S_y_corr2_avg);
    se->registerHisto(h.N_x_corr2_avg);
    se->registerHisto(h.N_y_corr2_avg);

    se->registerHisto(h.S_xx_corr_avg);
    se->registerHisto(h.S_yy_corr_avg);
    se->registerHisto(h.S_xy_corr_avg);

    se->registerHisto(h.N_xx_corr_avg);
    se->registerHisto(h.N_yy_corr_avg);
    se->registerHisto(h.N_xy_corr_avg);

    se->registerHisto(h.NS_xx_corr_avg);
    se->registerHisto(h.NS_yy_corr_avg);
    se->registerHisto(h.NS_xy_corr_avg);

    se->registerHisto(h.Psi_S_corr2);
    se->registerHisto(h.Psi_N_corr2);
    se->registerHisto(h.Psi_NS_corr2);

    m_flattening_hists.push_back(h);
  }
}

int QVecCalib::InitRun(PHCompositeNode *topNode)
{
  RunHeader* run_header = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!run_header)
  {
    std::cout << PHWHERE << "RunHeader Node missing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_runnumber = run_header->get_RunNumber();

  EpdGeom* epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  if (!epdgeom)
  {
    std::cout << PHWHERE << "TOWERGEOM_EPD Node missing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_trig_cache.assign(m_harmonics.size(), std::vector<std::pair<double, double>>(QVecShared::SEPD_CHANNELS));

  for (int channel = 0; channel < QVecShared::SEPD_CHANNELS; ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(channel);
    double phi = epdgeom->get_phi(key);

    for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
    {
      int n = m_harmonics[h_idx];
      m_trig_cache[h_idx][channel] = {std::cos(n * phi), std::sin(n * phi)};
    }
  }

  std::cout << "QVecCalib::InitRun - Trigonometry cache initialized for "
            << QVecShared::SEPD_CHANNELS << " channels." << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

void QVecCalib::process_averages(double cent, const QVecShared::QVec& q_S, const QVecShared::QVec& q_N, const AverageHists& h)
{
  double psi_S = std::atan2(q_S.y, q_S.x);
  double psi_N = std::atan2(q_N.y, q_N.x);
  double psi_NS = std::atan2(q_S.y + q_N.y, q_S.x + q_N.x);

  h.S_x_avg->Fill(cent, q_S.x);
  h.S_y_avg->Fill(cent, q_S.y);
  h.N_x_avg->Fill(cent, q_N.x);
  h.N_y_avg->Fill(cent, q_N.y);

  h.Psi_S->Fill(cent, psi_S);
  h.Psi_N->Fill(cent, psi_N);
  h.Psi_NS->Fill(cent, psi_NS);
}

void QVecCalib::process_recentering(double cent, size_t h_idx, const QVecShared::QVec& q_S, const QVecShared::QVec& q_N, const RecenterHists& h)
{
  size_t cent_bin = static_cast<size_t>(hCentrality->FindBin(cent) - 1);

  const auto& S = m_correction_data[cent_bin][h_idx][(size_t) QVecShared::Subdetector::S];
  const auto& N = m_correction_data[cent_bin][h_idx][(size_t) QVecShared::Subdetector::N];

  double Q_S_x_avg = S.avg_Q.x;
  double Q_S_y_avg = S.avg_Q.y;
  double Q_N_x_avg = N.avg_Q.x;
  double Q_N_y_avg = N.avg_Q.y;

  QVecShared::QVec q_S_corr = {q_S.x - Q_S_x_avg, q_S.y - Q_S_y_avg};
  QVecShared::QVec q_N_corr = {q_N.x - Q_N_x_avg, q_N.y - Q_N_y_avg};

  //  Construct Combined Recentered Vector
  //  We use the sum of the individually recentered vectors
  QVecShared::QVec q_NS_corr = {q_S_corr.x + q_N_corr.x, q_S_corr.y + q_N_corr.y};

  double psi_S_corr = std::atan2(q_S_corr.y, q_S_corr.x);
  double psi_N_corr = std::atan2(q_N_corr.y, q_N_corr.x);
  double psi_NS_corr = std::atan2(q_NS_corr.y, q_NS_corr.x);

  h.S_x_corr_avg->Fill(cent, q_S_corr.x);
  h.S_y_corr_avg->Fill(cent, q_S_corr.y);
  h.N_x_corr_avg->Fill(cent, q_N_corr.x);
  h.N_y_corr_avg->Fill(cent, q_N_corr.y);

  h.S_xx_avg->Fill(cent, q_S_corr.x * q_S_corr.x);
  h.S_yy_avg->Fill(cent, q_S_corr.y * q_S_corr.y);
  h.S_xy_avg->Fill(cent, q_S_corr.x * q_S_corr.y);
  h.N_xx_avg->Fill(cent, q_N_corr.x * q_N_corr.x);
  h.N_yy_avg->Fill(cent, q_N_corr.y * q_N_corr.y);
  h.N_xy_avg->Fill(cent, q_N_corr.x * q_N_corr.y);

  h.NS_xx_avg->Fill(cent, q_NS_corr.x * q_NS_corr.x);
  h.NS_yy_avg->Fill(cent, q_NS_corr.y * q_NS_corr.y);
  h.NS_xy_avg->Fill(cent, q_NS_corr.x * q_NS_corr.y);

  h.Psi_S_corr->Fill(cent, psi_S_corr);
  h.Psi_N_corr->Fill(cent, psi_N_corr);
  h.Psi_NS_corr->Fill(cent, psi_NS_corr);
}

void QVecCalib::process_flattening(double cent, size_t h_idx, const QVecShared::QVec& q_S, const QVecShared::QVec& q_N, const FlatteningHists& h)
{
  size_t cent_bin = static_cast<size_t>(hCentrality->FindBin(cent) - 1);

  const auto& S = m_correction_data[cent_bin][h_idx][(size_t) QVecShared::Subdetector::S];
  const auto& N = m_correction_data[cent_bin][h_idx][(size_t) QVecShared::Subdetector::N];
  const auto& NS = m_correction_data[cent_bin][h_idx][(size_t) QVecShared::Subdetector::NS];

  double Q_S_x_avg = S.avg_Q.x;
  double Q_S_y_avg = S.avg_Q.y;
  double Q_N_x_avg = N.avg_Q.x;
  double Q_N_y_avg = N.avg_Q.y;

  QVecShared::QVec q_S_corr = {q_S.x - Q_S_x_avg, q_S.y - Q_S_y_avg};
  QVecShared::QVec q_N_corr = {q_N.x - Q_N_x_avg, q_N.y - Q_N_y_avg};

  // Construct Combined Recentered Vector
  QVecShared::QVec q_NS_corr = {q_S_corr.x + q_N_corr.x, q_S_corr.y + q_N_corr.y};

  const auto& X_S = S.X_matrix;
  const auto& X_N = N.X_matrix;
  const auto& X_NS = NS.X_matrix;

  double Q_S_x_corr2 = X_S[0][0] * q_S_corr.x + X_S[0][1] * q_S_corr.y;
  double Q_S_y_corr2 = X_S[1][0] * q_S_corr.x + X_S[1][1] * q_S_corr.y;
  double Q_N_x_corr2 = X_N[0][0] * q_N_corr.x + X_N[0][1] * q_N_corr.y;
  double Q_N_y_corr2 = X_N[1][0] * q_N_corr.x + X_N[1][1] * q_N_corr.y;

  double Q_NS_x_corr2 = X_NS[0][0] * q_NS_corr.x + X_NS[0][1] * q_NS_corr.y;
  double Q_NS_y_corr2 = X_NS[1][0] * q_NS_corr.x + X_NS[1][1] * q_NS_corr.y;

  QVecShared::QVec q_S_corr2 = {Q_S_x_corr2, Q_S_y_corr2};
  QVecShared::QVec q_N_corr2 = {Q_N_x_corr2, Q_N_y_corr2};
  QVecShared::QVec q_NS_corr2 = {Q_NS_x_corr2, Q_NS_y_corr2};

  double psi_S = std::atan2(q_S_corr2.y, q_S_corr2.x);
  double psi_N = std::atan2(q_N_corr2.y, q_N_corr2.x);
  double psi_NS = std::atan2(q_NS_corr2.y, q_NS_corr2.x);

  h.S_x_corr2_avg->Fill(cent, q_S_corr2.x);
  h.S_y_corr2_avg->Fill(cent, q_S_corr2.y);
  h.N_x_corr2_avg->Fill(cent, q_N_corr2.x);
  h.N_y_corr2_avg->Fill(cent, q_N_corr2.y);

  h.S_xx_corr_avg->Fill(cent, q_S_corr2.x * q_S_corr2.x);
  h.S_yy_corr_avg->Fill(cent, q_S_corr2.y * q_S_corr2.y);
  h.S_xy_corr_avg->Fill(cent, q_S_corr2.x * q_S_corr2.y);
  h.N_xx_corr_avg->Fill(cent, q_N_corr2.x * q_N_corr2.x);
  h.N_yy_corr_avg->Fill(cent, q_N_corr2.y * q_N_corr2.y);
  h.N_xy_corr_avg->Fill(cent, q_N_corr2.x * q_N_corr2.y);

  h.NS_xx_corr_avg->Fill(cent, q_NS_corr2.x * q_NS_corr2.x);
  h.NS_yy_corr_avg->Fill(cent, q_NS_corr2.y * q_NS_corr2.y);
  h.NS_xy_corr_avg->Fill(cent, q_NS_corr2.x * q_NS_corr2.y);

  h.Psi_S_corr2->Fill(cent, psi_S);
  h.Psi_N_corr2->Fill(cent, psi_N);
  h.Psi_NS_corr2->Fill(cent, psi_NS);
}

bool QVecCalib::process_sEPD()
{
  double sepd_total_charge_south = 0;
  double sepd_total_charge_north = 0;

  // Loop over all sEPD Channels
  for (int channel = 0; channel < QVecShared::SEPD_CHANNELS; ++channel)
  {
    double charge = m_evtdata->get_sepd_charge(channel);

    // Skip Bad Channels
    if (m_bad_channels.contains(channel) || charge <= 0)
    {
      continue;
    }

    unsigned int key = TowerInfoDefs::encode_epd(channel);
    unsigned int arm = TowerInfoDefs::get_epd_arm(key);

    // arm = 0: South
    // arm = 1: North
    if (arm == 0)
    {
      sepd_total_charge_south += charge;
    }
    else
    {
      sepd_total_charge_north += charge;
    }

    // Compute Raw Q vectors for each harmonic and respective arm
    for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
    {
      // Optimized lookup instead of std::cos/std::sin calls
      const auto& [cached_cos, cached_sin] = m_trig_cache[h_idx][channel];

      m_q_vectors[h_idx][arm].x += charge * cached_cos;
      m_q_vectors[h_idx][arm].y += charge * cached_sin;
    }
  }

  // Skip Events with Zero sEPD Total Charge in either arm
  if (sepd_total_charge_south == 0 || sepd_total_charge_north == 0)
  {
    return false;
  }

  // Normalize the Q-vectors by total charge
  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    for (auto det : m_subdetectors)
    {
      size_t det_idx = (det == QVecShared::Subdetector::S) ? 0 : 1;
      double sepd_total_charge = (det_idx == 0) ? sepd_total_charge_south : sepd_total_charge_north;
      m_q_vectors[h_idx][det_idx].x /= sepd_total_charge;
      m_q_vectors[h_idx][det_idx].y /= sepd_total_charge;
    }
  }

  return true;
}

bool QVecCalib::process_event_check()
{
  double cent = m_evtdata->get_event_centrality();
  int cent_bin = hSEPD_Charge_Min->FindBin(cent);

  double sepd_totalcharge = m_evtdata->get_sepd_totalcharge();

  double sepd_totalcharge_min = hSEPD_Charge_Min->GetBinContent(cent_bin);
  double sepd_totalcharge_max = hSEPD_Charge_Max->GetBinContent(cent_bin);

  return sepd_totalcharge > sepd_totalcharge_min && sepd_totalcharge < sepd_totalcharge_max;
}

//____________________________________________________________________________..
int QVecCalib::process_event([[maybe_unused]] PHCompositeNode *topNode)
{
  m_evtdata = findNode::getClass<EventPlaneData>(topNode, "EventPlaneData");
  if (!m_evtdata)
  {
    std::cout << PHWHERE << "EventPlaneData Node missing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int event_id = m_evtdata->get_event_id();

  if (Verbosity() && m_event % PROGRESS_REPORT_INTERVAL == 0)
  {
    std::cout << "Progress: " << m_event << ", Global: " << event_id << std::endl;
  }
  ++m_event;

  double cent = m_evtdata->get_event_centrality();

  bool isGood = process_event_check();

  // Skip Events with non correlation between centrality and sEPD
  if (!isGood)
  {
    ++m_event_counters.bad_centrality_sepd_correlation;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  isGood = process_sEPD();

  // Skip Events with Zero sEPD Total Charge in either arm
  if (!isGood)
  {
    ++m_event_counters.zero_sepd_total_charge;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  hCentrality->Fill(cent);

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    const auto& q_S = m_q_vectors[h_idx][0];  // 0 for South
    const auto& q_N = m_q_vectors[h_idx][1];  // 1 for North

    // --- First Pass: Derive 1st Order ---
    if (m_pass == Pass::ComputeRecentering)
    {
      process_averages(cent, q_S, q_N, m_average_hists[h_idx]);
    }

    // --- Second Pass: Apply 1st Order, Derive 2nd Order ---
    else if (m_pass == Pass::ApplyRecentering)
    {
      process_recentering(cent, h_idx, q_S, q_N, m_recenter_hists[h_idx]);
    }

    // --- Third Pass: Apply 2nd Order, Validate ---
    else if (m_pass == Pass::ApplyFlattening)
    {
      process_flattening(cent, h_idx, q_S, q_N, m_flattening_hists[h_idx]);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QVecCalib::ResetEvent([[maybe_unused]] PHCompositeNode *topNode)
{
  m_q_vectors = {};

  return Fun4AllReturnCodes::EVENT_OK;
}

void QVecCalib::compute_averages(size_t cent_bin, int h_idx)
{
  int n = m_harmonics[h_idx];

  std::string S_x_avg_name = QVecShared::get_hist_name("S", "x", n);
  std::string S_y_avg_name = QVecShared::get_hist_name("S", "y", n);
  std::string N_x_avg_name = QVecShared::get_hist_name("N", "x", n);
  std::string N_y_avg_name = QVecShared::get_hist_name("N", "y", n);

  int bin = static_cast<int>(cent_bin + 1);

  double Q_S_x_avg = m_profiles[S_x_avg_name]->GetBinContent(bin);
  double Q_S_y_avg = m_profiles[S_y_avg_name]->GetBinContent(bin);
  double Q_N_x_avg = m_profiles[N_x_avg_name]->GetBinContent(bin);
  double Q_N_y_avg = m_profiles[N_y_avg_name]->GetBinContent(bin);

  m_correction_data[cent_bin][h_idx][(size_t)QVecShared::Subdetector::S].avg_Q = {Q_S_x_avg, Q_S_y_avg};
  m_correction_data[cent_bin][h_idx][(size_t)QVecShared::Subdetector::N].avg_Q = {Q_N_x_avg, Q_N_y_avg};

  std::cout << std::format(
      "Centrality Bin: {}, "
      "Harmonic: {}, "
      "Q_S_x_avg: {:13.10f}, "
      "Q_S_y_avg: {:13.10f}, "
      "Q_N_x_avg: {:13.10f}, "
      "Q_N_y_avg: {:13.10f}",
      cent_bin,
      n,
      Q_S_x_avg,
      Q_S_y_avg,
      Q_N_x_avg,
      Q_N_y_avg) << std::endl;
}

void QVecCalib::compute_recentering(size_t cent_bin, int h_idx)
{
  int n = m_harmonics[h_idx];

  std::string S_x_corr_avg_name = QVecShared::get_hist_name("S", "x", n, "_corr");
  std::string S_y_corr_avg_name = QVecShared::get_hist_name("S", "y", n, "_corr");
  std::string N_x_corr_avg_name = QVecShared::get_hist_name("N", "x", n, "_corr");
  std::string N_y_corr_avg_name = QVecShared::get_hist_name("N", "y", n, "_corr");

  int bin = static_cast<int>(cent_bin + 1);

  double Q_S_x_corr_avg = m_profiles[S_x_corr_avg_name]->GetBinContent(bin);
  double Q_S_y_corr_avg = m_profiles[S_y_corr_avg_name]->GetBinContent(bin);
  double Q_N_x_corr_avg = m_profiles[N_x_corr_avg_name]->GetBinContent(bin);
  double Q_N_y_corr_avg = m_profiles[N_y_corr_avg_name]->GetBinContent(bin);

  // -- Compute 2nd Order Correction --
  std::string S_xx_avg_name = QVecShared::get_hist_name("S", "xx", n);
  std::string S_yy_avg_name = QVecShared::get_hist_name("S", "yy", n);
  std::string S_xy_avg_name = QVecShared::get_hist_name("S", "xy", n);
  std::string N_xx_avg_name = QVecShared::get_hist_name("N", "xx", n);
  std::string N_yy_avg_name = QVecShared::get_hist_name("N", "yy", n);
  std::string N_xy_avg_name = QVecShared::get_hist_name("N", "xy", n);

  double Q_S_xx_avg = m_profiles[S_xx_avg_name]->GetBinContent(bin);
  double Q_S_yy_avg = m_profiles[S_yy_avg_name]->GetBinContent(bin);
  double Q_S_xy_avg = m_profiles[S_xy_avg_name]->GetBinContent(bin);
  double Q_N_xx_avg = m_profiles[N_xx_avg_name]->GetBinContent(bin);
  double Q_N_yy_avg = m_profiles[N_yy_avg_name]->GetBinContent(bin);
  double Q_N_xy_avg = m_profiles[N_xy_avg_name]->GetBinContent(bin);

  // -- Compute NS Matrix --
  std::string NS_xx_avg_name = QVecShared::get_hist_name("NS", "xx", n);
  std::string NS_yy_avg_name = QVecShared::get_hist_name("NS", "yy", n);
  std::string NS_xy_avg_name = QVecShared::get_hist_name("NS", "xy", n);

  double Q_NS_xx_avg = m_profiles[NS_xx_avg_name]->GetBinContent(bin);
  double Q_NS_yy_avg = m_profiles[NS_yy_avg_name]->GetBinContent(bin);
  double Q_NS_xy_avg = m_profiles[NS_xy_avg_name]->GetBinContent(bin);

  m_correction_data[cent_bin][h_idx][(size_t)QVecShared::Subdetector::NS].X_matrix = calculate_flattening_matrix(Q_NS_xx_avg, Q_NS_yy_avg, Q_NS_xy_avg, n, cent_bin, "NS");

  for (size_t det_idx = 0; det_idx < 2; ++det_idx)
  {
    double xx = (det_idx == 0) ? Q_S_xx_avg : Q_N_xx_avg;
    double yy = (det_idx == 0) ? Q_S_yy_avg : Q_N_yy_avg;
    double xy = (det_idx == 0) ? Q_S_xy_avg : Q_N_xy_avg;

    std::string label = (det_idx == 0) ? "S" : "N";

    m_correction_data[cent_bin][h_idx][det_idx].X_matrix = calculate_flattening_matrix(xx, yy, xy, n, cent_bin, label);
  }

  std::cout << std::format(
      "Centrality Bin: {}, "
      "Harmonic: {}, "
      "Q_S_x_corr_avg: {:13.10f}, "
      "Q_S_y_corr_avg: {:13.10f}, "
      "Q_N_x_corr_avg: {:13.10f}, "
      "Q_N_y_corr_avg: {:13.10f}, "
      "Q_S_xx_avg / Q_S_yy_avg: {:13.10f}, "
      "Q_N_xx_avg / Q_N_yy_avg: {:13.10f}, "
      "Q_NS_xx_avg / Q_NS_yy_avg: {:13.10f}, "
      "Q_S_xy_avg: {:13.10f}, "
      "Q_N_xy_avg: {:13.10f}, "
      "Q_NS_xy_avg: {:13.10f}",
      cent_bin,
      n,
      Q_S_x_corr_avg,
      Q_S_y_corr_avg,
      Q_N_x_corr_avg,
      Q_N_y_corr_avg,
      Q_S_xx_avg / Q_S_yy_avg,
      Q_N_xx_avg / Q_N_yy_avg,
      Q_NS_xx_avg / Q_NS_yy_avg,
      Q_S_xy_avg,
      Q_N_xy_avg,
      Q_NS_xy_avg) << std::endl;
}

void QVecCalib::print_flattening(size_t cent_bin, int n) const
{
  std::string S_x_corr2_avg_name = QVecShared::get_hist_name("S", "x", n, "_corr2");
  std::string S_y_corr2_avg_name = QVecShared::get_hist_name("S", "y", n, "_corr2");
  std::string N_x_corr2_avg_name = QVecShared::get_hist_name("N", "x", n, "_corr2");
  std::string N_y_corr2_avg_name = QVecShared::get_hist_name("N", "y", n, "_corr2");

  std::string S_xx_corr_avg_name = QVecShared::get_hist_name("S", "xx", n, "_corr");
  std::string S_yy_corr_avg_name = QVecShared::get_hist_name("S", "yy", n, "_corr");
  std::string S_xy_corr_avg_name = QVecShared::get_hist_name("S", "xy", n, "_corr");
  std::string N_xx_corr_avg_name = QVecShared::get_hist_name("N", "xx", n, "_corr");
  std::string N_yy_corr_avg_name = QVecShared::get_hist_name("N", "yy", n, "_corr");
  std::string N_xy_corr_avg_name = QVecShared::get_hist_name("N", "xy", n, "_corr");

  std::string NS_xx_corr_avg_name = QVecShared::get_hist_name("NS", "xx", n, "_corr");
  std::string NS_yy_corr_avg_name = QVecShared::get_hist_name("NS", "yy", n, "_corr");
  std::string NS_xy_corr_avg_name = QVecShared::get_hist_name("NS", "xy", n, "_corr");

  int bin = static_cast<int>(cent_bin + 1);

  double Q_S_x_corr2_avg = m_profiles.at(S_x_corr2_avg_name)->GetBinContent(bin);
  double Q_S_y_corr2_avg = m_profiles.at(S_y_corr2_avg_name)->GetBinContent(bin);
  double Q_N_x_corr2_avg = m_profiles.at(N_x_corr2_avg_name)->GetBinContent(bin);
  double Q_N_y_corr2_avg = m_profiles.at(N_y_corr2_avg_name)->GetBinContent(bin);

  double Q_S_xx_corr_avg = m_profiles.at(S_xx_corr_avg_name)->GetBinContent(bin);
  double Q_S_yy_corr_avg = m_profiles.at(S_yy_corr_avg_name)->GetBinContent(bin);
  double Q_S_xy_corr_avg = m_profiles.at(S_xy_corr_avg_name)->GetBinContent(bin);
  double Q_N_xx_corr_avg = m_profiles.at(N_xx_corr_avg_name)->GetBinContent(bin);
  double Q_N_yy_corr_avg = m_profiles.at(N_yy_corr_avg_name)->GetBinContent(bin);
  double Q_N_xy_corr_avg = m_profiles.at(N_xy_corr_avg_name)->GetBinContent(bin);

  double Q_NS_xx_corr_avg = m_profiles.at(NS_xx_corr_avg_name)->GetBinContent(bin);
  double Q_NS_yy_corr_avg = m_profiles.at(NS_yy_corr_avg_name)->GetBinContent(bin);
  double Q_NS_xy_corr_avg = m_profiles.at(NS_xy_corr_avg_name)->GetBinContent(bin);

  std::cout << std::format(
      "Centrality Bin: {}, "
      "Harmonic: {}, "
      "Q_S_x_corr2_avg: {:13.10f}, "
      "Q_S_y_corr2_avg: {:13.10f}, "
      "Q_N_x_corr2_avg: {:13.10f}, "
      "Q_N_y_corr2_avg: {:13.10f}, "
      "Q_S_xx_corr_avg / Q_S_yy_corr_avg: {:13.10f}, "
      "Q_N_xx_corr_avg / Q_N_yy_corr_avg: {:13.10f}, "
      "Q_NS_xx_corr_avg / Q_NS_yy_corr_avg: {:13.10f}, "
      "Q_S_xy_corr_avg: {:13.10f}, "
      "Q_N_xy_corr_avg: {:13.10f}, "
      "Q_NS_xy_corr_avg: {:13.10f}",
      cent_bin,
      n,
      Q_S_x_corr2_avg,
      Q_S_y_corr2_avg,
      Q_N_x_corr2_avg,
      Q_N_y_corr2_avg,
      Q_S_xx_corr_avg / Q_S_yy_corr_avg,
      Q_N_xx_corr_avg / Q_N_yy_corr_avg,
      Q_NS_xx_corr_avg / Q_NS_yy_corr_avg,
      Q_S_xy_corr_avg,
      Q_N_xy_corr_avg,
      Q_NS_xy_corr_avg) << std::endl;
}

void QVecCalib::write_cdb()
{
  std::error_code ec;
  if (std::filesystem::create_directories(m_cdb_output_dir, ec))
  {
    std::cout << "Success: Directory " << m_cdb_output_dir << " created" << std::endl;
  }
  else if (ec)
  {
    throw std::runtime_error(std::format("Failed to create directory {}: {}", m_cdb_output_dir, ec.message()));
  }
  else
  {
    std::cout << "Info: Directory " << m_cdb_output_dir << " already exists." << std::endl;
  }

  write_cdb_BadTowers();
  write_cdb_EventPlane();
}

void QVecCalib::write_cdb_BadTowers()
{
  std::cout << "Writing Bad Towers CDB" << std::endl;

  std::string payload = "SEPD_HotMap";
  std::string fieldname_status = "status";
  std::string fieldname_sigma = "SEPD_sigma";
  std::string output_file = std::format("{}/{}-{}-{}.root", m_cdb_output_dir, payload, m_dst_tag, m_runnumber);

  CDBTTree cdbttree(output_file);

  for (int channel = 0; channel < QVecShared::SEPD_CHANNELS; ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(channel);
    int status = hSEPD_Bad_Channels->GetBinContent(channel+1);

    float sigma = 0;

    // Hot
    if (status == static_cast<int>(QVecShared::ChannelStatus::Hot))
    {
      sigma = SIGMA_HOT;
    }

    // Cold
    else if (status == static_cast<int>(QVecShared::ChannelStatus::Cold))
    {
      sigma = SIGMA_COLD;
    }

    cdbttree.SetIntValue(key, fieldname_status, status);
    cdbttree.SetFloatValue(key, fieldname_sigma, sigma);
  }

  std::cout << "Saving CDB: " << payload << " to " << output_file << std::endl;

  cdbttree.Commit();
  cdbttree.WriteCDBTTree();
}

void QVecCalib::write_cdb_EventPlane()
{
  std::cout << "Writing Event Plane CDB" << std::endl;

  std::string payload = "SEPD_EventPlaneCalib";
  std::string output_file = std::format("{}/{}-{}-{}.root", m_cdb_output_dir, payload, m_dst_tag, m_runnumber);

  CDBTTree cdbttree(output_file);

  using SD = QVecShared::Subdetector;

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];

    // Define lambdas to generate field names consistently
    auto field = [&](const std::string& det, const std::string& var)
    {
      return std::format("Q_{}_{}_{}_avg", det, var, n);
    };

    for (size_t cent_bin = 0; cent_bin < m_cent_bins; ++cent_bin)
    {
      int key = static_cast<int>(cent_bin);

      // Iterate through all subdetectors (S, N, NS) using the Enum Count
      for (size_t d = 0; d < static_cast<size_t>(SD::Count); ++d)
      {
        auto det_enum = static_cast<SD>(d);

        // Map enum to the string labels used in the CDB field names
        std::string det_label;
        switch (det_enum)
        {
        case SD::S:
          det_label = "S";
          break;
        case SD::N:
          det_label = "N";
          break;
        case SD::NS:
          det_label = "NS";
          break;
        default:
          continue;
        }

        const auto& data = m_correction_data[cent_bin][h_idx][d];
        // 1st Order Moments (Recentering) - Skip for NS as it is a combined vector
        if (det_enum != SD::NS)
        {
          cdbttree.SetDoubleValue(key, field(det_label, "x"), data.avg_Q.x);
          cdbttree.SetDoubleValue(key, field(det_label, "y"), data.avg_Q.y);
        }

        // 2nd Order Moments (Flattening) - Applicable to S, N, and NS
        cdbttree.SetDoubleValue(key, field(det_label, "xx"), data.avg_Q_xx);
        cdbttree.SetDoubleValue(key, field(det_label, "yy"), data.avg_Q_yy);
        cdbttree.SetDoubleValue(key, field(det_label, "xy"), data.avg_Q_xy);
      }
    }
  }

  std::cout << "Saving CDB: " << payload << " to " << output_file << std::endl;

  cdbttree.Commit();
  cdbttree.WriteCDBTTree();
}

//____________________________________________________________________________..
int QVecCalib::End([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "QVecCalib::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  std::cout << "\n--- Event Counter Summary ---" << std::endl;
  std::cout << "Bad Centrality/sEPD corr: " << m_event_counters.bad_centrality_sepd_correlation << std::endl;
  std::cout << "Zero sEPD Charge:         " << m_event_counters.zero_sepd_total_charge << std::endl;
  std::cout << "Total Events Seen:        " << m_event << std::endl;
  std::cout << "-----------------------------\n" << std::endl;

  for (size_t cent_bin = 0; cent_bin < m_cent_bins; ++cent_bin)
  {
    for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
    {
      int n = m_harmonics[h_idx];

      if (m_pass == Pass::ComputeRecentering)
      {
        compute_averages(cent_bin, h_idx);
      }

      else if (m_pass == Pass::ApplyRecentering)
      {
        compute_recentering(cent_bin, h_idx);
      }

      else if (m_pass == Pass::ApplyFlattening)
      {
        print_flattening(cent_bin, n);
      }
    }
  }

  if (m_pass == Pass::ApplyFlattening)
  {
    write_cdb();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
