#include "QVecCalib.h"

// ====================================================================
// sPHENIX Includes
// ====================================================================
#include <calobase/TowerInfoDefs.h>

// ====================================================================
// Standard C++ Includes
// ====================================================================
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <numbers>

// ====================================================================
// ROOT Includes
// ====================================================================
#include <TTree.h>

std::unique_ptr<TChain> QVecCalib::setupTChain(const std::string& input_filepath, const std::string& tree_name_in_file)
{
  // 1. Pre-check: Does the file exist at all? (C++17 filesystem or traditional fstream)
  if (!std::filesystem::exists(input_filepath))
  {
    std::cout << "Error: Input file does not exist: " << input_filepath << std::endl;
    return nullptr;  // Return a null unique_ptr indicating failure
  }

  // 2. Open the file to check for the TTree directly
  // Use TFile::Open and unique_ptr for robust file handling (RAII)
  std::unique_ptr<TFile> file_checker(TFile::Open(input_filepath.c_str(), "READ"));

  if (!file_checker || file_checker->IsZombie())
  {
    std::cout << "Error: Could not open file " << input_filepath << " to check for TTree." << std::endl;
    return nullptr;
  }

  // Check if the TTree exists in the file
  // Get() returns a TObject*, which can be cast to TTree*.
  // If the object doesn't exist or isn't a TTree, Get() returns nullptr.
  TTree* tree_obj = dynamic_cast<TTree*>(file_checker->Get(tree_name_in_file.c_str()));
  if (!tree_obj)
  {
    std::cout << "Error: TTree '" << tree_name_in_file << "' not found in file " << input_filepath << std::endl;
    return nullptr;
  }
  // File will be automatically closed by file_checker's unique_ptr destructor

  // 3. If everything checks out, create and configure the TChain
  std::unique_ptr<TChain> chain = std::make_unique<TChain>(tree_name_in_file.c_str());
  if (!chain)
  {  // Check if make_unique failed (e.g. out of memory)
    std::cout << "Error: Could not create TChain object." << std::endl;
    return nullptr;
  }

  chain->Add(input_filepath.c_str());

  // 4. Verify TChain's state (optional but good final check)
  // GetEntries() will be -1 if no valid trees were added.
  if (chain->GetEntries() == 0)
  {
    std::cout << "Warning: TChain has 0 entries after adding file. This might indicate a problem." << std::endl;
    // Depending on your logic, you might return nullptr here too.
  }
  else
  {
    std::cout << "Successfully set up TChain for tree '" << tree_name_in_file
              << "' from file '" << input_filepath << "'. Entries: " << chain->GetEntries() << std::endl;
  }

  return chain;  // Return the successfully created and configured TChain
}

void QVecCalib::setup_chain()
{
  std::cout << "Processing... setup_chain" << std::endl;

  m_chain = setupTChain(m_input_file, "T");

  if (m_chain == nullptr)
  {
    throw std::runtime_error(std::format("Error in TChain Setup from file: {}", m_input_file));
  }

  // Setup branches
  m_chain->SetBranchStatus("*", false);
  m_chain->SetBranchStatus("event_id", true);
  m_chain->SetBranchStatus("event_centrality", true);
  m_chain->SetBranchStatus("sepd_totalcharge", true);
  m_chain->SetBranchStatus("sepd_channel", true);
  m_chain->SetBranchStatus("sepd_charge", true);
  m_chain->SetBranchStatus("sepd_phi", true);

  m_chain->SetBranchAddress("event_id", &m_event_data.event_id);
  m_chain->SetBranchAddress("event_centrality", &m_event_data.event_centrality);
  m_chain->SetBranchAddress("sepd_totalcharge", &m_event_data.sepd_totalcharge);
  m_chain->SetBranchAddress("sepd_channel", &m_event_data.sepd_channel);
  m_chain->SetBranchAddress("sepd_charge", &m_event_data.sepd_charge);
  m_chain->SetBranchAddress("sepd_phi", &m_event_data.sepd_phi);

  std::cout << "Finished... setup_chain" << std::endl;
}

void QVecCalib::process_QA_hist()
{
  TH1::AddDirectory(kFALSE);
  auto file = std::unique_ptr<TFile>(TFile::Open(m_input_hist.c_str()));

  // Check if the file was opened successfully.
  if (!file || file->IsZombie())
  {
    throw std::runtime_error(std::format("Could not open file '{}'", m_input_hist));
  }

  // Get List of Bad Channels
  process_bad_channels(file.get());

  // Get sEPD Total Charge Bounds as function of centrality
  process_sEPD_event_thresholds(file.get());
}

void QVecCalib::process_sEPD_event_thresholds(TFile* file)
{
  std::string sepd_totalcharge_centrality = "h2SEPD_totalcharge_centrality";

  auto* hist = file->Get(sepd_totalcharge_centrality.c_str());

  // Check if the hist is stored in the file
  if (hist == nullptr)
  {
    throw std::runtime_error(std::format("Cannot find hist: {}", sepd_totalcharge_centrality));
  }

  m_hists2D["h2SEPD_Charge"] = std::unique_ptr<TH2>(static_cast<TH2*>(hist->Clone("h2SEPD_Charge")));
  m_hists2D["h2SEPD_Chargev2"] = std::unique_ptr<TH2>(static_cast<TH2*>(hist->Clone("h2SEPD_Chargev2")));

  auto* h2SEPD_Charge = m_hists2D["h2SEPD_Charge"].get();
  auto* h2SEPD_Chargev2 = m_hists2D["h2SEPD_Chargev2"].get();

  auto* h2SEPD_Charge_py = h2SEPD_Charge->ProfileY("h2SEPD_Charge_py", 1, -1, "s");

  int binsx = h2SEPD_Charge->GetNbinsX();
  int binsy = h2SEPD_Charge->GetNbinsY();
  int ymin = h2SEPD_Charge->GetYaxis()->GetXmin();
  int ymax = h2SEPD_Charge->GetYaxis()->GetXmax();

  m_profiles["hSEPD_Charge_Min"] = std::make_unique<TProfile>("hSEPD_Charge_Min", "; Centrality [%]; sEPD Total Charge", binsy, ymin, ymax);
  m_profiles["hSEPD_Charge_Max"] = std::make_unique<TProfile>("hSEPD_Charge_Max", "; Centrality [%]; sEPD Total Charge", binsy, ymin, ymax);

  auto* hSEPD_Charge_Min = m_profiles["hSEPD_Charge_Min"].get();
  auto* hSEPD_Charge_Max = m_profiles["hSEPD_Charge_Max"].get();

  for (int y = 1; y <= binsy; ++y)
  {
    double cent = h2SEPD_Charge_py->GetBinCenter(y);
    double mean = h2SEPD_Charge_py->GetBinContent(y);
    double sigma = h2SEPD_Charge_py->GetBinError(y);
    double charge_low  = mean - m_sEPD_sigma_threshold * sigma;
    double charge_high = mean + m_sEPD_sigma_threshold * sigma;

    hSEPD_Charge_Min->Fill(cent, charge_low);
    hSEPD_Charge_Max->Fill(cent, charge_high);

    if (sigma == 0)
    {
      continue;
    }

    for (int x = 1; x <= binsx; ++x)
    {
      double charge = h2SEPD_Charge->GetXaxis()->GetBinCenter(x);
      double zscore = (charge - mean) / sigma;

      if (std::fabs(zscore) > m_sEPD_sigma_threshold)
      {
        h2SEPD_Chargev2->SetBinContent(x, y, 0);
      }
    }
  }
}

void QVecCalib::process_bad_channels(TFile* file)
{
  std::string sepd_charge_hist = "hSEPD_Charge";

  auto* hist = file->Get(sepd_charge_hist.c_str());

  // Check if the hist is stored in the file
  if (hist == nullptr)
  {
    throw std::runtime_error(std::format("Cannot find hist: {}", sepd_charge_hist));
  }

  auto* hSEPD_Charge = dynamic_cast<TH1*>(hist);

  int rbins = 16;
  int bins_charge = 40;

  m_hists2D["h2SEPD_South_Charge_rbin"] = std::make_unique<TH2F>("h2SEPD_South_Charge_rbin",
                                                                 "sEPD South; r_{bin}; Avg Charge",
                                                                 rbins, -0.5, rbins - 0.5,
                                                                 bins_charge, 0, bins_charge);

  m_hists2D["h2SEPD_North_Charge_rbin"] = std::make_unique<TH2F>("h2SEPD_North_Charge_rbin",
                                                                 "sEPD North; r_{bin}; Avg Charge",
                                                                 rbins, -0.5, rbins - 0.5,
                                                                 bins_charge, 0, bins_charge);

  m_hists2D["h2SEPD_South_Charge_rbinv2"] = std::make_unique<TH2F>("h2SEPD_South_Charge_rbinv2",
                                                                   "sEPD South; r_{bin}; Avg Charge",
                                                                   rbins, -0.5, rbins - 0.5,
                                                                   bins_charge, 0, bins_charge);

  m_hists2D["h2SEPD_North_Charge_rbinv2"] = std::make_unique<TH2F>("h2SEPD_North_Charge_rbinv2",
                                                                   "sEPD North; r_{bin}; Avg Charge",
                                                                   rbins, -0.5, rbins - 0.5,
                                                                   bins_charge, 0, bins_charge);

  m_profiles["h_sEPD_Bad_Channels"] = std::make_unique<TProfile>("h_sEPD_Bad_Channels", "sEPD Bad Channels; Channel; Status", QVecShared::sepd_channels, -0.5, QVecShared::sepd_channels-0.5);

  auto* h2S = m_hists2D["h2SEPD_South_Charge_rbin"].get();
  auto* h2N = m_hists2D["h2SEPD_North_Charge_rbin"].get();

  auto* h2Sv2 = m_hists2D["h2SEPD_South_Charge_rbinv2"].get();
  auto* h2Nv2 = m_hists2D["h2SEPD_North_Charge_rbinv2"].get();

  auto* hBad = m_profiles["h_sEPD_Bad_Channels"].get();

  for (int channel = 0; channel < QVecShared::sepd_channels; ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(static_cast<unsigned int>(channel));
    int rbin = static_cast<int>(TowerInfoDefs::get_epd_rbin(key));
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

  for (int channel = 0; channel < QVecShared::sepd_channels; ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(static_cast<unsigned int>(channel));
    int rbin = static_cast<int>(TowerInfoDefs::get_epd_rbin(key));
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

    if (charge < m_sEPD_min_avg_charge_threshold || std::fabs(zscore) > m_sEPD_sigma_threshold)
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
      std::cout << std::format("{:4} Channel: {:3d}, arm: {}, rbin: {:2d}, Mean: {:5.2f}, Charge: {:5.2f}, Z-Score: {:5.2f}\n", type, channel, arm, rbin, mean_charge, charge, zscore);
    }
    else
    {
      h2->Fill(rbin, charge);
    }
  }

  std::cout << std::format("Total Bad Channels: {}, Dead: {}, Hot: {}, Cold: {}\n", m_bad_channels.size(), ctr_dead, ctr_hot, ctr_cold);

  std::cout << "Finished processing Hot sEPD channels" << std::endl;
}

void QVecCalib::init_hists()
{
  unsigned int bins_psi = 126;
  double psi_low = -std::numbers::pi;
  double psi_high = std::numbers::pi;

  m_hists1D["h_Cent"] = std::make_unique<TH1F>("h_Cent", "", m_cent_bins, m_cent_low, m_cent_high);

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

    m_hists2D[name_S] = std::make_unique<TH2F>(name_S.c_str(), title_S.c_str(), m_cent_bins, m_cent_low, m_cent_high, bins_psi, psi_low, psi_high);
    m_hists2D[name_N] = std::make_unique<TH2F>(name_N.c_str(), title_N.c_str(), m_cent_bins, m_cent_low, m_cent_high, bins_psi, psi_low, psi_high);
    m_hists2D[name_NS] = std::make_unique<TH2F>(name_NS.c_str(), title_NS.c_str(), m_cent_bins, m_cent_low, m_cent_high, bins_psi, psi_low, psi_high);

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
        m_profiles[q_avg_sq_cross_name] = std::make_unique<TProfile>(q_avg_sq_cross_name.c_str(), q_avg_sq_cross_title.c_str(),
                                                                     m_cent_bins, m_cent_low, m_cent_high);
      }

      for (auto comp : m_components)
      {
        std::string comp_str = (comp == QVecShared::QComponent::X) ? "x" : "y";
        std::string name = QVecShared::get_hist_name(det_str, comp_str, n, pass_suffix);

        auto add_profile = [&](const std::string& prof_name, std::string_view label_suffix = "")
        {
          std::string title = std::format("sEPD {}; Centrality [%]; <Q_{{{},{}}}{}>", det_name, n, comp_str, label_suffix);
          m_profiles[prof_name] = std::make_unique<TProfile>(prof_name.c_str(), title.c_str(), m_cent_bins, m_cent_low, m_cent_high);
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
        m_profiles[name] = std::make_unique<TProfile>(name.c_str(), title.c_str(), m_cent_bins, m_cent_low, m_cent_high);
      }

      // Initialize Validation Profiles (Flattened NS)
      if (m_pass == Pass::ApplyFlattening)
      {
        for (const auto* comp : {"xx", "yy", "xy"})
        {
            std::string name = QVecShared::get_hist_name(det_str, comp, n, "_corr");
            std::string title = std::format("sEPD NS Corrected; Centrality [%]; <Q_{{{},{}}}^{{2}}>", n, comp);
            m_profiles[name] = std::make_unique<TProfile>(name.c_str(), title.c_str(), m_cent_bins, m_cent_low, m_cent_high);
        }
      }
    }
  }
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
  size_t cent_bin = static_cast<size_t>(m_hists1D["h_Cent"]->FindBin(cent) - 1);

  double Q_S_x_avg = m_correction_data[cent_bin][h_idx][0].avg_Q.x;
  double Q_S_y_avg = m_correction_data[cent_bin][h_idx][0].avg_Q.y;
  double Q_N_x_avg = m_correction_data[cent_bin][h_idx][1].avg_Q.x;
  double Q_N_y_avg = m_correction_data[cent_bin][h_idx][1].avg_Q.y;

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
  size_t cent_bin = static_cast<size_t>(m_hists1D["h_Cent"]->FindBin(cent) - 1);

  double Q_S_x_avg = m_correction_data[cent_bin][h_idx][0].avg_Q.x;
  double Q_S_y_avg = m_correction_data[cent_bin][h_idx][0].avg_Q.y;
  double Q_N_x_avg = m_correction_data[cent_bin][h_idx][1].avg_Q.x;
  double Q_N_y_avg = m_correction_data[cent_bin][h_idx][1].avg_Q.y;

  QVecShared::QVec q_S_corr = {q_S.x - Q_S_x_avg, q_S.y - Q_S_y_avg};
  QVecShared::QVec q_N_corr = {q_N.x - Q_N_x_avg, q_N.y - Q_N_y_avg};

  // Construct Combined Recentered Vector
  QVecShared::QVec q_NS_corr = {q_S_corr.x + q_N_corr.x, q_S_corr.y + q_N_corr.y};

  const auto& X_S = m_correction_data[cent_bin][h_idx][0].X_matrix;
  const auto& X_N = m_correction_data[cent_bin][h_idx][1].X_matrix;
  const auto& X_NS = m_correction_data[cent_bin][h_idx][IDX_NS].X_matrix;

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

  m_correction_data[cent_bin][h_idx][0].avg_Q = {Q_S_x_avg, Q_S_y_avg};
  m_correction_data[cent_bin][h_idx][1].avg_Q = {Q_N_x_avg, Q_N_y_avg};

  std::cout << std::format(
      "Centrality Bin: {}, "
      "Harmonic: {}, "
      "Q_S_x_avg: {:13.10f}, "
      "Q_S_y_avg: {:13.10f}, "
      "Q_N_x_avg: {:13.10f}, "
      "Q_N_y_avg: {:13.10f}\n",
      cent_bin,
      n,
      Q_S_x_avg,
      Q_S_y_avg,
      Q_N_x_avg,
      Q_N_y_avg);
}

std::array<std::array<double, 2>, 2> QVecCalib::calculate_flattening_matrix(double xx, double yy, double xy, int n, int cent_bin, const std::string& det_label)
{
  double D_arg = (xx * yy) - (xy * xy);
  if (D_arg <= 0)
  {
    throw std::runtime_error(std::format(
        "Invalid D-term ({}) for n={}, cent={}, det={}", D_arg, n, cent_bin, det_label));
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

  m_correction_data[cent_bin][h_idx][IDX_NS].X_matrix = calculate_flattening_matrix(Q_NS_xx_avg, Q_NS_yy_avg, Q_NS_xy_avg, n, cent_bin, "NS");

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
      "Q_NS_xy_avg: {:13.10f}\n",
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
      Q_NS_xy_avg);
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
      "Q_NS_xy_corr_avg: {:13.10f}\n",
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
      Q_NS_xy_corr_avg);
}

std::vector<QVecCalib::AverageHists> QVecCalib::prepare_average_hists()
{
  std::vector<AverageHists> hists_cache;
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

    h.S_x_avg = m_profiles.at(S_x_avg_name).get();
    h.S_y_avg = m_profiles.at(S_y_avg_name).get();
    h.N_x_avg = m_profiles.at(N_x_avg_name).get();
    h.N_y_avg = m_profiles.at(N_y_avg_name).get();

    h.Psi_S = m_hists2D.at(psi_S_name).get();
    h.Psi_N = m_hists2D.at(psi_N_name).get();
    h.Psi_NS = m_hists2D.at(psi_NS_name).get();

    hists_cache.push_back(h);
  }

  return hists_cache;
}

bool QVecCalib::process_sEPD()
{
  size_t nChannels = m_event_data.sepd_channel->size();

  double sepd_total_charge_south = 0;
  double sepd_total_charge_north = 0;

  // Loop over all sEPD Channels
  for (size_t idx = 0; idx < nChannels; ++idx)
  {
    int channel = m_event_data.sepd_channel->at(idx);
    double charge = m_event_data.sepd_charge->at(idx);
    double phi = m_event_data.sepd_phi->at(idx);

    // Skip Bad Channels
    if (m_bad_channels.contains(channel))
    {
      continue;
    }

    unsigned int key = TowerInfoDefs::encode_epd(static_cast<unsigned int>(channel));
    unsigned int arm = TowerInfoDefs::get_epd_arm(key);

    // arm = 0: South
    // arm = 1: North
    double& sepd_total_charge = (arm == 0) ? sepd_total_charge_south : sepd_total_charge_north;

    // Compute total charge for the respective sEPD arm
    sepd_total_charge += charge;

    // Compute Raw Q vectors for each harmonic and respective arm
    for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
    {
      int n = m_harmonics[h_idx];
      m_event_data.q_vectors[h_idx][arm].x += charge * std::cos(n * phi);
      m_event_data.q_vectors[h_idx][arm].y += charge * std::sin(n * phi);
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
      m_event_data.q_vectors[h_idx][det_idx].x /= sepd_total_charge;
      m_event_data.q_vectors[h_idx][det_idx].y /= sepd_total_charge;
    }
  }

  return true;
}

std::vector<QVecCalib::RecenterHists> QVecCalib::prepare_recenter_hists()
{
  std::vector<RecenterHists> hists_cache;
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

    h.S_x_corr_avg = m_profiles.at(S_x_corr_avg_name).get();
    h.S_y_corr_avg = m_profiles.at(S_y_corr_avg_name).get();
    h.N_x_corr_avg = m_profiles.at(N_x_corr_avg_name).get();
    h.N_y_corr_avg = m_profiles.at(N_y_corr_avg_name).get();

    h.S_xx_avg = m_profiles.at(S_xx_avg_name).get();
    h.S_yy_avg = m_profiles.at(S_yy_avg_name).get();
    h.S_xy_avg = m_profiles.at(S_xy_avg_name).get();
    h.N_xx_avg = m_profiles.at(N_xx_avg_name).get();
    h.N_yy_avg = m_profiles.at(N_yy_avg_name).get();
    h.N_xy_avg = m_profiles.at(N_xy_avg_name).get();

    h.NS_xx_avg = m_profiles.at(NS_xx_avg_name).get();
    h.NS_yy_avg = m_profiles.at(NS_yy_avg_name).get();
    h.NS_xy_avg = m_profiles.at(NS_xy_avg_name).get();

    h.Psi_S_corr = m_hists2D.at(psi_S_name).get();
    h.Psi_N_corr = m_hists2D.at(psi_N_name).get();
    h.Psi_NS_corr = m_hists2D.at(psi_NS_name).get();

    hists_cache.push_back(h);
  }

  return hists_cache;
}

std::vector<QVecCalib::FlatteningHists> QVecCalib::prepare_flattening_hists()
{
  std::vector<FlatteningHists> hists_cache;
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

    h.S_x_corr2_avg = m_profiles.at(S_x_corr2_avg_name).get();
    h.S_y_corr2_avg = m_profiles.at(S_y_corr2_avg_name).get();
    h.N_x_corr2_avg = m_profiles.at(N_x_corr2_avg_name).get();
    h.N_y_corr2_avg = m_profiles.at(N_y_corr2_avg_name).get();

    h.S_xx_corr_avg = m_profiles.at(S_xx_corr_avg_name).get();
    h.S_yy_corr_avg = m_profiles.at(S_yy_corr_avg_name).get();
    h.S_xy_corr_avg = m_profiles.at(S_xy_corr_avg_name).get();

    h.N_xx_corr_avg = m_profiles.at(N_xx_corr_avg_name).get();
    h.N_yy_corr_avg = m_profiles.at(N_yy_corr_avg_name).get();
    h.N_xy_corr_avg = m_profiles.at(N_xy_corr_avg_name).get();

    h.NS_xx_corr_avg = m_profiles.at(NS_xx_corr_avg_name).get();
    h.NS_yy_corr_avg = m_profiles.at(NS_yy_corr_avg_name).get();
    h.NS_xy_corr_avg = m_profiles.at(NS_xy_corr_avg_name).get();

    h.Psi_S_corr2 = m_hists2D.at(psi_S_name).get();
    h.Psi_N_corr2 = m_hists2D.at(psi_N_name).get();
    h.Psi_NS_corr2 = m_hists2D.at(psi_NS_name).get();

    hists_cache.push_back(h);
  }

  return hists_cache;
}

bool QVecCalib::process_event_check()
{
  auto* hSEPD_Charge_Min = m_profiles["hSEPD_Charge_Min"].get();
  auto* hSEPD_Charge_Max = m_profiles["hSEPD_Charge_Max"].get();

  double cent = m_event_data.event_centrality;
  int cent_bin = hSEPD_Charge_Min->FindBin(cent);

  double sepd_totalcharge = m_event_data.sepd_totalcharge;

  double sepd_totalcharge_min = hSEPD_Charge_Min->GetBinContent(cent_bin);
  double sepd_totalcharge_max = hSEPD_Charge_Max->GetBinContent(cent_bin);

  return sepd_totalcharge > sepd_totalcharge_min && sepd_totalcharge < sepd_totalcharge_max;
}

void QVecCalib::run_event_loop()
{
  std::cout << std::format("Pass: {}\n", static_cast<uint8_t>(m_pass));

  long long n_entries = m_chain->GetEntries();
  if (m_events_to_process > 0)
  {
    n_entries = std::min(m_events_to_process, n_entries);
  }

  std::vector<AverageHists> average_hists;
  std::vector<RecenterHists> recenter_hists;
  std::vector<FlatteningHists> flattening_hists;

  if (m_pass == Pass::ComputeRecentering)
  {
    average_hists = prepare_average_hists();
  }
  else if (m_pass == Pass::ApplyRecentering)
  {
    recenter_hists = prepare_recenter_hists();
  }
  else if (m_pass == Pass::ApplyFlattening)
  {
    flattening_hists = prepare_flattening_hists();
  }

  std::map<std::string, int> ctr;
  // Event Loop
  for (long long i = 0; i < n_entries; ++i)
  {
    // Load Event Data from TChain
    m_chain->GetEntry(i);
    m_event_data.reset();

    if (i % PROGRESS_REPORT_INTERVAL == 0)
    {
      std::cout << std::format("Processing {}/{}: {:.2f} %", i, n_entries, static_cast<double>(i) / n_entries * 100.) << std::endl;
    }

    double cent = m_event_data.event_centrality;

    int cent_bin_int = m_hists1D["h_Cent"]->FindBin(cent) - 1;

    // ensure centrality is valid
    if (cent_bin_int < 0 || static_cast<size_t>(cent_bin_int) >= m_cent_bins)
    {
      std::cout << std::format("Weird Centrality: {}, Skipping Event: {}\n", cent, m_event_data.event_id);
      ++ctr["invalid_cent_bin"];
      continue;
    }

    bool isGood = process_event_check();

    // Skip Events with non correlation between centrality and sEPD
    if (!isGood)
    {
      ++ctr["bad_centrality_sepd_correlation"];
      continue;
    }

    isGood = process_sEPD();

    // Skip Events with Zero sEPD Total Charge in either arm
    if (!isGood)
    {
      ++ctr["zero_sepd_total_charge"];
      continue;
    }

    m_hists1D["h_Cent"]->Fill(cent);

    for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
    {
      const auto& q_S = m_event_data.q_vectors[h_idx][0];  // 0 for South
      const auto& q_N = m_event_data.q_vectors[h_idx][1];  // 1 for North

      // --- First Pass: Derive 1st Order ---
      if (m_pass == Pass::ComputeRecentering)
      {
        process_averages(cent, q_S, q_N, average_hists[h_idx]);
      }

      // --- Second Pass: Apply 1st Order, Derive 2nd Order ---
      else if (m_pass == Pass::ApplyRecentering)
      {
        process_recentering(cent, h_idx, q_S, q_N, recenter_hists[h_idx]);
      }

      // --- Third Pass: Apply 2nd Order, Validate ---
      else if (m_pass == Pass::ApplyFlattening)
      {
        process_flattening(cent, h_idx, q_S, q_N, flattening_hists[h_idx]);
      }
    }
  }

  std::cout << "Skipped Event Types\n";
  for (const auto& [name, events] : ctr)
  {
    std::cout << std::format("{}: {}, {:.2f} %\n", name, events, static_cast<double>(events) / n_entries * 100.);
  }

  // ---------------

  std::cout << std::format("{:#<20}\n", "");
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

  std::cout << "Event loop finished." << std::endl;
}

template <typename T>
std::unique_ptr<T> QVecCalib::load_and_clone(TFile* file, const std::string& name) {
  auto* obj = dynamic_cast<T*>(file->Get(name.c_str()));
  if (!obj)
  {
    throw std::runtime_error(std::format("Could not find histogram '{}' in file '{}'", name, file->GetName()));
  }
  return std::unique_ptr<T>(static_cast<T*>(obj->Clone()));
}

void QVecCalib::load_correction_data()
{
  TH1::AddDirectory(kFALSE);

  auto file = std::unique_ptr<TFile>(TFile::Open(m_input_Q_calib.c_str()));

  // Check if the file was opened successfully.
  if (!file || file->IsZombie())
  {
    throw std::runtime_error(std::format("Could not open file '{}'", m_input_Q_calib));
  }

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];

    std::string S_x_avg_name = QVecShared::get_hist_name("S", "x", n);
    std::string S_y_avg_name = QVecShared::get_hist_name("S", "y", n);
    std::string N_x_avg_name = QVecShared::get_hist_name("N", "x", n);
    std::string N_y_avg_name = QVecShared::get_hist_name("N", "y", n);

    m_profiles[S_x_avg_name] = load_and_clone<TProfile>(file.get(), S_x_avg_name);
    m_profiles[S_y_avg_name] = load_and_clone<TProfile>(file.get(), S_y_avg_name);
    m_profiles[N_x_avg_name] = load_and_clone<TProfile>(file.get(), N_x_avg_name);
    m_profiles[N_y_avg_name] = load_and_clone<TProfile>(file.get(), N_y_avg_name);

    std::string S_xx_avg_name = QVecShared::get_hist_name("S", "xx", n);
    std::string S_yy_avg_name = QVecShared::get_hist_name("S", "yy", n);
    std::string S_xy_avg_name = QVecShared::get_hist_name("S", "xy", n);
    std::string N_xx_avg_name = QVecShared::get_hist_name("N", "xx", n);
    std::string N_yy_avg_name = QVecShared::get_hist_name("N", "yy", n);
    std::string N_xy_avg_name = QVecShared::get_hist_name("N", "xy", n);

    std::string NS_xx_avg_name = QVecShared::get_hist_name("NS", "xx", n);
    std::string NS_yy_avg_name = QVecShared::get_hist_name("NS", "yy", n);
    std::string NS_xy_avg_name = QVecShared::get_hist_name("NS", "xy", n);

    if(m_pass == Pass::ApplyFlattening)
    {
      m_profiles[S_xx_avg_name] = load_and_clone<TProfile>(file.get(), S_xx_avg_name);
      m_profiles[S_yy_avg_name] = load_and_clone<TProfile>(file.get(), S_yy_avg_name);
      m_profiles[S_xy_avg_name] = load_and_clone<TProfile>(file.get(), S_xy_avg_name);

      m_profiles[N_xx_avg_name] = load_and_clone<TProfile>(file.get(), N_xx_avg_name);
      m_profiles[N_yy_avg_name] = load_and_clone<TProfile>(file.get(), N_yy_avg_name);
      m_profiles[N_xy_avg_name] = load_and_clone<TProfile>(file.get(), N_xy_avg_name);

      m_profiles[NS_xx_avg_name] = load_and_clone<TProfile>(file.get(), NS_xx_avg_name);
      m_profiles[NS_yy_avg_name] = load_and_clone<TProfile>(file.get(), NS_yy_avg_name);
      m_profiles[NS_xy_avg_name] = load_and_clone<TProfile>(file.get(), NS_xy_avg_name);
    }

    size_t south_idx = static_cast<size_t>(QVecShared::Subdetector::S);
    size_t north_idx = static_cast<size_t>(QVecShared::Subdetector::N);

    for (size_t cent_bin = 0; cent_bin < m_cent_bins; ++cent_bin)
    {
      int bin = static_cast<int>(cent_bin) + 1;

      double Q_S_x_avg = m_profiles[S_x_avg_name]->GetBinContent(bin);
      double Q_S_y_avg = m_profiles[S_y_avg_name]->GetBinContent(bin);
      double Q_N_x_avg = m_profiles[N_x_avg_name]->GetBinContent(bin);
      double Q_N_y_avg = m_profiles[N_y_avg_name]->GetBinContent(bin);

      // Recentering Params
      m_correction_data[cent_bin][h_idx][south_idx].avg_Q = {Q_S_x_avg, Q_S_y_avg};
      m_correction_data[cent_bin][h_idx][north_idx].avg_Q = {Q_N_x_avg, Q_N_y_avg};

      if (m_pass == Pass::ApplyFlattening)
      {
        double Q_S_xx_avg = m_profiles[S_xx_avg_name]->GetBinContent(bin);
        double Q_S_yy_avg = m_profiles[S_yy_avg_name]->GetBinContent(bin);
        double Q_S_xy_avg = m_profiles[S_xy_avg_name]->GetBinContent(bin);

        double Q_N_xx_avg = m_profiles[N_xx_avg_name]->GetBinContent(bin);
        double Q_N_yy_avg = m_profiles[N_yy_avg_name]->GetBinContent(bin);
        double Q_N_xy_avg = m_profiles[N_xy_avg_name]->GetBinContent(bin);

        double Q_NS_xx_avg = m_profiles[NS_xx_avg_name]->GetBinContent(bin);
        double Q_NS_yy_avg = m_profiles[NS_yy_avg_name]->GetBinContent(bin);
        double Q_NS_xy_avg = m_profiles[NS_xy_avg_name]->GetBinContent(bin);

        m_correction_data[cent_bin][h_idx][IDX_NS].X_matrix = calculate_flattening_matrix(Q_NS_xx_avg, Q_NS_yy_avg, Q_NS_xy_avg, n, cent_bin, "NS");

        // Flattening Params
        for (size_t det_idx = 0; det_idx < 2; ++det_idx)
        {
          double xx = (det_idx == 0) ? Q_S_xx_avg : Q_N_xx_avg;
          double yy = (det_idx == 0) ? Q_S_yy_avg : Q_N_yy_avg;
          double xy = (det_idx == 0) ? Q_S_xy_avg : Q_N_xy_avg;

          std::string label = (det_idx == 0) ? "S" : "N";

          m_correction_data[cent_bin][h_idx][det_idx].X_matrix = calculate_flattening_matrix(xx, yy, xy, n, cent_bin, label);
        }
      }
    }
  }
}

void QVecCalib::process_events()
{
  if (m_pass == Pass::ApplyRecentering || m_pass == Pass::ApplyFlattening)
  {
    load_correction_data();
  }

  run_event_loop();
}

void QVecCalib::save_results() const
{
  std::filesystem::create_directories(m_output_dir);

  std::filesystem::path input_path(m_input_file);
  std::string output_stem = input_path.stem().string();
  std::string output_filename = std::format("{}/Q-vec-corr_Pass-{}_{}.root", m_output_dir, static_cast<int>(m_pass), output_stem);

  auto output_file = std::make_unique<TFile>(output_filename.c_str(), "RECREATE");

  for (const auto& [name, hist] : m_hists1D)
  {
    std::cout << std::format("Saving 1D: {}\n", name);
    hist->Write();
  }
  for (const auto& [name, hist] : m_hists2D)
  {
    std::cout << std::format("Saving 2D: {}\n", name);
    hist->Write();
  }
  for (const auto& [name, hist] : m_profiles)
  {
    std::cout << std::format("Saving Profile: {}\n", name);
    hist->Write();
  }
  output_file->Close();

  std::cout << std::format("Results saved to: {}", output_filename) << std::endl;
}
