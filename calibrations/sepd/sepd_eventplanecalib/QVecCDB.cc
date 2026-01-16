#include "QVecCDB.h"

// ====================================================================
// sPHENIX Includes
// ====================================================================
#include <calobase/TowerInfoDefs.h>
#include <cdbobjects/CDBTTree.h>

// ====================================================================
// ROOT Includes
// ====================================================================
#include <TProfile.h>

// ====================================================================
// Standard C++ Includes
// ====================================================================
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <numbers>

template <typename T>
std::unique_ptr<T> QVecCDB::load_and_clone(const std::string& name) {
  auto* obj = dynamic_cast<T*>(m_tfile->Get(name.c_str()));
  if (!obj)
  {
    throw std::runtime_error(std::format("Could not find histogram '{}' in file '{}'", name, m_tfile->GetName()));
  }
  return std::unique_ptr<T>(static_cast<T*>(obj->Clone()));
}

QVecShared::CorrectionMoments& QVecCDB::getData(size_t h_idx, size_t cent_bin, QVecShared::Subdetector sub) {
    return m_correction_data[h_idx][cent_bin][static_cast<size_t>(sub)];
}

void QVecCDB::load_data()
{
  m_tfile = std::unique_ptr<TFile>(TFile::Open(m_input_file.c_str()));

  // Check if the file was opened successfully.
  if (!m_tfile || m_tfile->IsZombie())
  {
    throw std::runtime_error(std::format("Could not open file '{}'", m_input_file));
  }

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    load_correction_data(h_idx);
  }
}

void QVecCDB::load_correction_data(size_t h_idx)
{
  int n = m_harmonics[h_idx];

  // Load recentering terms (x, y)
  auto pS_x = load_and_clone<TProfile>(QVecShared::get_hist_name("S", "x", n));
  auto pS_y = load_and_clone<TProfile>(QVecShared::get_hist_name("S", "y", n));
  auto pN_x = load_and_clone<TProfile>(QVecShared::get_hist_name("N", "x", n));
  auto pN_y = load_and_clone<TProfile>(QVecShared::get_hist_name("N", "y", n));

  // Load flattening terms (xx, yy, xy)
  auto pS_xx = load_and_clone<TProfile>(QVecShared::get_hist_name("S", "xx", n));
  auto pS_yy = load_and_clone<TProfile>(QVecShared::get_hist_name("S", "yy", n));
  auto pS_xy = load_and_clone<TProfile>(QVecShared::get_hist_name("S", "xy", n));
  auto pN_xx = load_and_clone<TProfile>(QVecShared::get_hist_name("N", "xx", n));
  auto pN_yy = load_and_clone<TProfile>(QVecShared::get_hist_name("N", "yy", n));
  auto pN_xy = load_and_clone<TProfile>(QVecShared::get_hist_name("N", "xy", n));

  for (size_t cent_bin = 0; cent_bin < m_cent_bins; ++cent_bin)
  {
    int bin = static_cast<int>(cent_bin) + 1; // ROOT bins start at 1

    // South
    auto& dataS = getData(h_idx, cent_bin, QVecShared::Subdetector::S);
    dataS.avg_Q = {pS_x->GetBinContent(bin), pS_y->GetBinContent(bin)};
    dataS.avg_Q_xx = pS_xx->GetBinContent(bin);
    dataS.avg_Q_yy = pS_yy->GetBinContent(bin);
    dataS.avg_Q_xy = pS_xy->GetBinContent(bin);

    // North
    auto& dataN = getData(h_idx, cent_bin, QVecShared::Subdetector::N);
    dataN.avg_Q = {pN_x->GetBinContent(bin), pN_y->GetBinContent(bin)};
    dataN.avg_Q_xx = pN_xx->GetBinContent(bin);
    dataN.avg_Q_yy = pN_yy->GetBinContent(bin);
    dataN.avg_Q_xy = pN_xy->GetBinContent(bin);
  }
}

void QVecCDB::write_cdb()
{
  std::string output_dir = std::format("{}/{}", m_output_dir, m_runnumber);

  if (std::filesystem::create_directories(output_dir))
  {
    std::cout << std::format("Success: Directory {} created.\n", output_dir);
  }
  else
  {
    std::cout << std::format("Info: Directory {} already exists.\n", output_dir);
  }

  write_cdb_EventPlane(output_dir);
  write_cdb_BadTowers(output_dir);
}

void QVecCDB::write_cdb_BadTowers(const std::string &output_dir)
{
  std::string payload = "SEPD_HotMap";
  std::string fieldname_status = "status";
  std::string fieldname_sigma = "SEPD_sigma";
  std::string output_file = std::format("{}/{}-{}-{}.root", output_dir, payload, m_cdb_tag, m_runnumber);

  std::unique_ptr<CDBTTree> cdbttree = std::make_unique<CDBTTree>(output_file);

  auto h_sEPD_Bad_Channels = load_and_clone<TProfile>("h_sEPD_Bad_Channels");

  for (int channel = 0; channel < h_sEPD_Bad_Channels->GetNbinsX(); ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(channel);
    int status = h_sEPD_Bad_Channels->GetBinContent(channel+1);

    float sigma = 0;

    // Hot
    if (status == 2)
    {
      sigma = SIGMA_HOT;
    }

    // Cold
    else if (status == 3)
    {
      sigma = SIGMA_COLD;
    }

    cdbttree->SetIntValue(key, fieldname_status, status);
    cdbttree->SetFloatValue(key, fieldname_sigma, sigma);
  }

  std::cout << std::format("Saving CDB: {} to {}\n", payload, output_file);

  cdbttree->Commit();
  cdbttree->WriteCDBTTree();
}

void QVecCDB::write_cdb_EventPlane(const std::string &output_dir)
{
  std::string payload = "SEPD_EventPlaneCalib";
  std::string output_file = std::format("{}/{}-{}-{}.root", output_dir, payload, m_cdb_tag, m_runnumber);

  std::unique_ptr<CDBTTree> cdbttree = std::make_unique<CDBTTree>(output_file);

  for (size_t h_idx = 0; h_idx < m_harmonics.size(); ++h_idx)
  {
    int n = m_harmonics[h_idx];

    // Define lambdas to generate field names consistently
    auto field = [&](const char* det, const char* var) {
      return std::format("Q_{}_{}_{}_avg", det, var, n);
    };

    for (size_t cent_bin = 0; cent_bin < m_cent_bins; ++cent_bin)
    {
      int key = static_cast<int>(cent_bin);

      // Access data references to clean up the calls
      const auto& S = getData(h_idx, cent_bin, QVecShared::Subdetector::S);
      const auto& N = getData(h_idx, cent_bin, QVecShared::Subdetector::N);

      // South
      cdbttree->SetDoubleValue(key, field("S", "x"), S.avg_Q.x);
      cdbttree->SetDoubleValue(key, field("S", "y"), S.avg_Q.y);
      cdbttree->SetDoubleValue(key, field("S", "xx"), S.avg_Q_xx);
      cdbttree->SetDoubleValue(key, field("S", "yy"), S.avg_Q_yy);
      cdbttree->SetDoubleValue(key, field("S", "xy"), S.avg_Q_xy);

      // North
      cdbttree->SetDoubleValue(key, field("N", "x"), N.avg_Q.x);
      cdbttree->SetDoubleValue(key, field("N", "y"), N.avg_Q.y);
      cdbttree->SetDoubleValue(key, field("N", "xx"), N.avg_Q_xx);
      cdbttree->SetDoubleValue(key, field("N", "yy"), N.avg_Q_yy);
      cdbttree->SetDoubleValue(key, field("N", "xy"), N.avg_Q_xy);
    }
  }

  std::cout << std::format("Saving CDB: {} to {}\n", payload, output_file);

  cdbttree->Commit();
  cdbttree->WriteCDBTTree();
}
