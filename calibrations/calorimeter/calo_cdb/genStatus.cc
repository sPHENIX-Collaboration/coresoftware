#include "genStatus.h"

#include "geometry_constants.h"

// -- sPHENIX includes --
#include <calobase/TowerInfoDefs.h>

#include <cdbobjects/CDBTTree.h>

#include <emcnoisytowerfinder/emcNoisyTowerFinder.h>

// -- root includes --
#include <TFile.h>

// c++ includes --
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

void GenStatus::setRunDataset(const std::string &input)
{
  std::string basename = std::filesystem::path(input).filename().stem().string();
  m_run = basename.substr(0, basename.find('_'));
  m_dataset = basename.substr(basename.find('_') + 1, basename.size() - basename.find('_'));
}

int GenStatus::readHists(const std::string &input)
{
  // Create an input stream
  std::ifstream file(input);

  // Check if the file was successfully opened
  if (!file.is_open())
  {
    std::cout << "Failed to open file list: " << input << std::endl;
    return 1;
  }

  std::cout << "Reading Hists" << std::endl;
  std::cout << "======================================" << std::endl;

  h_CaloValid_cemc_etaphi_badChi2 = std::make_unique<TProfile2D>("cemc_etaphi_badChi2", "", cemc_bins_eta, 0, cemc_bins_eta, cemc_bins_phi, 0, cemc_bins_phi);
  h_CaloValid_ihcal_etaphi_badChi2 = std::make_unique<TProfile2D>("ihcal_etaphi_badChi2", "", hcal_bins_eta, 0, hcal_bins_eta, hcal_bins_phi, 0, hcal_bins_phi);
  h_CaloValid_ohcal_etaphi_badChi2 = std::make_unique<TProfile2D>("ohcal_etaphi_badChi2", "", hcal_bins_eta, 0, hcal_bins_eta, hcal_bins_phi, 0, hcal_bins_phi);

  h_CaloValid_cemc_etaphi_time_raw = std::make_unique<TProfile2D>("cemc_etaphi_time_raw", "", cemc_bins_eta, 0, cemc_bins_eta, cemc_bins_phi, 0, cemc_bins_phi);
  h_CaloValid_ihcal_etaphi_time_raw = std::make_unique<TProfile2D>("ihcal_etaphi_time_raw", "", hcal_bins_eta, 0, hcal_bins_eta, hcal_bins_phi, 0, hcal_bins_phi);
  h_CaloValid_ohcal_etaphi_time_raw = std::make_unique<TProfile2D>("ohcal_etaphi_time_raw", "", hcal_bins_eta, 0, hcal_bins_eta, hcal_bins_phi, 0, hcal_bins_phi);

  std::string line;
  std::map<std::string, int> ctr;
  while (std::getline(file, line))
  {
    ++ctr["total_files"];

    std::cout << "Reading File: " << line << std::endl;
    std::unique_ptr<TFile> tf(TFile::Open(line.c_str()));
    if (!tf || tf->IsZombie())
    {
      std::cout << "Error: Could not open ROOT file: " << line << std::endl;
      continue;  // Indicate an error
    }

    ++ctr["successfully_opened_files"];

    auto *h = dynamic_cast<TProfile2D *>(tf->Get("h_CaloValid_cemc_etaphi_badChi2"));

    if (h)
    {
      h_CaloValid_cemc_etaphi_badChi2->Add(h);
      ++ctr["h_CaloValid_cemc_etaphi_badChi2"];
    }

    h = dynamic_cast<TProfile2D *>(tf->Get("h_CaloValid_ihcal_etaphi_badChi2"));

    if (h)
    {
      h_CaloValid_ihcal_etaphi_badChi2->Add(h);
      ++ctr["h_CaloValid_ihcal_etaphi_badChi2"];
    }

    h = dynamic_cast<TProfile2D *>(tf->Get("h_CaloValid_ohcal_etaphi_badChi2"));

    if (h)
    {
      h_CaloValid_ohcal_etaphi_badChi2->Add(h);
      ++ctr["h_CaloValid_ohcal_etaphi_badChi2"];
    }

    h = dynamic_cast<TProfile2D *>(tf->Get("h_CaloValid_cemc_etaphi_time_raw"));

    if (h)
    {
      h_CaloValid_cemc_etaphi_time_raw->Add(h);
      ++ctr["h_CaloValid_cemc_etaphi_time_raw"];
    }

    h = dynamic_cast<TProfile2D *>(tf->Get("h_CaloValid_ihcal_etaphi_time_raw"));

    if (h)
    {
      h_CaloValid_ihcal_etaphi_time_raw->Add(h);
      ++ctr["h_CaloValid_ihcal_etaphi_time_raw"];
    }

    h = dynamic_cast<TProfile2D *>(tf->Get("h_CaloValid_ohcal_etaphi_time_raw"));

    if (h)
    {
      h_CaloValid_ohcal_etaphi_time_raw->Add(h);
      ++ctr["h_CaloValid_ohcal_etaphi_time_raw"];
    }

    tf->Close();
  }

  std::cout << "===============================" << std::endl;
  std::cout << "Stats" << std::endl;
  std::cout << "Successfully opened files: " << ctr["successfully_opened_files"] << ", " << ctr["successfully_opened_files"] * 100. / ctr["total_files"] << " %" << std::endl;
  for (const auto &[name, value] : ctr)
  {
    if (name.starts_with("h_CaloValid"))
    {
      std::cout << "Hist: " << name << ", Found: " << value << ", " << value * 100. / ctr["successfully_opened_files"] << " %" << std::endl;
    }
  }
  std::cout << "===============================" << std::endl;

  // Close the file
  file.close();

  return 0;
}

void GenStatus::histToCaloCDBTree(const std::string &outputfile, const std::string &fieldName, int icalo, TProfile2D *hist)
{
  unsigned int neta;
  unsigned int nphi;

  if (icalo != 0 && icalo != 1)
  {
    return;
  }

  if (icalo == 0)
  {
    neta = CaloGeometry::CEMC_ETA_BINS;
    nphi = CaloGeometry::CEMC_PHI_BINS;
  }
  if (icalo == 1)
  {
    neta = CaloGeometry::HCAL_ETA_BINS;
    nphi = CaloGeometry::HCAL_PHI_BINS;
  }

  std::unique_ptr<CDBTTree> cdbttree = std::make_unique<CDBTTree>(outputfile);

  double mean = 0;
  int count = 0;

  for (unsigned int ie = 0; ie < neta; ie++)
  {
    for (unsigned int ip = 0; ip < nphi; ip++)
    {
      unsigned int key;
      if (icalo == 0)
      {
        key = TowerInfoDefs::encode_emcal(ie, ip);
      }
      if (icalo == 1)
      {
        key = TowerInfoDefs::encode_hcal(ie, ip);
      }
      float val = hist->GetBinContent(ie + 1, ip + 1);
      cdbttree->SetFloatValue(key, fieldName, val);
      mean += val;
      count++;
    }
  }

  std::cout << "Writing " << outputfile << "   with mean=" << mean / count << std::endl;
  cdbttree->Commit();
  cdbttree->WriteCDBTTree();
}

void GenStatus::analyze(const std::string &outputDir)
{
  std::string detector = "CEMC";
  // fracBadChi2
  std::string payloadName = outputDir + "/" + detector + "_hotTowers_fracBadChi2" + "_" + m_dataset + "_" + m_run + ".root";
  if (h_CaloValid_cemc_etaphi_badChi2)
  {
    histToCaloCDBTree(payloadName, "fraction", 0, h_CaloValid_cemc_etaphi_badChi2.get());
  }
  // time
  payloadName = outputDir + "/" + detector + "_meanTime" + "_" + m_dataset + "_" + m_run + ".root";
  if (h_CaloValid_cemc_etaphi_time_raw)
  {
    histToCaloCDBTree(payloadName, "time", 0, h_CaloValid_cemc_etaphi_time_raw.get());
  }

  detector = "HCALIN";
  // fracBadChi2
  payloadName = outputDir + "/" + detector + "_hotTowers_fracBadChi2" + "_" + m_dataset + "_" + m_run + ".root";
  if (h_CaloValid_ihcal_etaphi_badChi2)
  {
    histToCaloCDBTree(payloadName, "fraction", 1, h_CaloValid_ihcal_etaphi_badChi2.get());
  }
  // time
  payloadName = outputDir + "/" + detector + "_meanTime" + "_" + m_dataset + "_" + m_run + ".root";
  if (h_CaloValid_ihcal_etaphi_time_raw)
  {
    histToCaloCDBTree(payloadName, "time", 1, h_CaloValid_ihcal_etaphi_time_raw.get());
  }

  detector = "HCALOUT";
  // fracBadChi2
  payloadName = outputDir + "/" + detector + "_hotTowers_fracBadChi2" + "_" + m_dataset + "_" + m_run + ".root";
  if (h_CaloValid_ohcal_etaphi_badChi2)
  {
    histToCaloCDBTree(payloadName, "fraction", 1, h_CaloValid_ohcal_etaphi_badChi2.get());
  }
  // time
  payloadName = outputDir + "/" + detector + "_meanTime" + "_" + m_dataset + "_" + m_run + ".root";
  if (h_CaloValid_ohcal_etaphi_time_raw)
  {
    histToCaloCDBTree(payloadName, "time", 1, h_CaloValid_ohcal_etaphi_time_raw.get());
  }
}

void GenStatus::process(const std::string &input, const std::string &output)
{
  std::cout << "#############################" << std::endl;
  std::cout << "Run Parameters" << std::endl;
  std::cout << "input: " << input << std::endl;
  std::cout << "output: " << output << std::endl;
  std::cout << "#############################" << std::endl;

  setRunDataset(input);

  std::cout << "Processing: Run: " << m_run << ", Dataset: " << m_dataset << std::endl;

  std::stringstream datasetDir;
  datasetDir.str("");
  datasetDir << output << "/" << m_run << "_" << m_dataset;

  std::string outputDir = datasetDir.str();

  std::string hotMapFile = "EMCalHotMap_" + m_dataset + "_" + m_run + ".root";
  std::string hotMapOutput = datasetDir.str() + "/" + hotMapFile;

  datasetDir << "/QA";

  // create output & QA directory
  std::filesystem::create_directories(datasetDir.str());

  std::string s_input = input;
  std::string hotMapOutputQA = datasetDir.str() + "/" + hotMapFile;

  // merges individal qa into one per run
  readHists(input);
  analyze(outputDir);

  std::unique_ptr<emcNoisyTowerFinder> calo = std::make_unique<emcNoisyTowerFinder>();
  calo->FindHot(s_input, hotMapOutput, "h_CaloValid_cemc_etaphi");

  if (std::filesystem::exists(hotMapOutput))
  {
    std::filesystem::rename(hotMapOutput, hotMapOutputQA);
  }
  else
  {
    std::cout << "ERROR: EMCal Hot Map FAILED to Create." << std::endl;
  }
}
