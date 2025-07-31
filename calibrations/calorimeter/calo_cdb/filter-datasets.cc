
#include "filter-datasets.h"

// -- My Utils --
#include "myUtils.h"

// c++ includes --
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

FilterDatasets::FilterDatasets(Bool_t debug)
  : m_debug(debug)
{
}

void FilterDatasets::readRunInfo(const std::string &line)
{
  std::vector<std::string> tokens = myUtils::split(line, ',');

  if (tokens.size() < 2)
  {
    std::cout << "ERROR: line does not contain both RUN and DATASET information! Skipping: " << line << std::endl;
    return;
  }

  m_runInfo.emplace_back(tokens[0], tokens[1]);
}

std::string FilterDatasets::getCalibration(const std::string &pl_type, uint64_t iov)
{
  if (!uti)
  {
    uti = std::make_unique<CDBUtils>();
  }
  return uti->getUrl(pl_type, iov);
}

int FilterDatasets::setGlobalTag(const std::string &tagname)
{
  if (!uti)
  {
    uti = std::make_unique<CDBUtils>();
  }
  int iret = uti->setGlobalTag(tagname);
  return iret;
}

void FilterDatasets::analyze(const std::string &input, const std::string &outputDir)
{
  std::string filename = outputDir + "/" + std::filesystem::path(input).filename().stem().string() + "-process.csv";

  std::ofstream file(filename);

  if (!file.is_open())
  {
    std::cout << "ERROR: Could not open file " << filename << " for writing." << std::endl;
    return;  // Exit if file cannot be opened
  }

  // write the header for the CSV
  file << "runnumber" << std::endl;

  // loop over each run to check if any is missing the latest calibration
  // if a run is missing the latest calibration then write it to the CSV for processing
  for (const auto &run_dataset : m_runInfo)
  {
    std::string run = run_dataset.first;
    std::string dataset = run_dataset.second;
    ++m_ctr["ctr_run"];
    std::cout << "Run: " << run << ", Dataset: " << dataset << ", Processing: "
              << m_ctr["ctr_run"] << ", " << m_ctr["ctr_run"] * 100. / static_cast<Double_t>(m_runInfo.size()) << " %" << std::endl;

    Bool_t keep = false;

    for (const auto &cdbName : m_cdbName)
    {
      if (m_debug)
      {
        std::cout << "Attempt to get Calibration" << std::endl;
      }
      std::string cdb = getCalibration(cdbName, std::stoul(run));
      if (m_debug)
      {
        std::cout << "Get Calibration Calls: " << ++m_ctr["getCalib_calls"] << std::endl;
      }

      std::stringstream ss;
      if (cdbName == "CEMC_BadTowerMap")
      {
        ss << "EMCalHotMap_" << dataset << "_" << run << "cdb.root";
      }
      else
      {
        ss << cdbName << "_" << dataset << "_" << run << ".root";
      }

      std::string suffix = ss.str();

      if (!cdb.ends_with(suffix))
      {
        if (cdb.starts_with("DataBaseException"))
        {
          ++m_ctr["cdb_missing_" + cdbName];
        }
        else
        {
          ++m_ctr["cdb_outdated_" + cdbName];
        }
        keep = true;
      }
    }

    if (keep)
    {
      file << run << std::endl;
    }
  }

  std::cout << "===============================" << std::endl;
  std::cout << "Stats" << std::endl;
  std::cout << "Total Runs: " << m_runInfo.size() << std::endl;
  for (const auto &[name, value] : m_ctr)
  {
    if (name.starts_with("cdb"))
    {
      std::cout << "CDB: " << name << ", Counts: " << value << ", " << value * 100. / static_cast<Double_t>(m_runInfo.size()) << " %" << std::endl;
    }
  }
  std::cout << "===============================" << std::endl;

  file.close();
}

void FilterDatasets::process(const std::string &input, const std::string &output)
{
  std::cout << "#############################" << std::endl;
  std::cout << "Run Parameters" << std::endl;
  std::cout << "input: " << input << std::endl;
  std::cout << "output: " << output << std::endl;
  std::cout << "Debug: " << ((m_debug) ? "True" : "False") << std::endl;
  std::cout << "#############################" << std::endl;

  setGlobalTag("newcdbtag");

  std::filesystem::path input_filepath_obj(input);
  if (!myUtils::readCSV(input_filepath_obj, [this](const std::string &line)
                        { this->readRunInfo(line); }))
  {
    return;
  }

  analyze(input, output);
}
