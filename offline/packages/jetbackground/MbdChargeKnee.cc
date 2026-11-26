#include "MbdChargeKnee.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <TAxis.h>
#include <TH1.h>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

float MbdChargeKnee::find(const float user, const std::string &tag, const std::string &path, const TH1 *table, const int run,
                          std::string &source, bool &dataset_default)
{
  dataset_default = false;
  if (user > 0)
  {
    source = "set_mbd_knee";
    return user;
  }

  // per-run knee payload
  std::string missing;
  if (!path.empty() || !tag.empty())
  {
    std::string url = path;
    if (path.empty())
    {
      url = CDBInterface::instance()->getUrl(tag);
      source = "CDB tag " + tag + " (" + url + ")";
    }
    else
    {
      if (std::filesystem::is_directory(path))
      {
        url = path + "/mbdq_knee_" + std::to_string(run) + ".root";
      }
      source = url;
    }
    // CDBTTree exits on a file it cannot open
    if (!url.empty() && std::filesystem::exists(url))
    {
      CDBTTree cdb(url);
      cdb.LoadCalibrations();
      const float knee = cdb.GetSingleFloatValue("mbdq_knee");
      if (knee > 0)
      {
        return knee;
      }
    }
    missing = "no knee in " + source + "; ";
  }

  if (!table)
  {
    source = missing + "no per-run table in the calibration payload";
    return -1;
  }

  // the calibration payload's per-run table
  source = missing + "per-run table of the calibration payload";
  if (run >= table->GetXaxis()->GetXmin() && run < table->GetXaxis()->GetXmax())
  {
    const auto knee = static_cast<float>(table->GetBinContent(table->GetXaxis()->FindFixBin(run)));
    if (knee > 0)
    {
      return knee;
    }
  }

  // dataset default: the median knee of the table
  std::vector<double> knees;
  for (int i = 1; i <= table->GetNbinsX(); ++i)
  {
    if (table->GetBinContent(i) > 0)
    {
      knees.push_back(table->GetBinContent(i));
    }
  }
  if (knees.empty())
  {
    source = missing + "empty per-run table in the calibration payload";
    return -1;
  }
  std::nth_element(knees.begin(), knees.begin() + (knees.size() / 2), knees.end());
  source = missing + "dataset default: median knee of the " + std::to_string(knees.size()) + " runs of the calibration payload";
  dataset_default = true;
  return static_cast<float>(knees[knees.size() / 2]);
}
