// $Id: $

/*!
 * \file QAHistManagerDef.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "QAHistManagerDef.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

#include <TAxis.h>
#include <TH1.h>

#include <cassert>
#include <cmath>
#include <iosfwd>  // for std

namespace QAHistManagerDef
{
  //! Get a pointer to the default hist manager for QA modules
  Fun4AllHistoManager *
  getHistoManager()
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    Fun4AllHistoManager *hm = se->getHistoManager(HistoManagerName);

    if (! hm)
    {
      //        std::cout
      //            << "QAHistManagerDef::get_HistoManager - Making Fun4AllHistoManager EMCalAna_HISTOS"
      //            << std::endl;
      hm = new Fun4AllHistoManager(HistoManagerName);
      se->registerHistoManager(hm);
    }

    assert(hm);

    return hm;
  }
  std::vector<std::string> tokenize(const std::string& str, const char* delimiter)
  {
    std::vector<std::string> tokens;
  size_t start = 0;
  size_t end = str.find(delimiter);

  while (end != std::string::npos)
  {
    tokens.push_back(str.substr(start, end - start));
    start = end + 1;
    end = str.find(delimiter, start);
  }
  tokens.push_back(str.substr(start));

  return tokens;
  }
  //! Save hist to root files
  void saveQARootFile(const std::string &file_name)
  {
    // add provenance info
    std::string build = "";
    const std::string offlinemain = getenv("OFFLINE_MAIN");
    auto tokens = tokenize(offlinemain,"/");
    for(const auto& token : tokens)
    {
      if(token.find("new") != std::string::npos)
    {
      build = "new";
    }
      else if (token.find("ana") != std::string::npos)
      {
        build = tokens.back();
      }
    }
    auto rc = recoConsts::instance();
    std::string dbtag = rc->get_StringFlag("CDB_GLOBALTAG");
    std::string info = "Build: " + build + " , dbtag: " + dbtag;
    TH1* h = new TH1I("h_QAHistManagerDef_ProductionInfo","",10,0,10);
    h->SetTitle(info.c_str());
    getHistoManager()->registerHisto(h);
    

    // dump histos to file
    getHistoManager()->dumpHistos(file_name);
  }

  //! utility function to
  void useLogBins(TAxis *axis)
  {
    assert(axis);
    assert(axis->GetXmin() > 0);
    assert(axis->GetXmax() > 0);

    const int bins = axis->GetNbins();

    Axis_t from = log10(axis->GetXmin());
    Axis_t to = log10(axis->GetXmax());
    Axis_t width = (to - from) / bins;
    std::vector<Axis_t> new_bins(bins + 1);

    for (int i = 0; i <= bins; i++)
    {
      new_bins[i] = pow(10, from + i * width);
    }
    axis->Set(bins, new_bins.data());
  }
}  // namespace QAHistManagerDef
