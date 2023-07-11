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

#include <TAxis.h>

#include <cassert>
#include <cmath>
#include <iosfwd>  // for std
#include <vector>

namespace QAHistManagerDef
{
  //! Get a pointer to the default hist manager for QA modules
  Fun4AllHistoManager *
  getHistoManager()
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    Fun4AllHistoManager *hm = se->getHistoManager(HistoManagerName);

    if (not hm)
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

  //! Save hist to root files
  void saveQARootFile(const std::string &file_name)
  {
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
