// $Id: $                                                                                             
 
/*!
 * \file QAHistManagerDef.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "QAHistManagerDef.h"
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <cassert>

using namespace std;

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
//        cout
//            << "QAHistManagerDef::get_HistoManager - Making Fun4AllHistoManager EMCalAna_HISTOS"
//            << endl;
        hm = new Fun4AllHistoManager(HistoManagerName);
        se->registerHistoManager(hm);
      }

    assert(hm);

    return hm;
  }

  //! Save hist to root files
  void
  saveQARootFile(const std::string file_name)
  {
    getHistoManager()->dumpHistos(file_name);
  }
}


