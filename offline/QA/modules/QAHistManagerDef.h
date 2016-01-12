// $Id: $                                                                                             

/*!
 * \file QAHistManagerDef.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef QAHISTMANAGERDEF_H_
#define QAHISTMANAGERDEF_H_

#include <string>
class Fun4AllHistoManager;

namespace QAHistManagerDef
{

  //! Get a pointer to the default hist manager for QA modules
  Fun4AllHistoManager *
  getHistoManager();

  //! Save hist to root files. It will overwrite the old file if exist
  void
  saveQARootFile(const std::string file_name);

  //! default name for QA histogram manager
  static const std::string HistoManagerName = "QA_HISTOS";
}

#endif /* QAHISTMANAGERDEF_H_ */
