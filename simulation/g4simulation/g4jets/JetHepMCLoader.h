// $Id: $

/*!
 * \file JetHepMCLoader.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef JETHEPMCLOADER_H_
#define JETHEPMCLOADER_H_

#include <fun4all/SubsysReco.h>

/*!
 * \brief JetHepMCLoader loads special jet objects encoded in HepMC records to DST Jet nodes. Example use are loading sHijing HIJFRG jets
 */
class JetHepMCLoader : public SubsysReco
{
 public:
  JetHepMCLoader(const std::string &name = "JetHepMCLoader");
  virtual ~JetHepMCLoader();
};

#endif /* JETHEPMCLOADER_H_ */
