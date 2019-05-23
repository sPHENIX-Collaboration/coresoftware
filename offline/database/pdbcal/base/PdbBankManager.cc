//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/database/pdbcal/base/PdbBankManager.cc,v 1.2 2004/08/13 17:46:00 purschke Exp $
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PdbBankManager
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PdbBankManager.h"

#include <iostream>

PdbBankManager *PdbBankManager::__instance = nullptr;

PdbBankManager::~PdbBankManager()
{
  __instance = nullptr;
}

PdbBankManager *PdbBankManager::instance()
{
  if (!__instance)
  {
    std::cout << __FILE__ << "  " << __LINE__ << " No instance of PdbBankManager available" << std::endl;
  }

  return __instance;
}
