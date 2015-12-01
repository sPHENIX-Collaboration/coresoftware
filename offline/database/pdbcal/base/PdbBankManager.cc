//-----------------------------------------------------------------------------
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PdbBankManager
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PdbBankManager.h"

PdbBankManager *PdbBankManager::__instance = NULL;


PdbBankManager::PdbBankManager()
{
}

PdbBankManager::~PdbBankManager()
{
  __instance = NULL;
}

PdbBankManager *PdbBankManager::instance()
{
  if ( ! __instance )
    {
      std::cout << __FILE__ << "  " << __LINE__ << 
	" No instance of PdbBankManager available" << std::endl;
    }
  
  return __instance;
}

