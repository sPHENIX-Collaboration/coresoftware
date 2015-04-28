//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/framework/phool/PHIOManager.C,v 1.3 2014/01/12 04:29:18 pinkenbu Exp $
//
//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PHIOManager
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------

#include "PHIOManager.h"

PHIOManager::PHIOManager():
  filename("no file attached"),
  eventNumber(0)
{}

PHString
PHIOManager::getFilename() const
{
  if (filename.getString())
    {
      return filename;
    }
  else
    {
      PHString dummy("empty filename");
      return dummy;
    }
}
