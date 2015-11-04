//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/offline/database/pdbcal/base/PdbCalHeader.cc,v 1.4 2009/08/19 17:07:22 pinkenbu Exp $
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PdbCalHeader
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PdbCalHeader.h"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

PdbCalHeader::PdbCalHeader()
{
  setDescription("Unknown");
  setUserName("Unknown");
}

PdbCalHeader::PdbCalHeader(const PdbCalHeader& otherheader) :
  bankID(42), insertTime(otherheader.insertTime),
  startValTime(otherheader.startValTime),
  endValTime(otherheader.endValTime)
{
  setDescription(otherheader.description);
  setUserName(otherheader.userName);
}


PdbCalHeader::PdbCalHeader(PdbBankID id, const char *descr, PHTimeStamp &tStart, PHTimeStamp &tEnd)
{
  bankID       = id;
  startValTime = tStart;
  endValTime   = tEnd;

  insertTime.setToSystemTime();
  setDescription(descr);
  setUserName("Unknown");
}

void 
PdbCalHeader::print() const
{
  bankID.print();
  cout << "Insert   : " << insertTime   << endl;
  cout << "StartVal : " << startValTime << endl;
  cout << "EndVal   : " << endValTime   << endl;
  cout << "Calibration Description : " << description << endl;
  cout << "User Name = " << userName << endl;
}
