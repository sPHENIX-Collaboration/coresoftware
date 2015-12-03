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
#include <PdbCalHeader.h>
#include <PHString.h>

#include <cstdlib>
#include <iostream>

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

void
PdbCalHeader::setDescription(const char *val)
{
  // this construct protects against not zero terminated strings
  // it copies only a fixed number of chars into the theName array
  // (sizeof(theName)-1 leaves space for a zero terminated char
  // at the end. strlen has problems when used with non zero
  // terminated strings and might crash
  // this code might still crash in the cout of name, but that is
  // already in the exit() part
  // Flawfinder: ignore signals flawfinder to ignore this line, strncpy does not zero terminate the string
  // if the length is exceeded, the termination is done in the line afterwards where the last
  // character is set to \0
  strncpy(description, val, (sizeof(description) - 1)); /* Flawfinder: ignore */
  description[(sizeof(description)-1)] = '\0';
  if (strncmp(description, val, sizeof(description)))
    {
      cout << "Description exceeds maximum length of "
	   << (sizeof(description) - 1)
	   << " characters or is not zero terminated" << endl;
      cout << "Max lenght description: " << description << endl;
      cout << "There is no point in continuing, fix your code and try again" << endl;
      cout << "Name used (code might crash now when printing out not zero terminated string): " << val << endl;
      exit(1);
    }
  return;
}

void
PdbCalHeader::setUserName(const char *val)
{
  // this construct protects against not zero terminated strings
  // it copies only a fixed number of chars into the theName array
  // (sizeof(theName)-1 leaves space for a zero terminated char
  // at the end. strlen has problems when used with non zero
  // terminated strings and might crash
  // this code might still crash in the cout of name, but that is
  // already in the exit() part
  // Flawfinder: ignore signals flawfinder to ignore this line, strncpy does not zero terminate the string
  // if the length is exceeded, the termination is done in the line afterwards where the last
  // character is set to \0
  strncpy(userName, val, (sizeof(userName) - 1)); /* Flawfinder: ignore */
  userName[(sizeof(userName)-1)] = '\0';
  if (strncmp(userName, val, sizeof(userName)))
    {
      cout << "UserName exceeds maximum length of "
	   << (sizeof(userName) - 1)
	   << " characters or is not zero terminated" << endl;
      cout << "Max lenght userName: " << userName << endl;
      cout << "There is no point in continuing, fix your code and try again" << endl;
      cout << "Name used (code might crash now when printing out not zero terminated string): " << val << endl;
      exit(1);
    }
  return;
}
