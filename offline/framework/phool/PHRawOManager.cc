//-----------------------------------------------------------------------------
//
//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PHRawOManager
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PHRawOManager.h"
#include "PHCompositeNode.h"
#include "PHRawDataNode.h"

#include <Event/ogzBuffer.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

PHRawOManager::PHRawOManager(const string &newFile, const int run, const int bufl, const int evtl, const int complvl):
  filedesc(-1),
  runNumber(0),
  bufferSize(0),
  eventLength(0)
{
  if (!setFile(newFile, run, bufl, evtl, complvl))
    {
      filename         = "file open failed";
      filedesc         = -1;
      memBuffer        = NULL;
      fileBuffer       = NULL;
      compressionLevel = 0;
    }
}

PHRawOManager::PHRawOManager():
  filedesc(-1),
  memBuffer(NULL),
  fileBuffer(NULL),
  runNumber(0),
  bufferSize(0),
  eventLength(0),
  compressionLevel(0)
{}

PHRawOManager::~PHRawOManager()
{
  closeFile();
}

void
PHRawOManager::closeFile()
{
  if (fileBuffer)
    {
      delete fileBuffer;
    }
  if (memBuffer)
    {
      delete memBuffer;
    }
  fileBuffer = 0;
  memBuffer = 0;
  if (filedesc >= 0)
    {
      close(filedesc);
      filedesc = -1;
    }
}

PHBoolean
PHRawOManager::setFile(const string &setFile, const int setRun, const int setBufl, const int setEvtl, const int complvl)
{
  filename         = setFile;
  runNumber        = setRun;
  bufferSize       = setBufl;
  compressionLevel = complvl;

  if (setEvtl == -1)
    {
      eventLength = bufferSize / 4;
    }
  else
    {
      eventLength = setEvtl;
    }

  if (filedesc >= 0)
    {
      closeFile();  // close the file if it is open (originally unprotected close)
    }

  // open file
  filedesc =  open(filename.c_str(), 
		   O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE ,
		   S_IRWXU | S_IROTH | S_IRGRP );

  if (filedesc < 0)
    {

      PHMessage("PHRawOManager::setFile", PHError, "could not open file");
      return False;
    }
  memBuffer  = new PHDWORD[bufferSize];
  if (compressionLevel == 0)
    {
      int status = 0;
      fileBuffer = new oBuffer(filedesc, memBuffer, bufferSize, status, runNumber);
      if (status > 0)
        {
          PHMessage("PHRawOManager::setFile", PHError, "could not open file");
          return False;
        }
    }
  else if (1 <= compressionLevel && compressionLevel <= 9)
    {
      int status = 0;
      fileBuffer = new ogzBuffer(filedesc, memBuffer, bufferSize,  status, runNumber, compressionLevel);
      if (status > 0)
        {
          PHMessage("PHRawOManager::setFile", PHError, "could not open file");
          return False;
        }
    }
  else
    {
      PHMessage("PHRawOManager::setFile", PHWarning, "illegal compression level, no compression done");
    }

  return True;
}

PHBoolean
PHRawOManager::write(PHCompositeNode* topNode)
{
  //
  // The write function of the PHCompositeNode topNode will
  // recursively call the write functions of its subnodes.
  // The PHRawDataNodes among them call the write function found
  // below.
  //
  if (filedesc >= 0 && fileBuffer)
    {
      fileBuffer->nextEvent(eventLength, DATAEVENT);
      topNode->write(this);
      eventNumber++;
      return True;
    }
  return False;
}

PHBoolean
PHRawOManager::write(PHRawDataNode* node)
{
  if (filedesc >= 0 && fileBuffer)
    {
      int bytesAddedToBuffer = fileBuffer->addUnstructPacketData(node->getData(), node->getLength(), node->getID(),
								 node->getWordLength(), node->getHitFormat());
      if (bytesAddedToBuffer <= 0)
        {
          PHMessage("PHRawOManager::write", PHError, "Zero bytes added to buffer");
          return False;
        }
      return True;
    }
  return False;
}

void
PHRawOManager::print() const
{
  cout << "File attached : " << filename    << endl;
  cout << "Buffer size   : " << bufferSize  << endl;
  cout << "Event Length  : " << eventLength << endl;
  cout << "Run number    : " << runNumber   << endl;
  cout << "Events written: " << eventNumber << endl;
}
