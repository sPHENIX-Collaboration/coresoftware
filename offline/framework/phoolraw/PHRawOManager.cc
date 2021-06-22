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

#include "PHRawDataNode.h"

#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

#include <Event/EventTypes.h>
#include <Event/oBuffer.h>
#include <Event/ospBuffer.h>
#include <Event/phenixTypes.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>  // for close
#include <iostream>

using namespace std;

PHRawOManager::PHRawOManager(const string& newFile, const int run, const int bufl, const int evtl, const int complvl)
  : filedesc(-1)
  , runNumber(0)
  , bufferSize(0)
  , eventLength(0)
{
  if (!setFile(newFile, run, bufl, evtl, complvl))
  {
    filename = "file open failed";
    filedesc = -1;
    memBuffer = nullptr;
    fileBuffer = nullptr;
    compressionLevel = 0;
  }
}

PHRawOManager::PHRawOManager()
  : filedesc(-1)
  , memBuffer(nullptr)
  , fileBuffer(nullptr)
  , runNumber(0)
  , bufferSize(0)
  , eventLength(0)
  , compressionLevel(0)
{
}

PHRawOManager::~PHRawOManager()
{
  closeFile();
}

void PHRawOManager::closeFile()
{
  if (fileBuffer)
  {
    delete fileBuffer;
  }
  if (memBuffer)
  {
    delete memBuffer;
  }
  fileBuffer = nullptr;
  memBuffer = nullptr;
  if (filedesc >= 0)
  {
    close(filedesc);
    filedesc = -1;
  }
}

bool PHRawOManager::setFile(const string& setFile, const int setRun, const int setBufl, const int setEvtl, const int complvl)
{
  filename = setFile;
  runNumber = setRun;
  bufferSize = setBufl;
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
  filedesc = open(filename.c_str(),
                  O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                  S_IRWXU | S_IROTH | S_IRGRP);

  if (filedesc < 0)
  {
    PHMessage("PHRawOManager::setFile", PHError, "could not open file");
    return false;
  }
  memBuffer = new PHDWORD[bufferSize];
  fileBuffer = new ospBuffer(filedesc, memBuffer, bufferSize, runNumber);

  return true;
}

bool PHRawOManager::write(PHCompositeNode* topNode)
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
    return true;
  }
  return false;
}

bool PHRawOManager::write(PHRawDataNode* node)
{
  if (filedesc >= 0 && fileBuffer)
  {
    int bytesAddedToBuffer = fileBuffer->addUnstructPacketData(node->getData(), node->getLength(), node->getID(),
                                                               node->getWordLength(), node->getHitFormat());
    if (bytesAddedToBuffer <= 0)
    {
      PHMessage("PHRawOManager::write", PHError, "Zero bytes added to buffer");
      return false;
    }
    return true;
  }
  return false;
}

void PHRawOManager::print() const
{
  cout << "File attached : " << filename << endl;
  cout << "Buffer size   : " << bufferSize << endl;
  cout << "Event Length  : " << eventLength << endl;
  cout << "Run number    : " << runNumber << endl;
  cout << "Events written: " << eventNumber << endl;
}
