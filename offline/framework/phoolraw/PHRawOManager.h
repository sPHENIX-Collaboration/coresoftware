// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHOOLRAW_PHRAWOMANAGER_H
#define PHOOLRAW_PHRAWOMANAGER_H

//-----------------------------------------------------------------------------
//
//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Declaration of class PHRawOManager
//
//  Purpose: Write a node-tree into a PRDF file
//
//  Description:
//     - The default constructor sets all pointers to NULL. SetFile() would
//       have to be called in addition.
//     - The second constructor calls the setFile method to which it passes
//       on its arguments. In order of appearance:
//        1. The filename, type std::string, no default
//        2. The runnumber, int, defaults to 0
//        3. The buffersize, int, defaults to 100000
//        4. The eventlength, int, defaults to -1, in which case a quarter
//           of the selected buffersize will be chosen.
//        5. The compression-level. Values can range between 1 and 9, as in
//           gzip, the underlying compression program. Zero means no
//           compression at all, default is 3.
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------

#include <phool/PHIOManager.h>

#include <Event/phenixTypes.h>

#include <string>

class PHCompositeNode;
class PHRawDataNode;
class oBuffer;

class PHRawOManager : public PHIOManager
{
 public:
  PHRawOManager();
  PHRawOManager(const std::string &, const int run = 0, const int bufl = 100000, const int evtl = -1, const int complvl = 3);
  virtual ~PHRawOManager();

  bool setFile(const std::string &, const int setRun, const int setBufl, const int setEvtl, const int complvl);
  virtual void closeFile();
  virtual bool write(PHCompositeNode *);
  bool write(PHRawDataNode *);

  virtual void print() const;

 private:
  int filedesc;
  PHDWORD *memBuffer;
  oBuffer *fileBuffer;
  int runNumber;
  int bufferSize;
  int eventLength;
  int compressionLevel;
};

#endif
