#ifndef PHOOL_PHNODEIOMANAGER_H
#define PHOOL_PHNODEIOMANAGER_H

//  Declaration of class PHNodeIOManager
//  Purpose: manages file IO for PHIODataNodes
//  Author: Matthias Messer

#include "PHIOManager.h"

#include "phool.h"

#include <cstddef>
#include <map>
#include <string>

class PHCompositeNode;
class TBranch;
class TFile;
class TObject;
class TTree;

class PHNodeIOManager : public PHIOManager
{
 public:
  PHNodeIOManager();
  PHNodeIOManager(const std::string &, const PHAccessType = PHReadOnly);
  PHNodeIOManager(const std::string &, const std::string &, const PHAccessType = PHReadOnly);
  PHNodeIOManager(const std::string &, const PHAccessType, const PHTreeType);
  virtual ~PHNodeIOManager();

 public:
  virtual void closeFile();
  virtual bool write(PHCompositeNode *);
  virtual void print() const;

  bool setFile(const std::string &, const std::string &, const PHAccessType = PHReadOnly);
  PHCompositeNode *read(PHCompositeNode * = nullptr, size_t = 0);
  bool read(size_t requestedEvent);
  int readSpecific(size_t requestedEvent, const char *objectName);
  void selectObjectToRead(const char *objectName, bool readit);
  bool isSelected(const char *objectName);
  int isFunctional() const { return isFunctionalFlag; }
  bool SetCompressionLevel(const int level);
  double GetBytesWritten();
  std::map<std::string, TBranch *> *GetBranchMap();

 public:
  bool write(TObject **, const std::string &);

 private:
  int FillBranchMap();
  PHCompositeNode *reconstructNodeTree(PHCompositeNode *);
  bool readEventFromFile(size_t requestedEvent);
  std::string getBranchClassName(TBranch *);

  TFile *file;
  TTree *tree;
  std::string TreeName;
  int bufSize;
  int split;
  int accessMode;
  int CompressionLevel;
  std::map<std::string, TBranch *> fBranches;
  std::map<std::string, bool> objectToRead;

  int isFunctionalFlag;  // flag to tell if that object initialized properly
};

#endif
