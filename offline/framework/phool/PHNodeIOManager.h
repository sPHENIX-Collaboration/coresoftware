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
  PHNodeIOManager() {}
  PHNodeIOManager(const std::string &, const PHAccessType = PHReadOnly);
  PHNodeIOManager(const std::string &, const std::string &, const PHAccessType = PHReadOnly);
  PHNodeIOManager(const std::string &, const PHAccessType, const PHTreeType);
  ~PHNodeIOManager() override;

  void closeFile() override;
  bool write(PHCompositeNode *) override;
  void print() const override;

  bool setFile(const std::string &, const std::string &, const PHAccessType = PHReadOnly);
  PHCompositeNode *read(PHCompositeNode * = nullptr, size_t = 0);
  bool read(size_t requestedEvent);
  int readSpecific(size_t requestedEvent, const std::string &objectName);
  void selectObjectToRead(const std::string &objectName, bool readit);
  bool isSelected(const std::string &objectName);
  int isFunctional() const { return isFunctionalFlag; }
  bool SetCompressionLevel(const int level);
  double GetBytesWritten();
  std::map<std::string, TBranch *> *GetBranchMap();

  bool write(TObject **, const std::string &, int buffersize, int splitlevel);
  bool NodeExist(const std::string &nodename);

 private:
  int FillBranchMap();
  PHCompositeNode *reconstructNodeTree(PHCompositeNode *);
  bool readEventFromFile(size_t requestedEvent);
  std::string getBranchClassName(TBranch *);

  TFile *file = nullptr;
  TTree *tree = nullptr;
  std::string TreeName = "T";
  int accessMode = PHReadOnly;
  int CompressionLevel = 3;
  std::map<std::string, TBranch *> fBranches;
  std::map<std::string, bool> objectToRead;

  int isFunctionalFlag = 0;  // flag to tell if that object initialized properly
};

#endif
