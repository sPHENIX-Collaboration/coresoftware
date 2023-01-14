// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLDSTOUTPUTMANAGER_H
#define FUN4ALL_FUN4ALLDSTOUTPUTMANAGER_H

#include "Fun4AllOutputManager.h"

#include <set>
#include <string>

class PHNodeIOManager;
class PHCompositeNode;

class Fun4AllDstOutputManager : public Fun4AllOutputManager
{
 public:
  Fun4AllDstOutputManager(const std::string &myname = "DSTOUT", const std::string &filename = "dstout.root");
  ~Fun4AllDstOutputManager() override;
  // Fun4AllDstOutputManager contains pointer to memory
  // copy ctor and = operator  need explicit implementation, do just delete it here
  Fun4AllDstOutputManager(const Fun4AllDstOutputManager &) = delete;
  Fun4AllDstOutputManager &operator=(Fun4AllDstOutputManager const &) = delete;
  int AddNode(const std::string &nodename) override;
  int AddRunNode(const std::string &nodename) override;
  int StripNode(const std::string &nodename) override;
  int StripRunNode(const std::string &nodename) override;
  void SaveRunNode(const int i) override { m_SaveRunNodeFlag = i; }
  int outfileopen(const std::string &fname) override;

  void Print(const std::string &what = "ALL") const override;

  int Write(PHCompositeNode *startNode) override;
  int WriteNode(PHCompositeNode *thisNode) override;

 private:
  PHNodeIOManager *dstOut = nullptr;
  int m_SaveRunNodeFlag = 1;
  std::set<std::string> savenodes;
  std::set<std::string> saverunnodes;
  std::set<std::string> stripnodes;
  std::set<std::string> striprunnodes;
};

#endif
