#include "Fun4AllNoSyncDstInputManager.h"

#include <cstdlib>
#include <iostream>

Fun4AllNoSyncDstInputManager::Fun4AllNoSyncDstInputManager(const std::string &name,
                                                           const std::string &nodename,
                                                           const std::string &topnodename)
  : Fun4AllDstInputManager(name, nodename, topnodename)
{
  return;
}

int Fun4AllNoSyncDstInputManager::NoRunTTree()
{
  if (!IsOpen())
  {
    ReadRunTTree(0);
  }
  else
  {
    std::cout << Name()
              << "NoRunTTree() has to be done before opening a file" << std::endl;
    exit(1);
  }
  return 0;
}
