#include "Fun4AllRunNodeInputManager.h"

#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"

#include <ffaobjects/RunHeader.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIntegrate.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE, PHReadOnly, PHRunTree
#include <phool/phooldefs.h>

#include <TSystem.h>

#include <cstdlib>
#include <iostream>

Fun4AllRunNodeInputManager::Fun4AllRunNodeInputManager(const std::string &name,
                                                       const std::string &nodename,
                                                       const std::string &topnodename)
  : Fun4AllDstInputManager(name, nodename, topnodename)
{
  return;
}

int Fun4AllRunNodeInputManager::fileopen(const std::string &filenam)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  if (IsOpen())
  {
    std::cout << "Closing currently open file "
              << FileName()
              << " and opening " << filenam << std::endl;
    fileclose();
  }
  FileName(filenam);
  FROG frog;
  fullfilename = frog.location(FileName());
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": opening file " << fullfilename << std::endl;
  }
  // sanity check - the IManager must be nullptr when this method is executed
  // if not something is very very wrong and we must not continue
  if (IManager())
  {
    std::cout << PHWHERE << " IManager pointer is not nullptr but " << IManager()
              << std::endl;
    std::cout << "Send mail to off-l with this printout and the macro you used"
              << std::endl;
    std::cout << "Trying to execute IManager->print() to display more info"
              << std::endl;
    std::cout << "Code will probably segfault now" << std::endl;
    IManager()->print();
    std::cout << "Have someone look into this problem - Exiting now" << std::endl;
    exit(1);
  }
  IManager(new PHNodeIOManager(fullfilename, PHReadOnly, PHRunTree));
  if (IManager()->isFunctional())
  {
    runNode(se->getNode(RunNodeName(), TopNodeName()));
    IManager()->read(runNode());
    // get the current run number
    RunHeader *runheader = findNode::getClass<RunHeader>(runNode(), "RunHeader");
    if (runheader)
    {
      SetRunNumber(runheader->get_RunNumber());
    }
    // delete our internal copy of the runnode when opening subsequent files
    if (runNodeCopy())
    {
      std::cout << PHWHERE
                << " The impossible happened, we have a valid copy of the run node "
                << runNodeCopy()->getName() << " which should be a nullptr"
                << std::endl;
      gSystem->Exit(1);
    }
    runNodeCopy(new PHCompositeNode("RUNNODECOPY"));
    if (!runNodeSum())
    {
      runNodeSum(new PHCompositeNode("RUNNODESUM"));
    }
    PHNodeIOManager *tmpIman = new PHNodeIOManager(fullfilename, PHReadOnly, PHRunTree);
    tmpIman->read(runNodeCopy());
    delete tmpIman;

    PHNodeIntegrate integrate;
    integrate.RunNode(runNode());
    integrate.RunSumNode(runNodeSum());
    // run recursively over internal run node copy and integrate objects
    PHNodeIterator mainIter(runNodeCopy());
    mainIter.forEach(integrate);
    // we do not need to keep the internal copy, keeping it would crate
    // problems in case a subsequent file does not contain all the
    // runwise objects from the previous file. Keeping this copy would then
    // integrate the missing object again with the old copy
    delete runNodeCopy();
    runNodeCopy(nullptr);
  }
  // DLW: move the delete outside the if block to cover the case where isFunctional() fails
  delete IManager();
  return 0;
}

int Fun4AllRunNodeInputManager::run(const int /*nevents*/)
{
  if (!IsOpen())
  {
    if (FileListEmpty())
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": No Input file open" << std::endl;
      }
      return -1;
    }
    else
    {
      if (OpenNextFile())
      {
        std::cout << Name() << ": No Input file from filelist opened" << std::endl;
        return -1;
      }
    }
  }
  return 0;
}
