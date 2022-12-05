#include "Fun4AllDstInputManager.h"

#include "Fun4AllReturnCodes.h"
#include "Fun4AllServer.h"

#include <ffaobjects/RunHeader.h>
#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIntegrate.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE, PHReadOnly, PHRunTree
#include <phool/phooldefs.h>

#include <TSystem.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/algorithm/string.hpp>
#pragma GCC diagnostic pop

#include <cassert>
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair
#include <vector>    // for vector

class TBranch;

Fun4AllDstInputManager::Fun4AllDstInputManager(const std::string &name, const std::string &nodename, const std::string &topnodename)
  : Fun4AllInputManager(name, nodename, topnodename)
{
  return;
}

Fun4AllDstInputManager::~Fun4AllDstInputManager()
{
  delete IManager;
  delete runNodeSum;
  return;
}

int Fun4AllDstInputManager::fileopen(const std::string &filenam)
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
  if (IManager)
  {
    std::cout << PHWHERE << " IManager pointer is not nullptr but " << IManager
              << std::endl;
    std::cout << "Send mail to off-l with this printout and the macro you used"
              << std::endl;
    std::cout << "Trying to execute IManager->print() to display more info"
              << std::endl;
    std::cout << "Code will probably segfault now" << std::endl;
    IManager->print();
    std::cout << "Have someone look into this problem - Exiting now" << std::endl;
    exit(1);
  }
  // first read the runnode if not disabled
  if (m_ReadRunTTree)
  {
    IManager = new PHNodeIOManager(fullfilename, PHReadOnly, PHRunTree);
    if (IManager->isFunctional())
    {
      runNode = se->getNode(RunNode, TopNodeName());
      IManager->read(runNode);
      // get the current run number
      RunHeader *runheader = findNode::getClass<RunHeader>(runNode, "RunHeader");
      if (runheader)
      {
        SetRunNumber(runheader->get_RunNumber());
      }
      // delete our internal copy of the runnode when opening subsequent files
      if (runNodeCopy)
      {
        std::cout << PHWHERE
                  << " The impossible happened, we have a valid copy of the run node "
                  << runNodeCopy->getName() << " which should be a nullptr"
                  << std::endl;
        gSystem->Exit(1);
      }
      runNodeCopy = new PHCompositeNode("RUNNODECOPY");
      if (!runNodeSum)
      {
        runNodeSum = new PHCompositeNode("RUNNODESUM");
      }
      PHNodeIOManager *tmpIman = new PHNodeIOManager(fullfilename, PHReadOnly, PHRunTree);
      tmpIman->read(runNodeCopy);
      delete tmpIman;

      PHNodeIntegrate integrate;
      integrate.RunNode(runNode);
      integrate.RunSumNode(runNodeSum);
      // run recursively over internal run node copy and integrate objects
      PHNodeIterator mainIter(runNodeCopy);
      mainIter.forEach(integrate);
      // we do not need to keep the internal copy, keeping it would crate
      // problems in case a subsequent file does not contain all the
      // runwise objects from the previous file. Keeping this copy would then
      // integrate the missing object again with the old copy
      delete runNodeCopy;
      runNodeCopy = nullptr;
    }
    // DLW: move the delete outside the if block to cover the case where isFunctional() fails
    delete IManager;
  }
  // now open the dst node
  dstNode = se->getNode(InputNode(), TopNodeName());
  IManager = new PHNodeIOManager(fullfilename, PHReadOnly);
  if (IManager->isFunctional())
  {
    IsOpen(1);
    events_thisfile = 0;
    setBranches();                // set branch selections
    AddToFileOpened(FileName());  // add file to the list of files which were opened
                                  // check if our input file has a sync object or not
    if (IManager->NodeExist(syncdefs::SYNCNODENAME))
    {
      m_HaveSyncObject = 1;
    }
    else
    {
      m_HaveSyncObject = -1;
    }

    return 0;
  }
  else
  {
    std::cout << PHWHERE << ": " << Name() << " Could not open file "
              << FileName() << std::endl;
    delete IManager;
    IManager = nullptr;
    return -1;
  }
}

int Fun4AllDstInputManager::run(const int nevents)
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
  if (Verbosity() > 3)
  {
    std::cout << "Getting Event from " << Name() << std::endl;
  }
readagain:
  PHCompositeNode *dummy;
  int ncount = 0;
  dummy = IManager->read(dstNode);
  while (dummy)
  {
    ncount++;
    if (nevents > 0 && ncount >= nevents)
    {
      break;
    }
    dummy = IManager->read(dstNode);
  }
  if (!dummy)
  {
    fileclose();
    if (!OpenNextFile())
    {
      goto readagain;
    }
    return -1;
  }
  events_total += ncount;
  events_thisfile += ncount;
  // check if the local SubsysReco discards this event
  if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
  {
    goto readagain;
  }
  syncobject = findNode::getClass<SyncObject>(dstNode, syncdefs::SYNCNODENAME);
  return 0;
}

int Fun4AllDstInputManager::fileclose()
{
  if (!IsOpen())
  {
    std::cout << Name() << ": fileclose: No Input file open" << std::endl;
    return -1;
  }
  delete IManager;
  IManager = nullptr;
  IsOpen(0);
  UpdateFileList();
  m_HaveSyncObject = 0;
  return 0;
}

int Fun4AllDstInputManager::GetSyncObject(SyncObject **mastersync)
{
  // here we copy the sync object from the current file to the
  // location pointed to by mastersync. If mastersync is a 0 pointer
  // the syncobject is cloned. If mastersync allready exists the content
  // of syncobject is copied
  if (!(*mastersync))
  {
    if (syncobject)
    {
      *mastersync = dynamic_cast<SyncObject *>(syncobject->CloneMe());
      assert(*mastersync);
    }
  }
  else
  {
    *(*mastersync) = *syncobject;  // copy syncobject content
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

int Fun4AllDstInputManager::SyncIt(const SyncObject *mastersync)
{
  if (!mastersync)
  {
    std::cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << std::endl;
    std::cout << "Most likely your first file does not contain a SyncObject and the file" << std::endl;
    std::cout << "opened by the Fun4AllDstInputManager with Name " << Name() << " has one" << std::endl;
    std::cout << "Change your macro and use the file opened by this input manager as first input" << std::endl;
    std::cout << "and you will be okay. Fun4All will not process the current configuration" << std::endl
              << std::endl;
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  if (!syncobject)
  {
    std::cout << Name() << " no sync object found in this manager but synchronization needed" << std::endl;
    std::cout << "  Check if you have used an empty listfile. If this is not the case - please ask for help" << std::endl;
    std::cout << "This may be a really bad internal problem and cannot continue, exiting now " << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  int iret = syncobject->Different(mastersync);
  if (iret)  // what to do if files are not in sync
  {
    if (mastersync->EventNumber() == -999999)  // first file does not contain sync object
    {
      std::cout << PHWHERE << " Mastersync not filled, your first file does not contain a SyncObject" << std::endl;
      std::cout << "This Event will not be processed further" << std::endl;
    }
    else  // okay try to resync here
    {
      if (Verbosity() > 3)
      {
        std::cout << "Need to Resync, mastersync evt no: " << mastersync->EventNumber()
                  << ", this Event no: " << syncobject->EventNumber() << std::endl;
        std::cout << "mastersync evt counter: " << mastersync->EventCounter()
                  << ", this Event counter: " << syncobject->EventCounter() << std::endl;
        std::cout << "mastersync run number: " << mastersync->RunNumber()
                  << ", this run number: " << syncobject->RunNumber() << std::endl;
      }

      while (syncobject->RunNumber() < mastersync->RunNumber())
      {
        events_skipped_during_sync++;
        if (Verbosity() > 2)
        {
          std::cout << Name() << " Run Number: " << syncobject->RunNumber()
                    << ", master: " << mastersync->RunNumber()
                    << std::endl;
        }
        iret = ReadNextEventSyncObject();
        if (iret)
        {
          return iret;
        }
      }
      int igood = 0;
      if (syncobject->RunNumber() == mastersync->RunNumber())
      {
        igood = 1;
      }
      // only run up the Segment Number if run numbers are identical
      while (syncobject->SegmentNumber() < mastersync->SegmentNumber() && igood)
      {
        events_skipped_during_sync++;
        if (Verbosity() > 2)
        {
          std::cout << Name() << " Segment Number: " << syncobject->SegmentNumber()
                    << ", master: " << mastersync->SegmentNumber()
                    << std::endl;
        }
        iret = ReadNextEventSyncObject();
        if (iret)
        {
          return iret;
        }
      }
      // only run up the Event Counter if run number and segment number are identical
      if (syncobject->SegmentNumber() == mastersync->SegmentNumber() && syncobject->RunNumber() == mastersync->RunNumber())
      {
        igood = 1;
      }
      else
      {
        igood = 0;
      }
      while (syncobject->EventCounter() < mastersync->EventCounter() && igood)
      {
        events_skipped_during_sync++;
        if (Verbosity() > 2)
        {
          std::cout << Name()
                    << ", EventCounter: " << syncobject->EventCounter()
                    << ", master: " << mastersync->EventCounter()
                    << std::endl;
        }
        iret = ReadNextEventSyncObject();
        if (iret)
        {
          return iret;
        }
      }
      // Since up to here we only read the sync object we need to push
      // the current event back inot the root file (subtract one from the
      // local root file event counter) so we can read the full event
      // if it syncs, if it does not sync we also read one event too many
      // (otherwise we cannot determine that we are "too far")
      // and also have to push this one back
      PushBackEvents(1);
      if (syncobject->RunNumber() > mastersync->RunNumber() ||        // check if run number too large
          syncobject->EventCounter() > mastersync->EventCounter() ||  // check if event counter too large
          syncobject->SegmentNumber() > mastersync->SegmentNumber())  // check segment number too large
      {
        // the event from first file which determines the mastersync
        // and which we are trying to find on this file does not exist on this file
        // so: return failure. This will cause the input managers to read
        // the next event from the input files file
        return Fun4AllReturnCodes::SYNC_FAIL;
      }
      // Here the event counter and segment number and run number do agree - we found the right match
      // now read the full event (previously we only read the sync object)
      PHCompositeNode *dummy;
      dummy = IManager->read(dstNode);
      if (!dummy)
      {
        std::cout << PHWHERE << " " << Name() << " Could not read full Event" << std::endl;
        std::cout << "PLEASE NOTIFY PHENIX-OFF-L and post the macro you used" << std::endl;
        fileclose();
        return Fun4AllReturnCodes::SYNC_FAIL;
      }
      iret = syncobject->Different(mastersync);  // final check if they really agree
      if (iret)                                  // if not things are severely wrong
      {
        std::cout << PHWHERE << " MasterSync and SyncObject of " << Name() << " are different" << std::endl;
        std::cout << "This Event will not be processed further, here is some debugging info:" << std::endl;
        std::cout << "PLEASE NOTIFY PHENIX-OFF-L and post the macro you used" << std::endl;
        std::cout << "MasterSync->identify:" << std::endl;
        mastersync->identify();
        std::cout << Name() << ": SyncObject->identify:" << std::endl;
        syncobject->identify();
        return Fun4AllReturnCodes::SYNC_FAIL;
      }
      else if (Verbosity() > 3)
      {
        std::cout << PHWHERE << " Resynchronization successfull for " << Name() << std::endl;
        std::cout << "MasterSync->identify:" << std::endl;
        mastersync->identify();
        std::cout << Name() << ": SyncObject->identify:" << std::endl;
        syncobject->identify();
      }
    }
  }
  //	else
  //		{
  //			std::cout << "No Need to Resync" << std::endl;
  //		}
  return Fun4AllReturnCodes::SYNC_OK;
}

int Fun4AllDstInputManager::ReadNextEventSyncObject()
{
readnextsync:
  static int readfull = 0;  // for some reason all the input managers need to see the same (I think, should look at this at some point)
  if (!IManager)            // in case the old file was exhausted and there is no new file opened
  {
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  if (syncbranchname.empty())
  {
    readfull = 1;  // we need to read a full event to set the root branches to phool nodes right when a new file has been opened
    std::map<std::string, TBranch *>::const_iterator bIter;
    for (bIter = IManager->GetBranchMap()->begin(); bIter != IManager->GetBranchMap()->end(); ++bIter)
    {
      if (Verbosity() > 2)
      {
        std::cout << Name() << ": branch: " << bIter->first << std::endl;
      }
      std::string delimeters = phooldefs::branchpathdelim;  // + phooldefs::legacypathdelims;
      std::vector<std::string> splitvec;
      boost::split(splitvec, bIter->first, boost::is_any_of(delimeters));
      for (auto & ia : splitvec)  // -1 so we skip the node name
      {
        if (ia == syncdefs::SYNCNODENAME)
        {
          syncbranchname = bIter->first;
          break;
        }
      }
      if (!syncbranchname.empty())
      {
        break;
      }
    }
    if (syncbranchname.empty())
    {
      std::cout << PHWHERE << "Could not locate Sync Branch" << std::endl;
      std::cout << "Please check for it in the following list of branch names and" << std::endl;
      std::cout << "PLEASE NOTIFY PHENIX-OFF-L and post the macro you used" << std::endl;
      for (bIter = IManager->GetBranchMap()->begin(); bIter != IManager->GetBranchMap()->end(); ++bIter)
      {
        std::cout << bIter->first << std::endl;
      }
      return Fun4AllReturnCodes::SYNC_FAIL;
    }
  }
  size_t EventOnDst = 0;
  int itest = 0;
  if (!readfull)
  {
    // if all files are exhausted, the IManager is deleted and set to nullptr
    // so check if IManager is valid before getting a new event
    if (IManager)
    {
      EventOnDst = IManager->getEventNumber();  // this returns the next number of the event
      itest = IManager->readSpecific(EventOnDst, syncbranchname.c_str());
    }
    else
    {
      if (Verbosity() > 2)
      {
        std::cout << Name() << ": File exhausted while resyncing" << std::endl;
      }
      return Fun4AllReturnCodes::SYNC_FAIL;
    }
  }
  else
  {
    if (IManager->read(dstNode))
    {
      itest = 1;
    }
    else
    {
      itest = 0;
    }
  }
  if (!itest)
  {
    if (Verbosity() > 2)
    {
      std::cout << Name() << ": File exhausted while resyncing" << std::endl;
    }
    fileclose();
    if (OpenNextFile())
    {
      return Fun4AllReturnCodes::SYNC_FAIL;
    }
    syncbranchname.clear();  // clear the sync branch name, who knows - it might be different on new file
    goto readnextsync;
  }
  if (!readfull)
  {
    EventOnDst++;
    IManager->setEventNumber(EventOnDst);  // update event number in phool io manager
  }
  else
  {
    readfull = 0;
  }
  return 0;
}

int Fun4AllDstInputManager::BranchSelect(const std::string &branch, const int iflag)
{
  if (IsOpen())
  {
    std::cout << "BranchSelect(\"" << branch << "\", " << iflag
              << ") : Input branches can only selected for reading before fileopen is called proceeding without input branch selection" << std::endl;
    return -1;
  }
  // if iflag > 0 the branch is set to read
  // if iflag = 0, the branch is set to NOT read
  // if iflag < 0 the branchname is erased from our internal branch read map
  // this does not have any effect on phool yet
  if (iflag < 0)
  {
    std::map<const std::string, int>::iterator branchiter;
    branchiter = branchread.find(branch);
    if (branchiter != branchread.end())
    {
      branchread.erase(branchiter);
    }
    return 0;
  }
  int readit = 0;
  if (iflag > 0)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Setting Root Tree Branch: " << branch << " to read" << std::endl;
    }
    readit = 1;
  }
  else
  {
    if (Verbosity() > 1)
    {
      std::cout << "Setting Root Tree Branch: " << branch << " to NOT read" << std::endl;
    }
  }
  branchread[branch] = readit;
  return 0;
}

int Fun4AllDstInputManager::setBranches()
{
  if (IManager)
  {
    if (!branchread.empty())
    {
      std::map<const std::string, int>::const_iterator branchiter;
      for (branchiter = branchread.begin(); branchiter != branchread.end(); ++branchiter)
      {
        IManager->selectObjectToRead(branchiter->first.c_str(), branchiter->second);
        if (Verbosity() > 0)
        {
          std::cout << branchiter->first << " set to " << branchiter->second << std::endl;
        }
      }
      // protection against switching off the sync variables
      // only implemented in the Sync Manager
      setSyncBranches(IManager);
    }
  }
  else
  {
    std::cout << PHWHERE << " " << Name() << ": You can only call this function after a file has been opened" << std::endl;
    std::cout << "Do not worry, the branches will be set as soon as you open a file" << std::endl;
    return -1;
  }
  return 0;
}

int Fun4AllDstInputManager::setSyncBranches(PHNodeIOManager *IMan)
{
  // protection against switching off the sync variables
  for (auto & i : syncdefs::SYNCVARS)
  {
    IMan->selectObjectToRead(i, true);
  }
  return 0;
}

void Fun4AllDstInputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "BRANCH")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of selected branches in Fun4AllDstInputManager " << Name() << ":" << std::endl;

    std::map<const std::string, int>::const_iterator iter;
    for (iter = branchread.begin(); iter != branchread.end(); ++iter)
    {
      std::cout << iter->first << " is switched ";
      if (iter->second)
      {
        std::cout << "ON";
      }
      else
      {
        std::cout << "OFF";
      }
      std::cout << std::endl;
    }
  }
  if ((what == "ALL" || what == "PHOOL") && IManager)
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "PHNodeIOManager print in Fun4AllDstInputManager " << Name() << ":" << std::endl;
    IManager->print();
  }
  Fun4AllInputManager::Print(what);
  return;
}

int Fun4AllDstInputManager::PushBackEvents(const int i)
{
  if (IManager)
  {
    unsigned EventOnDst = IManager->getEventNumber();
    EventOnDst -= static_cast<unsigned>(i);
    IManager->setEventNumber(EventOnDst);
    return 0;
  }
  std::cout << PHWHERE << Name() << ": could not push back events, Imanager is NULL"
            << " probably the dst is not open yet (you need to call fileopen or run 1 event for lists)" << std::endl;
  return -1;
}

int Fun4AllDstInputManager::HasSyncObject() const
{
  if (m_HaveSyncObject)
  {
    return m_HaveSyncObject;
  }
  if (IsOpen())
  {
    std::cout << PHWHERE << "HasSyncObject() not initialized check the calling order" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  return 0;
}
