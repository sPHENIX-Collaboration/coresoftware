/*!
 * \file Fun4AllDstPileupInputManager.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "Fun4AllDstPileupInputManager.h"

#include "PHG4Hit.h"  // for PHG4Hit
#include "PHG4HitContainer.h"
#include "PHG4Hitv1.h"
#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Particlev3.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4VtxPoint.h"  // for PHG4VtxPoint
#include "PHG4VtxPointv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/EventHeaderv2.h>
#include <ffaobjects/RunHeader.h>
#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <frog/FROG.h>

#include <phhepmc/PHHepMCGenEvent.h>  // for PHHepMCGenEvent
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNode.h>        // for PHNode
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIntegrate.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHNodeOperation.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE, PHReadOnly, PHRunTree

#include <TSystem.h>

#include <HepMC/GenEvent.h>

#include <cassert>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <iterator>  // for reverse_iterator, operator!=
#include <utility>   // for pair

// convenient aliases for deep copying nodes
namespace
{
  using PHG4Particle_t = PHG4Particlev3;
  using PHG4VtxPoint_t = PHG4VtxPointv1;
  using PHG4Hit_t = PHG4Hitv1;

  //! utility class to find all PHG4Hit container nodes from the DST node
  class FindG4HitContainer : public PHNodeOperation
  {
   public:
    //! container map alias
    using ContainerMap = std::map<std::string, PHG4HitContainer *>;

    //! get container map
    const ContainerMap &containers() const
    {
      return m_containers;
    }

   protected:
    //! iterator action
    void perform(PHNode *node) override
    {
      // check type name. Only load PHIODataNode
      if (node->getType() != "PHIODataNode") return;

      // cast to IODataNode and check data
      auto ionode = static_cast<PHIODataNode<TObject> *>(node);
      auto data = dynamic_cast<PHG4HitContainer *>(ionode->getData());
      if (data)
      {
        m_containers.insert(std::make_pair(node->getName(), data));
      }
    }

   private:
    //! container map
    ContainerMap m_containers;
  };

}  // namespace

//_____________________________________________________________________________
Fun4AllDstPileupInputManager::Fun4AllDstPileupInputManager(const std::string &name, const std::string &nodename, const std::string &topnodename)
  : Fun4AllInputManager(name, nodename, topnodename)
{
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::fileopen(const std::string &filenam)
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
  m_fullfilename = frog.location(FileName());
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": opening file " << m_fullfilename << std::endl;
  }
  // sanity check - the IManager must be nullptr when this method is executed
  // if not something is very very wrong and we must not continue
  if (m_IManager)
  {
    std::cout << PHWHERE << " IManager pointer is not nullptr but " << m_IManager.get()
              << std::endl;
    std::cout << "Send mail to off-l with this printout and the macro you used"
              << std::endl;
    std::cout << "Trying to execute IManager->print() to display more info"
              << std::endl;
    std::cout << "Code will probably segfault now" << std::endl;
    m_IManager->print();
    std::cout << "Have someone look into this problem - Exiting now" << std::endl;
    exit(1);
  }
  // first read the runnode if not disabled
  if (m_ReadRunTTree)
  {
    m_IManager.reset(new PHNodeIOManager(m_fullfilename, PHReadOnly, PHRunTree));
    if (m_IManager->isFunctional())
    {
      m_runNode = se->getNode(m_RunNode, TopNodeName());
      m_IManager->read(m_runNode);

      // get the current run number
      RunHeader *runheader = findNode::getClass<RunHeader>(m_runNode, "RunHeader");
      if (runheader)
      {
        SetRunNumber(runheader->get_RunNumber());
      }
      // delete our internal copy of the runnode when opening subsequent files
      if (m_runNodeCopy)
      {
        std::cout << PHWHERE
                  << " The impossible happened, we have a valid copy of the run node "
                  << m_runNodeCopy->getName() << " which should be a nullptr"
                  << std::endl;
        gSystem->Exit(1);
      }
      m_runNodeCopy.reset(new PHCompositeNode("RUNNODECOPY"));
      if (!m_runNodeSum)
      {
        m_runNodeSum.reset(new PHCompositeNode("RUNNODESUM"));
      }

      {
        std::unique_ptr<PHNodeIOManager> tmpIman(new PHNodeIOManager(m_fullfilename, PHReadOnly, PHRunTree));
        tmpIman->read(m_runNodeCopy.get());
      }

      PHNodeIntegrate integrate;
      integrate.RunNode(m_runNode);
      integrate.RunSumNode(m_runNodeSum.get());
      // run recursively over internal run node copy and integrate objects
      PHNodeIterator mainIter(m_runNodeCopy.get());
      mainIter.forEach(integrate);
      // we do not need to keep the internal copy, keeping it would crate
      // problems in case a subsequent file does not contain all the
      // runwise objects from the previous file. Keeping this copy would then
      // integrate the missing object again with the old copy
      m_runNodeCopy.reset();
    }
    // DLW: move the delete outside the if block to cover the case where isFunctional() fails
    m_IManager.reset();
  }

  // create internal dst node
  if (!m_dstNodeInternal)
  {
    m_dstNodeInternal.reset(new PHCompositeNode("DST_INTERNAL"));
  }

  // update dst node from fun4all server
  m_dstNode = se->getNode(InputNode(), TopNodeName());

  // open file in both active and background input manager
  m_IManager.reset(new PHNodeIOManager(m_fullfilename, PHReadOnly));
  m_IManager_background.reset(new PHNodeIOManager(m_fullfilename, PHReadOnly));
  if (m_IManager->isFunctional())
  {
    IsOpen(1);
    m_ievent_thisfile = 0;
    setBranches();                // set branch selections
    AddToFileOpened(FileName());  // add file to the list of files which were opened
    return 0;
  }
  else
  {
    std::cout << PHWHERE << ": " << Name() << " Could not open file " << FileName() << std::endl;
    m_IManager.reset();
    m_IManager_background.reset();
    return -1;
  }
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::run(const int nevents)
{
  if (!IsOpen())
  {
    if (FileListEmpty())
    {
      if (Verbosity() > 0) std::cout << Name() << ": No Input file open" << std::endl;
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

  // read main event to dstNode
  auto dummy = m_IManager->read(m_dstNode);
  int ncount = 0;
  while (dummy)
  {
    ++ncount;
    if (nevents > 0 && ncount >= nevents) break;
    dummy = m_IManager->read(m_dstNode);
  }

  if (!dummy)
  {
    fileclose();
    if (!OpenNextFile())
      goto readagain;
    else
      return -1;
  }

  m_ievent_total += ncount;
  m_ievent_thisfile += ncount;

  // check events consistency
  if (m_ievent_thisfile != static_cast<int>(m_IManager->getEventNumber()))
  {
    std::cout << PHWHERE
              << " inconsistent event counters between inputmanager and ionode manager: "
              << " m_ievent_thisfile: " << m_ievent_thisfile
              << " m_IManager->getEventNumber(): " << m_IManager->getEventNumber()
              << std::endl;
    gSystem->Exit(1);
  }

  // get bunchCrossing  associated to this event
  const size_t ievent = m_ievent_total + m_event_offset - 1;
  if (ievent >= m_bunchCrossings.size()) return Fun4AllReturnCodes::ABORTRUN;

  const auto bunchCrossing = m_bunchCrossings[ievent];

  if (m_events_accepted > 0 )
  {
    // reject if event is so close from last included event that the last one would be counted as background to this one
    const double delta_t = m_time_between_crossings * (m_last_bunchCrossing - bunchCrossing);
    if( delta_t >= m_tmin )
    {
      // reject if event if too close to previous trigger
      std::cout << "Fun4AllDstPileupInputManager::run - skipped event " << m_ievent_thisfile - 1 << std::endl;
      goto readagain;
    }
  }

  // check if the local SubsysReco discards this event
  if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << "Fun4AllDstPileupInputManager::run - skipped event " << m_ievent_thisfile - 1 << std::endl;
    goto readagain;
  }

  // load relevant DST nodes to internal pointers
  std::cout << "Fun4AllDstPileupInputManager::run - loaded event " << m_ievent_thisfile - 1 << std::endl;

  // store bunch crossing both inside the event header and as the last bunch crossing
  load_nodes(m_dstNode);
  if (m_eventheader) m_eventheader->set_BunchCrossing(bunchCrossing);

  // store as last used bunch crossing to avoid overlap with next event
  m_last_bunchCrossing = bunchCrossing;

  // fetch past events falling into the TPC integration time
  int neventspast = 0;
  for (int ieventpast = ievent - 1; ieventpast >= 0; ieventpast--)
  {
    // check time
    const double delta_t = m_time_between_crossings * (m_bunchCrossings[ieventpast] - bunchCrossing);
    if (delta_t < m_tmin) break;

    // get relevant event in file index
    const int ievent_thisfile = m_ievent_thisfile + ieventpast - ievent - 1;

    // try read
    if (ievent_thisfile < 0 || !m_IManager_background->read(m_dstNodeInternal.get(), ievent_thisfile)) break;

    // merge
    copy_background_event(m_dstNodeInternal.get(), delta_t);

    // increment number of read events
    std::cout << "Fun4AllDstPileupInputManager::run - merged past event " << ievent_thisfile << std::endl;
    ++neventspast;
  }

  // fetch future events falling into the TPC integration time
  int neventsfuture = 0;
  for (size_t ieventfuture = ievent + 1; ieventfuture < m_bunchCrossings.size(); ++ieventfuture)
  {
    // check time
    const double delta_t = m_time_between_crossings * (m_bunchCrossings[ieventfuture] - bunchCrossing);
    if (delta_t > m_tmax) break;

    // get relevant event in file index
    const int ievent_thisfile = m_ievent_thisfile + ieventfuture - ievent - 1;

    // try read
    if (!m_IManager_background->read(m_dstNodeInternal.get(), ievent_thisfile)) break;

    // merge
    copy_background_event(m_dstNodeInternal.get(), delta_t);

    // increment number of read events
    std::cout << "Fun4AllDstPileupInputManager::run - merged future event " << ievent_thisfile << std::endl;
    ++neventsfuture;

    // store as last used bunch crossing to avoid overlap with next event
    m_last_bunchCrossing = m_bunchCrossings[ieventfuture];

  }

  // update event counter
  ++m_events_accepted;

  // jump event counter to the last background accepted event
  if( neventsfuture > 0 ) PushBackEvents( -neventsfuture );

  // update syncobject
  m_syncobject = findNode::getClass<SyncObject>(m_dstNode, "Sync");
  return 0;
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::fileclose()
{
  if (!IsOpen())
  {
    std::cout << Name() << ": fileclose: No Input file open" << std::endl;
    return -1;
  }
  m_IManager.reset();
  m_IManager_background.reset();
  IsOpen(0);
  UpdateFileList();
  return 0;
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::GetSyncObject(SyncObject **mastersync)
{
  // here we copy the sync object from the current file to the
  // location pointed to by mastersync. If mastersync is a 0 pointer
  // the syncobject is cloned. If mastersync allready exists the content
  // of syncobject is copied
  if (!(*mastersync))
  {
    if (m_syncobject)
    {
      *mastersync = dynamic_cast<SyncObject *>(m_syncobject->CloneMe());
      assert(*mastersync);
    }
  }
  else
  {
    *(*mastersync) = *m_syncobject;  // copy syncobject content
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::SyncIt(const SyncObject *mastersync)
{
  if (!mastersync)
  {
    std::cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << std::endl;
    std::cout << "Most likely your first file does not contain a SyncObject and the file" << std::endl;
    std::cout << "opened by the Fun4AllDstPileupInputManager with Name " << Name() << " has one" << std::endl;
    std::cout << "Change your macro and use the file opened by this input manager as first input" << std::endl;
    std::cout << "and you will be okay. Fun4All will not process the current configuration" << std::endl
              << std::endl;
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  int iret = m_syncobject->Different(mastersync);
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
        std::cout
            << "Need to Resync, mastersync evt no: " << mastersync->EventNumber()
            << ", this Event no: " << m_syncobject->EventNumber()
            << std::endl;
        std::cout
            << "mastersync evt counter: " << mastersync->EventCounter()
            << ", this Event counter: " << m_syncobject->EventCounter()
            << std::endl;
        std::cout
            << "mastersync run number: " << mastersync->RunNumber()
            << ", this run number: " << m_syncobject->RunNumber()
            << std::endl;
      }
      while (m_syncobject->RunNumber() < mastersync->RunNumber())
      {
        m_events_skipped_during_sync++;
        if (Verbosity() > 2)
        {
          std::cout
              << Name() << " Run Number: " << m_syncobject->RunNumber()
              << ", master: " << mastersync->RunNumber()
              << std::endl;
        }
        int iret = ReadNextEventSyncObject();
        if (iret) return iret;
      }
      bool igood(m_syncobject->RunNumber() == mastersync->RunNumber());

      // only run up the Segment Number if run numbers are identical
      while (igood && m_syncobject->SegmentNumber() < mastersync->SegmentNumber())
      {
        m_events_skipped_during_sync++;
        if (Verbosity() > 2)
        {
          std::cout << Name() << " Segment Number: " << m_syncobject->SegmentNumber()
                    << ", master: " << mastersync->SegmentNumber()
                    << std::endl;
        }
        int iret = ReadNextEventSyncObject();
        if (iret)
        {
          return iret;
        }
      }
      // only run up the Event Counter if run number and segment number are identical
      igood = (m_syncobject->SegmentNumber() == mastersync->SegmentNumber() && m_syncobject->RunNumber() == mastersync->RunNumber());
      while (igood && m_syncobject->EventCounter() < mastersync->EventCounter())
      {
        m_events_skipped_during_sync++;
        if (Verbosity() > 2)
        {
          std::cout << Name()
                    << ", EventCounter: " << m_syncobject->EventCounter()
                    << ", master: " << mastersync->EventCounter()
                    << std::endl;
        }
        int iret = ReadNextEventSyncObject();
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
      if (m_syncobject->RunNumber() > mastersync->RunNumber() ||        // check if run number too large
          m_syncobject->EventCounter() > mastersync->EventCounter() ||  // check if event counter too large
          m_syncobject->SegmentNumber() > mastersync->SegmentNumber())  // check segment number too large
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
      dummy = m_IManager->read(m_dstNode);
      if (!dummy)
      {
        std::cout << PHWHERE << " " << Name() << " Could not read full Event" << std::endl;
        std::cout << "PLEASE NOTIFY PHENIX-OFF-L and post the macro you used" << std::endl;
        fileclose();
        return Fun4AllReturnCodes::SYNC_FAIL;
      }
      int iret = m_syncobject->Different(mastersync);  // final check if they really agree
      if (iret)                                        // if not things are severely wrong
      {
        std::cout << PHWHERE << " MasterSync and SyncObject of " << Name() << " are different" << std::endl;
        std::cout << "This Event will not be processed further, here is some debugging info:" << std::endl;
        std::cout << "PLEASE NOTIFY PHENIX-OFF-L and post the macro you used" << std::endl;
        std::cout << "MasterSync->identify:" << std::endl;
        mastersync->identify();
        std::cout << Name() << ": SyncObject->identify:" << std::endl;
        m_syncobject->identify();
        return Fun4AllReturnCodes::SYNC_FAIL;
      }
      else if (Verbosity() > 3)
      {
        std::cout << PHWHERE << " Resynchronization successfull for " << Name() << std::endl;
        std::cout << "MasterSync->identify:" << std::endl;
        mastersync->identify();
        std::cout << Name() << ": SyncObject->identify:" << std::endl;
        m_syncobject->identify();
      }
    }
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::ReadNextEventSyncObject()
{
readnextsync:
  static int readfull = 0;
  if (!m_IManager)
  {
    // in case the old file was exhausted and there is no new file opened
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  if (m_syncbranchname.empty())
  {
    readfull = 1;  // we need to read a full events to set the root branches to phool nodes right for a new file
    for (auto bIter = m_IManager->GetBranchMap()->begin(); bIter != m_IManager->GetBranchMap()->end(); ++bIter)
    {
      if (Verbosity() > 2)
      {
        std::cout << Name() << ": branch: " << bIter->first << std::endl;
      }

      const auto pos = bIter->first.find("/Sync");
      if (pos != std::string::npos)
      {
        m_syncbranchname = bIter->first;
        break;
      }
    }
    if (m_syncbranchname.empty())
    {
      std::cout << PHWHERE << "Could not locate Sync Branch" << std::endl;
      std::cout << "Please check for it in the following list of branch names and" << std::endl;
      std::cout << "PLEASE NOTIFY PHENIX-OFF-L and post the macro you used" << std::endl;
      for (auto bIter = m_IManager->GetBranchMap()->begin(); bIter != m_IManager->GetBranchMap()->end(); ++bIter)
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
    if (m_IManager)
    {
      EventOnDst = m_IManager->getEventNumber();  // this returns the next number of the event
      itest = m_IManager->readSpecific(EventOnDst, m_syncbranchname.c_str());
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
    if (m_IManager->read(m_dstNode))
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
    m_syncbranchname.clear();  // clear the sync branch name, who knows - it might be different on new file
    goto readnextsync;
  }
  if (!readfull)
  {
    ++EventOnDst;
    ++m_ievent_thisfile;
    ++m_ievent_total;
    m_IManager->setEventNumber(EventOnDst);  // update event number in phool io manager
  }
  else
  {
    readfull = 0;
  }
  return 0;
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::BranchSelect(const std::string &branch, const int iflag)
{
  int myflag = iflag;
  // if iflag > 0 the branch is set to read
  // if iflag = 0, the branch is set to NOT read
  // if iflag < 0 the branchname is erased from our internal branch read map
  // this does not have any effect on phool yet
  if (myflag < 0)
  {
    std::map<const std::string, int>::iterator branchiter;
    branchiter = m_branchread.find(branch);
    if (branchiter != m_branchread.end())
    {
      m_branchread.erase(branchiter);
    }
    return 0;
  }

  if (myflag > 0)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Setting Root Tree Branch: " << branch << " to read" << std::endl;
    }
    myflag = 1;
  }
  else
  {
    if (Verbosity() > 1)
    {
      std::cout << "Setting Root Tree Branch: " << branch << " to NOT read" << std::endl;
    }
  }
  m_branchread[branch] = myflag;
  return 0;
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::setBranches()
{
  if (m_IManager)
  {
    if (!m_branchread.empty())
    {
      std::map<const std::string, int>::const_iterator branchiter;
      for (branchiter = m_branchread.begin(); branchiter != m_branchread.end(); ++branchiter)
      {
        m_IManager->selectObjectToRead(branchiter->first.c_str(), branchiter->second);
        if (Verbosity() > 0)
        {
          std::cout << branchiter->first << " set to " << branchiter->second << std::endl;
        }
      }
      // protection against switching off the sync variables
      // only implemented in the Sync Manager
      setSyncBranches(m_IManager.get());
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

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::setSyncBranches(PHNodeIOManager *IManager)
{
  // protection against switching off the sync variables
  for (int i = 0; i < syncdefs::NUM_SYNC_VARS; i++)
  {
    IManager->selectObjectToRead(syncdefs::SYNCVARS[i], 1);
  }
  return 0;
}

//_____________________________________________________________________________
void Fun4AllDstPileupInputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "BRANCH")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of selected branches in Fun4AllDstPileupInputManager " << Name() << ":" << std::endl;

    std::map<const std::string, int>::const_iterator iter;
    for (iter = m_branchread.begin(); iter != m_branchread.end(); ++iter)
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
  if ((what == "ALL" || what == "PHOOL") && m_IManager)
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "PHNodeIOManager print in Fun4AllDstPileupInputManager " << Name() << ":" << std::endl;
    m_IManager->print();
  }
  Fun4AllInputManager::Print(what);
  return;
}

//_____________________________________________________________________________
int Fun4AllDstPileupInputManager::PushBackEvents(const int i)
{
  if (m_IManager)
  {
    unsigned EventOnDst = m_IManager->getEventNumber();
    EventOnDst -= static_cast<unsigned>(i);
    m_ievent_thisfile -= i;
    m_ievent_total -= i;
    m_IManager->setEventNumber(EventOnDst);
    return 0;
  }
  std::cout << PHWHERE << Name() << ": could not push back events, Imanager is NULL"
            << " probably the dst is not open yet (you need to call fileopen or run 1 event for lists)" << std::endl;
  return -1;
}

//_____________________________________________________________________________
void Fun4AllDstPileupInputManager::load_nodes(PHCompositeNode *dstNode)
{
  // event header
  m_eventheader = findNode::getClass<EventHeader>(dstNode, "EventHeader");
  if (!m_eventheader)
  {
    std::cout << "Fun4AllDstPileupInputManager::load_nodes - creating EventHeader" << std::endl;
    m_eventheader = new EventHeaderv2();
    dstNode->addNode(new PHIODataNode<PHObject>(m_eventheader, "EventHeader", "PHObject"));
  }

  // hep mc
  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(dstNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "Fun4AllDstPileupInputManager::load_nodes - creating PHHepMCGenEventMap" << std::endl;
    m_geneventmap = new PHHepMCGenEventMap();
    dstNode->addNode(new PHIODataNode<PHObject>(m_geneventmap, "PHHepMCGenEventMap", "PHObject"));
  }

  // find all G4Hit containers under dstNode
  FindG4HitContainer nodeFinder;
  PHNodeIterator(dstNode).forEach(nodeFinder);
  m_g4hitscontainers = nodeFinder.containers();

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(dstNode, "G4TruthInfo");
  if (!m_g4truthinfo)
  {
    std::cout << "Fun4AllDstPileupInputManager::load_nodes - creating node G4TruthInfo" << std::endl;
    m_g4truthinfo = new PHG4TruthInfoContainer();
    dstNode->addNode(new PHIODataNode<PHObject>(m_g4truthinfo, "G4TruthInfo", "PHObject"));
  }
}

//_____________________________________________________________________________
void Fun4AllDstPileupInputManager::copy_background_event(PHCompositeNode *dstNode, double delta_t)
{
  // copy PHHepMCGenEventMap
  const auto map = findNode::getClass<PHHepMCGenEventMap>(dstNode, "PHHepMCGenEventMap");

  // keep track of new embed id, after insertion as background event
  int new_embed_id = 0;

  if (map && m_geneventmap)
  {
    if (map->size() != 1)
    {
      std::cout << "Fun4AllDstPileupInputManager::copy_background_event - cannot merge events that contain more than one PHHepMCGenEventMap" << std::endl;
      return;
    }

    // get event and insert in new map
    auto genevent = map->get_map().begin()->second;
    auto newevent = m_geneventmap->insert_background_event(genevent);

    /*
     * this hack prevents a crash when writting out
     * it boils down to root trying to write deleted items from the HepMC::GenEvent copy if the source has been deleted
     * it does not happen if the source gets written while the copy is deleted
     */
    newevent->getEvent()->swap(*genevent->getEvent());

    // shift vertex time and store new embed id
    newevent->moveVertex(0, 0, 0, delta_t);
    new_embed_id = newevent->get_embedding_id();
  }

  // copy truth container
  // keep track of the correspondance between source index and destination index for vertices, tracks and showers
  using ConversionMap = std::map<int, int>;
  ConversionMap vtxid_map;
  ConversionMap trkid_map;

  const auto container = findNode::getClass<PHG4TruthInfoContainer>(dstNode, "G4TruthInfo");
  if (container && m_g4truthinfo)
  {
    {
      // primary vertices
      auto key = m_g4truthinfo->maxvtxindex();
      const auto range = container->GetPrimaryVtxRange();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        // clone vertex, insert in map, and add index conversion
        const auto &sourceVertex = iter->second;
        auto newVertex = new PHG4VtxPoint_t(sourceVertex);
        newVertex->set_t(sourceVertex->get_t() + delta_t);
        m_g4truthinfo->AddVertex(++key, newVertex);
        vtxid_map.insert(std::make_pair(sourceVertex->get_id(), key));
      }
    }

    {
      // secondary vertices
      auto key = m_g4truthinfo->minvtxindex();
      const auto range = container->GetSecondaryVtxRange();

      // loop from last to first to preserve order with respect to the original event
      for (
          auto iter = std::reverse_iterator<PHG4TruthInfoContainer::ConstVtxIterator>(range.second);
          iter != std::reverse_iterator<PHG4TruthInfoContainer::ConstVtxIterator>(range.first);
          ++iter)
      {
        // clone vertex, shift time, insert in map, and add index conversion
        const auto &sourceVertex = iter->second;
        auto newVertex = new PHG4VtxPoint_t(sourceVertex);
        newVertex->set_t(sourceVertex->get_t() + delta_t);
        m_g4truthinfo->AddVertex(--key, newVertex);
        vtxid_map.insert(std::make_pair(sourceVertex->get_id(), key));
      }
    }

    {
      // primary particles
      auto key = m_g4truthinfo->maxtrkindex();
      const auto range = container->GetPrimaryParticleRange();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto &source = iter->second;
        auto dest = new PHG4Particle_t(source);
        m_g4truthinfo->AddParticle(++key, dest);
        dest->set_track_id(key);

        // set parent to zero
        dest->set_parent_id(0);

        // set primary to itself
        dest->set_primary_id(dest->get_track_id());

        // update vertex
        const auto keyiter = vtxid_map.find(source->get_vtx_id());
        if (keyiter != vtxid_map.end())
          dest->set_vtx_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupInputManager::copy_background_event - vertex id " << source->get_vtx_id() << " not found in map" << std::endl;

        // insert in map
        trkid_map.insert(std::make_pair(source->get_track_id(), dest->get_track_id()));
      }
    }

    {
      // secondary particles
      auto key = m_g4truthinfo->mintrkindex();
      const auto range = container->GetSecondaryParticleRange();

      /*
       * loop from last to first to preserve order with respect to the original event
       * also this ensures that for a given particle its parent has already been converted and thus found in the map
       */
      for (
          auto iter = std::reverse_iterator<PHG4TruthInfoContainer::ConstIterator>(range.second);
          iter != std::reverse_iterator<PHG4TruthInfoContainer::ConstIterator>(range.first);
          ++iter)
      {
        const auto &source = iter->second;
        auto dest = new PHG4Particle_t(source);
        m_g4truthinfo->AddParticle(--key, dest);
        dest->set_track_id(key);

        // update parent id
        auto keyiter = trkid_map.find(source->get_parent_id());
        if (keyiter != trkid_map.end())
          dest->set_parent_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupInputManager::copy_background_event - track id " << source->get_parent_id() << " not found in map" << std::endl;

        // update primary id
        keyiter = trkid_map.find(source->get_primary_id());
        if (keyiter != trkid_map.end())
          dest->set_primary_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupInputManager::copy_background_event - track id " << source->get_primary_id() << " not found in map" << std::endl;

        // update vertex
        keyiter = vtxid_map.find(source->get_vtx_id());
        if (keyiter != vtxid_map.end())
          dest->set_vtx_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupInputManager::copy_background_event - vertex id " << source->get_vtx_id() << " not found in map" << std::endl;

        // insert in map
        trkid_map.insert(std::make_pair(source->get_track_id(), dest->get_track_id()));
      }
    }

    // vertex embed flags
    /* embed flag is stored only for primary vertices, consistently with PHG4TruthEventAction */
    for (const auto &pair : vtxid_map)
    {
      if (pair.first > 0) m_g4truthinfo->AddEmbededVtxId(pair.second, new_embed_id);
    }

    // track embed flags
    /* embed flag is stored only for primary tracks, consistently with PHG4TruthEventAction */
    for (const auto &pair : trkid_map)
    {
      if (pair.first > 0) m_g4truthinfo->AddEmbededTrkId(pair.second, new_embed_id);
    }
  }

  // copy g4hits
  // loop over registered maps
  for (const auto &pair : m_g4hitscontainers)
  {
    // check destination node
    if (!pair.second)
    {
      std::cout << "Fun4AllDstPileupInputManager::copy_background_event - invalid destination container " << pair.first << std::endl;
      continue;
    }

    // find source node
    auto container = findNode::getClass<PHG4HitContainer>(dstNode, pair.first);
    if (!container)
    {
      std::cout << "Fun4AllDstPileupInputManager::copy_background_event - invalid source container " << pair.first << std::endl;
      continue;
    }

    {
      // hits
      const auto range = container->getHits();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        // clone hit
        const auto &sourceHit = iter->second;
        auto newHit = new PHG4Hit_t(sourceHit);

        // shift time
        newHit->set_t(0, sourceHit->get_t(0) + delta_t);
        newHit->set_t(1, sourceHit->get_t(1) + delta_t);

        // update track id
        const auto keyiter = trkid_map.find(sourceHit->get_trkid());
        if (keyiter != trkid_map.end())
          newHit->set_trkid(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupInputManager::copy_background_event - track id " << sourceHit->get_trkid() << " not found in map" << std::endl;

        /*
         * reset shower ids
         * it was decided that showers from the background events will not be copied to the merged event
         * as such we just reset the hits shower id
         */
        newHit->set_shower_id(INT_MIN);

        /*
         * this will generate a new key for the hit and assign it to the hit
         * this ensures that there is no conflict with the hits from the 'main' event
         */
        pair.second->AddHit(newHit->get_detid(), newHit);
      }
    }

    {
      // layers
      const auto range = container->getLayers();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        pair.second->AddLayer(*iter);
      }
    }
  }
}
