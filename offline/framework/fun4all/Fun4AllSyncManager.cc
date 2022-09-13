#include "Fun4AllSyncManager.h"

#include "Fun4AllInputManager.h"
#include "Fun4AllReturnCodes.h"  // for EVENT_OK, RESET_NODE_TREE
#include "Fun4AllServer.h"

#include <ffaobjects/SyncObject.h>

#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <cstdlib>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <list>      // for list<>::const_iterator, _List_con...
#include <string>
#include <utility>  // for pair
#include <vector>

Fun4AllSyncManager::Fun4AllSyncManager(const std::string &name)
  : Fun4AllBase(name)
{
  return;
}

Fun4AllSyncManager::~Fun4AllSyncManager()
{
  delete m_MasterSync;
  while (m_InManager.begin() != m_InManager.end())
  {
    if (Verbosity())
    {
      m_InManager.back()->Verbosity(Verbosity());
    }
    delete m_InManager.back();
    m_InManager.pop_back();
  }
  return;
}

int Fun4AllSyncManager::registerInputManager(Fun4AllInputManager *InputManager)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (inman->Name() == InputManager->Name())
    {
      std::cout << "InputManager " << InputManager->Name() << " allready in list" << std::endl;
      return -1;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "Registering InputManager " << InputManager->Name() << std::endl;
  }
  m_InManager.push_back(InputManager);
  m_iretInManager.push_back(0);
  InputManager->setSyncManager(this);
  return 0;
}

Fun4AllInputManager *
Fun4AllSyncManager::getInputManager(const std::string &name)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (name == inman->Name())
    {
      return inman;
    }
  }
  std::cout << Name() << ": Could not find InputManager" << name << std::endl;
  return nullptr;
}

//_________________________________________________________________
int Fun4AllSyncManager::run(const int nevnts)
{
  int iret = 0;
  int icnt = 0;
  int iretsync = 0;
  int resetnodetree = 0;
  // on read errors (assumed to be caused that a file is empty and we need to open the next one)
  // we have to go through this 3 times
  // 1st pass: The error is detected (for all input mgrs), for the failed input manager(s) a fileclose(), fileopen() is executed, if this fails control goes back to Fun4All since we are done. For input managers without errors, the event is pushed back, so it is read again in the next pass
  // 2nd pass: Events are read from all input managers. If no errors all input managers are pushed back
  // The reason for this is that we opened a new files with maybe different content and we need to clean
  // the node tree so we do not propagate old objects from the previous event
  // The node tree reset is done by the Fun4AllServer so we give control back and signal via resetnodetree return code
  // 3rd pass: read from every input manager and go back to Fun4All

  while (!iret)
  {
    unsigned iman = 0;
    int ifirst = 0;
    int hassync = 0;
    for (auto & iter : m_InManager)
    {
      m_iretInManager[iman] = iter->run(1);
      iret += m_iretInManager[iman];
      // one can run DSTs without sync object via the DST input manager
      // this only poses a problem if one runs two of them and expects the syncing to work
      // or mix DSTs with sync object and without
      if (!hassync && iter->HasSyncObject())  // only update if hassync is 0 and input mgr is non zero
      {
        hassync = iter->HasSyncObject();
      }
      else
      {
        if (iter->HasSyncObject())  // if zero (no syncing) no need to go further
        {
          if (hassync != iter->HasSyncObject())  // we have sync and no sync mixed
          {
            PrintSyncProblem();
            gSystem->Exit(1);
          }
          else if (hassync < 0)  // we have more than one nosync input
          {
            PrintSyncProblem();
            gSystem->Exit(1);
          }
        }
      }
      if (!ifirst)
      {
        if (!m_iretInManager[iman])
        {
          if (!(iter->GetSyncObject(&m_MasterSync)))  // NoSync managers return non zero
          {
            ifirst = 1;
          }
        }
      }
      else
      {
        iretsync = CheckSync(iman);
        if (iretsync)
        {
          break;
        }
      }
      iman++;
    }

    // check event reading, syncronisation
    if (iret || iretsync)
    {
      // tell the server to reset the node tree
      resetnodetree = Fun4AllReturnCodes::RESET_NODE_TREE;

      // if there was an io error (file exhausted) we nee to push back
      // the events from files which are not exhausted yet into the root files
      // here we check the return codes for each input manager and if the
      // read was successful (iret = 0) we push the event back
      if (iret)
      {
        unsigned inputmgr_cnt = 0;
        std::vector<Fun4AllInputManager *>::const_iterator InIter;
        // set in the macro for Sync Manager. Permanently enabled (default when using syncman->Repeat(),
        // m_Repeat = -1 so the m_Repeat--; is not called, this is used when you give it a positive number of
        // repetitions
        if (m_Repeat)
        {
          for (InIter = m_InManager.begin(); InIter != m_InManager.end(); ++InIter)
          {
            if (m_iretInManager[inputmgr_cnt] == Fun4AllReturnCodes::EVENT_OK)
            {
              (*InIter)->PushBackEvents(1);
            }
            else
            {
              if ((*InIter)->IsOpen())
              {
                (*InIter)->fileclose();
              }
              int ireset = (*InIter)->ResetFileList();
              if (ireset)
              {
                std::cout << "Resetting input manager " << (*InIter)->Name() << " failed during Repeat" << std::endl;
                exit(1);
              }
              inputmgr_cnt++;
            }
          }
          if (m_Repeat > 0)
          {
            m_Repeat--;
          }
          iret = 0;
          continue;  // got back and run again
        }
        // push back events where the Imanager did not report an error
        InIter = m_InManager.begin();
        for (int iter : m_iretInManager)
        {
          if (Verbosity() > 0)
          {
            std::cout << (*InIter)->Name() << ": return code: " << iter << std::endl;
          }
          if (!iter)
          {
            (*InIter)->PushBackEvents(1);
            if (Verbosity() > 0)
            {
              std::cout << (*InIter)->Name() << ": push evts: " << iter << std::endl;
            }
          }
          ++InIter;
        }
        goto readerror;
      }
      else
      {
        // just read the next event and hope it syncs
        // this won't update the event counter
        for (unsigned nman = 0; nman < iman; nman++)
        {
          m_InManager[nman]->NoSyncPushBackEvents(1);
        }
        continue;
      }
    }
    if (!resetnodetree)
    {
      m_EventsTotal++;
    }
    // this check is meaningless nowadays since we call this method with 1 event every time
    // so we can just break here but maybe this changes in the future
    if (nevnts > 0 && ++icnt >= nevnts)
    {
      break;
    }
  }

readerror:
  if (!iret)
  {
    if (!resetnodetree)  // all syncing is done and no read errors --> we have a good event in memory
    {
      m_CurrentRun = 0;  // reset current run to 0
      for (auto & iter : m_InManager)
      {
        int runno = iter->RunNumber();
        if (Verbosity() > 2)
        {
          std::cout << Name() << " input mgr " << iter->Name() << " run: " << runno << std::endl;
        }
        if (runno != 0)
        {
          if (m_CurrentRun == 0)
          {
            m_CurrentRun = runno;
            continue;
          }
          else
          {
            if (m_CurrentRun != runno)
            {
              std::cout << "Mixing run numbers (except runnumber=0 which means no valid runnumber) is not supported" << std::endl;
              std::cout << "Here are the list of input managers and runnumbers:" << std::endl;
              for (Fun4AllInputManager *inman : m_InManager)
              {
                std::cout << inman->Name() << " runno: " << inman->RunNumber() << std::endl;
              }
              std::cout << "Exiting now" << std::endl;
              exit(1);
            }
          }
        }
      }
    }
    return resetnodetree;
  }
  return iret;
}

//_________________________________________________________________
int Fun4AllSyncManager::skip(const int nevnts)
{
  if (!m_InManager.empty())
  {
    int Npushback = -nevnts;
    // the PushBackEvents(const int nevents) "pushes back" events events into the root file
    // (technically it just decrements the local counter in the PHNodeIOManager)
    // giving it a negative argument will skip events
    // this is much faster than actually reading the events in
    int iret = m_InManager[0]->PushBackEvents(Npushback);
    for (unsigned int i = 1; i < m_InManager.size(); ++i)
    {
      iret += m_InManager[i]->SkipForThisManager(nevnts);
    }
    if (!iret)
    {
      return 0;
    }
    else
    {
      std::cout << PHWHERE << " Error during skipping events" << std::endl;
      return iret;
    }
  }
  std::cout << PHWHERE << " Cannot skip events: No Input Managers registered?" << std::endl;
  Print("INPUTMANAGER");
  std::cout << "If there are Input Managers in this list, send mail with this" << std::endl;
  std::cout << "error message to off-l" << std::endl;
  std::cout << "and include the macro you used" << std::endl;
  return -1;
}

//_________________________________________________________________
int Fun4AllSyncManager::fileopen(const std::string &managername, const std::string &filename)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (managername == inman->Name())
    {
      int iret = inman->fileopen(filename);
      return iret;
    }
  }
  std::cout << "No Input Manager " << managername << " registered" << std::endl;
  return -1;
}

int Fun4AllSyncManager::BranchSelect(const std::string &managername, const std::string &branch, const int iflag)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (managername == inman->Name())
    {
      int iret = inman->BranchSelect(branch, iflag);
      return iret;
    }
  }
  std::cout << "No Input Manager " << managername << " registered" << std::endl;
  return -1;
}

int Fun4AllSyncManager::BranchSelect(const std::string &branch, const int iflag)
{
  int iret = 0;
  for (Fun4AllInputManager *inman : m_InManager)
  {
    iret += inman->BranchSelect(branch, iflag);
  }
  return iret;
}

int Fun4AllSyncManager::setBranches(const std::string &managername)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (managername == inman->Name())
    {
      int iret = inman->setBranches();
      return iret;
    }
  }
  std::cout << "No Input Manager " << managername << " registered" << std::endl;
  return -1;
}

int Fun4AllSyncManager::setBranches()
{
  int iret = 0;
  for (Fun4AllInputManager *inman : m_InManager)
  {
    iret += inman->setBranches();
  }
  return iret;
}

int Fun4AllSyncManager::fileclose(const std::string &managername)
{
  int foundIt = 0;
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (managername == inman->Name() || managername.empty())
    {
      inman->fileclose();
      foundIt = 1;
    }
  }
  if (foundIt)
  {
    return 0;
  }
  std::cout << "No Input Manager " << managername << " registered" << std::endl;
  return -1;
}

void Fun4AllSyncManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "INPUTMANAGER")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of InputManagers in Fun4AllSyncManager "
              << Name() << ":" << std::endl;

    for (Fun4AllInputManager *inman : m_InManager)
    {
      std::cout << inman->Name() << std::endl;
    }
    std::cout << std::endl;
  }
  return;
}

int Fun4AllSyncManager::CheckSync(const unsigned i)
{
  int iret;
  iret = m_InManager[i]->SyncIt(m_MasterSync);
  return iret;
}

void Fun4AllSyncManager::GetInputFullFileList(std::vector<std::string> &fnames) const
{
  for (Fun4AllInputManager *InMan : m_InManager)
  {
    std::pair<std::list<std::string>::const_iterator, std::list<std::string>::const_iterator> beginend = InMan->FileOpenListBeginEnd();
    for (auto iter = beginend.first; iter != beginend.second; ++iter)
    {
      fnames.push_back(*iter);
    }
  }
  return;
}

void Fun4AllSyncManager::PushBackInputMgrsEvents(const int i)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    inman->PushBackEvents(i);
  }
  return;
}

int Fun4AllSyncManager::ResetEvent()
{
  int iret = 0;
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Resetting Event for Input Manager " << inman->Name() << std::endl;
    }
    iret += inman->ResetEvent();
  }
  return iret;
}

void Fun4AllSyncManager::CurrentEvent(const int evt)
{
  m_CurrentEvent = evt;
  Fun4AllServer *se = Fun4AllServer::instance();
  se->EventNumber(evt);
  return;
}

void Fun4AllSyncManager::PrintSyncProblem() const
{
  std::cout << "Bad use of Fun4AllDstInputManager for file(s) which do not have a synchronization object" << std::endl;
  std::cout << "This works for single streams but if you run with multiple input streams this might lead to event mixing" << std::endl;
  std::cout << "If you insist to run this (you take full responsibility), change the following in your macro: " << std::endl;
  for (auto iter : m_InManager)
  {
    if (iter->HasSyncObject() < 0)
    {
      std::cout << "File " << iter->FileName() << " does not contain a sync object" << std::endl;
      std::cout << "Change its Fun4AllDstInputManager with name " << iter->Name() << " from Fun4AllDstInputManager to Fun4AllNoSyncDstInputManager" << std::endl;
    }
  }
  return;
}
