#include "Fun4AllSyncManager.h"

#include "Fun4AllInputManager.h"
#include "Fun4AllReturnCodes.h"  // for EVENT_OK, RESET_NODE_TREE
#include "Fun4AllServer.h"

#include <ffaobjects/SyncObject.h>

#include <phool/phool.h>  // for PHWHERE

#include <cstdlib>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <list>      // for list<>::const_iterator, _List_con...
#include <string>
#include <utility>  // for pair
#include <vector>

using namespace std;

Fun4AllSyncManager::Fun4AllSyncManager(const string &name)
  : Fun4AllBase(name)
  , m_PrdfSegment(0)
  , m_PrdfEvents(0)
  , m_EventsTotal(0)
  , m_CurrentRun(0)
  , m_CurrentEvent(0)
  , m_Repeat(0)
  , m_MasterSync(nullptr)
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
      cout << "InputManager " << InputManager->Name() << " allready in list" << endl;
      return -1;
    }
  }

  if (Verbosity() > 0)
  {
    cout << "Registering InputManager " << InputManager->Name() << endl;
  }
  m_InManager.push_back(InputManager);
  m_iretInManager.push_back(0);
  InputManager->setSyncManager(this);
  return 0;
}

Fun4AllInputManager *
Fun4AllSyncManager::getInputManager(const string &name)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (name == inman->Name())
    {
      return inman;
    }
  }
  cout << Name() << ": Could not find InputManager" << name << endl;
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
    for (vector<Fun4AllInputManager *>::iterator iter = m_InManager.begin(); iter != m_InManager.end(); ++iter)
    {
      m_iretInManager[iman] = (*iter)->run(1);
      iret += m_iretInManager[iman];
      if (!ifirst)
      {
        if (!m_iretInManager[iman])
        {
          if (!((*iter)->GetSyncObject(&m_MasterSync)))  // NoSync managers return non zero
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
        vector<Fun4AllInputManager *>::const_iterator InIter;
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
                cout << "Resetting input manager " << (*InIter)->Name() << " failed during Repeat" << endl;
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
        for (vector<int>::const_iterator iter = m_iretInManager.begin(); iter != m_iretInManager.end(); ++iter)
        {
          if (Verbosity() > 0)
          {
            cout << (*InIter)->Name() << ": return code: " << *iter << endl;
          }
          if (!(*iter))
          {
            (*InIter)->PushBackEvents(1);
            if (Verbosity() > 0)
            {
              cout << (*InIter)->Name() << ": push evts: " << *iter << endl;
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
      for (vector<Fun4AllInputManager *>::iterator iter = m_InManager.begin(); iter != m_InManager.end(); ++iter)
      {
        int runno = (*iter)->RunNumber();
        if (Verbosity() > 2)
        {
          cout << Name() << " input mgr " << (*iter)->Name() << " run: " << runno << endl;
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
              cout << "Mixing run numbers (except runnumber=0 which means no valid runnumber) is not supported" << endl;
              cout << "Here are the list of input managers and runnumbers:" << endl;
              for (Fun4AllInputManager *inman : m_InManager)
              {
                cout << inman->Name() << " runno: " << inman->RunNumber() << endl;
              }
              cout << "Exiting now" << endl;
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
      cout << PHWHERE << " Error during skipping events" << endl;
      return iret;
    }
  }
  cout << PHWHERE << " Cannot skip events: No Input Managers registered?" << endl;
  Print("INPUTMANAGER");
  cout << "If there are Input Managers in this list, send mail with this" << endl;
  cout << "error message to off-l" << endl;
  cout << "and include the macro you used" << endl;
  return -1;
}

//_________________________________________________________________
int Fun4AllSyncManager::fileopen(const string &managername, const string &filename)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (managername == inman->Name())
    {
      int iret = inman->fileopen(filename);
      return iret;
    }
  }
  cout << "No Input Manager " << managername << " registered" << endl;
  return -1;
}

int Fun4AllSyncManager::BranchSelect(const string &managername, const string &branch, const int iflag)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (managername == inman->Name())
    {
      int iret = inman->BranchSelect(branch, iflag);
      return iret;
    }
  }
  cout << "No Input Manager " << managername << " registered" << endl;
  return -1;
}

int Fun4AllSyncManager::BranchSelect(const string &branch, const int iflag)
{
  int iret = 0;
  for (Fun4AllInputManager *inman : m_InManager)
  {
    iret += inman->BranchSelect(branch, iflag);
  }
  return iret;
}

int Fun4AllSyncManager::setBranches(const string &managername)
{
  for (Fun4AllInputManager *inman : m_InManager)
  {
    if (managername == inman->Name())
    {
      int iret = inman->setBranches();
      return iret;
    }
  }
  cout << "No Input Manager " << managername << " registered" << endl;
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

int Fun4AllSyncManager::fileclose(const string &managername)
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
  cout << "No Input Manager " << managername << " registered" << endl;
  return -1;
}

void Fun4AllSyncManager::Print(const string &what) const
{
  if (what == "ALL" || what == "INPUTMANAGER")
  {
    // loop over the map and print out the content (name and location in memory)
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of InputManagers in Fun4AllSyncManager "
         << Name() << ":" << endl;

    for (Fun4AllInputManager *inman : m_InManager)
    {
      cout << inman->Name() << endl;
    }
    cout << endl;
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
      cout << "Resetting Event for Input Manager " << inman->Name() << endl;
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
