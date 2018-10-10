#include "Fun4AllSyncManager.h"
#include "Fun4AllInputManager.h"
#include "Fun4AllHistoBinDefs.h"
#include "Fun4AllServer.h"
#include <phool/recoConsts.h>

#include <ffaobjects/RunHeader.h>
#include <ffaobjects/SyncObject.h>

#include <phool/PHIODataNode.h>
#include <phool/PHTypedNodeIterator.h>

#include <boost/foreach.hpp>

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace std;

Fun4AllSyncManager::Fun4AllSyncManager(const string &name):
  Fun4AllBase(name),
  prdf_segment(0),
  prdf_events(0),
  events_total(0),
  currentrun(0),
  currentevent(0),
  repeat(0),
  MasterSync(NULL)
{
  return;
}

Fun4AllSyncManager::~Fun4AllSyncManager()
{
  delete MasterSync;
  while (InManager.begin() != InManager.end())
    {
      if (Verbosity())
        {
          InManager.back()->Verbosity(Verbosity());
        }
      delete InManager.back();
      InManager.pop_back();
    }
  return ;
}

int
Fun4AllSyncManager::registerInputManager(Fun4AllInputManager *InputManager)
{
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
    {
      if ( inman->Name() == InputManager->Name() )
	{
          cout << "InputManager " << InputManager->Name() << " allready in list" << endl;
          return -1;
        }
    }

  if (Verbosity() > 0)
    {
      cout << "Registering InputManager " << InputManager->Name() << endl;
    }
  InManager.push_back(InputManager);
  iretInManager.push_back(0);
  InputManager->setSyncManager(this);
  return 0;
}

Fun4AllInputManager *
Fun4AllSyncManager::getInputManager(const string &name)
{
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
    {
      if (name == inman->Name())
	{
	  return inman;
	}
    }
  cout << Name() << ": Could not find InputManager" << name << endl;
  return NULL;
}

//_________________________________________________________________
int
Fun4AllSyncManager::run(const int nevnts)
{
  vector<Fun4AllInputManager *>::iterator iter;
  int iret = 0;
  int icnt = 0;
  int iretsync = 0;
  int resetnodetree = 0;
  while (! iret)
    {
      unsigned iman = 0;
      int ifirst = 0;
      for (iter = InManager.begin(); iter != InManager.end(); ++iter)
        {
          iretInManager[iman] = (*iter)->run(1);
          iret += iretInManager[iman];
          if (!ifirst)
            {
              if (!iretInManager[iman])
                {
                  if (!((*iter)->GetSyncObject(&MasterSync))) // NoSync managers return non zero
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
          resetnodetree = 1;

          // if there was an io error (file exhausted) we nee to push back
          // the events from files which are not exhausted yet into the root files
          // here we check the return codes for each input manager and if the
          // read was successful (iret = 0) we push the event back
          if (iret)
            {
              vector<Fun4AllInputManager *>::const_iterator InIter;
              if (repeat)
                {
                  for (InIter = InManager.begin(); InIter != InManager.end(); ++InIter)
                    {
                      if ((*InIter)->isOpen())
                        {
                          (*InIter)->fileclose();
                        }
                      int ireset = (*InIter)->ResetFileList();
                      if (ireset)
                        {

                          cout << "Resetting input manager " << (*InIter)->Name() << " failed during Repeat" << endl;
                          exit(1);
                        }
                    }
                  if (repeat > 0)
                    {
                      repeat--;
                    }
                  iret = 0;
                  continue;
                }
              vector<int>::const_iterator iter;
              // push back events where the Imanager did not report an error
              InIter = InManager.begin();
              for (iter = iretInManager.begin(); iter != iretInManager.end(); ++iter)
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
                          cout << (*InIter)->Name() << ": push evts: " << *iter  << endl;
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
                  InManager[nman]->NoSyncPushBackEvents(1);
                }
              continue;
            }
        }

      events_total++;
      if (nevnts > 0 && ++icnt >= nevnts)
        {
          break;
        }
    }

 readerror:
  if (!iret)
    {
      if (!resetnodetree) // all syncing is done and no read errors --> we have a good event in memory
	{
	  currentrun = 0; // reset current run to 0
	  for (iter = InManager.begin(); iter != InManager.end(); ++iter)
	    {
	      int runno = (*iter)->RunNumber();
	      if (Verbosity() > 2)
		{
		  cout << Name() << " input mgr " << (*iter)->Name() << " run: " << runno << endl;
		}
	      if (runno != 0)
		{
		  if (currentrun == 0)
		    {
		      currentrun = runno;
		      continue;

		    }
		  else
		    {
		      if (currentrun != runno)
			{
			  cout << "Mixing run numbers (except runnumber=0 which means no valid runnumber) is not supported" << endl;
			  cout << "Here are the list of input managers and runnumbers:" << endl;
			  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
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
int
Fun4AllSyncManager::skip(const int nevnts)
{
  if (!InManager.empty())
    {
      int Npushback = -nevnts;
      // the PushBackEvents(const int nevents) "pushes back" events events into the root file
      // (technically it just decrements the local counter in the PHNodeIOManager)
      // giving it a negative argument will skip events
      // this is much faster than actually reading the events in
      int iret = InManager[0]->PushBackEvents(Npushback);
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
int
Fun4AllSyncManager::fileopen(const string &managername, const string &filename)
{
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
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

int
Fun4AllSyncManager::BranchSelect(const string &managername, const string &branch, const int iflag)
{
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
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

int
Fun4AllSyncManager::BranchSelect(const string &branch, const int iflag)
{
  int iret = 0;
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
    {
      iret += inman->BranchSelect(branch, iflag);
    }
  return iret;
}

int
Fun4AllSyncManager::setBranches(const string &managername)
{
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
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

int
Fun4AllSyncManager::setBranches()
{
  int iret = 0;
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
    {
      iret += inman->setBranches();
    }
  return iret;
}

int Fun4AllSyncManager::fileclose(const string &managername)
{
  int foundIt = 0;
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
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
      cout << "--------------------------------------" << endl << endl;
      cout << "List of InputManagers in Fun4AllSyncManager "
	   << Name() << ":" << endl;

      BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
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
  iret = InManager[i]->SyncIt(MasterSync);
  return iret;
}

void
Fun4AllSyncManager::GetInputFullFileList(std::vector<std::string> &fnames) const
{
  list<string>::const_iterator listiter;
  vector<Fun4AllInputManager *>::const_iterator iter;
  for (iter = InManager.begin(); iter != InManager.end(); ++iter)
    {
      list<string> fl = (*iter)->GetFileOpenedList();
      for (listiter = fl.begin(); listiter != fl.end(); ++listiter)
        {
          fnames.push_back(*listiter);
        }
    }
  return;
}

void
Fun4AllSyncManager::PushBackInputMgrsEvents(const int i)
{
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
    {
      inman->PushBackEvents(i);
    }
    return;
}

int
Fun4AllSyncManager::ResetEvent()
{
  int iret = 0;
  BOOST_FOREACH(Fun4AllInputManager * inman, InManager)
    {
      if (Verbosity() > 0)
        {
          cout << "Resetting Event for Input Manager " << inman->Name() << endl;
        }
      iret += inman->ResetEvent();
    }
  return iret;
}

void
Fun4AllSyncManager::CurrentEvent(const int evt)
{
  currentevent = evt;
  Fun4AllServer *se = Fun4AllServer::instance();
  se->EventNumber(evt);
  return;
}
