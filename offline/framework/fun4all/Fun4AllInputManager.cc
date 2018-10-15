#include "Fun4AllInputManager.h"
#include <phool/phool.h>
#include "Fun4AllServer.h"
#include "SubsysReco.h"

#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>

using namespace std;

Fun4AllInputManager::Fun4AllInputManager(const string &name, const string &nodename, const string &topnodename)
  : Fun4AllBase(name)
  , m_InputNode(nodename)
  , topNodeName(topnodename)
  , mySyncManager(NULL)
  , repeat(0)
  , myrunnumber(0)
  , initrun(0)
{
  return;
}

Fun4AllInputManager::~Fun4AllInputManager()
{
  while (Subsystems.begin() != Subsystems.end())
  {
    if (Verbosity())
    {
      Subsystems.back()->Verbosity(Verbosity());
    }
    delete Subsystems.back();
    Subsystems.pop_back();
  }
}

int Fun4AllInputManager::AddFile(const string &filename)
{
  if (Verbosity() > 0)
  {
    cout << "Adding " << filename << " to list of input files for "
         << Name() << endl;
  }
  filelist.push_back(filename);
  filelist_copy.push_back(filename);
  return 0;
}

int Fun4AllInputManager::AddListFile(const string &filename, const int do_it)
{
  // checking filesize to see if we have a text file
  if (boost::filesystem::exists(filename.c_str()))
  {
    if (boost::filesystem::is_regular_file(filename.c_str()))
    {
      uintmax_t fsize = boost::filesystem::file_size(filename.c_str());
      if (fsize > 1000000 && !do_it)
      {
        cout << "size of " << filename
             << " is suspiciously large for a text file: "
             << fsize << " bytes" << endl;
        cout << "if you really want to use " << filename
             << " as list file (it will be used as a text file containing a list of input files), use AddListFile(\""
             << filename << "\",1)" << endl;
        return -1;
      }
    }
    else
    {
      cout << filename << " is not a regular file" << endl;
      return -1;
    }
  }
  else
  {
    cout << PHWHERE << "Could not open " << filename << endl;
    return -1;
  }
  ifstream infile;
  infile.open(filename.c_str(), ios_base::in);
  if (!infile)
  {
    cout << PHWHERE << "Could not open " << filename << endl;
    return -1;
  }
  string FullLine;
  getline(infile, FullLine);
  while (!infile.eof())
  {
    if (FullLine.size() && FullLine[0] != '#')  // remove comments
    {
      AddFile(FullLine);
    }
    else if (FullLine.size())
    {
      if (Verbosity() > 0)
      {
        cout << "Found Comment: " << FullLine << endl;
      }
    }

    getline(infile, FullLine);
  }
  infile.close();
  return 0;
}

void Fun4AllInputManager::Print(const string &what) const
{
  if (what == "ALL" || what == "FILELIST")
  {
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of input files in Fun4AllInputManager " << Name() << ":" << endl;

    list<string>::const_iterator iter;
    for (iter = filelist.begin(); iter != filelist.end(); ++iter)
    {
      cout << *iter << endl;
    }
  }
  if (what == "ALL" || what == "SUBSYSTEMS")
  {
    // loop over the map and print out the content (name and location in memory)
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of SubsysRecos in Fun4AllInputManager " << Name() << ":" << endl;

    vector<SubsysReco *>::const_iterator miter;
    for (miter = Subsystems.begin(); miter != Subsystems.end(); ++miter)
    {
      cout << (*miter)->Name() << endl;
    }
    cout << endl;
  }
  return;
}

int Fun4AllInputManager::registerSubsystem(SubsysReco *subsystem)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  int iret = subsystem->Init(se->topNode(topNodeName));
  if (iret)
  {
    cout << PHWHERE << " Error initializing subsystem "
         << subsystem->Name() << ", return code: " << iret << endl;
    return iret;
  }
  if (Verbosity() > 0)
  {
    cout << "Registering Subsystem " << subsystem->Name() << endl;
  }
  Subsystems.push_back(subsystem);
  return 0;
}

int Fun4AllInputManager::RejectEvent()
{
  if (!Subsystems.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    vector<SubsysReco *>::iterator iter;
    for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
    {
      if (!initrun)
      {
        (*iter)->InitRun(se->topNode(topNodeName));
        initrun = 1;
      }
      if (Verbosity() > 0)
      {
        cout << Name() << ": Fun4AllInpuManager::EventReject processing " << (*iter)->Name() << endl;
      }
      if ((*iter)->process_event(se->topNode(topNodeName)) != Fun4AllReturnCodes::EVENT_OK)
      {
        return Fun4AllReturnCodes::DISCARDEVENT;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int Fun4AllInputManager::ResetFileList()
{
  if (filelist_copy.empty())
  {
    cout << Name() << ": ResetFileList can only be used with filelists" << endl;
    return -1;
  }
  filelist.clear();
  filelist = filelist_copy;
  return 0;
}
