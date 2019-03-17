#include "Fun4AllInputManager.h"

#include "Fun4AllServer.h"
#include "SubsysReco.h"

#include <phool/phool.h>

#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>

using namespace std;

Fun4AllInputManager::Fun4AllInputManager(const string &name, const string &nodename, const string &topnodename)
  : Fun4AllBase(name)
  , m_MySyncManager(nullptr)
  , m_IsOpen(0)
  , m_InputNode(nodename)
  , m_TopNodeName(topnodename)
  , repeat(0)
  , myrunnumber(0)
  , initrun(0)
{
  return;
}

Fun4AllInputManager::~Fun4AllInputManager()
{
  while (m_SubsystemsVector.begin() != m_SubsystemsVector.end())
  {
    if (Verbosity())
    {
      m_SubsystemsVector.back()->Verbosity(Verbosity());
    }
    delete m_SubsystemsVector.back();
    m_SubsystemsVector.pop_back();
  }
}

int Fun4AllInputManager::AddFile(const string &filename)
{
  if (Verbosity() > 0)
  {
    cout << "Adding " << filename << " to list of input files for "
         << Name() << endl;
  }
  m_FileList.push_back(filename);
  m_FileListCopy.push_back(filename);
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

    for (string file: m_FileList)
    {
      cout << file << endl;
    }
  }
  if (what == "ALL" || what == "SUBSYSTEMS")
  {
    // loop over the map and print out the content (name and location in memory)
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of SubsysRecos in Fun4AllInputManager " << Name() << ":" << endl;

      for (SubsysReco *subsys: m_SubsystemsVector)
    {
      cout << subsys->Name() << endl;
    }
    cout << endl;
  }
  return;
}

int Fun4AllInputManager::registerSubsystem(SubsysReco *subsystem)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  int iret = subsystem->Init(se->topNode(m_TopNodeName));
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
  m_SubsystemsVector.push_back(subsystem);
  return 0;
}

int Fun4AllInputManager::RejectEvent()
{
  if (!m_SubsystemsVector.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    for (SubsysReco *subsys:  m_SubsystemsVector)
    {
      if (!initrun)
      {
        subsys->InitRun(se->topNode(m_TopNodeName));
        initrun = 1;
      }
      if (Verbosity() > 0)
      {
        cout << Name() << ": Fun4AllInpuManager::EventReject processing " << subsys->Name() << endl;
      }
      if (subsys->process_event(se->topNode(m_TopNodeName)) != Fun4AllReturnCodes::EVENT_OK)
      {
        return Fun4AllReturnCodes::DISCARDEVENT;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int Fun4AllInputManager::ResetFileList()
{
  if (m_FileListCopy.empty())
  {
    cout << Name() << ": ResetFileList can only be used with filelists" << endl;
    return -1;
  }
  m_FileList.clear();
  m_FileList = m_FileListCopy;
  return 0;
}

void Fun4AllInputManager::UpdateFileList()
{
  if (!m_FileList.empty())
  {
    if (repeat)
    {
      m_FileList.push_back(*(m_FileList.begin()));
      if (repeat > 0)
      {
        repeat--;
      }
    }
    m_FileList.pop_front();
  }
  return;
}

int Fun4AllInputManager::OpenNextFile()
{
  while (!m_FileList.empty())
    {
      list<string>::const_iterator iter = m_FileList.begin();
      if (Verbosity())
        {
          cout << PHWHERE << " opening next file: " << *iter << endl;
        }
      if (fileopen(*iter))
        {
          cout << PHWHERE << " could not open file: " << *iter << endl;
          m_FileList.pop_front();
        }
      else
        {
          return 0;
        }

    }
  return -1;
}
