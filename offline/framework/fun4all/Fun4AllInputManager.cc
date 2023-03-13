#include "Fun4AllInputManager.h"

#include "Fun4AllServer.h"
#include "SubsysReco.h"

#include <phool/phool.h>

#include <boost/filesystem.hpp>

#include <cstdint>  // for uintmax_t
#include <fstream>
#include <iostream>

Fun4AllInputManager::Fun4AllInputManager(const std::string &name, const std::string &nodename, const std::string &topnodename)
  : Fun4AllBase(name)
  , m_InputNode(nodename)
  , m_TopNodeName(topnodename)
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

int Fun4AllInputManager::AddFile(const std::string &filename)
{
  if (Verbosity() > 0)
  {
    std::cout << "Adding " << filename << " to list of input files for "
              << Name() << std::endl;
  }
  m_FileList.push_back(filename);
  m_FileListCopy.push_back(filename);
  return 0;
}

int Fun4AllInputManager::AddListFile(const std::string &filename, const int do_it)
{
  // checking filesize to see if we have a text file
  if (boost::filesystem::exists(filename.c_str()))
  {
    if (boost::filesystem::is_regular_file(filename.c_str()))
    {
      uintmax_t fsize = boost::filesystem::file_size(filename.c_str());
      if (fsize > 1000000 && !do_it)
      {
        std::cout << "size of " << filename
                  << " is suspiciously large for a text file: "
                  << fsize << " bytes" << std::endl;
        std::cout << "if you really want to use " << filename
                  << " as list file (it will be used as a text file containing a list of input files), use AddListFile(\""
                  << filename << "\",1)" << std::endl;
        return -1;
      }
    }
    else
    {
      std::cout << filename << " is not a regular file" << std::endl;
      return -1;
    }
  }
  else
  {
    std::cout << PHWHERE << "Could not open " << filename << std::endl;
    return -1;
  }
  std::ifstream infile;
  infile.open(filename, std::ios_base::in);
  if (!infile)
  {
    std::cout << PHWHERE << "Could not open " << filename << std::endl;
    return -1;
  }
  std::string FullLine;
  int nfiles = 0;
  getline(infile, FullLine);
  while (!infile.eof())
  {
    if (FullLine.size() && FullLine[0] != '#')  // remove comments
    {
      AddFile(FullLine);
      nfiles++;
    }
    else if (FullLine.size())
    {
      if (Verbosity() > 0)
      {
        std::cout << "Found Comment: " << FullLine << std::endl;
      }
    }
    getline(infile, FullLine);
  }
  infile.close();
  if (nfiles == 0)
  {
    std::cout << Name() << " listfile " << filename << " does not contain filenames "
              << "if this is the only list you load into this Input Manager your code will exit very soon" << std::endl;
  }
  return 0;
}

void Fun4AllInputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "FILELIST")
  {
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of input files in Fun4AllInputManager " << Name() << ":" << std::endl;

    for (const std::string &file : m_FileList)
    {
      std::cout << file << std::endl;
    }
  }
  if (what == "ALL" || what == "SUBSYSTEMS")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of SubsysRecos in Fun4AllInputManager " << Name() << ":" << std::endl;

    for (SubsysReco *subsys : m_SubsystemsVector)
    {
      std::cout << subsys->Name() << std::endl;
    }
    std::cout << std::endl;
  }
  return;
}

int Fun4AllInputManager::registerSubsystem(SubsysReco *subsystem)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  int iret = subsystem->Init(se->topNode(m_TopNodeName));
  if (iret)
  {
    std::cout << PHWHERE << " Error initializing subsystem "
              << subsystem->Name() << ", return code: " << iret << std::endl;
    return iret;
  }
  if (Verbosity() > 0)
  {
    std::cout << "Registering Subsystem " << subsystem->Name() << std::endl;
  }
  m_SubsystemsVector.push_back(subsystem);
  return 0;
}

int Fun4AllInputManager::RejectEvent()
{
  if (!m_SubsystemsVector.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    for (SubsysReco *subsys : m_SubsystemsVector)
    {
      if (!m_InitRun)
      {
        subsys->InitRun(se->topNode(m_TopNodeName));
        m_InitRun = 1;
      }
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": Fun4AllInpuManager::EventReject processing " << subsys->Name() << std::endl;
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
    std::cout << Name() << ": ResetFileList can only be used with filelists" << std::endl;
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
    if (m_Repeat)
    {
      m_FileList.push_back(*(m_FileList.begin()));
      if (m_Repeat > 0)
      {
        m_Repeat--;
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
    std::list<std::string>::const_iterator iter = m_FileList.begin();
    if (Verbosity())
    {
      std::cout << PHWHERE << " opening next file: " << *iter << std::endl;
    }
    if (fileopen(*iter))
    {
      std::cout << PHWHERE << " could not open file: " << *iter << std::endl;
      m_FileList.pop_front();
    }
    else
    {
      return 0;
    }
  }
  return -1;
}
