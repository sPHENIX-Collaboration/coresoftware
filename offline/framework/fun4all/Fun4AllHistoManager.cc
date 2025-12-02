#include "Fun4AllHistoManager.h"

#include "Fun4AllOutputManager.h"
#include "TDirectoryHelper.h"

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TFile.h>
#include <TH1.h>
#include <THnBase.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TTree.h>

#include <filesystem>
#include <format>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>  // for pair

Fun4AllHistoManager::Fun4AllHistoManager(const std::string &name)
  : Fun4AllBase(name)
{
  return;
}

Fun4AllHistoManager::~Fun4AllHistoManager()
{
  while (Histo.begin() != Histo.end())
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << "deleting " << Histo.begin()->first
		<< ", histo name " << Histo.begin()->second->GetName()
		<< std::endl;
    }
    delete Histo.begin()->second;
    Histo.erase(Histo.begin());
  }
  return;
}

int Fun4AllHistoManager::RunAfterClosing()
{
  unsigned int iret = 0;

  if (!m_RunAfterClosingScript.empty())
  {
    if (!std::filesystem::exists(m_RunAfterClosingScript))
    {
      std::cout << PHWHERE << "RunAfterClosing() closing script " << m_RunAfterClosingScript << " not found" << std::endl;
      return -1;
    }
    if (!((std::filesystem::status(m_RunAfterClosingScript).permissions() & std::filesystem::perms::owner_exec) == std::filesystem::perms::owner_exec))
    {
      std::cout << PHWHERE << "RunAfterClosing() closing script " << m_RunAfterClosingScript << " is not owner executable" << std::endl;
      return -1;
    }
    std::string fullcmd = m_RunAfterClosingScript + " " + m_LastClosedFileName + " " + m_ClosingArgs;
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << " running " << fullcmd << std::endl;
    }
    iret = gSystem->Exec(fullcmd.c_str());
  }
  if (iret)
  {
    iret = iret >> 8U;
  }
  return iret;
}

int Fun4AllHistoManager::dumpHistos(const std::string &filename, const std::string &openmode)
{
  int iret = 0;
  if (!filename.empty())
  {
    m_outfilename = filename;
  }
  else
  {
    if (m_outfilename.empty())
    {
      recoConsts *rc = recoConsts::instance();
      int runnumber = 0;
      if (rc->FlagExist("RUNNUMBER"))
      {
        runnumber = rc->get_IntFlag("RUNNUMBER");
      }
      // this will set the filename to the name of the manager
      // add the runnumber in the std 10 digit format and
      // end it with a .root extension
      m_outfilename = Name() + std::format("-{:08}.root", runnumber);
    }
  }
  std::string theoutfile = m_outfilename;
  if (ApplyFileRule())
  {
    
  std::filesystem::path p = m_outfilename;

  recoConsts *rc = recoConsts::instance();
  int runnumber = 0;
  if (rc->FlagExist("RUNNUMBER"))
  {
    runnumber = rc->get_IntFlag("RUNNUMBER");
  }
    std::string fullpath = ".";
    if (p.has_parent_path())
    {
      fullpath = p.parent_path();
    }
    std::string runseg = std::format("-{:08}-{:05}",runnumber,m_CurrentSegment);
    theoutfile = fullpath + std::string("/") + std::string(p.stem()) + runseg + std::string(p.extension());
    m_CurrentSegment++;
  }
  m_LastClosedFileName = theoutfile;
  std::filesystem::path pout = theoutfile;
  theoutfile = theoutfile + std::string("?reproducible=") + std::string(pout.filename());
  std::cout << "Fun4AllHistoManager::dump() Writing root file: " <<  m_LastClosedFileName << std::endl;

  const int compress = 505;
  std::string creator = "Created by " + Name();
  TFile hfile(theoutfile.c_str(), openmode.c_str(), creator.c_str());
  if (!hfile.IsOpen())
  {
    std::cout << PHWHERE << " Could not open output file" << theoutfile.c_str() << std::endl;
    return -1;
  }
  hfile.SetCompressionSettings(compress);
  std::map<const std::string, TNamed *>::const_iterator hiter;
  for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
  {
    const std::string &hname = hiter->first;
    const TNamed *hptr = hiter->second;
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << " Saving histo "
                << hname
                << std::endl;
    }

    //  Decode the string to see if it wants a directory
    std::string::size_type pos = hname.find_last_of('/');
    std::string dirname;
    if (pos != std::string::npos)  // string::npos is the result if search unsuccessful
    {
      dirname = hname.substr(0, pos);
    }
    else
    {
      dirname = "";
    }

    if (Verbosity())
    {
      std::cout << " Histogram named " << hptr->GetName();
      std::cout << " key " << hname;
      if (!dirname.empty())
      {
        std::cout << " being saved to directory " << dirname;
      }
      std::cout << std::endl;
    }

    if (!dirname.empty())
    {
      TDirectoryHelper::mkdir(&hfile, dirname);
      hfile.cd(dirname.c_str());
    }

    if (hptr)
    {
      int byteswritten = hptr->Write();
      if (!byteswritten)
      {
        std::cout << PHWHERE << "Error saving histogram "
                  << hptr->GetName()
                  << std::endl;
        iret = -2;
      }
    }
    else
    {
      std::cout << PHWHERE << "dumpHistos : histogram "
                << hname << " is a null pointer! Won't be saved."
                << std::endl;
    }
  }
  hfile.Close();
  RunAfterClosing();
  return iret;
}

bool Fun4AllHistoManager::registerHisto(TNamed *h1d, const int replace)
{
  return registerHisto(h1d->GetName(), h1d, replace);
}

bool Fun4AllHistoManager::registerHisto(const std::string &hname, TNamed *h1d, const int replace)
{
  std::map<const std::string, TNamed *>::const_iterator histoiter = Histo.find(hname);
  if (histoiter != Histo.end() && replace == 0)
  {
    std::cout << "Histogram " << hname << " already registered, I won't overwrite it" << std::endl;
    std::cout << "Use a different name and try again" << std::endl;
    return false;
  }

  std::string::size_type pos = hname.find_last_of('/');
  std::string histoname = hname;
  if (pos != std::string::npos)  // okay someone wants damn TDirectories
  {
    histoname = hname.substr(pos + 1);
  }
  if (Verbosity() > 1)
  {
    if (histoname != h1d->GetName())
    {
      std::cout << PHWHERE << "Histogram " << h1d->GetName()
                << " at " << h1d << " renamed to " << histoname << std::endl;
    }
  }
  // this one did some very ugly mutilation to a const char *
  // using a string seems to avoid the damage
  h1d->SetName(histoname.c_str());
  Histo[hname] = h1d;

  // reset directory for TTree
  if (h1d->InheritsFrom("TTree"))
  {
    static_cast<TTree *>(h1d)->SetDirectory(nullptr); // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  }

  // For histograms, enforce error calculation and propagation
  if (h1d->InheritsFrom("TH1"))
  {
    static_cast<TH1 *>(h1d)->Sumw2();// NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  }

  return true;
}

int Fun4AllHistoManager::isHistoRegistered(const std::string &name) const
{
  std::map<const std::string, TNamed *>::const_iterator histoiter = Histo.find(name);
  if (histoiter != Histo.end())
  {
    return 1;
  }
  return 0;
}

TNamed *
Fun4AllHistoManager::getHisto(const unsigned int ihisto) const
{
  std::map<const std::string, TNamed *>::const_iterator histoiter = Histo.begin();
  unsigned int size = Histo.size();
  if (Verbosity() > 3)
  {
    std::cout << "Map contains " << size << " Elements" << std::endl;
  }
  if (ihisto < size)
  {
    for (unsigned int i = 0; i < ihisto; i++)
    {
      ++histoiter;
    }
    return histoiter->second;
  }

  std::cout << "Fun4AllHistoManager::getHisto: ERROR Invalid histogram number: "
            << ihisto << ", maximum number is " << size << std::endl;

  return nullptr;
}

std::string
Fun4AllHistoManager::getHistoName(const unsigned int ihisto) const
{
  std::map<const std::string, TNamed *>::const_iterator histoiter = Histo.begin();
  unsigned int size = Histo.size();
  if (Verbosity() > 3)
  {
    std::cout << "Map contains " << size << " Elements" << std::endl;
  }
  if (ihisto < size)
  {
    for (unsigned int i = 0; i < ihisto; i++)
    {
      ++histoiter;
    }
    return histoiter->first;
  }

  std::cout << "Fun4AllHistoManager::getHisto: ERROR Invalid histogram number: "
            << ihisto << ", maximum number is " << size << std::endl;

  return "";
}

TNamed *
Fun4AllHistoManager::getHisto(const std::string &hname) const
{
  std::map<const std::string, TNamed *>::const_iterator histoiter = Histo.find(hname);
  if (histoiter != Histo.end())
  {
    return histoiter->second;
  }
  std::cout << "Fun4AllHistoManager::getHisto: ERROR Unknown Histogram " << hname
            << ", The following are implemented: " << std::endl;
  Print("ALL");
  return nullptr;
}

void Fun4AllHistoManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "HISTOS")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Histos in Fun4AllHistoManager "
              << Name() << ":" << std::endl;

    for (const auto &hiter : Histo)
    {
      std::cout << hiter.first << " is " << hiter.second << std::endl;
    }
    std::cout << std::endl;
  }
  if (what == "ALL" || what == "CONFIG")
  {
    std::cout << Name() << std::endl;
    std::cout << "Output Filename: " << m_outfilename << std::endl;
    std::cout << "use filerule: " << m_UseFileRuleFlag << std::endl;
    std::cout << "Event Rollover at " << GetEventNumberRollover() << std::endl;
    std::cout << "current last event " << LastEventNumber() << std::endl;
    std::cout << "current segment number: " << m_CurrentSegment << std::endl;
  }
  return;
}

void Fun4AllHistoManager::Reset()
{
  std::map<const std::string, TNamed *>::const_iterator hiter;
  for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
  {
    TNamed *h = hiter->second;
    if (h->InheritsFrom("TH1"))
    {
      (dynamic_cast<TH1 *>(h))->Reset();
    }
    else if (h->InheritsFrom("THnBase"))
    {
      (dynamic_cast<THnBase *>(h))->Reset();
    }
  }
  return;
}

// this method figures out the last event number to be saved before rolling over
// an integer div of the current event by the number of events gives the first event we can expect
// in this process (this is not needed), then adding the number of events we want gives us the last event
// since want ranges like 1-99999, 100,000 - 199,999 we need to subtract 1
// from the calculated range.
// This is just run for the first event - later we just add the number of events to the last event number
void Fun4AllHistoManager::InitializeLastEvent(int eventnumber)
{
  if (GetEventNumberRollover() == 0 || m_LastEventInitializedFlag || eventnumber < 0)
  {
    return;
  }
  m_LastEventInitializedFlag = true;
  unsigned int firstevent = eventnumber / GetEventNumberRollover();
  unsigned int newlastevent = firstevent * GetEventNumberRollover() + GetEventNumberRollover() - 1;
  if (Verbosity() > 1)
  {
    std::cout << "event number: " << eventnumber << ", rollover: " << GetEventNumberRollover() << ", multiple: "
              << eventnumber / GetEventNumberRollover() << ", new last event number "
              << newlastevent << std::endl;
  }
  SetLastEventNumber(firstevent * GetEventNumberRollover() + GetEventNumberRollover() - 1);
  return;
}

void Fun4AllHistoManager::CopyRolloverSetting(const Fun4AllOutputManager *outman)
{
  SetEventNumberRollover(outman->GetEventNumberRollover());
  UseFileRule(outman->ApplyFileRule());
  StartSegment(outman->Segment());
  SetClosingScript(outman->GetClosingScript());
  SetClosingScriptArgs(outman->GetClosingScript());
}

bool Fun4AllHistoManager::isEmpty() const
{
  bool thisempty = true;
  for (const auto &hiter : Histo)
  {
    TH1 *h1 = dynamic_cast<TH1 *> (hiter.second);
    if (h1)
    {
      if (h1->GetEntries() > 0 || h1->Integral() > 0)
      {
	std::cout << hiter.first << " has " << h1->GetEntries() << " entries"
		  << " with integral " << h1->Integral() << std::endl;
	thisempty = false;
	break;
      }
    }
    else
    {
      THnBase *hn = dynamic_cast<THnBase *> (hiter.second);
      if (hn)
      {
	if (hn->GetEntries() > 0)
	{
	  std::cout << hiter.first << " has " << hn->GetEntries() << " entries" << std::endl;
	  thisempty = false;
	  break;
	}
      }
      else
      {
	std::cout << hiter.first << " is not a TH1 or THnBase" << std::endl;
      }
    }
  }
  return thisempty;
}
