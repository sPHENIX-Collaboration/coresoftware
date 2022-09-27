#include "Fun4AllHistoManager.h"

#include "TDirectoryHelper.h"

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>
#include <TTree.h>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 20, 0)
#define HAS_THNSPARSE 1
#include <THnSparse.h>
#endif

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
    delete Histo.begin()->second;
    Histo.erase(Histo.begin());
  }
  return;
}

int Fun4AllHistoManager::dumpHistos(const std::string &filename, const std::string &openmode)
{
  int iret = 0;
  if (!filename.empty())
  {
    outfilename = filename;
  }
  else
  {
    if (outfilename.empty())
    {
      recoConsts *rc = recoConsts::instance();
      std::ostringstream filnam;
      int runnumber = -1;
      if (rc->FlagExist("RUNNUMBER"))
      {
        runnumber = rc->get_IntFlag("RUNNUMBER");
      }
      // this will set the filename to the name of the manager
      // add the runnumber in the std 10 digit format and
      // end it with a .root extension
      filnam << Name() << "-"
             << std::setfill('0') << std::setw(10)
             << runnumber << ".root";
      outfilename = filnam.str();
    }
  }
  std::cout << "Fun4AllHistoManager::dumpHistos() Writing root file: " << outfilename << std::endl;

  const int compress = 9;
  std::ostringstream creator;
  creator << "Created by " << Name();
  TFile hfile(outfilename.c_str(), openmode.c_str(), creator.str().c_str(), compress);
  if (!hfile.IsOpen())
  {
    std::cout << PHWHERE << " Could not open output file" << outfilename << std::endl;
    return -1;
  }

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
      if (dirname.size())
      {
        std::cout << " being saved to directory " << dirname;
      }
      std::cout << std::endl;
    }

    if (dirname.size())
    {
      TDirectoryHelper::mkdir(&hfile, dirname.c_str());
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
    static_cast<TTree *>(h1d)->SetDirectory(nullptr);

  // For histograms, enforce error calculation and propagation
  if (h1d->InheritsFrom("TH1"))
    static_cast<TH1 *>(h1d)->Sumw2();

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
  else
  {
    std::cout << "Fun4AllHistoManager::getHisto: ERROR Invalid histogram number: "
              << ihisto << ", maximum number is " << size << std::endl;
  }
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
  else
  {
    std::cout << "Fun4AllHistoManager::getHisto: ERROR Invalid histogram number: "
              << ihisto << ", maximum number is " << size << std::endl;
  }
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

    std::map<const std::string, TNamed *>::const_iterator hiter;
    for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
    {
      std::cout << hiter->first << " is " << hiter->second << std::endl;
    }
    std::cout << std::endl;
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
      (dynamic_cast<TH1 *>(h))->Reset();
#if HAS_THNSPARSE
    else if (h->InheritsFrom("THnSparse"))
      (dynamic_cast<THnSparse *>(h))->Reset();
#endif
  }
  return;
}
