#include "Fun4AllHistoManager.h"
#include "TDirectoryHelper.h"

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TFile.h>
#include <TH1.h>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,20,0)
#define HAS_THNSPARSE 1
#include <THnSparse.h>
#endif

#include <TNamed.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

Fun4AllHistoManager::Fun4AllHistoManager(const string &name): Fun4AllBase(name)
{
  return ;
}

Fun4AllHistoManager::~Fun4AllHistoManager()
{
  while(Histo.begin() != Histo.end())
    {
      delete Histo.begin()->second;
      Histo.erase(Histo.begin());
    }
  return ;
}

int
Fun4AllHistoManager::dumpHistos(const string &filename, const string &openmode)
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
          ostringstream filnam;
          int runnumber = -1;
          if (rc->FlagExist("RUNNUMBER"))
            {
              runnumber = rc->get_IntFlag("RUNNUMBER");
            }
          // this will set the filename to the name of the manager
          // add the runnumber in the std 10 digit format and
          // end it with a .root extension
          filnam << Name() << "-"
		 << setfill('0') << setw(10)
		 << runnumber << ".root";
          outfilename = filnam.str();
        }
    }
  cout << "Fun4AllHistoManager::dumpHistos() Writing root file: " << outfilename << endl;

  const int compress = 9;
  ostringstream creator;
  creator << "Created by " << Name();
  TFile hfile(outfilename.c_str(), openmode.c_str(), creator.str().c_str(), compress);
  if (!hfile.IsOpen())
    {
      cout << PHWHERE << " Could not open output file" << outfilename << endl;
      return -1;
    }

  map<const string, TNamed *>::const_iterator hiter;
  for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
    {
      const std::string & hname = hiter->first;
      const TNamed*       hptr  = hiter->second;
      if ( Verbosity() > 0 )
        {
          std::cout << PHWHERE << " Saving histo "
		    << hname
		    << std::endl;
        }

      //  Decode the string to see if it wants a directory
      string::size_type pos = hname.find_last_of('/');
      string dirname;
      if ( pos != string::npos ) // string::npos is the result if search unsuccessful
        {
          dirname = hname.substr(0, pos);
        }
      else
        {
          dirname = "";
        }

      if (Verbosity())
        {
          cout << " Histogram named " << hptr->GetName();
	  cout << " key " << hname;
	  if (dirname.size())
	    {
	      cout << " being saved to directory " << dirname;
	    }
	  cout << endl;
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
	      cout << PHWHERE << "Error saving histogram " 
		   << hptr->GetName()
		   << endl;
	      iret = -2;
	    }
        }
      else
        {
          cout << PHWHERE << "dumpHistos : histogram "
	       << hname << " is a null pointer! Won't be saved."
	       << std::endl;
        }
    }
  hfile.Close();
  return iret;
}

bool
Fun4AllHistoManager::registerHisto(TNamed *h1d, const int replace)
{
  return registerHisto(h1d->GetName(), h1d, replace);
}

bool
Fun4AllHistoManager::registerHisto(const string &hname, TNamed *h1d, const int replace)
{
  map<const string, TNamed *>::const_iterator histoiter = Histo.find(hname);
  if (histoiter != Histo.end() && replace == 0)
    {
      cerr << "Histogram " << hname << " already registered, I won't overwrite it" << endl;
      cerr << "Use a different name and try again" << endl;
      return false;
    }

  string::size_type pos = hname.find_last_of('/');
  string histoname = hname;
  if ( pos != string::npos ) // okay someone wants damn TDirectories
  {
    histoname = hname.substr(pos + 1);
  }
  if (Verbosity() > 1)
  {
    if (histoname != h1d->GetName())
    {
      cout << PHWHERE << "Histogram " << h1d->GetName()
        << " at " << h1d << " renamed to " << histoname << endl;
    }
  }
  // this one did some very ugly mutilation to a const char *
  // using a string seems to avoid the damage
  h1d->SetName(histoname.c_str());
  Histo[hname] = h1d;

  if (h1d->InheritsFrom("TTree"))
    static_cast<TTree*>(h1d)->SetDirectory(0);

  return true;
}

int
Fun4AllHistoManager::isHistoRegistered(const std::string &name) const
{

  map<const string, TNamed *>::const_iterator histoiter = Histo.find(name);
  if (histoiter != Histo.end())
    {
      return 1;
    }
  return 0;
}

TNamed *
Fun4AllHistoManager::getHisto(const unsigned int ihisto) const
{
  map<const string, TNamed *>::const_iterator histoiter = Histo.begin();
  unsigned int size = Histo.size();
  if (Verbosity() > 3)
    {
      cout << "Map contains " << size << " Elements" << endl;
    }
  if (ihisto < size)
    {
      for (unsigned int i = 0;i < ihisto;i++)
	{
	  ++histoiter;
	}
      return histoiter->second;
    }
  else
    {
      cout << "Fun4AllHistoManager::getHisto: ERROR Invalid histogram number: "
	   << ihisto << ", maximum number is " << size << endl;
    }
  return NULL;
}

const char *
Fun4AllHistoManager::getHistoName(const unsigned int ihisto) const
{
  map<const string, TNamed *>::const_iterator histoiter = Histo.begin();
  unsigned int size = Histo.size();
  if (Verbosity() > 3)
    {
      cout << "Map contains " << size << " Elements" << endl;
    }
  if (ihisto < size)
    {
      for (unsigned int i = 0;i < ihisto;i++)
	{
	  ++histoiter;
	}
      return histoiter->first.c_str();
    }
  else
    {
      cout << "Fun4AllHistoManager::getHisto: ERROR Invalid histogram number: "
	   << ihisto << ", maximum number is " << size << endl;
    }
  return NULL;
}

TNamed *
Fun4AllHistoManager::getHisto(const string &hname) const
{
  map<const string, TNamed *>::const_iterator histoiter = Histo.find(hname);
  if (histoiter != Histo.end())
    {
      return histoiter->second;
    }
  cout << "Fun4AllHistoManager::getHisto: ERROR Unknown Histogram " << hname
       << ", The following are implemented: " << endl;
  Print("ALL");
  return NULL;
}



void
Fun4AllHistoManager::Print(const string &what) const
{
  if (what == "ALL" || what == "HISTOS")
    {
      // loop over the map and print out the content (name and location in memory)
      cout << "--------------------------------------" << endl << endl;
      cout << "List of Histos in Fun4AllHistoManager "
	   << Name() << ":" << endl;

      map<const string, TNamed *>::const_iterator hiter;
      for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
	{
	  cout << hiter->first << " is " << hiter->second << endl;
	}
      cout << endl;
    }
  return ;
}

void
Fun4AllHistoManager::Reset()
{
  map<const string, TNamed *>::const_iterator hiter;
  for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
    {
      TNamed* h = hiter->second;
      if (h->InheritsFrom("TH1"))
        (dynamic_cast<TH1*>(h))->Reset();
#if HAS_THNSPARSE
      else if (h->InheritsFrom("THnSparse"))
        (dynamic_cast<THnSparse*>(h))->Reset();
#endif

    }
  return ;
}
