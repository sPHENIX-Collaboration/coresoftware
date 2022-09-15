#include "TDirectoryHelper.h"

#include <TCollection.h>  // for TIter
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TList.h>  // for TList
#include <TObject.h>
#include <TROOT.h>

#include <algorithm>  // for max
#include <cassert>
#include <cstddef>  // for size_t
#include <iostream>
#include <memory>  // for allocator_traits<>::value_type
#include <string>
#include <vector>

//_____________________________________________________________________________
void TDirectoryHelper::copyToFile(TDirectory* src, TFile* dest)
{
  TDirectory* save = gDirectory;

  // We basically have two cases here to consider, depending
  // on whether (1) or not (2) a TFile was opened prior to calling
  // the PhotonHistogrammer ctor.
  // 1. Nothing special to do, all the TDirectory structure
  //    is already attached to the correct file.
  // 2. We have to "duplicate" the TDirectory structure into
  //    the newly opened file and save the HistogramCollection(s) there.

  if (!dest || !src)
  {
    return;
  }

  if (!dest->IsWritable())
  {
    std::cout << "TDirectoryHelper::copyToFile : destination file is not "
                 " writeable"
              << std::endl;
    return;
  }

  duplicateDir(dest, src);

  //       list<TDirectory*> mothers;

  //       TDirectory* mother = fDir;

  //       while ( (mother = dynamic_cast<TDirectory*>(mother->GetMother()) ) )
  // 	{
  // 	  std::string motherName = mother->GetName();
  // 	  if (motherName != "Rint" )
  // 	    {
  // 	      mothers.push_front(mother);
  // 	    }
  // 	}

  //       TDirectory* currentDir = save;

  //       list<TDirectory*>::const_iterator it;

  //       for ( it = mothers.begin(); it != mothers.end() ; it++ )
  // 	{

  // 	  TDirectory* dir;
  // 	  if ( (dir=(TDirectory*)currentDir->FindObject((*it)->GetName()) ))
  // 	    {
  // 	      currentDir = dir;
  // 	    }
  // 	  else
  // 	    {
  // 	      currentDir = currentDir->mkdir((*it)->GetName(),
  // 					     (*it)->GetTitle());
  // 	    }
  // 	}

  //       TDirectoryHelper::duplicateDir
  // 	((TDirectory*)save->FindObject(mothers.back()->GetName()),fDir);

  //     }

  save->cd();
}

//_____________________________________________________________________________
void TDirectoryHelper::duplicateDir(TDirectory* dest, TDirectory* source)
{
  dest->cd();

  TDirectory* newdir;

  newdir = static_cast<TDirectory*>(gDirectory->FindObject(source->GetName()));

  if (!newdir)
  {
    newdir = dest->mkdir(source->GetName(), source->GetTitle());
  }

  newdir->cd();

  TIter next(source->GetList());
  TObject* obj;

  while ((obj = next()))
  {
    TDirectory* dir = dynamic_cast<TDirectory*>(obj);
    if (dir)
    {
      duplicateDir(newdir, dir);
    }
    else
    {
      obj->Write();
    }
  }
}

//_____________________________________________________________________________
bool TDirectoryHelper::mkpath(TDirectory* dir, const std::string& pathin)
{
  static std::vector<std::string> paths;

  splitPath(pathin, paths);

  TDirectory* currentdir = dir;

  for (auto & path : paths)
  {
    currentdir->cd();

    currentdir = dynamic_cast<TDirectory*>(gDirectory->Get(path.c_str()));
    if (!currentdir)
    {
      currentdir = gDirectory->mkdir(path.c_str());
      assert(currentdir != nullptr);
    }
  }

  return true;
}

//_____________________________________________________________________________
TDirectory*
TDirectoryHelper::mkdir(TDirectory* topDir,
                        const std::string& path,
                        std::vector<std::string>* titles)
{
  TDirectory* save = gDirectory;

  TDirectory* dir = topDir;
  TDirectory* tdir = dir;

  if (topDir == nullptr)
  {
    gROOT->cd();
    tdir = gDirectory;
  }

  dir = tdir;

  dir->cd();
  std::vector<std::string> paths;

  splitPath(path, paths);

  for (size_t i = 0; i < paths.size(); i++)
  {
    TDirectory* subdir = static_cast<TDirectory*>(dir->FindObject(paths[i].c_str()));
    if (subdir == nullptr)
    {
      if (titles && i < titles->size())
      {
        dir = dir->mkdir(paths[i].c_str(), (*titles)[i].c_str());
      }
      else
      {
        dir = dir->mkdir(paths[i].c_str());
      }
    }
    else
    {
      dir = subdir;
    }
    dir->cd();
  }

  save->cd();

  return dir;
}

//_____________________________________________________________________________
bool TDirectoryHelper::pathIsInDir(const std::string& path, TDirectory* dir)
{
  // This is to avoid annoying ROOT message when a directory does not exist
  // in Cd(), so we provide this small method to check whereas
  // a path exists under a directory, but without issuing error message
  // in case of failure (just returning false in this case).

  TDirectory* dirsave = gDirectory;

  static std::vector<std::string> paths;

  paths.clear();
  splitPath(path, paths);

  bool ok = true;

  TDirectory* cdir = dir;

  for (size_t i = 0; i < paths.size() && ok; i++)
  {
    cdir->cd();

    cdir = dynamic_cast<TDirectory*>(cdir->Get(paths[i].c_str()));
    if (!cdir)
    {
      ok = false;
    }
  }

  dirsave->cd();

  return ok;
}

//_____________________________________________________________________________
TH1* TDirectoryHelper::getHisto(TDirectory* dir, const std::string& histoname,
                                const std::string& where)
{
  // Try to find histogram named histoname into directory dir, under
  // path=where (where e.g. = "/Cut##/OK/C#/V#").

  TH1* rv = nullptr;

  bool ok = pathIsInDir(where, dir);

  if (ok)
  {
    // Path is in dir, we can safely (i.e. without getting ROOT error message
    // on stdout) cd into it.
    //    dir->cd();
    ok = dir->cd(where.c_str());
    assert(ok == true);
    TObject* obj = gDirectory->Get(histoname.c_str());
    if (obj)
    {
      rv = dynamic_cast<TH1*>(obj);
      if (!rv)
      {
        std::cout << "GetHisto : object " << histoname << " is not a TH1" << std::endl;
      }
    }
  }
  return rv;
}

//_____________________________________________________________________________
void TDirectoryHelper::splitPath(const std::string& path,
                                 std::vector<std::string>& paths)
{
  // Given a path e.g. /Cut##/OK/C#/V#, will return
  // a vector of string with Cut#,

  paths.clear();

  std::string str = path;

  if (str.empty())
  {
    return;
  }

  std::vector<size_t> slashes_pos;

  if (str[0] != '/')
  {
    str.insert(str.begin(), '/');
  }

  if (str[str.size() - 1] != '/')
  {
    str.push_back('/');
  }

  for (size_t i = 0; i < str.size(); i++)
  {
    if (str[i] == '/')
    {
      slashes_pos.push_back(i);
    }
  }

  if (not slashes_pos.empty())
  {
    for (size_t i = 0; i < slashes_pos.size() - 1; i++)
    {
      paths.push_back(str.substr(slashes_pos[i] + 1,
                                 slashes_pos[i + 1] - slashes_pos[i] - 1));
    }
  }
}
