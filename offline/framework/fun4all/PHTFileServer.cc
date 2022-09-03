//////////////////////////////////////////////////////////////////
/*! 
  \file PHTFileServer.cc
  \brief TFile clean handling
  \author  Hugo Pereira
  \version $Revision: 1.8 $
  \date    $Date: 2011/08/11 14:25:49 $
*/
//////////////////////////////////////////////////////////////////

#include "PHTFileServer.h"

#include <TObject.h>  // for TObject, TObject::kWriteDelete

#include <iostream>  // for operator<<, basic_ostream, ostringstream, endl
#include <sstream>
#include <utility>  // for pair, make_pair

//_________________________________________________
PHTFileServer::SafeTFile::TFileMap PHTFileServer::SafeTFile::_map;

//_________________________________________________
PHTFileServer::~PHTFileServer()
{
  if (!SafeTFile::file_map().empty()) close();
}

//_________________________________________________
void PHTFileServer::open(const std::string& filename, const std::string& type)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
  {
    std::ostringstream what;
    what << "PHTFileServer::open - file " << filename << " already opened.";
    std::cout << (what.str()) << std::endl;

    // increment counter; change TDirectory
    iter->second->counter()++;
    iter->second->cd();
  }
  else
  {
    std::ostringstream what;
    what << "PHTFileServer::open - opening file " << filename << " (" << type << ")";
    std::cout << (what.str()) << std::endl;

    // create new SafeTFile; insert in map; change TDirectory
    SafeTFile* file(new SafeTFile(filename, type));
    if (!file->IsOpen()) std::cout << ("PHTFileServer::open - error opening TFile") << std::endl;
    SafeTFile::file_map().insert(make_pair(filename, file));
    file->cd();
  }
}

//_________________________________________________
bool PHTFileServer::flush(const std::string& filename)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
    iter->second->Flush();
  else
  {
    std::ostringstream what;
    what << "PHTFileServer::flush - file " << filename << " not found";
    std::cout << (what.str()) << std::endl;
    return false;
  }

  return true;
}

//_________________________________________________
bool PHTFileServer::cd(const std::string& filename)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
    iter->second->cd();
  else
  {
    std::ostringstream what;
    what << "PHTFileServer::flush - file " << filename << " not found";
    std::cout << (what.str()) << std::endl;
    return false;
  }

  return true;
}

//_________________________________________________
bool PHTFileServer::write(const std::string& filename)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
  {
    if (iter->second->counter() > 1)
    {
      iter->second->counter()--;
      std::ostringstream what;
      what << "PHTFileServer::write - file " << filename << " still in use.";
      std::cout << (what.str()) << std::endl;
    }
    else if (iter->second->counter() == 1)
    {
      iter->second->Write();
      iter->second->counter()--;
      std::ostringstream what;
      what << "PHTFileServer::write - writing file " << filename << ".";
      std::cout << (what.str()) << std::endl;
    }
    else
    {
      iter->second->Write();
      std::ostringstream what;
      what << "PHTFileServer::write - warning: too many calls for file " << filename << ".";
      std::cout << (what.str()) << std::endl;
    }
  }
  else
  {
    std::ostringstream what;
    what << "PHTFileServer::write - file " << filename << " not found";
    std::cout << (what.str()) << std::endl;
    return false;
  }

  return true;
}

//__________________________________________
void PHTFileServer::close()
{
  // close
  //  MUTOO::TRACE( "PHTFileServer::close" );
  for (auto & iter : SafeTFile::file_map())
  {
    if (iter.second->IsOpen())
    {
      if (iter.second->counter())
      {
        std::ostringstream what;
        what << "PHTFileServer::close - file " << iter.first << " forced write with kWriteDelete.";
        std::cout << (what.str()) << std::endl;
        iter.second->Write("0", TObject::kWriteDelete);
      }

      // close TFile
      std::ostringstream what;
      what << "PHTFileServer::close - closing " << iter.first << ".";
      iter.second->Close();
      std::cout << (what.str()) << std::endl;
    }
  }

  // clear file map
  SafeTFile::file_map().clear();
}

//__________________________________________________________________________________
PHTFileServer::SafeTFile::~SafeTFile()
{
  // see if TFile is still open
  if (IsOpen())
  {
    // check if TFile needs writing first
    if (_counter)
    {
      std::ostringstream what;
      what << "PHTFileServer::SafeTFile::~SafeTFile - file " << _filename << " forced write with kWriteDelete.";
      std::cout << (what.str()) << std::endl;
      Write("0", TObject::kWriteDelete);
    }

    std::ostringstream what;
    what << "PHTFileServer::SafeTFile::~SafeTFile - closing " << _filename << ".";
    std::cout << (what.str()) << std::endl;
    Close();
  }

  /* 
  remove this filename from the make to make sure that PHTFileServer
  does not try to write/close this TFile during the destructor
  */
  _map.erase(_filename);
}
