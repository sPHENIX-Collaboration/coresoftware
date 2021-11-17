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

using namespace std;

//_________________________________________________
PHTFileServer::SafeTFile::TFileMap PHTFileServer::SafeTFile::_map;

//_________________________________________________
PHTFileServer::~PHTFileServer(void)
{
  if (!SafeTFile::file_map().empty()) close();
}

//_________________________________________________
void PHTFileServer::open(const string& filename, const string& type)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
  {
    ostringstream what;
    what << "PHTFileServer::open - file " << filename << " already opened.";
    cout << (what.str()) << endl;

    // increment counter; change TDirectory
    iter->second->counter()++;
    iter->second->cd();
  }
  else
  {
    ostringstream what;
    what << "PHTFileServer::open - opening file " << filename << " (" << type << ")";
    cout << (what.str()) << endl;

    // create new SafeTFile; insert in map; change TDirectory
    SafeTFile* file(new SafeTFile(filename, type));
    if (!file->IsOpen()) cout << ("PHTFileServer::open - error opening TFile") << endl;
    SafeTFile::file_map().insert(make_pair(filename, file));
    file->cd();
  }
}

//_________________________________________________
bool PHTFileServer::flush(const string& filename)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
    iter->second->Flush();
  else
  {
    ostringstream what;
    what << "PHTFileServer::flush - file " << filename << " not found";
    cout << (what.str()) << endl;
    return false;
  }

  return true;
}

//_________________________________________________
bool PHTFileServer::cd(const string& filename)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
    iter->second->cd();
  else
  {
    ostringstream what;
    what << "PHTFileServer::flush - file " << filename << " not found";
    cout << (what.str()) << endl;
    return false;
  }

  return true;
}

//_________________________________________________
bool PHTFileServer::write(const string& filename)
{
  SafeTFile::TFileMap::iterator iter(SafeTFile::file_map().find(filename));
  if (iter != SafeTFile::file_map().end())
  {
    if (iter->second->counter() > 1)
    {
      iter->second->counter()--;
      ostringstream what;
      what << "PHTFileServer::write - file " << filename << " still in use.";
      cout << (what.str()) << endl;
    }
    else if (iter->second->counter() == 1)
    {
      iter->second->Write();
      iter->second->counter()--;
      ostringstream what;
      what << "PHTFileServer::write - writing file " << filename << ".";
      cout << (what.str()) << endl;
    }
    else
    {
      iter->second->Write();
      ostringstream what;
      what << "PHTFileServer::write - warning: too many calls for file " << filename << ".";
      cout << (what.str()) << endl;
    }
  }
  else
  {
    ostringstream what;
    what << "PHTFileServer::write - file " << filename << " not found";
    cout << (what.str()) << endl;
    return false;
  }

  return true;
}

//__________________________________________
void PHTFileServer::close(void)
{
  // close
  //  MUTOO::TRACE( "PHTFileServer::close" );
  for (SafeTFile::TFileMap::iterator iter = SafeTFile::file_map().begin(); iter != SafeTFile::file_map().end(); ++iter)
  {
    if (iter->second->IsOpen())
    {
      if (iter->second->counter())
      {
        ostringstream what;
        what << "PHTFileServer::close - file " << iter->first << " forced write with kWriteDelete.";
        cout << (what.str()) << endl;
        iter->second->Write("0", TObject::kWriteDelete);
      }

      // close TFile
      ostringstream what;
      what << "PHTFileServer::close - closing " << iter->first << ".";
      iter->second->Close();
      cout << (what.str()) << endl;
    }
  }

  // clear file map
  SafeTFile::file_map().clear();
}

//__________________________________________________________________________________
PHTFileServer::SafeTFile::~SafeTFile(void)
{
  // see if TFile is still open
  if (IsOpen())
  {
    // check if TFile needs writing first
    if (_counter)
    {
      ostringstream what;
      what << "PHTFileServer::SafeTFile::~SafeTFile - file " << _filename << " forced write with kWriteDelete.";
      cout << (what.str()) << endl;
      Write("0", TObject::kWriteDelete);
    }

    ostringstream what;
    what << "PHTFileServer::SafeTFile::~SafeTFile - closing " << _filename << ".";
    cout << (what.str()) << endl;
    Close();
  }

  /* 
  remove this filename from the make to make sure that PHTFileServer
  does not try to write/close this TFile during the destructor
  */
  _map.erase(_filename);
}
