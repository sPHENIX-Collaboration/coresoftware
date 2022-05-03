// Tell emacs that this is a C++ source
//  -*- C++ -*-.
//////////////////////////////////////////////////////////////////
/*! 
  \file PHTFileServer.h
  \brief TFile clean handling
  \author  Hugo Pereira
  \version $Revision: 1.5 $
  \date    $Date: 2007/10/27 14:27:57 $
*/
//////////////////////////////////////////////////////////////////

#ifndef FUN4ALL_PHTFILESERVER_H
#define FUN4ALL_PHTFILESERVER_H

#include <TFile.h>

#include <map>
#include <sstream>
#include <string>

/*! 
  \brief TFile clean handling. It allow independant classes to access
  the same TFile and write ntuple to it. TFiles get written only when as many
  write request are achieved as open request. 
  It get closed when the server is deleted
*/
class PHTFileServer
{
 public:
  //! return reference to class singleton
  static PHTFileServer& get(void)
  {
    static PHTFileServer singleton;
    return singleton;
  }

  //! destructor. All non close TFiles are closed, with a warning.
  virtual ~PHTFileServer();

  /*! \brief 
    open a SafeTFile. If filename is not found in the map, create a new TFile
    and append to the map; increment counter otherwise
  */
  void open(const std::string& filename, const std::string& type = "RECREATE");

  //! flush TFile matching filename
  bool flush(const std::string& filename);

  //! change to directory of TFile matching filename
  bool cd(const std::string& filename);

  /*! \brief
    if TFile is found in map and counter is 0, close the TFile, 
    decrement counter otherwise
  */
  bool write(const std::string& filename);

  //! close all TFiles
  void close(void);

 private:
  //! constructor
  PHTFileServer(void)
  {
  }

  //! local class to store TFile and counter
  class SafeTFile : public TFile
  {
   public:
    //! constructor
    SafeTFile(const std::string& filename, const std::string& type = "RECREATE")
      : TFile(filename.c_str(), type.c_str())
      , _filename(filename)
      , _counter(1)
    {
    }

    //! destructor
    ~SafeTFile(void) override;

    //! get reference to counter
    int& counter()
    {
      return _counter;
    }

    //! get const reference to counter
    const int& counter() const
    {
      return _counter;
    }

    //! shortcut for SafeTFile map
    typedef std::map<std::string, SafeTFile*> TFileMap;

    static TFileMap& file_map(void)
    {
      return _map;
    }

   private:
    //! filename (for debugging)
    std::string _filename;

    //! call counter
    int _counter;

    //! filename/SafeTFile map
    static TFileMap _map;
  };
};

#endif
