#include "Fun4AllHepMCOutputManager.h"

#include "PHHepMCGenEvent.h"

#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>

#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

#include <iostream>
#include <string>
#include <vector>

#include <TString.h>
#include <TPRegexp.h>

#include <ostream>
#include <fstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;

static boost::iostreams::filtering_streambuf<boost::iostreams::output> zoutbuffer;

Fun4AllHepMCOutputManager::Fun4AllHepMCOutputManager(const string &myname,
						     const string &filename):
  Fun4AllOutputManager( myname ),
  outfilename(filename),
  comment_written(0),
  filestream(NULL),
  zipstream(NULL)
{
  TString tstr(filename);
  TPRegexp bzip_ext(".bz2$");
  TPRegexp gzip_ext(".gz$");
  
  if (tstr.Contains(bzip_ext)) {
    // use boost iosteam library to compress to bz2 file on the fly
    filestream = new ofstream(filename.c_str(), std::ios::out | std::ios::binary);
    zoutbuffer.push(boost::iostreams::bzip2_compressor(9));
    zoutbuffer.push(*filestream);
    zipstream = new ostream(&zoutbuffer);
    ascii_out = new HepMC::IO_GenEvent(*zipstream);
    
  } else if (tstr.Contains(gzip_ext)) {
    // use boost iosream to compress to gzip file on the fly
    filestream = new ofstream(filename.c_str(), std::ios::out | std::ios::binary);
    zoutbuffer.push(boost::iostreams::gzip_compressor(9));
    zoutbuffer.push(*filestream);
    zipstream = new ostream(&zoutbuffer);
    ascii_out = new HepMC::IO_GenEvent(*zipstream);

  } else {
    // produces normal ascii hepmc file
    ascii_out = new HepMC::IO_GenEvent(filename,std::ios::out);
  }

  if (ascii_out->rdstate())
    {
      cout << "error opening " << outfilename << " exiting " << endl;
      exit(1);
    }
  return ;
}

Fun4AllHepMCOutputManager::~Fun4AllHepMCOutputManager()
{
  if (ascii_out)
    {
      if (!comment_written)
        {
          if (comment.size())
            {
              ascii_out->write_comment(comment);
            }
        }
      ascii_out->clear();
    }

  delete ascii_out;

  if (zoutbuffer.size() > 0) zoutbuffer.reset();
  
  if (zipstream) delete zipstream;
  if (filestream) delete filestream;
  
  return ;
}


void
Fun4AllHepMCOutputManager::Print(const string &what) const
{
  cout << ThisName << " writes " << outfilename << endl;
  if (comment.size())
    {
      cout << "comment : " << comment << endl;
    }
  // base class print method
  Fun4AllOutputManager::Print( what );

  return ;
}

int
Fun4AllHepMCOutputManager::Write(PHCompositeNode* topNode)
{
  if (!comment_written)
    {
      if (comment.size())
	{
      ascii_out->write_comment(comment);
	}
      comment_written = 1;
    }
PHHepMCGenEvent *genevt = findNode::getClass<PHHepMCGenEvent>(topNode,"PHHepMCGenEvent");
HepMC::GenEvent *evt = genevt->getEvent();
  if (!evt)
    {
      cout << PHWHERE << "0 HepMC Pointer" << endl;
      return -1;
    }
  nEvents++;
  ascii_out->write_event(evt);
  return 0;
}

int
Fun4AllHepMCOutputManager::AddComment(const std::string &text)
{
  comment = text;
  comment_written = 0;
  return 0;
}
