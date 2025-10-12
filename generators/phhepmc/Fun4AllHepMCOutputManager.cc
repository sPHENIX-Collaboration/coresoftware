#include "Fun4AllHepMCOutputManager.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"

#include <fun4all/Fun4AllOutputManager.h>  // for Fun4AllOutp...
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <HepMC/IO_GenEvent.h>

#include <TPRegexp.h>
#include <TString.h>

#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <cassert>
#include <cstdlib>  // for exit
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>  // for swap

namespace HepMC
{
  class GenEvent;
}

namespace
{
  boost::iostreams::filtering_streambuf<boost::iostreams::output> zoutbuffer;
}

Fun4AllHepMCOutputManager::Fun4AllHepMCOutputManager(const std::string &myname,
                                                     const std::string &filename)
  : Fun4AllOutputManager(myname)
  , outfilename(filename)
  , comment_written(0)
  , filestream(nullptr)
  , zipstream(nullptr)
  , _embedding_id(0)
{
  TString tstr(filename);
  TPRegexp bzip_ext(".bz2$");
  TPRegexp gzip_ext(".gz$");

  if (tstr.Contains(bzip_ext))
  {
    // use boost iosteam library to compress to bz2 file on the fly
    filestream = new std::ofstream(filename.c_str(), std::ios::out | std::ios::binary);
    zoutbuffer.push(boost::iostreams::bzip2_compressor(9));
    zoutbuffer.push(*filestream);
    zipstream = new std::ostream(&zoutbuffer);
    ascii_out = new HepMC::IO_GenEvent(*zipstream);
  }
  else if (tstr.Contains(gzip_ext))
  {
    // use boost iosream to compress to gzip file on the fly
    filestream = new std::ofstream(filename.c_str(), std::ios::out | std::ios::binary);
    zoutbuffer.push(boost::iostreams::gzip_compressor(9));
    zoutbuffer.push(*filestream);
    zipstream = new std::ostream(&zoutbuffer);
    ascii_out = new HepMC::IO_GenEvent(*zipstream);
  }
  else
  {
    // produces normal ascii hepmc file
    ascii_out = new HepMC::IO_GenEvent(filename, std::ios::out);
  }

  if (ascii_out->rdstate())
  {
    std::cout << "error opening " << outfilename << " exiting " << std::endl;
    exit(1);
  }
  return;
}

Fun4AllHepMCOutputManager::~Fun4AllHepMCOutputManager()
{
  try
  {
    if (ascii_out)
    {
      if (!comment_written && !comment.empty())
      {
        ascii_out->write_comment(comment);
      }
      ascii_out->clear();
      delete ascii_out;
      ascii_out = nullptr;
    }

    if (!zoutbuffer.empty())
    {
      zoutbuffer.reset();
    }

    delete zipstream;
    zipstream = nullptr;

    delete filestream;
    filestream = nullptr;
  }
  catch (const std::exception &e)
  {
    std::cout << "Exception caught in ~Fun4AllHepMCOutputManager: "
              << e.what() << std::endl;
  }
  catch (...)
  {
    std::cout << "Unknown exception caught in ~Fun4AllHepMCOutputManager"
              << std::endl;
  }

  return;
}

void Fun4AllHepMCOutputManager::Print(const std::string &what) const
{
  std::cout << Name() << " writes " << outfilename << std::endl;
  if (!comment.empty())
  {
    std::cout << "comment : " << comment << std::endl;
  }
  // base class print method
  Fun4AllOutputManager::Print(what);

  return;
}

int Fun4AllHepMCOutputManager::Write(PHCompositeNode *topNode)
{
  if (!comment_written)
  {
    if (!comment.empty())
    {
      ascii_out->write_comment(comment);
    }
    comment_written = 1;
  }

  PHHepMCGenEventMap *geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  if (!geneventmap)
  {
    std::cout << "Fun4AllHepMCOutputManager::Write - Fatal Error - missing source node PHHepMCGenEventMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  assert(geneventmap);

  PHHepMCGenEvent *genevt = geneventmap->get(_embedding_id);
  if (!genevt)
  {
    std::cout << "Fun4AllHepMCOutputManager::Write - Warning - missing sub-event with embedding ID" << _embedding_id << " on node PHHepMCGenEventMap" << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  assert(genevt);

  HepMC::GenEvent *evt = genevt->getEvent();
  if (!evt)
  {
    std::cout << PHWHERE << "0 HepMC Pointer" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  assert(evt);

  IncrementEvents(1);
  ascii_out->write_event(evt);
  return Fun4AllReturnCodes::EVENT_OK;
}

int Fun4AllHepMCOutputManager::AddComment(const std::string &text)
{
  comment = text;
  comment_written = 0;
  return 0;
}
