#include "MetaDataTreeWriter.h"

#include "MetaDataTree.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/phool.h>  // PHWHERE

#include <TDirectory.h>
#include <TTree.h>

#include <iostream>

MetaDataTreeWriter::MetaDataTreeWriter(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_FileName(filename)
{
}

int MetaDataTreeWriter::Init(PHCompositeNode * /*topNode*/)
{
  if (m_FileName.empty())
  {
    std::cout << PHWHERE << " " << Name() << ": no output file given (constructor or set_filename)" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // PHTFileServer keeps one TFile per name and counts its users: if an analysis module already
  // opened this file, this only adds us as a user (the "RECREATE" is then ignored), so we do not
  // overwrite its file. open() also cd's into the file; TContext puts gDirectory back afterwards.
  TDirectory::TContext ctx;
  PHTFileServer::open(m_FileName, "RECREATE");
  m_Opened = true;
  return Fun4AllReturnCodes::EVENT_OK;
}

int MetaDataTreeWriter::End(PHCompositeNode * /*topNode*/)
{
  if (!m_Opened)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  {
    // make the tree inside the output file (gDirectory after cd), restoring gDirectory afterwards
    TDirectory::TContext ctx;
    if (PHTFileServer::cd(m_FileName))
    {
      TTree *tree = MetaDataTree::Make(gDirectory, m_TreeName);
      if (!tree)
      {
        std::cout << PHWHERE << " " << Name() << ": could not make the provenance tree, " << m_FileName
                  << " will not have one" << std::endl;
      }
      else if (Verbosity() > 0)
      {
        std::cout << Name() << ": provenance tree " << m_TreeName << " made in " << m_FileName << std::endl;
      }
    }
  }
  // No tree->Write() here: the tree belongs to the file, and PHTFileServer::write() calls
  // TFile::Write() (which writes the tree) once the last user of the file is done. If other
  // modules still use the file, this only takes us off its user count, and the tree is written
  // with their last write(). Writing it ourselves as well would store it twice (provenance;1, ;2).
  PHTFileServer::write(m_FileName);
  m_Opened = false;
  return Fun4AllReturnCodes::EVENT_OK;
}
