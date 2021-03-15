/*!
 * \file Fun4AllSingleDstPileupInputManager.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "Fun4AllSingleDstPileupInputManager.h"
#include "Fun4AllDstPileupMerger.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <ffaobjects/RunHeader.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIntegrate.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE, PHReadOnly, PHRunTree

#include <TSystem.h>

#include <gsl/gsl_randist.h>

#include <cassert>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair

//_____________________________________________________________________________
Fun4AllSingleDstPileupInputManager::Fun4AllSingleDstPileupInputManager(const std::string &name, const std::string &nodename, const std::string &topnodename)
  : Fun4AllInputManager(name, nodename, topnodename)
{
  // initialize random generator
  const uint seed = PHRandomSeed();
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );
}

//_____________________________________________________________________________
int Fun4AllSingleDstPileupInputManager::fileopen(const std::string &filenam)
{
  /*
  this is largely copied from fun4all/Fun4AllDstInputManager::fileopen
  with additional code to handle the background IManager
  */

  auto se = Fun4AllServer::instance();
  if (IsOpen())
  {
    std::cout << "Closing currently open file "
              << FileName()
              << " and opening " << filenam << std::endl;
    fileclose();
  }
  FileName(filenam);
  FROG frog;
  m_fullfilename = frog.location(FileName());
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": opening file " << m_fullfilename << std::endl;
  }
  // sanity check - the IManager must be nullptr when this method is executed
  // if not something is very very wrong and we must not continue
  assert( !m_IManager );

  // first read the runnode if not disabled
  if (m_ReadRunTTree)
  {
    m_IManager.reset(new PHNodeIOManager(m_fullfilename, PHReadOnly, PHRunTree));
    if (m_IManager->isFunctional())
    {
      m_runNode = se->getNode(m_RunNode, TopNodeName());
      m_IManager->read(m_runNode);

      // get the current run number
      auto runheader = findNode::getClass<RunHeader>(m_runNode, "RunHeader");
      if (runheader)
      {
        SetRunNumber(runheader->get_RunNumber());
      }
      // delete our internal copy of the runnode when opening subsequent files
      assert( !m_runNodeCopy );
      m_runNodeCopy.reset(new PHCompositeNode("RUNNODECOPY"));
      if (!m_runNodeSum)
      {
        m_runNodeSum.reset(new PHCompositeNode("RUNNODESUM"));
      }

      {
        // read run node using temporary node iomanager
        PHNodeIOManager(m_fullfilename, PHReadOnly, PHRunTree).read(m_runNodeCopy.get());
      }

      PHNodeIntegrate integrate;
      integrate.RunNode(m_runNode);
      integrate.RunSumNode(m_runNodeSum.get());
      // run recursively over internal run node copy and integrate objects
      PHNodeIterator mainIter(m_runNodeCopy.get());
      mainIter.forEach(integrate);
      // we do not need to keep the internal copy, keeping it would crate
      // problems in case a subsequent file does not contain all the
      // runwise objects from the previous file. Keeping this copy would then
      // integrate the missing object again with the old copy
      m_runNodeCopy.reset();
    }
    // DLW: move the delete outside the if block to cover the case where isFunctional() fails
    m_IManager.reset();
  }

  // create internal dst node
  if (!m_dstNodeInternal)
  {
    m_dstNodeInternal.reset(new PHCompositeNode("DST_INTERNAL"));
  }

  // update dst node from fun4all server
  m_dstNode = se->getNode(InputNode(), TopNodeName());

  // open file in both active and background input manager
  m_IManager.reset(new PHNodeIOManager(m_fullfilename, PHReadOnly));
  m_IManager_background.reset(new PHNodeIOManager(m_fullfilename, PHReadOnly));
  if (m_IManager->isFunctional())
  {
    IsOpen(1);
    m_ievent_thisfile = 0;
    setBranches();                // set branch selections
    AddToFileOpened(FileName());  // add file to the list of files which were opened
    return 0;
  }
  else
  {
    std::cout << PHWHERE << ": " << Name() << " Could not open file " << FileName() << std::endl;
    m_IManager.reset();
    m_IManager_background.reset();
    return -1;
  }
}

//_____________________________________________________________________________
int Fun4AllSingleDstPileupInputManager::run(const int nevents)
{

  if (!IsOpen())
  {
    if (FileListEmpty())
    {
      if (Verbosity() > 0) std::cout << Name() << ": No Input file open" << std::endl;
      return -1;
    }
    else
    {
      if (OpenNextFile())
      {
        std::cout << Name() << ": No Input file from filelist opened" << std::endl;
        return -1;
      }
    }
  }

  if (Verbosity() > 3)
  {
    std::cout << "Getting Event from " << Name() << std::endl;
  }

readagain:

  // read main event to dstNode
  auto dummy = m_IManager->read(m_dstNode);
  int ncount = 0;
  while (dummy)
  {
    ++ncount;
    if (nevents > 0 && ncount >= nevents) break;
    dummy = m_IManager->read(m_dstNode);
  }

  if (!dummy)
  {
    fileclose();
    if (!OpenNextFile())
      goto readagain;
    else
      return -1;
  }

  m_ievent_total += ncount;
  m_ievent_thisfile += ncount;

  // check events consistency
  if (m_ievent_thisfile != static_cast<int>(m_IManager->getEventNumber()))
  {
    std::cout << PHWHERE
              << " inconsistent event counters between inputmanager and ionode manager: "
              << " m_ievent_thisfile: " << m_ievent_thisfile
              << " m_IManager->getEventNumber(): " << m_IManager->getEventNumber()
              << std::endl;
    gSystem->Exit(1);
  }

  // check if the local SubsysReco discards this event
  if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << "Fun4AllSingleDstPileupInputManager::run - skipped event " << m_ievent_thisfile - 1 << std::endl;
    goto readagain;
  }

  // load relevant DST nodes to internal pointers
  if (Verbosity() > 0)
  {
    std::cout << "Fun4AllSingleDstPileupInputManager::run - loaded event " << m_ievent_thisfile - 1 << std::endl;
  }

  Fun4AllDstPileupMerger merger;
  merger.load_nodes(m_dstNode);

  // generate background collisions
  const double mu = m_collision_rate*m_time_between_crossings*1e-9;

  const int min_crossing = m_tmin/m_time_between_crossings;
  const int max_crossing = m_tmax/m_time_between_crossings;
  int ievent_thisfile = m_ievent_thisfile;
  int neventsbackground = 0;
  for( int icrossing = min_crossing; icrossing <= max_crossing; ++icrossing )
  {
    const double crossing_time = m_time_between_crossings * icrossing;
    const int ncollisions = gsl_ran_poisson(m_rng.get(), mu);
    for (int icollision = 0; icollision < ncollisions; ++icollision)
    {

      // try read
      if(!m_IManager_background->read(m_dstNodeInternal.get(), ievent_thisfile) ) break;

      // merge

      if (Verbosity() > 0)
      {
	std::cout << "Fun4AllSingleDstPileupInputManager::run - merged background event " << ievent_thisfile << " time: " << crossing_time << std::endl;
      }
      merger.copy_background_event(m_dstNodeInternal.get(), crossing_time);

      ++neventsbackground;
      ++ievent_thisfile;
    }
  }

  // jump event counter to the last background accepted event
  if( neventsbackground > 0 ) PushBackEvents( -neventsbackground );
  return 0;
}

//_____________________________________________________________________________
int Fun4AllSingleDstPileupInputManager::fileclose()
{
  if (!IsOpen())
  {
    std::cout << Name() << ": fileclose: No Input file open" << std::endl;
    return -1;
  }
  m_IManager.reset();
  m_IManager_background.reset();
  IsOpen(0);
  UpdateFileList();
  return 0;
}

//_____________________________________________________________________________
int Fun4AllSingleDstPileupInputManager::BranchSelect(const std::string &branch, const int iflag)
{
  int myflag = iflag;
  // if iflag > 0 the branch is set to read
  // if iflag = 0, the branch is set to NOT read
  // if iflag < 0 the branchname is erased from our internal branch read map
  // this does not have any effect on phool yet
  if (myflag < 0)
  {
    std::map<const std::string, int>::iterator branchiter;
    branchiter = m_branchread.find(branch);
    if (branchiter != m_branchread.end())
    {
      m_branchread.erase(branchiter);
    }
    return 0;
  }

  if (myflag > 0)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Setting Root Tree Branch: " << branch << " to read" << std::endl;
    }
    myflag = 1;
  }
  else
  {
    if (Verbosity() > 1)
    {
      std::cout << "Setting Root Tree Branch: " << branch << " to NOT read" << std::endl;
    }
  }
  m_branchread[branch] = myflag;
  return 0;
}

//_____________________________________________________________________________
int Fun4AllSingleDstPileupInputManager::setBranches()
{
  if (m_IManager)
  {
    if (!m_branchread.empty())
    {
      std::map<const std::string, int>::const_iterator branchiter;
      for (branchiter = m_branchread.begin(); branchiter != m_branchread.end(); ++branchiter)
      {
        m_IManager->selectObjectToRead(branchiter->first.c_str(), branchiter->second);
        if( m_IManager_background )
        {
          m_IManager_background->selectObjectToRead(branchiter->first.c_str(), branchiter->second);
        }
        if (Verbosity() > 0)
        {
          std::cout << branchiter->first << " set to " << branchiter->second << std::endl;
        }
      }
    }
  }
  else
  {
    std::cout << PHWHERE << " " << Name() << ": You can only call this function after a file has been opened" << std::endl;
    std::cout << "Do not worry, the branches will be set as soon as you open a file" << std::endl;
    return -1;
  }
  return 0;
}

//_____________________________________________________________________________
void Fun4AllSingleDstPileupInputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "BRANCH")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of selected branches in Fun4AllSingleDstPileupInputManager " << Name() << ":" << std::endl;

    std::map<const std::string, int>::const_iterator iter;
    for (iter = m_branchread.begin(); iter != m_branchread.end(); ++iter)
    {
      std::cout << iter->first << " is switched ";
      if (iter->second)
      {
        std::cout << "ON";
      }
      else
      {
        std::cout << "OFF";
      }
      std::cout << std::endl;
    }
  }
  if ((what == "ALL" || what == "PHOOL") && m_IManager)
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "PHNodeIOManager print in Fun4AllSingleDstPileupInputManager " << Name() << ":" << std::endl;
    m_IManager->print();
  }
  Fun4AllInputManager::Print(what);
  return;
}

//_____________________________________________________________________________
int Fun4AllSingleDstPileupInputManager::PushBackEvents(const int i)
{
  if (m_IManager)
  {
    unsigned EventOnDst = m_IManager->getEventNumber();
    EventOnDst -= static_cast<unsigned>(i);
    m_ievent_thisfile -= i;
    m_ievent_total -= i;
    m_IManager->setEventNumber(EventOnDst);
    return 0;
  }
  std::cout << PHWHERE << Name() << ": could not push back events, Imanager is NULL"
            << " probably the dst is not open yet (you need to call fileopen or run 1 event for lists)" << std::endl;
  return -1;
}

