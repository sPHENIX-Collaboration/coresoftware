// $Id: $

/*!
 * \file PHG4InttDeadMapLoader.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4InttDeadMapLoader.h"

#include "InttDeadMap.h"  // for InttDeadMap
#include "InttDeadMapv1.h"

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>

// boost headers
#include <boost/tokenizer.hpp>
// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp>  // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 4 && BOOST_VERSION == 105700)
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma message "ignoring bogus gcc warning in boost header lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#pragma GCC diagnostic warning "-Wuninitialized"
#else
#include <boost/lexical_cast.hpp>
#endif

#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>  // for pair

PHG4InttDeadMapLoader::PHG4InttDeadMapLoader(const std::string &detector)
  : SubsysReco("PHG4InttDeadMapLoader_" + detector)
  , m_detector(detector)
{
}

int PHG4InttDeadMapLoader::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator topiter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(topiter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in RawTowerCalibration::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(runNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_detector);
    runNode->addNode(DetNode);
  }

  // Be careful as a previous calibrator may have been registered for this detector
  std::string deadMapName = "DEADMAP_" + m_detector;
  InttDeadMap *deadmap = findNode::getClass<InttDeadMapv1>(DetNode, deadMapName);
  if (!deadmap)
  {
    deadmap = new InttDeadMapv1();
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(deadmap, deadMapName, "PHObject");
    DetNode->addNode(towerNode);
  }

  assert(deadmap);

  for (const auto &pathiter : m_deadMapPathMap)
  {
    const unsigned int ilayer = pathiter.first;
    const std::string &deadMapPath = pathiter.second;

    int counter = 0;

    PHParameters deadMapParam(m_detector);
    deadMapParam.ReadFromFile(m_detector, "xml", 0, 0, deadMapPath);

    const auto in_par_ranges = deadMapParam.get_all_int_params();

    for (auto iter = in_par_ranges.first; iter != in_par_ranges.second; ++iter)
    {
      const std::string &deadChanName = iter->first;

      if (Verbosity())
      {
        std::cout << "HG4InttDeadMapLoader::InitRun - deadMapParam[" << deadChanName << "] = " << iter->second << ": ";
      }

      boost::char_separator<char> sep("_");
      boost::tokenizer<boost::char_separator<char> > tok(deadChanName, sep);
      boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;

      for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
      {
        if (*tokeniter == "INTT")
        {
          // Form("INTT_%d_%d_%d_%d", ladder_phi, ladder_z, strip_z, strip_phi)

          ++tokeniter;
          assert(tokeniter != tok.end());
          int ladder_phi = boost::lexical_cast<int>(*tokeniter);

          ++tokeniter;
          assert(tokeniter != tok.end());
          int ladder_z = boost::lexical_cast<int>(*tokeniter);

          ++tokeniter;
          assert(tokeniter != tok.end());
          int strip_z = boost::lexical_cast<int>(*tokeniter);

          ++tokeniter;
          assert(tokeniter != tok.end());
          int strip_phi = boost::lexical_cast<int>(*tokeniter);

          deadmap->addDeadChannelIntt(ilayer, ladder_phi, ladder_z, strip_z, strip_phi);
          ++counter;

          if (Verbosity())
          {
            std::cout << "add Intt dead channel ladder_phi" << ladder_phi << " ladder_z" << ladder_z
                      << " strip_z" << strip_z << " strip_phi" << strip_phi;
          }
        }  // if (*tokeniter == "INTT")
        else
        {
          if (Verbosity())
          {
            std::cout << "skip " << deadChanName;
          }
        }

      }  //     for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)

      if (Verbosity())
      {
        std::cout << std::endl;
      }

    }  //  for (const auto iter = in_par_ranges.first; iter != in_par_ranges.second; ++iter)

    std::cout << "PHG4InttDeadMapLoader::" << m_detector << "::InitRun - loading " << counter << " dead channel for layer "
              << ilayer << " from " << deadMapPath << ". Total dead chan = " << deadmap->size() << std::endl;
  }

  if (Verbosity())
  {
    std::cout << "PHG4InttDeadMapLoader::" << m_detector << "::InitRun - loading dead map completed : ";
    deadmap->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
