// $Id: $

/*!
 * \file RawTowerDeadMapLoader.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "RawTowerDeadMapLoader.h"

#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerDeadMapv1.h>
#include <calobase/RawTowerDefs.h>

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

// boost headers
#include <boost/token_iterator.hpp>  // for token_iterator
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
#include <utility>

RawTowerDeadMapLoader::RawTowerDeadMapLoader(const std::string &detector)
  : SubsysReco("RawTowerDeadMapLoader_" + detector)
  , m_detector(detector)
{
}

int RawTowerDeadMapLoader::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find Run node in RawTowerCalibration::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(runNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst(
      "PHCompositeNode", m_detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_detector);
    runNode->addNode(DetNode);
  }

  // Be careful as a previous calibrator may have been registered for this detector
  std::string deadMapName = "DEADMAP_" + m_detector;
  RawTowerDeadMap *m_deadmap = findNode::getClass<RawTowerDeadMapv1>(DetNode, deadMapName);
  if (!m_deadmap)
  {
    const RawTowerDefs::CalorimeterId caloid = RawTowerDefs::convert_name_to_caloid(m_detector);

    m_deadmap = new RawTowerDeadMapv1(caloid);
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(
        m_deadmap, deadMapName, "PHObject");
    DetNode->addNode(towerNode);
  }

  assert(m_deadmap);

  std::cout << "RawTowerDeadMapLoader::" << m_detector << "::InitRun - loading dead map from " << m_deadMapPath << std::endl;

  PHParameters deadMapParam(m_detector);
  deadMapParam.ReadFromFile(m_detector, "xml", 0, 0, m_deadMapPath);

  const auto in_par_ranges = deadMapParam.get_all_int_params();

  for (auto iterA = in_par_ranges.first; iterA != in_par_ranges.second; ++iterA)
  {
    const std::string &deadChanName = iterA->first;

    if (Verbosity())
    {
      std::cout << "deadMapParam[" << deadChanName << "] = " << iterA->second << ": ";
    }

    boost::char_separator<char> sep("_");
    boost::tokenizer<boost::char_separator<char> > tok(deadChanName, sep);
    boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;

    for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
    {
      if (*tokeniter == "deadtower")
      {
        // Form("deadtower_eta_%d_phi_%d", eta, phi)

        ++tokeniter;
        assert(tokeniter != tok.end());

        if (*tokeniter == "eta")
        {
          ++tokeniter;
          assert(tokeniter != tok.end());
          const int eta = boost::lexical_cast<int>(*tokeniter);

          ++tokeniter;
          assert(tokeniter != tok.end());
          assert(*tokeniter == "phi");

          ++tokeniter;
          assert(tokeniter != tok.end());
          const int phi = boost::lexical_cast<int>(*tokeniter);

          m_deadmap->addDeadTower(eta, phi);

          if (Verbosity())
          {
            std::cout << "add dead channel eta" << eta << " phi" << phi;
          }
        }  // if (*tokeniter == "eta")
        else
        {
          if (Verbosity())
          {
            std::cout << "skip " << deadChanName;
          }
        }

      }  //      if (*tokeniter == "deadtower")
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

  }  //  for (const auto iterA = in_par_ranges.first; iterA != in_par_ranges.second; ++iterA)

  if (Verbosity())
  {
    std::cout << "RawTowerDeadMapLoader::" << m_detector << "::InitRun - loading dead map completed : ";
    m_deadmap->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
