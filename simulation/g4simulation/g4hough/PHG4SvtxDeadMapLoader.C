// $Id: $

/*!
 * \file PHG4SvtxDeadMapLoader.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4SvtxDeadMapLoader.h"

#include <calobase/SvtxDeadMapv1.h>
#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

// boost headers
#include <boost/foreach.hpp>
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

using namespace std;

PHG4SvtxDeadMapLoader::PHG4SvtxDeadMapLoader(const std::string &detector)
  : SubsysReco("PHG4SvtxDeadMapLoader_" + detector)
  , m_detector(detector)
  , m_deadmap(nullptr)
{
}

PHG4SvtxDeadMapLoader::~PHG4SvtxDeadMapLoader()
{
}

int PHG4SvtxDeadMapLoader::InitRun(PHCompositeNode *topNode)
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
  string deadMapName = "DEADMAP_" + m_detector;
  SvtxDeadMap *m_deadmap = findNode::getClass<SvtxDeadMapv1>(DetNode, deadMapName);
  if (!m_deadmap)
  {
    m_deadmap = new SvtxDeadMapv1();
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(
        m_deadmap, deadMapName, "PHObject");
    DetNode->addNode(towerNode);
  }

  assert(m_deadmap);

  for (const auto pathiter : m_deadMapPathMap)
  {
    const unsigned int ilayer = pathiter.first;
    const string &deadMapPath = pathiter.second;

    int counter = 0;

    PHParameters deadMapParam(m_detector);
    deadMapParam.ReadFromFile(m_detector, "xml", 0, 0, deadMapPath);

    const auto in_par_ranges = deadMapParam.get_all_int_params();

    for (auto iter = in_par_ranges.first; iter != in_par_ranges.second; ++iter)
    {
      const string &deadChanName = iter->first;

      if (Verbosity())
      {
        cout << "deadMapParam[" << deadChanName << "] = " << iter->second << ": ";
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
            const unsigned int eta = boost::lexical_cast<unsigned int>(*tokeniter);

            ++tokeniter;
            assert(tokeniter != tok.end());
            assert(*tokeniter == "phi");

            ++tokeniter;
            assert(tokeniter != tok.end());
            const unsigned int phi = boost::lexical_cast<unsigned int>(*tokeniter);

            m_deadmap->addDeadChannel(ilayer, eta, phi);
            ++counter;

            if (Verbosity())
            {
              cout << "add dead channel eta" << eta << " phi" << phi;
            }
          }  // if (*tokeniter == "eta")
          else
          {
            if (Verbosity())
            {
              cout << "skip " << deadChanName;
            }
          }

        }  //      if (*tokeniter == "deadtower")
        else
        {
          if (Verbosity())
          {
            cout << "skip " << deadChanName;
          }
        }

      }  //     for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)

      if (Verbosity())
      {
        cout << endl;
      }

    }  //  for (const auto iter = in_par_ranges.first; iter != in_par_ranges.second; ++iter)

    cout << "PHG4SvtxDeadMapLoader::" << m_detector << "::InitRun - loading " << counter << " dead channel for layer "
         << ilayer << " from " << deadMapPath << ". Total dead chan = " << m_deadmap->size() << endl;
  }

  if (Verbosity())
  {
    cout << "PHG4SvtxDeadMapLoader::" << m_detector << "::InitRun - loading dead map completed : ";
    m_deadmap->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
