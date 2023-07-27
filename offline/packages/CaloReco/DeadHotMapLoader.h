// $Id: $

/*!
 * \file DeadHotMapLoader.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef CALORECO_DEADHOTMAPLOADER_H
#define CALORECO_DEADHOTMAPLOADER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class CDBTTree;

/*!
 * \brief DeadHotMapLoader loads dead map at inti run
 */
class DeadHotMapLoader : public SubsysReco
{
 public:
  explicit DeadHotMapLoader(const std::string& detector);

  ~DeadHotMapLoader() override {}

  int InitRun(PHCompositeNode* topNode) override;

  const std::string& detector() const
  {
    return m_detector;
  }

  void detector(const std::string& detector)
  {
    m_detector = detector;
  }

 private:
  std::string m_detector;
  CDBTTree *m_CDBTTree = nullptr;
};

#endif
