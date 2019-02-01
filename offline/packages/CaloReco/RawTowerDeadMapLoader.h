// $Id: $

/*!
 * \file RawTowerDeadMapLoader.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef CALORECO_RAWTOWERDEADMAPLOADER_H
#define CALORECO_RAWTOWERDEADMAPLOADER_H

#include <fun4all/SubsysReco.h>

#include <string>

class RawTowerDeadMap;


/*!
 * \brief RawTowerDeadMapLoader loads dead map at inti run
 */
class RawTowerDeadMapLoader : public SubsysReco
{
 public:
  RawTowerDeadMapLoader(const std::string& detector);

  virtual ~RawTowerDeadMapLoader() {}

  virtual int InitRun(PHCompositeNode *topNode);

  const std::string& deadMapPath() const
  {
    return m_deadMapPath;
  }

  void deadMapPath(const std::string& deadMapPath)
  {
    m_deadMapPath = deadMapPath;
  }

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
  std::string m_deadMapPath;

  RawTowerDeadMap * m_deadmap;
};

#endif
