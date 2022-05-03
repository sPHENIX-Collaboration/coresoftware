// $Id: $

/*!
 * \file PHG4InttDeadMapLoader.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4INTT_PHG4INTTDEADMAPLOADER_H
#define G4INTT_PHG4INTTDEADMAPLOADER_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;

/*!
 * \brief PHG4InttDeadMapLoader loads dead map at inti run
 */
class PHG4InttDeadMapLoader : public SubsysReco
{
 public:
  explicit PHG4InttDeadMapLoader(const std::string& detector = "SILICON_TRACKER");

  ~PHG4InttDeadMapLoader() override {}

  int InitRun(PHCompositeNode* topNode) override;

  void deadMapPath(unsigned int layer, const std::string& deadMapPath)
  {
    m_deadMapPathMap[layer] = deadMapPath;
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
  std::map<unsigned int, std::string> m_deadMapPathMap;

  std::string m_detector;
};

#endif /* G4INTT_PHG4INTTDeadMapLoader_H */
