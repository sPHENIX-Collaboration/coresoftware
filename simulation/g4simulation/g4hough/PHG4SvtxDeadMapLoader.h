// $Id: $

/*!
 * \file PHG4SvtxDeadMapLoader.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4CEMC_PHG4SvtxDeadMapLoader_H_
#define SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4CEMC_PHG4SvtxDeadMapLoader_H_

#include <fun4all/SubsysReco.h>
#include <string>
#include <map>

class SvtxDeadMap;


/*!
 * \brief PHG4SvtxDeadMapLoader loads dead map at inti run
 */
class PHG4SvtxDeadMapLoader : public SubsysReco
{
 public:
  explicit PHG4SvtxDeadMapLoader(const std::string& detector = "SILICON_TRACKER");

  virtual ~PHG4SvtxDeadMapLoader();

  virtual int InitRun(PHCompositeNode *topNode);

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
  SvtxDeadMap * m_deadmap;
};

#endif /* SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4CEMC_PHG4SvtxDeadMapLoader_H_ */
