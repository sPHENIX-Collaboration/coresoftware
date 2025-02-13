// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALORECO_CALOGEOMMAPPING_H
#define CALORECO_CALOGEOMMAPPING_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawTowerGeomContainer;

class CaloGeomMapping : public SubsysReco
{
 public:
  CaloGeomMapping(const std::string &name = "CaloGeomMapping");

  ~CaloGeomMapping() override = default;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  // Create tower geometry mapping node
  void CreateGeomNode(PHCompositeNode *topNode);

  void set_detector_name(const std::string &name)
  {
    m_Detector = name;
  }
  const std::string &get_detector_name()
  {
    return m_Detector;
  }

 protected:
  std::string m_Detector;  // CEMC, HCALIN or HCALOUT
  std::string m_TowerGeomNodeName;
  RawTowerGeomContainer *m_RawTowerGeomContainer {nullptr};
};

#endif  // CALOGEOMMAPPING_H
