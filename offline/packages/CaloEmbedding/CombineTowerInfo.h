#ifndef COMBINETOWERINFO_H
#define COMBINETOWERINFO_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TowerInfoContainer;

class CombineTowerInfo : public SubsysReco
{
 public:
  explicit CombineTowerInfo(const std::string& name = "CombineTowerInfo");
  ~CombineTowerInfo() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  void set_inputNodeA(const std::string& name) { m_inputNodeA = name; }
  void set_inputNodeB(const std::string& name) { m_inputNodeB = name; }
  void set_outputNode(const std::string& name) { m_outputNode = name; }
  void set_detector(const std::string& name) { m_detector = name; }

 private:
  void CreateNodes(PHCompositeNode* topNode);

  std::string m_inputNodeA;
  std::string m_inputNodeB;
  std::string m_outputNode;
  std::string m_detector;

  TowerInfoContainer* m_towersA{nullptr};
  TowerInfoContainer* m_towersB{nullptr};
  TowerInfoContainer* m_towersOut{nullptr};
};

#endif

