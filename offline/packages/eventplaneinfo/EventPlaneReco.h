// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANERECO_H
#define EVENTPLANERECO_H

//===========================================================
/// \author Ejiro Umaka
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>  // for string
#include <vector>  // for vector

class PHCompositeNode;

class EventPlaneReco : public SubsysReco
{
 public:
  EventPlaneReco(const std::string &name = "EventPlaneReco");
  ~EventPlaneReco() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void ResetMe();


 private:
  int CreateNodes(PHCompositeNode *topNode);
  double GetPsiInRange(const double Qx, const double Qy, const unsigned int order) const;
    
  unsigned int m_MaxOrder = 0;
    
  std::vector<double> south_norm;
  std::vector<std::vector<double>> south_qraw;
  std::vector<double> south_rawpsi;
  
  std::vector<double> north_norm;
  std::vector<std::vector<double>> north_qraw;
  std::vector<double> north_rawpsi;
  
};

#endif  // EVENTPLANERECO_H
