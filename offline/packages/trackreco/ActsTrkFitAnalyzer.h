#ifndef ACTSTRKFITANALYZER_H
#define ACTSTRKFITANALYZER_H

#include <fun4all/SubsysReco.h>

class TTree;
class TFile;

#include <string>

class ActsTrkFitAnalyzer : public SubsysReco
{

 public:
  ActsTrkFitAnalyzer(const std::string& name = "ActsTrkFitAnalyzer.root");
  ~ActsTrkFitAnalyzer();
  
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 private:

  void initializeTree();

  TFile *m_trackFile;
  TTree *m_trackTree;
  int m_event;
  int m_traj;
  

};

#endif
