#ifndef PHG4FullProjSpacalCellReco_H
#define PHG4FullProjSpacalCellReco_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>
#include <vector>

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4FullProjSpacalCellReco : public SubsysReco
{
 public:

  PHG4FullProjSpacalCellReco(const std::string &name = "HCALCELLRECO");

  virtual ~PHG4FullProjSpacalCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d) {detector = d;}

  void checkenergy(const int i=1) {chkenergyconservation = i;}

 protected:

  int CheckEnergy(PHCompositeNode *topNode);


  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;

  PHTimeServer::timer _timer;
  int chkenergyconservation;
  std::map<unsigned int, PHG4CylinderCell *> celllist;
};

#endif
