#ifndef EVENTDISPLAY_TRACKEREVENTDISPLAY_H
#define EVENTDISPLAY_TRACKEREVENTDISPLAY_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <set>
#include <fstream>

class PHCompositeNode;
class TrkrCluster;

class TrackerEventDisplay : public SubsysReco
{
 public:
  TrackerEventDisplay(const std::string &name = "TRACKEREVENTDISPLAY",
                const std::string &filename = "./TrackerCosmicsDisplay_25926", //Don't include .root
                const std::string &runnumber = "25926",
                const std::string &date = "2023-08-23");
  ~TrackerEventDisplay() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_first_event(int value) { _ievent = value; }
  void makeHitsJson(bool value) { _hit = value; }
  void makeClustersJson(bool value) { _cluster = value; }

 private:
  std::ofstream outdata;

  bool _hit;
  bool _cluster;
  unsigned int _ievent;
  std::string _filename;
  std::string _runnumber;
  std::string _date;

  double AdcClockPeriod = 53.0; // ns

  // output subroutines
  void makeJsonFile(PHCompositeNode *topNode);
};

#endif  // EVENTDISPLAY_TRACKEREVENTDISPLAY_H
