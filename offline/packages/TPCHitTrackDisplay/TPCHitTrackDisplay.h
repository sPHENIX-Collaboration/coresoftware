#ifndef TPCHITTRACKDISPLAY_TPCHITTRACKDISPLAY_H
#define TPCHITTRACKDISPLAY_TPCHITTRACKDISPLAY_H

#include <fun4all/SubsysReco.h>

#include <string>

// Forward declarations
class PHCompositeNode;

// Writes json file to be used to display an event with:
// https://www.sphenix.bnl.gov/edisplay/
class TPCHitTrackDisplay : public SubsysReco
{
 public:
  // Default constructor
  TPCHitTrackDisplay(const std::string &name = "TPCHitTrackDisplay" /*, bool &tpcRaw=True*/);

  // Process Event, called for each event
  int process_event(PHCompositeNode *) override;

  void set_pdgcode(const int thispdgcode) { _pdgcode = thispdgcode; }

  // set the ADC cut for displaying trackless clusters
  void setCutADC(float value) { m_cut_ADC = value; }

  // Boolean for whether or not to include clusters without an associted track above a certain ADC value
  void setIncludeTracklessClusters(float value) { m_trackless_clusters = value; }

 private:
  // User modules
  void SimulationOut(PHCompositeNode *);

  float m_cut_ADC;
  bool m_trackless_clusters;

  // Event counter
  int _event{0};
  int _pdgcode{0};
  std::string _fileName;

};

#endif  //* TPCHITTRACKDISPLAY_TPCHITTRACKDISPLAY_H *//
