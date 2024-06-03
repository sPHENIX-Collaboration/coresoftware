#ifndef MICROMEGAS_MicromegasRawDataTimingEvaluation_H
#define MICROMEGAS_MicromegasRawDataTimingEvaluation_H

/*!
 * \file MicromegasRawDataTimingEvaluation.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"

#include <fun4all/SubsysReco.h>
#include <fun4allraw/MicromegasBcoMatchingInformation.h>
#include <phool/PHObject.h>

#include <TTree.h>

#include <list>
#include <map>
#include <memory>
#include <string>
#include <utility>

class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TProfile;

class Packet;

/// micromegas raw data decoder
class MicromegasRawDataTimingEvaluation : public SubsysReco
{
 public:
  /// constructor
  MicromegasRawDataTimingEvaluation(const std::string& name = "MicromegasRawDataTimingEvaluation");

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// output file name for evaluation histograms
  void set_evaluation_outputfile(const std::string& outputfile) { m_evaluation_filename = outputfile; }

  /**
   * waveform is similar to sample except that there is only one of which per waveform,
   * and that it stores the max adc and corresponding sample_id
   */
  class Waveform
  {
   public:
    /// packet
    unsigned int packet_id = 0;

    /// ll1 bco
    uint64_t gtm_bco = 0;

    /// fee bco
    unsigned int fee_bco = 0;

    /// fee bco predicted (from gtm)
    unsigned int fee_bco_predicted = 0;

    /// fee
    unsigned short fee_id = 0;

    /// channel id
    unsigned short channel = 0;

    using List = std::vector<Waveform>;
  };

  class Container : public PHObject
  {
   public:
    void Reset();

    Waveform::List waveforms;
    ClassDef(Container, 1)
  };

 private:

  //! evaluation output filename
  std::string m_evaluation_filename = "MicromegasRawDataTimingEvaluation.root";
  std::unique_ptr<TFile> m_evaluation_file;

  //! tree
  TTree* m_evaluation_tree = nullptr;

  //! main branch
  Container* m_container = nullptr;

  //! map bco_information_t to packet id
  using bco_matching_information_map_t = std::map<unsigned int, MicromegasBcoMatchingInformation>;
  bco_matching_information_map_t m_bco_matching_information_map;

  //! map waveforms to bco
  /** this is used to count how many waveforms are found for a given lvl1 bco */
  using bco_map_t = std::map<uint64_t, unsigned int>;
  bco_map_t m_bco_map;

  //! keep track of total number of waveforms
  uint64_t m_waveform_count_total = 0;

  //! keep track of dropped waveforms
  uint64_t m_waveform_count_dropped = 0;

};

#endif
