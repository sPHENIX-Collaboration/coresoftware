#ifndef MICROMEGAS_MicromegasRawDataEvaluation_H
#define MICROMEGAS_MicromegasRawDataEvaluation_H

/*!
 * \file MicromegasRawDataEvaluation.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCalibrationData.h"
#include "MicromegasMapping.h"

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

#include <TTree.h>

#include <map>
#include <memory>
#include <string>
#include <utility>

class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TProfile;

/// micromegas raw data decoder
class MicromegasRawDataEvaluation : public SubsysReco
{
  public:

  /// constructor
  MicromegasRawDataEvaluation( const std::string &name = "MicromegasRawDataEvaluation" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// calibration file
  void set_calibration_file( const std::string& value ) { m_calibration_filename = value; }

  /// set number of RMS sigma used to defined static threshold on a given channel
  void set_n_sigma( double value ) { m_n_sigma = value; }

  /// set minimum ADC value, disregarding pedestal and RMS. 
  /** This removes faulty channels for which calibration has failed */
  void set_min_adc( double value ) { m_min_adc = value; }

  /// set min sample for noise estimation
  void set_sample_min( int value ) { m_sample_min = value; }

  /// set min sample for noise estimation
  void set_sample_max( int value ) { m_sample_max = value; }

  /// output file name for evaluation histograms
  void set_evaluation_outputfile(const std::string &outputfile) {m_evaluation_filename = outputfile;}

  class Sample
  {
    public:
    /// packet
    unsigned int packet_id = 0;

    /// ll1 bco
    uint64_t lvl1_bco = 0;

    /// fee bco
    unsigned int fee_bco = 0;

    /// checksum and checksum error
    unsigned int checksum = 0;
    unsigned int checksum_error = 0;

    /// fee
    unsigned short fee_id = 0;
    unsigned short layer = 0;
    unsigned short tile = 0;

    /// sampa channel and sampa address
    unsigned short sampa_address = 0;
    unsigned short sampa_channel = 0;

    /// channel id
    unsigned short channel = 0;

    /// physical strip id
    unsigned short strip = 0;

    unsigned short sample = 0;
    unsigned short adc = 0;
    
    double pedestal = 0;
    double rms = 0;

    using List = std::vector<Sample>;
  };

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
    uint64_t lvl1_bco = 0;

    /// fee bco
    unsigned int fee_bco = 0;

    /// checksum and checksum error
    unsigned int checksum = 0;
    unsigned int checksum_error = 0;

    /// fee
    unsigned short fee_id = 0;
    unsigned short layer = 0;
    unsigned short tile = 0;

    /// sampa channel and sampa address
    unsigned short sampa_address = 0;
    unsigned short sampa_channel = 0;

    /// channel id
    unsigned short channel = 0;

    /// physical strip id
    unsigned short strip = 0;

    unsigned short sample_max = 0;
    unsigned short adc_max = 0;
    
    double pedestal = 0;
    double rms = 0;

    bool is_signal = false;

    //! default constructor
    Waveform() = default;

    //! construct from sample
    Waveform( const Sample& sample )
    { copy_from( sample ); }

    //! copy from sample
    void copy_from( const Sample& );

    using List = std::vector<Waveform>;
  };


  class Container: public PHObject
  {
    public:
    void Reset();
    
    // number of taggers for each packet
    std::vector<int> n_tagger;
    
    // number of waveform for each packet
    std::vector<int> n_waveform;

    Waveform::List waveforms;
    Sample::List samples;
    
    // bco for this event
    std::vector<uint64_t> lvl1_bco_list;
    
    // lvl1 count for this event
    std::vector<uint32_t> lvl1_count_list;
    
    ClassDef(Container,1)
  };

  private:

  //! calibration filename
  std::string m_calibration_filename = "TPOT_Pedestal_000.root";

  //! calibration data
  MicromegasCalibrationData m_calibration_data;

  //! mapping
  MicromegasMapping m_mapping;

  /// number of RMS sigma used to define threshold
  double m_n_sigma = 5;

  //! minimum ADC value, disregarding pedestal and RMS. 
  /* This removes faulty channels for which calibration has failed */
  double m_min_adc = 50;
  
  /// min sample for signal
  int m_sample_min = 0;

  /// max sample for signal
  int m_sample_max = 100;

  //! evaluation output filename
  std::string m_evaluation_filename = "MicromegasRawDataEvaluation.root";
  std::unique_ptr<TFile> m_evaluation_file;

  //! tree
  TTree* m_evaluation_tree = nullptr;

  //! main branch
  Container* m_container = nullptr;
  
  //! map fee bco to lvl1 bco
  using bco_matching_pair_t = std::pair<unsigned int, uint64_t>;

  //! map fee_id to bco maps
  using fee_bco_matching_map_t = std::map<unsigned short, bco_matching_pair_t>;
  fee_bco_matching_map_t m_fee_bco_matching_map;
  
  /// map waveforms to bco
  /** this is used to count how many waveforms are found for a given lvl1 bco */
  using bco_map_t = std::map<uint64_t,unsigned int>;
  bco_map_t m_bco_map;

};

#endif
