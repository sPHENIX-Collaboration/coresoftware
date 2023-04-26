#ifndef MICROMEGAS_MICROMEGASRAWDATACALIBRATION_H
#define MICROMEGAS_MICROMEGASRAWDATACALIBRATION_H

/*!
 * \file MicromegasRawDataCalibration.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */
#include <fun4all/SubsysReco.h>

#include <map>
#include <memory>
#include <string>

class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TProfile;

/// micromegas raw data decoder
class MicromegasRawDataCalibration : public SubsysReco
{
  public:

  /// constructor
  MicromegasRawDataCalibration( const std::string &name = "MicromegasRawDataCalibration" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// set min sample for noise estimation
  void set_sample_min( int value ) { m_sample_min = value; }
  
  /// set min sample for noise estimation
  void set_sample_max( int value ) { m_sample_max = value; }

  /// set to true to store evaluation histograms and ntuples
  void set_calibration_file( const std::string& value ) { m_calibration_filename = value; }
  
  /// set to true to store evaluation histograms and ntuples
  void set_save_histograms( bool value ) { m_savehistograms = value; }

  /// output file name for evaluation histograms
  void set_histogram_outputfile(const std::string &outputfile) {m_histogramfilename = outputfile;}

  private:

  /// create evaluation histograms
  void create_histograms();
  
  /// min sample for noise estimation 
  int m_sample_min = 0;
  
  /// max sample for noise estimation
  int m_sample_max = 100;
  
  /// calibration output file
  std::string m_calibration_filename = "TPOT_Pedestal_000.root";
    
  /// map fee id to Profile histogram
  using profile_map_t = std::map<int, TProfile*>;
  profile_map_t m_profile_map;
  
  ///@name evaluation histograms
  //@{

  /// Output root histograms
  bool m_savehistograms = true;

  /// histogram output file name
  std::string m_histogramfilename = "MicromegasRawDataCalibration.root";
  std::unique_ptr<TFile> m_histogramfile;

  /// Fired FEE
  TH1* m_h_fee_id = nullptr;

  /// ADC distribution vs channel number
  TH2* m_h_adc_channel = nullptr;

  //@}

};

#endif
