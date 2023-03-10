#ifndef MICROMEGAS_MICROMEGASRAWDATADECODER_H
#define MICROMEGAS_MICROMEGASRAWDATADECODER_H

/*!
 * \file MicromegasRawDataDecoder.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class PHCompositeNode;
class TFile;
class TH1;
class TH2;

//! micromegas raw data decoder
class MicromegasRawDataDecoder : public SubsysReco
{
  public:

  //! constructor
  MicromegasRawDataDecoder( const std::string &name = "MicromegasRawDataDecoder" );

  //! global initialization
  int Init(PHCompositeNode*) override;

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// set to true to store evaluation histograms and ntuples
  void set_save_histograms( bool value ) { m_savehistograms = value; }

  /// output file name for evaluation histograms
  void set_histogram_outputfile(const std::string &outputfile) {m_histogramfilename = outputfile;}

  private:

  /// create evaluation histograms
  void create_histograms();

  static constexpr int m_nchannels_fee = 256;
  static constexpr int m_nfee = 8;
  static constexpr int m_nchannels_total = m_nfee*m_nchannels_fee;
  static constexpr int m_max_adc = 1024;

  //!@name evaluation histograms
  //@{

  /// Output root histograms
  bool m_savehistograms = true;

  /// histogram output file name
  std::string m_histogramfilename = "MicromegasRawDataDecoder.root";
  std::unique_ptr<TFile> m_histogramfile;

  //! Fired FEE
  TH1* m_h_fee_id = nullptr;

  //! ADC distribution vs channel number
  TH2* m_h_adc_channel = nullptr;

  //@}

};

#endif
