#ifndef MICROMEGAS_MICROMEGASRAWDATADECODER_H
#define MICROMEGAS_MICROMEGASRAWDATADECODER_H

/*!
 * \file MicromegasRawDataDecoder.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>

#include <string>        

class PHCompositeNode;
class TH1;
class TH2;

//! micromegas clusterizer
class MicromegasRawDataDecoder : public SubsysReco
{
 public:

  //! constructor
  MicromegasRawDataDecoder( const std::string &name = "MicromegasRawDataDecoder" );

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  private:  
  static constexpr int m_nchannels_fee = 256;
  static constexpr int m_nfee = 16;
  static constexpr int m_nchannels_total = m_nfee*m_nchannels_fee;
  static constexpr int m_max_adc = 1024;
  
  //!@name evaluation histograms
  //@{
  
  //! Fired FEE
  TH1* m_h_fee = nullptr;
  
  //! ADC distribution vs channel number
  TH2* m_h_adc = nullptr;
  
  //@}
  
};

#endif
