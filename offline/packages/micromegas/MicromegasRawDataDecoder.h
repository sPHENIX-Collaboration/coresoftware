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

};

#endif
