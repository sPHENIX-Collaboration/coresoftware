#ifndef MVTX_MVTXRAWDATADECODER_H
#define MVTX_MVTXRAWDATADECODER_H

/*!
 * \file MvtxRawDataDecoder.h
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include <fun4all/SubsysReco.h>
#include <trackbase/MvtxDefs.h>

#include <memory>
#include <string>
#include <map>

class PHCompositeNode;

/// mvtx raw data decoder
class MvtxRawDataDecoder : public SubsysReco
{
  public:

  /// constructor
  MvtxRawDataDecoder( const std::string &name = "MvtxRawDataDecoder" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;
  
  /// end of processing
  int End(PHCompositeNode*) override;

  private:
  /// keep track of number of hits per hitsetid
  using hitcountmap_t = std::map<TrkrDefs::hitsetkey,int>;
  hitcountmap_t m_hitcounts;
  
};

#endif
