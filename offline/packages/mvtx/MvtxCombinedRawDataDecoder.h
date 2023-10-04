#ifndef MVTX_MVTXCOMBINEDRAWDATADECODER_H
#define MVTX_MVTXCOMBINEDRAWDATADECODER_H

/*!
 * \file MvtxCombinedRawDataDecoder.h
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include <fun4all/SubsysReco.h>
#include <trackbase/MvtxDefs.h>

#include <memory>
#include <string>
#include <map>

class PHCompositeNode;

/// mvtx raw data decoder
class MvtxCombinedRawDataDecoder : public SubsysReco
{
  public:

  /// constructor
  MvtxCombinedRawDataDecoder( const std::string &name = "MvtxCombinedRawDataDecoder" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;
  
  /// end of processing
  int End(PHCompositeNode*) override;

  void useRawHitNodeName(const std::string &name){ m_MvtxRawNodeName = name; }

  private:
  /// keep track of number of hits per hitsetid
  using hitcountmap_t = std::map<TrkrDefs::hitsetkey,int>;
  hitcountmap_t m_hitcounts;

  std::string m_MvtxRawNodeName = "MVTXRAWHIT";
  
};

#endif
