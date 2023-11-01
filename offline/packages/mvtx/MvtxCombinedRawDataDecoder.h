#ifndef MVTX_MVTXCOMBINEDRAWDATADECODER_H
#define MVTX_MVTXCOMBINEDRAWDATADECODER_H

/*!
 * \file MvtxCombinedRawDataDecoder.h
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include <fun4all/SubsysReco.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/MvtxEventInfov1.h>
#include <trackbase/TrkrHitSetContainerv1.h>

#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitContainer.h>
#include <ffarawobjects/MvtxRawEvtHeaderv1.h>

#include <memory>
#include <string>
#include <map>
#include <vector>

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

   void useRawHitNodeName(const std::string &name){ m_MvtxRawHitNodeName = name; }

   void useRawEvtHeaderNodeName(const std::string &name){ m_MvtxRawEvtHeaderNodeName = name; }

   void writeMvtxEventHeader(bool write){ m_writeMvtxEventHeader = write; }

  protected:
    void removeDuplicates(std::vector<std::pair<uint64_t, uint32_t>> &v);

  private:
    TrkrHitSetContainer* hit_set_container = nullptr;
    MvtxEventInfov1* mvtx_event_header = nullptr;
    MvtxRawEvtHeaderv1* mvtx_raw_event_header = nullptr;
    MvtxRawHitContainer* mvtx_hit_container = nullptr;
    MvtxRawHit* mvtx_hit = nullptr;
 
    std::string m_MvtxRawHitNodeName = "MVTXRAWHIT";
    std::string m_MvtxRawEvtHeaderNodeName = "MVTXRAWEVTHEADER";

    bool m_writeMvtxEventHeader = true;
};

#endif
