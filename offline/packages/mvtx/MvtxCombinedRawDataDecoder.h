#ifndef MVTX_MVTXCOMBINEDRAWDATADECODER_H
#define MVTX_MVTXCOMBINEDRAWDATADECODER_H

/*!
 * \file MvtxCombinedRawDataDecoder.h
 * \author Cameron Dean <cameron.dean@cern.ch>
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>
#include <map>
#include <vector>

class MvtxEventInfo;
class MvtxRawEvtHeader;
class MvtxRawHitContainer;
class MvtxRawHit;
class PHCompositeNode;
class TrkrHitSetContainer;

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

  private:
    void removeDuplicates(std::vector<std::pair<uint64_t, uint32_t>> &v);
    TrkrHitSetContainer* hit_set_container = nullptr;
    MvtxEventInfo* mvtx_event_header = nullptr;
    MvtxRawEvtHeader* mvtx_raw_event_header = nullptr;
    MvtxRawHitContainer* mvtx_hit_container = nullptr;
    MvtxRawHit* mvtx_hit = nullptr;
 
    std::string m_MvtxRawHitNodeName = "MVTXRAWHIT";
    std::string m_MvtxRawEvtHeaderNodeName = "MVTXRAWEVTHEADER";

    bool m_writeMvtxEventHeader = true;
};

#endif
