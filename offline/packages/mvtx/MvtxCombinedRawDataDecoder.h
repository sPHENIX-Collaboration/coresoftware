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
#include <ffarawobjects/MvtxRawEvtHeader.h>

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
    //void removeDuplicates(std::vector<int> &v);
    void removeDuplicates(std::vector<std::pair<uint64_t, uint32_t>> &v);
    //std::vector<int> getAssociatedL1BCOs(int start);
    //bool isInRange(int min, int value, int max);

  private:
    TrkrHitSetContainer* hit_set_container = nullptr;
    MvtxEventInfov1* mvtx_event_header = nullptr;
    MvtxRawEvtHeader* mvtx_raw_event_header = nullptr;
    MvtxRawHitContainer* mvtx_hit_container = nullptr;
    MvtxRawHit* mvtx_hit = nullptr;
 
    /// keep track of number of hits per hitsetid
    //using hitcountmap_t = std::map<TrkrDefs::hitsetkey,int>;
    //hitcountmap_t m_hitcounts;

    std::string m_MvtxRawHitNodeName = "MVTXRAWHIT";
    std::string m_MvtxRawEvtHeaderNodeName = "MVTXRAWEVTHEADER";

    bool m_writeMvtxEventHeader = true;

    //std::vector <int> strobe_bco;
    //std::vector <int> l1_bco;
};

#endif
