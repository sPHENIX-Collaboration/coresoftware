#ifndef MVTX_MVTXCOMBINEDRAWDATADECODER_H
#define MVTX_MVTXCOMBINEDRAWDATADECODER_H

/*!
 * \file MvtxCombinedRawDataDecoder.h
 * \author Cameron Dean <cameron.dean@cern.ch>
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "MvtxPixelMask.h"

class MvtxEventInfo;
class MvtxRawEvtHeader;
class MvtxRawHitContainer;
class MvtxRawHit;
class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrHitSetContMvtxHelper;

/// mvtx raw data decoder
class MvtxCombinedRawDataDecoder : public SubsysReco
{
 public:
  /// constructor
  explicit MvtxCombinedRawDataDecoder(const std::string& name = "MvtxCombinedRawDataDecoder");

  /// global initialization
  int Init(PHCompositeNode* /*dummy*/) override;

  /// run initialization
  int InitRun(PHCompositeNode* /*dummy*/) override;

  /// event processing
  int process_event(PHCompositeNode* /*dummy*/) override;

  /// end of processing
  int End(PHCompositeNode* /*dummy*/) override;

  void useRawHitNodeName(const std::string& name) { m_MvtxRawHitNodeName = name; }

  void useRawEvtHeaderNodeName(const std::string& name) { m_MvtxRawEvtHeaderNodeName = name; }

  void doOfflineMasking(bool do_masking) { m_doOfflineMasking = do_masking; }

  void runMvtxTriggered(bool b = true) { m_mvtx_is_triggered = b; }

  void SetReadStrWidthFromDB(const bool val) { m_readStrWidthFromDB = val; }
  bool GetReadStrWidthFromDB() const { return m_readStrWidthFromDB; }
  void SetStrobeWidth(const float val) { m_strobeWidth = val; }
  float GetStrobeWidth() const { return m_strobeWidth; }

 private:
  // static void removeDuplicates(std::vector<std::pair<uint64_t, uint32_t>>& v);
  void CreateNodes(PHCompositeNode*);
  void GetNodes(PHCompositeNode*);

  uint64_t gl1rawhitbco = 0;

  TrkrHitSetContainer* hit_set_container = nullptr;
  TrkrHitSetContMvtxHelper* mvtx_hit_set_helper = nullptr;
  MvtxEventInfo* mvtx_event_header = nullptr;
  MvtxRawEvtHeader* mvtx_raw_event_header = nullptr;
  MvtxRawHitContainer* mvtx_raw_hit_container = nullptr;
  MvtxRawHit* mvtx_rawhit = nullptr;

  std::string m_MvtxRawHitNodeName = "MVTXRAWHIT";
  std::string m_MvtxRawEvtHeaderNodeName = "MVTXRAWEVTHEADER";

  bool m_readStrWidthFromDB = true;
  float m_strobeWidth = 89.;  //! microseconds

  // mask hot pixels
  bool m_doOfflineMasking{false};
  MvtxPixelMask* m_hot_pixel_mask{nullptr};

  bool m_mvtx_is_triggered{false};
};

#endif
