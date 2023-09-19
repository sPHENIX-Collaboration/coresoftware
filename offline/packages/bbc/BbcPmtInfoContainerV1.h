#ifndef __BBCPMTINFOCONTAINERV1_H__
#define __BBCPMTINFOCONTAINERV1_H__

#include <calobase/TowerInfoContainer.h>
#include "BbcPmtInfoV1.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>

class BbcPmtInfoContainerV1 : public TowerInfoContainer
{
public:

  BbcPmtInfoContainerV1();

  ~BbcPmtInfoContainerV1() override;

  void Reset() override;

  BbcPmtInfoV1 *get_tower_at_channel(int pos) override;

  BbcPmtInfoV1 *get_pmt(int ich) { return (BbcPmtInfoV1*) _clones->At(ich); }
  TClonesArray *getarray() const { return _clones; }

  size_t size() override { return _clones->GetEntries(); }

  //void AddBbcPmtInfo(const short ipmt, const Float_t q, const Float_t tt, const Float_t tq);

protected:
  TClonesArray *_clones = nullptr;
  DETECTOR _detector = DETECTOR_INVALID;

private:
  ClassDefOverride(BbcPmtInfoContainerV1, 1);
};

#endif
