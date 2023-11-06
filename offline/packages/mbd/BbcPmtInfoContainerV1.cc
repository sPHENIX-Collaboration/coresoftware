#include "BbcPmtInfoContainerV1.h"
#include "BbcPmtInfoV1.h"
#include "BbcDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <TClonesArray.h>

#include <cassert>


BbcPmtInfoContainerV1::BbcPmtInfoContainerV1()
  : _detector( DETECTOR::MBD )
{
  const int nchannels = BbcDefs::BBC_N_PMT;
  _clones = new TClonesArray("BbcPmtInfoV1", nchannels);
  _clones->SetOwner();
  _clones->SetName("BbcPmtInfoContainerV1");
  for (int ipmt = 0; ipmt < nchannels; ipmt++)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt( ipmt );
  }
}

BbcPmtInfoContainerV1::~BbcPmtInfoContainerV1()
{
  delete _clones;
}

void BbcPmtInfoContainerV1::Reset()
{
  // clear content of towers in the container for the next event
  //_clones->Clear();
  for (int ipmt = 0; ipmt < BbcDefs::BBC_N_PMT; ipmt++)
  {
    _clones->ConstructedAt( ipmt )->Clear();
  }
}

BbcPmtInfoV1* BbcPmtInfoContainerV1::get_tower_at_channel(int pos)
{
  return (BbcPmtInfoV1*)_clones->ConstructedAt( pos );
}

