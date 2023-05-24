#include "InttRawDataDecoder.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/TrkrHitSet.h>

InttRawDataDecoder::InttRawDataDecoder(std::string const& name):
	SubsysReco(name)
{
	//Do nothing
}

//Init
int InttRawDataDecoder::Init(PHCompositeNode* /*topNode*/)
{
	//Do nothing (yet)
	return Fun4AllReturnCodes::EVENT_OK;
}

//InitRun
int InttRawDataDecoder::InitRun(PHCompositeNode* /*topNode*/)
{
	//Do nothing (yet)
	return Fun4AllReturnCodes::EVENT_OK;
}

//process_event
int InttRawDataDecoder::process_event(PHCompositeNode* topNode)
{
	TrkrHitSetContainer trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
	if(!trkr_hit_set_container)return Fun4AllReturnCodes::DISCARDEVENT;

	Event evt = findNode::getClass<Event>(topNode, "PRDF");
	if(!evt)return Fun4AllReturnCodes::DISCARDEVENT;

	for(int pid = 3001; pid < 3009; ++pid)
	{
		Packet* p = evt->getPacket(pid);
		if(!p)continue;

		int N = p->iValue(0, "NR_HITS");
		if(!N)continue;

		bco_full = p->lValue(0, "BCO");

		int adc = 0;
		int amp = 0;
		int chp = 0;
		int chn = 0;
		int fee = 0;
		int bco = 0;

		TrkrDefs::hitsetkey hit_set_key = 0;
		TrkrDefs::hitkey hit_key = 0;
		TrkrHitSetContainer::Iterator hit_set_container_itr;
		TrkrHit* hit = nullptr;

		for(int n = 0; n < N; ++n)
		{
			adc = p->iValue(i, "ADC");
			amp = p->iValue(i, "AMPLITUE");
			chp = p->iValue(i, "CHIP_ID");
			chn = p->iValue(i, "CHANNEL_ID");
			fee = p->iValue(i, "FEE");
			bco = p->iValue(i, "FPHX_BCO");

			//make hit_set_key
			//...

			hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);

			//make hit_key
			//...

			hit = hit_set_container_itr->second->getHit(hit_key);
			if(hit)continue; //Not a typo; this is to avoid duplicate hits

			hit = new TrkrHit;
			//hit->setAdc(...)
			//...
			hit_setcontainer_itr->second->addHitSpecificKey(hit_key, hit);
		}
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

//End
int InttRawDataDecoder::End(PHCompositeNode* /*topNode*/)
{
	//Do nothing (yet)
	return Fun4AllReturnCodes::EVENT_OK;
}

