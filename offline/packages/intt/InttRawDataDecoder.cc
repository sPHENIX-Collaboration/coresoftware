#include "InttRawDataDecoder.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/InttDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

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
	TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
	if(!trkr_hit_set_container)return Fun4AllReturnCodes::DISCARDEVENT;

	Event* evt = findNode::getClass<Event>(topNode, "PRDF");
	if(!evt)return Fun4AllReturnCodes::DISCARDEVENT;

	int adc = 0;
	//int amp = 0;
	int bco = 0;
	int chp = 0;
	int chn = 0;
	int fee = 0;

	struct INTT_Felix::Ladder_s ldr_struct;
	int layer = 0;
	int ladder_z = 0;
	int ladder_phi = 0;
	int arm = 0;
	int strip_x = 0;
	int strip_y = 0;

	TrkrDefs::hitsetkey hit_set_key = 0;
	TrkrDefs::hitkey hit_key = 0;
	TrkrHitSetContainer::Iterator hit_set_container_itr;
	TrkrHit* hit = nullptr;

	for(int pid = 3001; pid < 3009; ++pid)
	{
		Packet* p = evt->getPacket(pid);
		if(!p)continue;

		int N = p->iValue(0, "NR_HITS");
		if(!N)continue;

		full_bco = p->lValue(0, "BCO");

		for(int n = 0; n < N; ++n)
		{
			adc = p->iValue(n, "ADC");
			//amp = p->iValue(n, "AMPLITUE");
			bco = p->iValue(n, "FPHX_BCO");
			chp = p->iValue(n, "CHIP_ID");
			chn = p->iValue(n, "CHANNEL_ID");
			fee = p->iValue(n, "FEE");

			INTT_Felix::FelixMap(pid - 3001, fee, ldr_struct);
			layer = 2 * ldr_struct.barrel + ldr_struct.ladder;
			ladder_phi = ldr_struct.ladder;				//        B  A  A  B
			ladder_z = arm * 2 + (chp % 13 < 5);			//South<- 1, 0, 2, 3 ->North
			arm = (pid - 3001) / 4;
			strip_x = 25 * arm - (2 * arm - 1) * chp % 13;
			strip_y = 128 * ((arm + chp / 13) % 2) + chn; //need to check this is the convention

			hit_key = InttDefs::genHitKey(strip_x, strip_y);
			hit_set_key = InttDefs::genHitSetKey(layer, ladder_z, ladder_phi, bco);

			hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
			hit = hit_set_container_itr->second->getHit(hit_key);
			if(hit)continue;

			hit = new TrkrHitv2;
			hit->setAdc(adc);
			hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
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

