#include "InttRawDataDecoder.h"
#include "InttMapping.h"

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

	struct Intt::RawData_s rawdata;
	struct Intt::Offline_s offline;

	int adc = 0;
	//int amp = 0;
	int bco = 0;

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
			rawdata.felix_server = pid - 3001;
			rawdata.felix_channel = p->iValue(n, "FEE");
			rawdata.chip = p->iValue(n, "CHIP_ID");
			rawdata.channel = p->iValue(n, "CHANNEL_ID");

			adc = p->iValue(n, "ADC");
			//amp = p->iValue(n, "AMPLITUE");
			bco = p->iValue(n, "FPHX_BCO");

			offline = Intt::ToOffline(rawdata);

			hit_key = InttDefs::genHitKey(offline.strip_x, offline.strip_y);
			hit_set_key = InttDefs::genHitSetKey(offline.layer, offline.ladder_z, offline.ladder_phi, bco);

			hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
			hit = hit_set_container_itr->second->getHit(hit_key);
			if(hit)continue;

			hit = new TrkrHitv2;
			hit->setAdc(adc);
			hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
		}

		delete p;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

//End
int InttRawDataDecoder::End(PHCompositeNode* /*topNode*/)
{
	//Do nothing (yet)
	return Fun4AllReturnCodes::EVENT_OK;
}

