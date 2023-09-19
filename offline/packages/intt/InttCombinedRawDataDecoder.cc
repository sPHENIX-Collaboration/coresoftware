#include "InttCombinedRawDataDecoder.h"

#include <trackbase/InttDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <ffarawobjects/InttRawHitContainer.h>
#include <ffarawobjects/InttRawHit.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TSystem.h>

InttCombinedRawDataDecoder::InttCombinedRawDataDecoder(std::string const& name)
	: SubsysReco(name)
{
	// Do nothing
}

int InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)
{
	if (!topNode)
	{
		std::cout << "InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)" << std::endl;
		std::cout << "\tCould not retrieve topNode; doing nothing" << std::endl;
		exit(1);
		gSystem->Exit(1);

		return 1;
	}

	PHNodeIterator dst_itr(topNode);
	PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
	if (!dst_node)
	{
		if (Verbosity()) std::cout << "InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)" << std::endl;
		if (Verbosity()) std::cout << "\tCould not retrieve dst_node; doing nothing" << std::endl;
		exit(1);
		gSystem->Exit(1);

		return 1;
	}

	PHNodeIterator trkr_itr(dst_node);
	PHCompositeNode* trkr_node = dynamic_cast<PHCompositeNode*>(trkr_itr.findFirst("PHCompositeNode", "TRKR"));
	if (!trkr_node)
	{
		trkr_node = new PHCompositeNode("TRKR");
		dst_node->addNode(trkr_node);
	}

	TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
	if (!trkr_hit_set_container)
	{
		if (Verbosity()) std::cout << "InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)" << std::endl;
		if (Verbosity()) std::cout << "\tMaking TrkrHitSetContainer" << std::endl;

		trkr_hit_set_container = new TrkrHitSetContainerv1;
		PHIODataNode<PHObject>* new_node = new PHIODataNode<PHObject>(trkr_hit_set_container, "TRKR_HITSET", "PHObject");
		trkr_node->addNode(new_node);
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)
{
	TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
	if (!trkr_hit_set_container)
	{
		std::cout << PHWHERE << std::endl;
		std::cout << "InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)" << std::endl;
		std::cout << "Could not get \"TRKR_HITSET\" from Node Tree" << std::endl;
		std::cout << "Exiting" << std::endl;
		gSystem->Exit(1);
		exit(1);

		return Fun4AllReturnCodes::DISCARDEVENT;
	}

  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode,m_InttRawNodeName);
	if(!inttcont)
	{
		std::cout << PHWHERE << std::endl;
		std::cout << "InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)" << std::endl;
		std::cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << std::endl;
		std::cout << "Exiting" << std::endl;

		gSystem->Exit(1);
		exit(1);

	}

//	int adc = 0;
	// int amp = 0;
//	int bco = 0;

	TrkrDefs::hitsetkey hit_set_key = 0;
	TrkrDefs::hitkey hit_key = 0;
	TrkrHitSetContainer::Iterator hit_set_container_itr;
	TrkrHit* hit = nullptr;

	for (unsigned int i=0; i<inttcont->get_nhits(); i++)
	{
          InttRawHit *intthit = inttcont->get_hit(i);
//	  uint64_t gtm_bco = intthit->get_bco();

	  int adc = intthit->get_adc();
	  //amp = intthit->get_amplitude();
	  //int bco = intthit->get_FPHX_BCO();
					
	  // hit_key = InttDefs::genHitKey(offline.strip_y, offline.strip_x); //col, row <trackbase/InttDefs.h>
	  // hit_set_key = InttDefs::genHitSetKey(offline.layer, offline.ladder_z, offline.ladder_phi, bco);
					
	  hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
	  hit = hit_set_container_itr->second->getHit(hit_key);
	  if(hit)continue;
					
	  hit = new TrkrHitv2;
	  hit->setAdc(adc);
	  hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
	}


	return Fun4AllReturnCodes::EVENT_OK;
}

/*
	for (int pktid = 3001; pktid <= 3007; pktid++)
	{
		for (PacketMap::PacketIterator pktiter = pktmap->begin_end(pktid).first;
				 pktiter != pktmap->begin_end(pktid).second;
				 ++pktiter)
		{
			for (PacketMap::BclkIterator bclkitr = pktmap->begin_end_bclk(pktid).first;
					 bclkitr != pktmap->begin_end_bclk(pktid).second;
					 ++bclkitr)
			{
				uint64_t bit40bclk = *bclkitr & 0xFFFFFFFFFF;
				std::cout << "going over packet " << std::endl;	//<< pktiter->getIdentifier()
				std::cout << "bclk: " << *bclkitr << std::endl;
				std::cout << "packet : " << (*pktiter)->getIdentifier() << std::endl;
				//<< " for bclk 0x" << std::hex << bclkitr << std::dec << std::endl;
				int num_hits = (*pktiter)->iValue(0, "NR_HITS");
				for (int j = 0; j < num_hits; j++)
				{
					int FEE = (*pktiter)->iValue(j, "FEE");
					uint64_t gtm_bco = (*pktiter)->lValue(j, "BCO");
					if (bit40bclk == gtm_bco)
					{
						std::cout << "Chosen hit: FEE " << FEE << " bclk 0x" << std::hex
											<< gtm_bco << std::dec << std::endl;

						rawdata = Intt::RawFromPacket(Intt::Packet_Id.find(pktid) != Intt::Packet_Id.end() ? Intt::Packet_Id.find(pktid)->second : -1, j, (*pktiter));
						adc = (*pktiter)->iValue(j, "ADC");
						//amp = p->iValue(j, "AMPLITUE");
						bco = (*pktiter)->iValue(j, "FPHX_BCO");
				
						offline = Intt::ToOffline(rawdata);
				
						hit_key = InttDefs::genHitKey(offline.strip_y, offline.strip_x); //col, row <trackbase/InttDefs.h>
						hit_set_key = InttDefs::genHitSetKey(offline.layer, offline.ladder_z, offline.ladder_phi, bco);
				
						hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
						hit = hit_set_container_itr->second->getHit(hit_key);
						if(hit)continue;
				
						hit = new TrkrHitv2;
						hit->setAdc(adc);
						hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
					}
				}
			}
		}
	}

		Packet* p = evt->getPacket(itr->first);
		if(!p)continue;

		int N = p->iValue(0, "NR_HITS");
		full_bco = p->lValue(0, "BCO");

		if(Verbosity() > 20)std::cout << N << std::endl;

		for(int n = 0; n < N; ++n)
		{
		rawdata = Intt::RawFromPacket(itr->second, n, p);

		adc = p->iValue(n, "ADC");
		//amp = p->iValue(n, "AMPLITUE");
		bco = p->iValue(n, "FPHX_BCO");

		offline = Intt::ToOffline(rawdata);

		hit_key = InttDefs::genHitKey(offline.strip_y, offline.strip_x); //col, row <trackbase/InttDefs.h>
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
		}
		if(Verbosity() > 20)
		{
		std::cout << std::endl;
		std::cout << "Identify():" << std::endl;
		trkr_hit_set_container->identify();
		std::cout << std::endl;
		}
	*/
