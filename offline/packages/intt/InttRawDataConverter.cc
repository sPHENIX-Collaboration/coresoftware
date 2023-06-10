#include "InttRawDataConverter.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

InttRawDataConverter::InttRawDataConverter(std::string const& name):
	SubsysReco(name)
{
	//Do nothing
}

int InttRawDataConverter::WriteOutputFile(std::string const& filename)
{
	if(!tree)
	{
		std::cout << "InttRawDataConverter::WriteOputputFile(std:: string const& filename)" << std::endl;
		std::cout << "Member \"tree\" is uninitialized" << std::endl;
		std::cout << "No output was written" << std::endl;
		return 1;
	}

	if(filename.empty())
	{
		std::cout << "InttRawDataConverter::WriteOputputFile(std:: string const& filename)" << std::endl;
		std::cout << "Argument \"filename\" is empty string" << std::endl;
		std::cout << "No output was written" << std::endl;

		return 1;
	}

	std::cout << "Writing to file:" << std::endl;
	std::cout << "\t" << filename << std::endl;

	TFile* file = TFile::Open(filename.c_str(), "RECREATE");
	tree->SetDirectory(file);
	tree->Write();
	file->Write();
	file->Close();

	return 0;
}

int InttRawDataConverter::Init(PHCompositeNode* /*topNode*/)
{
	tree = new TTree("prdf_tree", "prdf_tree");

	tree->Branch("num_hits", &num_hits);
	tree->Branch("gtm_bco", &gtm_bco);
	tree->Branch("flx_svr", &flx_svr);

	branches =
	{
		{"flx_chn",	nullptr},
		{"chp",		nullptr},
		{"chn",		nullptr},

		{"flx_bco",	nullptr},
		{"adc",		nullptr},
		{"amp",		nullptr},
	};

	for(Branches_t::const_iterator itr = branches.begin(); itr != branches.end(); ++itr)
	{
		tree->Branch(itr->first.c_str(), itr->second, Form("%s[num_hits]/I", itr->first.c_str()));
	}

	return 0;
}

int InttRawDataConverter::InitRun(PHCompositeNode* /*topNode*/)
{
	return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawDataConverter::process_event(PHCompositeNode* topNode)
{
	if(!tree)return Fun4AllReturnCodes::DISCARDEVENT;

	Event* evt = findNode::getClass<Event>(topNode, "PRDF");
	if(!evt)return Fun4AllReturnCodes::DISCARDEVENT;

	for(std::map<int, int>::const_iterator pkt_itr = Intt::Packet_Id.begin(); pkt_itr != Intt::Packet_Id.end(); ++pkt_itr)
	{
		Packet* pkt = evt->getPacket(pkt_itr->first);
		if(!pkt)continue;

		num_hits = pkt->iValue(0, "NR_HITS");
		if(!num_hits)
		{
			delete pkt;
			return Fun4AllReturnCodes::DISCARDEVENT;
		}

		if(Verbosity() > 20)std::cout << num_hits << std::endl;

		gtm_bco = pkt->lValue(0, "BCO");
		flx_svr = pkt_itr->second;

		for(Branches_t::iterator itr = branches.begin(); itr != branches.end(); ++itr)
		{
			itr->second = new Int_t[num_hits];
			tree->GetBranch(itr->first.c_str())->SetAddress(itr->second);
		}
	
		for(int n = 0; n < num_hits; ++n)
		{
			branches["flx_chn"][n] = pkt->iValue(n, "FEE");
			branches["chp"][n] = pkt->iValue(n, "CHIP_ID") % 26;
			branches["chn"][n] = pkt->iValue(n, "CHANNEL_ID");

			branches["flx_bco"][n] = pkt->iValue(n, "FPHX_BCO");
			branches["adc"][n] = pkt->iValue(n, "ADC");
			branches["amp"][n] = pkt->iValue(n, "AMPLITUDE");
		}
	
		tree->Fill();
		for(Branches_t::iterator itr = branches.begin(); itr != branches.end(); ++itr)
		{
			delete[] itr->second;
		}

		delete pkt;
	}	
	
	return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawDataConverter::End(PHCompositeNode* /*topNode*/)
{
	return Fun4AllReturnCodes::EVENT_OK;
}
