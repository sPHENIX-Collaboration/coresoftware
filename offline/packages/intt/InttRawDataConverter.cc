#include "InttRawDataConverter.h"

#include <Event/Event.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TString.h> // for Form
#include <TTree.h>

#include <iostream>  // for operator<<, endl, basic_ostream
#include <utility>   // for pair

class PHCompositeNode;

InttRawDataConverter::InttRawDataConverter(std::string const& name):
	SubsysReco(name)
{
	//Do nothing
}

int InttRawDataConverter::SetOutputFile(std::string const& filename)
{
	if(filename.empty())
	{
		std::cout << "InttRawDataConverter::SetOputputFile(std:: string const& filename)" << std::endl;
		std::cout << "Argument \"filename\" is empty string" << std::endl;
		std::cout << "No output was written" << std::endl;

		return 1;
	}

	std::cout << "Will write to file:" << std::endl;
	std::cout << "\t" << filename << std::endl;

	file = TFile::Open(filename.c_str(), "RECREATE");
	if(tree)tree->SetDirectory(file);

	return 0;
}

int InttRawDataConverter::WriteOutputFile()
{
	if(!file)
	{
		std::cout << "InttRawDataConverter::WriteOputputFile()" << std::endl;
		std::cout << "Member \"file\" is uninitialized" << std::endl;
		std::cout << "Did you call SetOutputFile()?" << std::endl;
		std::cout << "No output was written" << std::endl;
		return 1;
	}

	if(!tree)
	{
		std::cout << "InttRawDataConverter::WriteOputputFile()" << std::endl;
		std::cout << "Member \"tree\" is uninitialized" << std::endl;
		std::cout << "Did you call SetOutputFile()?" << std::endl;
		std::cout << "No output was written" << std::endl;
		return 1;
	}

	file->cd();
	tree->Write();
	file->Write();
	file->Close();
	delete tree;

	return 0;
}

int InttRawDataConverter::Init(PHCompositeNode* /*topNode*/)
{
	if(tree)delete tree;
	tree = new TTree("prdf_tree", "prdf_tree");
	if(file)tree->SetDirectory(file);

	tree->Branch("n_evt", &n_evt);
	tree->Branch("num_hits", &num_hits);
	tree->Branch("gtm_bco", &gtm_bco);
	tree->Branch("flx_svr", &flx_svr);

	branches =
	{
		{"flx_chn",	nullptr},
		{"lyr",		nullptr},
		{"ldr",		nullptr},
		{"arm",		nullptr},
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
	n_evt = -1;

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

		if(Verbosity() > 20)std::cout << num_hits << std::endl;

		++n_evt;// = pkt->lValue(0, "");
		gtm_bco = pkt->lValue(0, "BCO");
		flx_svr = pkt_itr->second;

		for(Branches_t::iterator itr = branches.begin(); itr != branches.end(); ++itr)
		{
			itr->second = new Int_t[num_hits];
			tree->GetBranch(itr->first.c_str())->SetAddress(itr->second);
		}
	
		for(int n = 0; n < num_hits; ++n)
		{
			raw = Intt::RawFromPacket(pkt_itr->second, n, pkt);
			branches["flx_chn"][n] = raw.felix_channel;

			onl = Intt::ToOnline(raw);
			branches["lyr"][n] = onl.lyr;
			branches["ldr"][n] = onl.ldr;
			branches["arm"][n] = onl.arm;
			branches["chp"][n] = onl.chp;
			branches["chn"][n] = onl.chn;

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
