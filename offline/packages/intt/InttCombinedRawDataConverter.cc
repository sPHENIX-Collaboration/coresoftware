#include "InttCombinedRawDataConverter.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4allraw/PacketMap.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

InttCombinedRawDataConverter::InttCombinedRawDataConverter(std::string const& name):
	SubsysReco(name)
{
	//Do nothing
}

int InttCombinedRawDataConverter::SetOutputFile(std::string const& filename)
{
	if(filename.empty())
	{
		std::cout << "InttCombinedRawDataConverter::SetOputputFile(std:: string const& filename)" << std::endl;
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

int InttCombinedRawDataConverter::WriteOutputFile()
{
	if(!file)
	{
		std::cout << "InttCombinedRawDataConverter::WriteOputputFile()" << std::endl;
		std::cout << "Member \"file\" is uninitialized" << std::endl;
		std::cout << "Did you call SetOutputFile()?" << std::endl;
		std::cout << "No output was written" << std::endl;
		return 1;
	}

	if(!tree)
	{
		std::cout << "InttCombinedRawDataConverter::WriteOputputFile()" << std::endl;
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

int InttCombinedRawDataConverter::Init(PHCompositeNode* /*topNode*/)
{
	tree = new TTree("prdf_tree", "prdf_tree");
	if(file)tree->SetDirectory(file);

	tree->Branch("n_evt", &n_evt);
	tree->Branch("num_hits", &num_hits);

	branches =
	{
		{"flx_svr",	new std::vector<Long64_t>()},
		{"flx_chn",	new std::vector<Long64_t>()},
		{"lyr",		new std::vector<Long64_t>()},
		{"ldr",		new std::vector<Long64_t>()},
		{"arm",		new std::vector<Long64_t>()},
		{"chp",		new std::vector<Long64_t>()},
		{"chn",		new std::vector<Long64_t>()},

		{"gtm_bco",	new std::vector<Long64_t>()},
		{"flx_bco",	new std::vector<Long64_t>()},
		{"adc",		new std::vector<Long64_t>()},
		{"amp",		new std::vector<Long64_t>()},
	};

	for(Branches_t::iterator itr = branches.begin(); itr != branches.end(); ++itr)
	{
		tree->Branch(itr->first.c_str(), &(itr->second));
	}

	if(Verbosity() > 20)std::cout << "int InttCombinedRawDataConverter::Init(PHCompositeNode* /*topNode*/)" << std::endl;
	if(Verbosity() > 20)std::cout << "\tDone";

	return 0;
}

int InttCombinedRawDataConverter::InitRun(PHCompositeNode* /*topNode*/)
{
	n_evt = 0;
	for(Branches_t::iterator itr = branches.begin(); itr != branches.end(); ++itr)itr->second->clear();

	if(Verbosity() > 20)std::cout << "int InttCombinedRawDataConverter::InitRun(PHCompositeNode* /*topNode*/)" << std::endl;
	if(Verbosity() > 20)std::cout << "\tDone";

	return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataConverter::process_event(PHCompositeNode* topNode)
{
	if(!tree)return Fun4AllReturnCodes::DISCARDEVENT;

	if(Verbosity() > 20)std::cout << "int InttCombinedRawDataConverter::process_event(PHCompositeNode* topNode)" << std::endl;
	if(Verbosity() > 20)std::cout << "\t" << n_evt << std::endl;

	PacketMap* pktmap = findNode::getClass<PacketMap>(topNode, m_EvtNodeName);
	if(!pktmap)
	{
		std::cout << "int InttCombinedRawDataConverter::process_event(PHCompositeNode* topNode)" << std::endl;
		std::cout << "\t" << PHWHERE << " could not find node " << m_EvtNodeName << std::endl;
		std::cout << "\tExiting" << std::endl; 

		gSystem->Exit(1);
		exit(1);
	}

	for(Branches_t::iterator itr = branches.begin(); itr != branches.end(); ++itr)itr->second->clear();

	PacketMap::PacketListRange pktrange = pktmap->first_last_packet();
	for(auto iter = pktrange.first; iter != pktrange.second; iter++)
	{
		for(auto bclkiter : iter->second.m_BeamClockSet)
		{
			for(auto pktiter : iter->second.m_PacketVector)
			{
				int num_hits = pktiter->iValue(0, "NR_HITS");
				for(int j = 0; j < num_hits; j++)
				{
					uint64_t gtm_bco = pktiter->lValue(j, "BCO");

					if((bclkiter & 0xFFFFFFFFFF) != gtm_bco)continue;

					if(Verbosity() > 40)std::cout << "\tfound match at bclk 0x" << std::hex << gtm_bco << std::dec << std::endl;

					raw = Intt::RawFromPacket(iter->first, j, pktiter);

					branches["flx_svr"]->push_back(raw.felix_channel);
					branches["flx_chn"]->push_back(raw.felix_channel);

					onl = Intt::ToOnline(raw);
					branches["lyr"]->push_back(onl.lyr);
					branches["ldr"]->push_back(onl.ldr);
					branches["arm"]->push_back(onl.arm);
					branches["chp"]->push_back(onl.chp);
					branches["chn"]->push_back(onl.chn);

					branches["gtm_bco"]->push_back(gtm_bco);
					branches["flx_bco"]->push_back(pktiter->iValue(j, "FPHX_BCO"));
					branches["adc"]->push_back(pktiter->iValue(j, "ADC"));
					branches["amp"]->push_back(pktiter->iValue(j, "AMPLITUDE"));
				}
			}
		}
	}

	num_hits = branches.begin()->second->size();
	tree->Fill();
	++n_evt;

	return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataConverter::End(PHCompositeNode* /*topNode*/)
{
	return Fun4AllReturnCodes::EVENT_OK;
}
