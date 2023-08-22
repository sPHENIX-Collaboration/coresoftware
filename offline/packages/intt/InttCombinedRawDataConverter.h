#ifndef INTT_COMBINED_RAW_DATA_CONVERTER_H
#define INTT_COMBINED_RAW_DATA_CONVERTER_H

#include "InttMapping.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <fun4all/SubsysReco.h>

class InttCombinedRawDataConverter : public SubsysReco
{
public:
	InttCombinedRawDataConverter(std::string const& name = "InttRawDataConverter");

	int SetOutputFile(std::string const&);
	int WriteOutputFile();

	int Init(PHCompositeNode*) override;
	int InitRun(PHCompositeNode*) override;
	int process_event(PHCompositeNode*) override;
	int End(PHCompositeNode*) override;

private:
	std::string m_EvtNodeName = "EVT";

	TFile* file = nullptr;
	TTree* tree = nullptr;

	Int_t n_evt = 0;
	Int_t num_hits = 0;

	Intt::RawData_s raw;
	Intt::Online_s onl;

	typedef std::map<std::string, std::vector<Long64_t>*> Branches_t;
	Branches_t branches;
};

#endif//INTT_COMBINED_RAW_DATA_CONVERTER_H
