#ifndef INTT_RAW_DATA_CONVERTER_H
#define INTT_RAW_DATA_CONVERTER_H

#include "InttMapping.h"

#include <fun4all/SubsysReco.h>

#include <Rtypes.h>

#include <string>
#include <map>

class PHCompositeNode;
class TFile;
class TTree;

class InttRawDataConverter : public SubsysReco
{
public:
	InttRawDataConverter(std::string const& name = "InttRawDataConverter");

	int SetOutputFile(std::string const&);
	int WriteOutputFile();

	int Init(PHCompositeNode*) override;
	int InitRun(PHCompositeNode*) override;
	int process_event(PHCompositeNode*) override;
	int End(PHCompositeNode*) override;

private:
	TFile* file = nullptr;
	TTree* tree = nullptr;

	Int_t n_evt = -1;
	Int_t num_hits = 0;
	Long64_t gtm_bco = 0;
	Int_t flx_svr = 0;

	InttNameSpace::RawData_s raw;
	InttNameSpace::Online_s onl;

	typedef std::map<std::string, Int_t*> Branches_t;
	Branches_t branches;
};

#endif//INTT_RAW_DATA_CONVERTER_H
