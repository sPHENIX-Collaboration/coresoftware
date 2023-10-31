#ifndef INTT_RAW_DATA_DECODER_H
#define INTT_RAW_DATA_DECODER_H

//////////////////////////////////////////////////////////////////////////////
//	Shamelessly imitated from the micromegas version:
//
//	coresoftware/offline/packages/micromegas/MicromegasRawDataDecoder.h
//
//	Thank you Hugo

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>

class PHCompositeNode;

class InttRawDataDecoder : public SubsysReco
{
public:
	InttRawDataDecoder(std::string const& name = "InttRawDataDecoder");

	int Init(PHCompositeNode*) override;
	int InitRun(PHCompositeNode*) override;
	int process_event(PHCompositeNode*) override;
	int End(PHCompositeNode*) override;

private:
	int64_t full_bco = 0;
};

#endif//INTT_RAW_DATA_DECODER_H
