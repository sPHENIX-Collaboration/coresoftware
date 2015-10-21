#ifndef __MOMENTUM_EVALUATOR__
#define __MOMENTUM_EVALUATOR__

#include <string>
#include <phool/PHCompositeNode.h>
#include <fun4all/SubsysReco.h>

class TNtuple;

class MomentumEvaluator : public SubsysReco
{
	public:
		MomentumEvaluator( std::string fname, float pt_s=0.1, float pz_s=0.2 );
		~MomentumEvaluator();

		int Init(PHCompositeNode *topNode);
		int process_event(PHCompositeNode *topNode);
		int End(PHCompositeNode *topNode);


	private:
		TNtuple* ntp_true;
		TNtuple* ntp_reco;
		float pt_search_scale;
		float pz_search_scale;
		unsigned int event_counter;
		std::string file_name;
};


#endif

