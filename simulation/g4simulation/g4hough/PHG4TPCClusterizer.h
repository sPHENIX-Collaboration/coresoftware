#ifndef __PHG4TPCCLUSTERIZER__
#define __PHG4TPCCLUSTERIZER__

#include <fun4all/SubsysReco.h>
#include <vector>

class PHG4TPCClusterizer : public SubsysReco
{
	public:
		PHG4TPCClusterizer(const char * name = "PHG4SvtxClusterizer", unsigned int phi_s=10,unsigned int z_s=5) : SubsysReco(name), _phi_span(phi_s), _z_span(z_s) {}
		~PHG4TPCClusterizer(){}

		//! module initialization
		int Init(PHCompositeNode *topNode){return 0;}

		//! run initialization
		int InitRun(PHCompositeNode *topNode);

		//! event processing
		int process_event(PHCompositeNode *topNode);

		//! end of process
		int End(PHCompositeNode *topNode){return 0;}


	private:
		// std::vector<std::vector<std::vector<float> > > amps;
		// std::vector<std::vector<std::vector<int> > > cellids;
		// std::vector<std::vector<int> > nhits;

		std::vector<int> nhits;
		std::vector<float> amps;
		std::vector<int> cellids;
		int nphibins;
		int nzbins;

		void reset();

		unsigned int _phi_span;
		unsigned int _z_span;
};


#endif
