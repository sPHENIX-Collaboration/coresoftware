/*
 * TPCFEETestRecov1.h
 *
 *  Created on: Sep 19, 2018
 *      Author: jinhuang
 */

#ifndef CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_
#define CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class Fun4AllHistoManager;
class TTree;

class TPCFEETestRecov1: public SubsysReco {
public:
	TPCFEETestRecov1(const std::string &outputfilename =
			"TPCFEETestRecov1.root");
	virtual ~TPCFEETestRecov1();

	int Init(PHCompositeNode *topNode);
	int InitRun(PHCompositeNode *topNode);
	int process_event(PHCompositeNode *topNode);
	int End(PHCompositeNode *topNode);

private:
#ifndef __CINT__

	Fun4AllHistoManager *getHistoManager();

	std::string m_outputFileName;

	TTree * m_T;

	int m_event;

#endif  // #ifndef __CINT__
};

#endif /* CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_ */
