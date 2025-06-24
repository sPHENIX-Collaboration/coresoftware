// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKFITTINGQA_H
#define TRACKFITTINGQA_H

#include <fun4all/SubsysReco.h>

#include <string>

class TH1;
class PHCompositeNode;

class TrackFittingQA : public SubsysReco {
public:
	TrackFittingQA(const std::string &name = "TrackFittingQA");
	~TrackFittingQA() override;

	/// sets the name of node to retrieve the track map from (default member value is "SvtxTrackMap")
	void set_track_map_name (std::string const& track_map_node_name) {m_track_map_node_name = track_map_node_name;}

	/// sets the name of node to retrieve the state map from (default member value is "SvtxAlignmentStateMap")
	void set_state_map_name (std::string const& state_map_node_name) {m_state_map_node_name = state_map_node_name;}

	int Init(PHCompositeNode*) override;
	int InitRun(PHCompositeNode*) override;
	
	int process_event(PHCompositeNode*) override;
	
	int ResetEvent(PHCompositeNode*) override;
	
	int EndRun(const int runnumber) override;
	int End(PHCompositeNode*) override;
	
	int Reset(PHCompositeNode*) override;
	void Print(std::string const& = "ALL") const override;

private:
	TH1* m_quality_hist[2]{};
	TH1* m_p_hist[2]{};
	TH1* m_pt_hist[2]{};
	TH1* m_eta_hist[2]{};
	TH1* m_phi_eta_hist[2]{};
	TH1* m_intt_states_hist[2]{};
	TH1* m_mvtx_states_hist[2]{};
	TH1* m_tpc_states_hist[2]{};
	TH1* m_tpot_states_hist[2]{};
	// ...
	// residual plots as a function of tpc sector (Mariia)

	std::string m_track_map_node_name = "SvtxTrackMap";
	std::string m_state_map_node_name = "SvtxAlignmentStateMap";

};

#endif // TRACKFITTINGQA_H
