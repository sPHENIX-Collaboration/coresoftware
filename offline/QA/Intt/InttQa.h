#pragma once

// We really just need timing and hitmap plots
// I don't know what the status of the Event Mixup QA is
//
// This will eventually contain everything

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class TH1;
class PHCompositeNode;
// class Gl1Packet; // For streaming timing QA; correct implementations are only in inttcalib
// class InttEventInfo; // EventMixUpQa (eventually)
class InttRawHitContainer;

class InttQa : public SubsysReco {
public:
	InttQa (std::string const& = "InttQa");
	// The HistoManager owns our histograms after we register them,
	// and we do not dynamically allocate anything else
	virtual ~InttQa() = default;

	int InitRun(PHCompositeNode*) override;
	int process_event(PHCompositeNode*) override;

private:
	static std::string constexpr m_prefix{"h_InttQa"};

	static int const n_felix_servers{8};
	static int const n_felix_channels{14};
	static int const n_chips{26};
	static int const n_channels{128};
	static int const n_bcos{128};
	static int const n_adcs{8};

	typedef std::array<TH1*, n_felix_servers> per_felix_server_array_t;
	typedef std::array<std::array<TH1*, n_felix_channels>, n_felix_servers> per_felix_channel_array_t;
	typedef std::array<std::array<std::array<TH1*, n_chips>, n_felix_channels>, n_felix_servers> per_chip_array_t;

	per_felix_server_array_t m_felix_server_hit_distribution{};
	per_felix_channel_array_t m_felix_channel_bco_distribution{};
	per_felix_channel_array_t m_felix_channel_hit_distribution{};
	per_chip_array_t m_chip_hit_distribution{};
	per_chip_array_t m_chip_adc_distribution{};

	static int const n_barrels{2};
	static int const n_ladder_z{4};
	static std::array<int, n_barrels> constexpr n_ladders{24, 32, };
	static std::array<int, n_ladder_z> constexpr n_strip_y{5, 8, 8, 5, };

	typedef std::array<TH1*, n_barrels> per_barrel_array_t;
	per_barrel_array_t m_barrel_hit_distribution{};

	// Gl1Packet* m_gl1_packet{};
	// std::vector<InttEventHeader*> m_intt_event_headers{};
	std::vector<InttRawHitContainer*> m_intt_raw_hit_containers{};

};
