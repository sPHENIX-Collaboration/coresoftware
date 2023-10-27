#ifndef INTT_FELIX_MAP_H
#define INTT_FELIX_MAP_H

namespace InttDefs { struct Online_s; }
namespace InttDefs { struct RawData_s; }

namespace InttFelix
{
	int RawDataToOnline(struct InttDefs::RawData_s const&, struct InttDefs::Online_s&);
	int OnlineToRawData(struct InttDefs::Online_s const&, struct InttDefs::RawData_s&);
};

#endif//INTT_FELIX_MAP_H
