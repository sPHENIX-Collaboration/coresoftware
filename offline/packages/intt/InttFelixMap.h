#ifndef INTT_FELIX_MAP_H
#define INTT_FELIX_MAP_H

namespace Intt { struct Online_s; }
namespace Intt { struct RawData_s; }

namespace InttFelix
{
	int RawDataToOnline(struct Intt::RawData_s const&, struct Intt::Online_s&);
	int OnlineToRawData(struct Intt::Online_s const&, struct Intt::RawData_s&);
};

#endif//INTT_FELIX_MAP_H
