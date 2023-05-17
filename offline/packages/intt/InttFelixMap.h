#ifndef FELIX_MAP_H
#define FELIX_MAP_H

#include <cstdlib>

namespace INTT_Felix
{
	struct Ladder_s
	{
		int barrel;
		int layer;
		int ladder;
	};

	int FelixMap(int const&, int const&, struct Ladder_s&);
};

#endif//FELIX_MAP_H
