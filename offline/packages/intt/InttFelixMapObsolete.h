#ifndef INTT_FELIXMAPOBSOLETE_H
#define INTT_FELIXMAPOBSOLETE_H

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

#endif//INTT_FELIXMAPOBSOLETE_H
