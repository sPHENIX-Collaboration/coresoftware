#include "InttFelixMapObsolete.h"

#include <cstdlib>

int INTT_Felix::FelixMap(int const& felix, int const& felix_channel, struct Ladder_s& ladder_struct)
{
	switch(felix)
	{
		case 0:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		case 1:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		case 2:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		case 3:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 12;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 12;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 13;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 13;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 14;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 14;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 15;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 15;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		case 4:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 0;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 1;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 2;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		case 5:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 3;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 5;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 4;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		case 6:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 6;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 7;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 8;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		case 7:
			switch(felix_channel)
			{
				case 0:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 1:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 9;
					return EXIT_SUCCESS;
				case 2:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 3:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 12;
					return EXIT_SUCCESS;
				case 4:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 12;
					return EXIT_SUCCESS;
				case 5:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 13;
					return EXIT_SUCCESS;
				case 6:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 13;
					return EXIT_SUCCESS;
				case 7:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 10;
					return EXIT_SUCCESS;
				case 8:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				case 9:
					ladder_struct.barrel = 0;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 11;
					return EXIT_SUCCESS;
				case 10:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 14;
					return EXIT_SUCCESS;
				case 11:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 14;
					return EXIT_SUCCESS;
				case 12:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 1;
					ladder_struct.ladder = 15;
					return EXIT_SUCCESS;
				case 13:
					ladder_struct.barrel = 1;
					ladder_struct.layer = 0;
					ladder_struct.ladder = 15;
					return EXIT_SUCCESS;
				default:
					ladder_struct.barrel = -1;
					ladder_struct.layer = -1;
					ladder_struct.ladder = -1;
					return EXIT_FAILURE;
			}
		default:
			ladder_struct.barrel = -1;
			ladder_struct.layer = -1;
			ladder_struct.ladder = -1;
			return EXIT_FAILURE;
	}
	return EXIT_FAILURE;
}
