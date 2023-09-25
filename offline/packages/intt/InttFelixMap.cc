#include "InttFelixMap.h"

#include "InttMapping.h"

int InttFelix::RawDataToOnline(struct Intt::RawData_s const& raw, struct Intt::Online_s& onl)
{
	switch(raw.felix_server) {
		case  0:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  3;onl.ldr =  1;onl.arm =  0;return 0;
			case  1:
			onl.lyr =  1;onl.ldr =  1;onl.arm =  0;return 0;
			case  2:
			onl.lyr =  2;onl.ldr =  1;onl.arm =  0;return 0;
			case  3:
			onl.lyr =  2;onl.ldr =  0;onl.arm =  0;return 0;
			case  4:
			onl.lyr =  3;onl.ldr =  0;onl.arm =  0;return 0;
			case  5:
			onl.lyr =  0;onl.ldr =  0;onl.arm =  0;return 0;
			case  6:
			onl.lyr =  1;onl.ldr =  0;onl.arm =  0;return 0;
			case  7:
			onl.lyr =  3;onl.ldr =  3;onl.arm =  0;return 0;
			case  8:
			onl.lyr =  0;onl.ldr =  2;onl.arm =  0;return 0;
			case  9:
			onl.lyr =  0;onl.ldr =  1;onl.arm =  0;return 0;
			case 10:
			onl.lyr =  2;onl.ldr =  2;onl.arm =  0;return 0;
			case 11:
			onl.lyr =  3;onl.ldr =  2;onl.arm =  0;return 0;
			case 12:
			onl.lyr =  1;onl.ldr =  2;onl.arm =  0;return 0;
			case 13:
			onl.lyr =  2;onl.ldr =  3;onl.arm =  0;return 0;
			default:
			return 1;
		}
		case  1:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  3;onl.ldr =  5;onl.arm =  0;return 0;
			case  1:
			onl.lyr =  1;onl.ldr =  4;onl.arm =  0;return 0;
			case  2:
			onl.lyr =  1;onl.ldr =  3;onl.arm =  0;return 0;
			case  3:
			onl.lyr =  2;onl.ldr =  4;onl.arm =  0;return 0;
			case  4:
			onl.lyr =  3;onl.ldr =  4;onl.arm =  0;return 0;
			case  5:
			onl.lyr =  0;onl.ldr =  3;onl.arm =  0;return 0;
			case  6:
			onl.lyr =  2;onl.ldr =  5;onl.arm =  0;return 0;
			case  7:
			onl.lyr =  3;onl.ldr =  7;onl.arm =  0;return 0;
			case  8:
			onl.lyr =  0;onl.ldr =  5;onl.arm =  0;return 0;
			case  9:
			onl.lyr =  0;onl.ldr =  4;onl.arm =  0;return 0;
			case 10:
			onl.lyr =  2;onl.ldr =  6;onl.arm =  0;return 0;
			case 11:
			onl.lyr =  3;onl.ldr =  6;onl.arm =  0;return 0;
			case 12:
			onl.lyr =  1;onl.ldr =  5;onl.arm =  0;return 0;
			case 13:
			onl.lyr =  2;onl.ldr =  7;onl.arm =  0;return 0;
			default:
			return 1;
		}
		case  2:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  1;onl.ldr =  6;onl.arm =  0;return 0;
			case  1:
			onl.lyr =  0;onl.ldr =  6;onl.arm =  0;return 0;
			case  2:
			onl.lyr =  1;onl.ldr =  7;onl.arm =  0;return 0;
			case  3:
			onl.lyr =  3;onl.ldr =  8;onl.arm =  0;return 0;
			case  4:
			onl.lyr =  2;onl.ldr =  8;onl.arm =  0;return 0;
			case  5:
			onl.lyr =  3;onl.ldr =  9;onl.arm =  0;return 0;
			case  6:
			onl.lyr =  2;onl.ldr =  9;onl.arm =  0;return 0;
			case  7:
			onl.lyr =  0;onl.ldr =  7;onl.arm =  0;return 0;
			case  8:
			onl.lyr =  1;onl.ldr =  8;onl.arm =  0;return 0;
			case  9:
			onl.lyr =  0;onl.ldr =  8;onl.arm =  0;return 0;
			case 10:
			onl.lyr =  3;onl.ldr = 10;onl.arm =  0;return 0;
			case 11:
			onl.lyr =  2;onl.ldr = 10;onl.arm =  0;return 0;
			case 12:
			onl.lyr =  3;onl.ldr = 11;onl.arm =  0;return 0;
			case 13:
			onl.lyr =  2;onl.ldr = 11;onl.arm =  0;return 0;
			default:
			return 1;
		}
		case  3:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  1;onl.ldr =  9;onl.arm =  0;return 0;
			case  1:
			onl.lyr =  0;onl.ldr =  9;onl.arm =  0;return 0;
			case  2:
			onl.lyr =  1;onl.ldr = 10;onl.arm =  0;return 0;
			case  3:
			onl.lyr =  3;onl.ldr = 12;onl.arm =  0;return 0;
			case  4:
			onl.lyr =  2;onl.ldr = 12;onl.arm =  0;return 0;
			case  5:
			onl.lyr =  3;onl.ldr = 13;onl.arm =  0;return 0;
			case  6:
			onl.lyr =  2;onl.ldr = 13;onl.arm =  0;return 0;
			case  7:
			onl.lyr =  0;onl.ldr = 10;onl.arm =  0;return 0;
			case  8:
			onl.lyr =  1;onl.ldr = 11;onl.arm =  0;return 0;
			case  9:
			onl.lyr =  0;onl.ldr = 11;onl.arm =  0;return 0;
			case 10:
			onl.lyr =  3;onl.ldr = 14;onl.arm =  0;return 0;
			case 11:
			onl.lyr =  2;onl.ldr = 14;onl.arm =  0;return 0;
			case 12:
			onl.lyr =  3;onl.ldr = 15;onl.arm =  0;return 0;
			case 13:
			onl.lyr =  2;onl.ldr = 15;onl.arm =  0;return 0;
			default:
			return 1;
		}
		case  4:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  3;onl.ldr =  1;onl.arm =  1;return 0;
			case  1:
			onl.lyr =  0;onl.ldr =  0;onl.arm =  1;return 0;
			case  2:
			onl.lyr =  1;onl.ldr =  0;onl.arm =  1;return 0;
			case  3:
			onl.lyr =  2;onl.ldr =  1;onl.arm =  1;return 0;
			case  4:
			onl.lyr =  1;onl.ldr =  1;onl.arm =  1;return 0;
			case  5:
			onl.lyr =  2;onl.ldr =  0;onl.arm =  1;return 0;
			case  6:
			onl.lyr =  3;onl.ldr =  0;onl.arm =  1;return 0;
			case  7:
			onl.lyr =  0;onl.ldr =  2;onl.arm =  1;return 0;
			case  8:
			onl.lyr =  1;onl.ldr =  2;onl.arm =  1;return 0;
			case  9:
			onl.lyr =  0;onl.ldr =  1;onl.arm =  1;return 0;
			case 10:
			onl.lyr =  2;onl.ldr =  3;onl.arm =  1;return 0;
			case 11:
			onl.lyr =  3;onl.ldr =  3;onl.arm =  1;return 0;
			case 12:
			onl.lyr =  2;onl.ldr =  2;onl.arm =  1;return 0;
			case 13:
			onl.lyr =  3;onl.ldr =  2;onl.arm =  1;return 0;
			default:
			return 1;
		}
		case  5:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  0;onl.ldr =  3;onl.arm =  1;return 0;
			case  1:
			onl.lyr =  1;onl.ldr =  4;onl.arm =  1;return 0;
			case  2:
			onl.lyr =  1;onl.ldr =  3;onl.arm =  1;return 0;
			case  3:
			onl.lyr =  2;onl.ldr =  4;onl.arm =  1;return 0;
			case  4:
			onl.lyr =  2;onl.ldr =  5;onl.arm =  1;return 0;
			case  5:
			onl.lyr =  3;onl.ldr =  4;onl.arm =  1;return 0;
			case  6:
			onl.lyr =  3;onl.ldr =  5;onl.arm =  1;return 0;
			case  7:
			onl.lyr =  3;onl.ldr =  7;onl.arm =  1;return 0;
			case  8:
			onl.lyr =  2;onl.ldr =  7;onl.arm =  1;return 0;
			case  9:
			onl.lyr =  2;onl.ldr =  6;onl.arm =  1;return 0;
			case 10:
			onl.lyr =  3;onl.ldr =  6;onl.arm =  1;return 0;
			case 11:
			onl.lyr =  0;onl.ldr =  5;onl.arm =  1;return 0;
			case 12:
			onl.lyr =  1;onl.ldr =  5;onl.arm =  1;return 0;
			case 13:
			onl.lyr =  0;onl.ldr =  4;onl.arm =  1;return 0;
			default:
			return 1;
		}
		case  6:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  1;onl.ldr =  6;onl.arm =  1;return 0;
			case  1:
			onl.lyr =  0;onl.ldr =  6;onl.arm =  1;return 0;
			case  2:
			onl.lyr =  1;onl.ldr =  7;onl.arm =  1;return 0;
			case  3:
			onl.lyr =  3;onl.ldr =  8;onl.arm =  1;return 0;
			case  4:
			onl.lyr =  2;onl.ldr =  8;onl.arm =  1;return 0;
			case  5:
			onl.lyr =  3;onl.ldr =  9;onl.arm =  1;return 0;
			case  6:
			onl.lyr =  2;onl.ldr =  9;onl.arm =  1;return 0;
			case  7:
			onl.lyr =  0;onl.ldr =  7;onl.arm =  1;return 0;
			case  8:
			onl.lyr =  1;onl.ldr =  8;onl.arm =  1;return 0;
			case  9:
			onl.lyr =  0;onl.ldr =  8;onl.arm =  1;return 0;
			case 10:
			onl.lyr =  3;onl.ldr = 10;onl.arm =  1;return 0;
			case 11:
			onl.lyr =  2;onl.ldr = 10;onl.arm =  1;return 0;
			case 12:
			onl.lyr =  3;onl.ldr = 11;onl.arm =  1;return 0;
			case 13:
			onl.lyr =  2;onl.ldr = 11;onl.arm =  1;return 0;
			default:
			return 1;
		}
		case  7:
		switch(raw.felix_channel) {
			case  0:
			onl.lyr =  1;onl.ldr =  9;onl.arm =  1;return 0;
			case  1:
			onl.lyr =  0;onl.ldr =  9;onl.arm =  1;return 0;
			case  2:
			onl.lyr =  1;onl.ldr = 10;onl.arm =  1;return 0;
			case  3:
			onl.lyr =  3;onl.ldr = 12;onl.arm =  1;return 0;
			case  4:
			onl.lyr =  2;onl.ldr = 12;onl.arm =  1;return 0;
			case  5:
			onl.lyr =  3;onl.ldr = 13;onl.arm =  1;return 0;
			case  6:
			onl.lyr =  2;onl.ldr = 13;onl.arm =  1;return 0;
			case  7:
			onl.lyr =  0;onl.ldr = 10;onl.arm =  1;return 0;
			case  8:
			onl.lyr =  1;onl.ldr = 11;onl.arm =  1;return 0;
			case  9:
			onl.lyr =  0;onl.ldr = 11;onl.arm =  1;return 0;
			case 10:
			onl.lyr =  3;onl.ldr = 14;onl.arm =  1;return 0;
			case 11:
			onl.lyr =  2;onl.ldr = 14;onl.arm =  1;return 0;
			case 12:
			onl.lyr =  3;onl.ldr = 15;onl.arm =  1;return 0;
			case 13:
			onl.lyr =  2;onl.ldr = 15;onl.arm =  1;return 0;
			default:
			return 1;
		}
		default:
		return 1;
	}

	return 1;
}

int InttFelix::OnlineToRawData(struct Intt::Online_s const& onl, struct Intt::RawData_s& raw)
{
	switch(onl.lyr) {
		case  0:
		switch(onl.ldr) {
			case  0:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  5;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  1;return 0;
				default:
				return 1;
			}
			case  1:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  9;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  9;return 0;
				default:
				return 1;
			}
			case  2:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  8;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  7;return 0;
				default:
				return 1;
			}
			case  3:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  5;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  0;return 0;
				default:
				return 1;
			}
			case  4:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  9;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel = 13;return 0;
				default:
				return 1;
			}
			case  5:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  8;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel = 11;return 0;
				default:
				return 1;
			}
			case  6:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  1;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  1;return 0;
				default:
				return 1;
			}
			case  7:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  7;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  7;return 0;
				default:
				return 1;
			}
			case  8:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  9;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  9;return 0;
				default:
				return 1;
			}
			case  9:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  1;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  1;return 0;
				default:
				return 1;
			}
			case 10:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  7;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  7;return 0;
				default:
				return 1;
			}
			case 11:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  9;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  9;return 0;
				default:
				return 1;
			}
			default:
			return 1;
		}
		case  1:
		switch(onl.ldr) {
			case  0:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  6;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  2;return 0;
				default:
				return 1;
			}
			case  1:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  1;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  4;return 0;
				default:
				return 1;
			}
			case  2:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel = 12;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  8;return 0;
				default:
				return 1;
			}
			case  3:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  2;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  2;return 0;
				default:
				return 1;
			}
			case  4:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  1;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  1;return 0;
				default:
				return 1;
			}
			case  5:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel = 12;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel = 12;return 0;
				default:
				return 1;
			}
			case  6:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  0;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  0;return 0;
				default:
				return 1;
			}
			case  7:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  2;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  2;return 0;
				default:
				return 1;
			}
			case  8:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  8;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  8;return 0;
				default:
				return 1;
			}
			case  9:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  0;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  0;return 0;
				default:
				return 1;
			}
			case 10:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  2;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  2;return 0;
				default:
				return 1;
			}
			case 11:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  8;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  8;return 0;
				default:
				return 1;
			}
			default:
			return 1;
		}
		case  2:
		switch(onl.ldr) {
			case  0:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  3;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  5;return 0;
				default:
				return 1;
			}
			case  1:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  2;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  3;return 0;
				default:
				return 1;
			}
			case  2:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel = 10;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel = 12;return 0;
				default:
				return 1;
			}
			case  3:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel = 13;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel = 10;return 0;
				default:
				return 1;
			}
			case  4:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  3;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  3;return 0;
				default:
				return 1;
			}
			case  5:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  6;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  4;return 0;
				default:
				return 1;
			}
			case  6:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel = 10;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  9;return 0;
				default:
				return 1;
			}
			case  7:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel = 13;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  8;return 0;
				default:
				return 1;
			}
			case  8:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  4;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  4;return 0;
				default:
				return 1;
			}
			case  9:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  6;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  6;return 0;
				default:
				return 1;
			}
			case 10:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel = 11;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel = 11;return 0;
				default:
				return 1;
			}
			case 11:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel = 13;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel = 13;return 0;
				default:
				return 1;
			}
			case 12:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  4;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  4;return 0;
				default:
				return 1;
			}
			case 13:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  6;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  6;return 0;
				default:
				return 1;
			}
			case 14:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel = 11;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel = 11;return 0;
				default:
				return 1;
			}
			case 15:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel = 13;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel = 13;return 0;
				default:
				return 1;
			}
			default:
			return 1;
		}
		case  3:
		switch(onl.ldr) {
			case  0:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  4;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  6;return 0;
				default:
				return 1;
			}
			case  1:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  0;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel =  0;return 0;
				default:
				return 1;
			}
			case  2:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel = 11;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel = 13;return 0;
				default:
				return 1;
			}
			case  3:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  0;raw.felix_channel =  7;return 0;
				case  1:
				raw.felix_server =  4;raw.felix_channel = 11;return 0;
				default:
				return 1;
			}
			case  4:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  4;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  5;return 0;
				default:
				return 1;
			}
			case  5:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  0;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  6;return 0;
				default:
				return 1;
			}
			case  6:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel = 11;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel = 10;return 0;
				default:
				return 1;
			}
			case  7:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  1;raw.felix_channel =  7;return 0;
				case  1:
				raw.felix_server =  5;raw.felix_channel =  7;return 0;
				default:
				return 1;
			}
			case  8:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  3;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  3;return 0;
				default:
				return 1;
			}
			case  9:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel =  5;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel =  5;return 0;
				default:
				return 1;
			}
			case 10:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel = 10;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel = 10;return 0;
				default:
				return 1;
			}
			case 11:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  2;raw.felix_channel = 12;return 0;
				case  1:
				raw.felix_server =  6;raw.felix_channel = 12;return 0;
				default:
				return 1;
			}
			case 12:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  3;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  3;return 0;
				default:
				return 1;
			}
			case 13:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel =  5;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel =  5;return 0;
				default:
				return 1;
			}
			case 14:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel = 10;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel = 10;return 0;
				default:
				return 1;
			}
			case 15:
			switch(onl.arm) {
				case  0:
				raw.felix_server =  3;raw.felix_channel = 12;return 0;
				case  1:
				raw.felix_server =  7;raw.felix_channel = 12;return 0;
				default:
				return 1;
			}
			default:
			return 1;
		}
		default:
		return 1;
	}

	return 1;
}
