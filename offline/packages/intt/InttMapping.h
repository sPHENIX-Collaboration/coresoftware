#ifndef INTT_MAPPING_H
#define INTT_MAPPING_H

#include <map>

class Packet;
class InttRawHit;

namespace Intt
{
	extern const std::map<int, int> Packet_Id;
	int FelixFromPacket(int);

	struct RawData_s
	{
		int felix_server = 0;
		int felix_channel = 0;
		int chip = 0;
		int channel = 0;
	};

	struct Online_s
	{
		int lyr = 0;
		int ldr = 0;
		int arm = 0;
		int chp = 0;
		int chn = 0;
	};

	struct Offline_s
	{
		int layer = 0;
		int ladder_phi = 0;
		int ladder_z = 0;
		int strip_x = 0;
		int strip_y = 0;
	};

	struct Intt::RawData_s RawFromPacket(int const, int const, Packet*);
	void RawFromHit(struct Intt::RawData_s&, InttRawHit*);

	//nontrivial
	struct Online_s ToOnline(struct Offline_s const&);
	struct Offline_s ToOffline(struct Online_s const&);

	struct RawData_s ToRawData(struct Online_s const&);
	struct Online_s ToOnline(struct RawData_s const&);

	//trivial
	struct RawData_s ToRawData(struct Offline_s const&);
	struct Offline_s ToOffline(struct RawData_s const&);

	//Eigen::Affine3d GetTransform(TTree*, struct Offline_s const&);
	//Eigen::Vector4d GetLocalPos(struct Offline_s const&);
};

bool operator==(struct Intt::RawData_s const&, struct Intt::RawData_s const&);
bool operator==(struct Intt::Online_s const&, struct Intt::Online_s const&);
bool operator==(struct Intt::Offline_s const&, struct Intt::Offline_s const&);

bool operator!=(struct Intt::RawData_s const&, struct Intt::RawData_s const&);
bool operator!=(struct Intt::Online_s const&, struct Intt::Online_s const&);
bool operator!=(struct Intt::Offline_s const&, struct Intt::Offline_s const&);

bool operator<(struct Intt::RawData_s const&, struct Intt::RawData_s const&);
bool operator<(struct Intt::Online_s const&, struct Intt::Online_s const&);
bool operator<(struct Intt::Offline_s const&, struct Intt::Offline_s const&);

bool operator>(struct Intt::RawData_s const&, struct Intt::RawData_s const&);
bool operator>(struct Intt::Online_s const&, struct Intt::Online_s const&);
bool operator>(struct Intt::Offline_s const&, struct Intt::Offline_s const&);

bool operator<=(struct Intt::RawData_s const&, struct Intt::RawData_s const&);
bool operator<=(struct Intt::Online_s const&, struct Intt::Online_s const&);
bool operator<=(struct Intt::Offline_s const&, struct Intt::Offline_s const&);

bool operator>=(struct Intt::RawData_s const&, struct Intt::RawData_s const&);
bool operator>=(struct Intt::Online_s const&, struct Intt::Online_s const&);
bool operator>=(struct Intt::Offline_s const&, struct Intt::Offline_s const&);

#endif//INTT_MAPPING_H

//      Ladder Structure                //
//======================================//
//      U14     U1      Type B  North   //
//      U15     U2        .       .     //
//      U16     U3        .       .     //
//      U17     U4        .       .     //
//      U18     U5      Type B    .     //
//------------------------------  .     //
//      U19     U6      Type A    .     //
//      U20     U7        .       .     //
//      U21     U8        .       .     //
//      U22     U9        .       .     //
//      U23     U10       .       .     //
//      U24     U11       .       .     //
//      U25     U12       .       .     //
//      U26     U13     Type A  North   //
//--------------------------------------//
//      U13     U26     Type A  South   //
//      U12     U25       .       .     //
//      U11     U24       .       .     //
//      U10     U23       .       .     //
//      U9      U22       .       .     //
//      U8      U21       .       .     //
//      U7      U20       .       .     //
//      U6      U19     Type A    .     //
//------------------------------  .     //
//      U5      U18     Type B    .     //
//      U4      U17       .       .     //
//      U3      U16       .       .     //
//      U2      U15       .       .     //
//      U1      U14     Type B  South   //
//======================================//

//|<-            South                        |                       North            ->|
//|         B      |            A             |             A           |       B        |
//+----------------+--------------------------+-------------------------+----------------+
//|         1      |            0             |             2           |       3        | ladder_z  (offline)
//|  0  1  2  3  4 |  0  1  2  3  4  5  6  7  |  0  1  2  3  4  5  6  7 |  0  1  2  3  4 | strip_x   (offline)
//+----------------+--------------------------+-------------------------+----------------+
//|  0  1  2  3  4 |  5  6  7  8  9 10 11 12  | 25 24 23 22 21 20 19 18 | 17 16 15 14 13 | chp       (online)
//| 13 14 15 16 17 | 18 19 20 21 22 23 24 25  | 12 11 10  9  8  7  6  5 |  4  3  2  1  0 | chp       (online)
