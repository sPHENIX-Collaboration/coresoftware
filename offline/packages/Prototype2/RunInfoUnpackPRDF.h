#ifndef __CaloUnpackPRDFF__
#define __CaloUnpackPRDFF__

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <string>
#include <map>
#include <utility>

class Event;
class Packet;
class RawTowerContainer;
class RawTower;

class RunInfoUnpackPRDF : public SubsysReco
{
public:
  RunInfoUnpackPRDF();

  int
  Init(PHCompositeNode *topNode);

  int
  InitRun(PHCompositeNode *topNode);

  int
  process_event(PHCompositeNode *topNode);

  int
  End(PHCompositeNode *topNode);

  void
  CreateNodeTree(PHCompositeNode *topNode);

  //! add stuff to be unpacked
  void
  add_channel(const std::string & name, //! name of the channel
      const int packet_id, //! packet id
      const unsigned int offset, //! offset in packet data
      const double calibration_const = +1 //! conversion constant from integer to meaningful value
      );

private:

  class channel_info
  {
  public:
    channel_info(int p, unsigned int o, double c) :
        packet_id(p), offset(o), calibration_const(c)
    {
    }

    int packet_id;
    unsigned offset;
    double calibration_const;
  };

  //! list of channel name -> channel info
  typedef std::map<std::string, channel_info> typ_channel_map;

  typ_channel_map channel_map;

  std::string runinfo_node_name;
};

#endif //**CaloUnpackPRDFF**//
