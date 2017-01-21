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


class GenericUnpackPRDF : public SubsysReco
{
 public:
  GenericUnpackPRDF(const std::string & detector);

  int Init(PHCompositeNode *topNode);

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);
  
  void CreateNodeTree(PHCompositeNode *topNode);

  //! add stuff to be unpacked
  void add_channel(
      const int packet_id, //! packet id
      const int channel, //! channel in packet
      const int tower_id //! output tower id
  );

 private:
  std::string _detector;

  //!packet_id, channel number to define a hbd_channel
  typedef std::pair<int, int> hbd_channel_typ;

  //! list of hbd_channel -> channel id which is also tower id
  typedef std::map<hbd_channel_typ, int> hbd_channel_map;

  hbd_channel_map _hbd_channel_map;

  Event* _event;

  //output -> Towers
  RawTowerContainer* _towers;
};


#endif //**CaloUnpackPRDFF**//
