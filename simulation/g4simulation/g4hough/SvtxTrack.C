#include "SvtxTrack.h"
#include <math.h>
#include <phool/phool.h>

ClassImp(SvtxTrack)

using namespace std;

SvtxTrack::SvtxTrack() : phi(0.),d(0.),kappa(0.),z0(0.),dzdl(0.), x(0.), y(0.), z(0.), covariance(6,6)
{
  Reset();
}

SvtxTrack::SvtxTrack(SvtxTrack *track) : covariance( *(track->getCovariance()) )
{
  phi   = track->phi   ;
  d     = track->d     ;
  kappa = track->kappa ;
  z0    = track->z0    ;
  dzdl  = track->dzdl  ;
  setDCA2Dsigma(track->getDCA2Dsigma());
  
  for(int i=0;i<100;i++)
  {
    clusterID[i] = track->getClusterID(i);
    for(int j=0;j<3;j++){
      position[i][j] = track->getHitPosition(i,j);
    }
  }
  
  trackID = track->getTrackID();
  momentum = track->getMomentum();
  for(int j=0;j<3;j++){
    mom3[j] = track->get3Momentum(j);
  }

  charge = track->getCharge();
  isprimary = track->getPrimary();
  ispositive = track->getPositive();
  quality = track->getQuality();
  
  DCA = track->getDCA();
  DCA2D = track->getDCA2D();

  chisq = track->getChisq();
  chisqv = track->getChisqv();
  ndf = track->getNDF();
  
  for(int i=0;i<9;i++){
    scatter[i] = track->getScatter(i);
  }

  for(int i=0;i<4;++i) {
    cal_dphi[i] = track->get_cal_dphi(i);
    cal_deta[i] = track->get_cal_deta(i);
    cal_energy_3x3[i] = track->get_cal_energy_3x3(i);
    cal_cluster_id[i] = track->get_cal_cluster_id(i);
    cal_cluster_e[i] = track->get_cal_cluster_e(i);
  }

  for (prop_map_t::const_iterator i = track->prop_map.begin();
      i != track->prop_map.end(); ++i)
    {
      PROPERTY prop_id = static_cast<PROPERTY>(i->first);
      set_property_nocheck(prop_id, track->get_property_nocheck(prop_id));
    }
}

SvtxTrack::SvtxTrack(const SvtxTrack& track) : covariance( *(track.getCovariance()) )
{
  phi   = track.phi   ;
  d     = track.d     ;
  kappa = track.kappa ;
  z0    = track.z0    ;
  dzdl  = track.dzdl  ;
  setDCA2Dsigma(track.getDCA2Dsigma());
  
  for(int i=0;i<100;i++)
  {
    clusterID[i] = track.getClusterID(i);
    for(int j=0;j<3;j++){
      position[i][j] = track.getHitPosition(i,j);
    }
  }
  
  trackID = track.getTrackID();
  momentum = track.getMomentum();
  for(int j=0;j<3;j++){
    mom3[j] = track.get3Momentum(j);
  }

  charge = track.getCharge();
  isprimary = track.getPrimary();
  ispositive = track.getPositive();
  quality = track.getQuality();
  
  DCA = track.getDCA();
  DCA2D = track.getDCA2D();

  chisq = track.getChisq();
  chisqv = track.getChisqv();
  ndf = track.getNDF();
  
  for(int i=0;i<9;i++){
    scatter[i] = track.getScatter(i);
  }

  for(int i=0;i<4;++i) {
    cal_dphi[i] = track.get_cal_dphi(i);
    cal_deta[i] = track.get_cal_deta(i);
    cal_energy_3x3[i] = track.get_cal_energy_3x3(i);
    cal_cluster_id[i] = track.get_cal_cluster_id(i);
    cal_cluster_e[i] = track.get_cal_cluster_e(i);
  }

  x = track.get_x();
  y = track.get_y();
  z = track.get_z();


  for (prop_map_t::const_iterator i = track.prop_map.begin();
      i != track.prop_map.end(); ++i)
    {
      PROPERTY prop_id = static_cast<PROPERTY>(i->first);
      set_property_nocheck(prop_id, track.get_property_nocheck(prop_id));
    }

}

void SvtxTrack::identify(ostream& os) const
{
  os << "SvtxTrack Object ";
  os << "id: " << getTrackID() << " ";

  if(getNhits() > 0)
    {
      os << "clusters: ";
      for(unsigned int i = 0; i < 100; i++)
	{
	  if(hasCluster(i))
	    {
	      os << getClusterID(i) << " ";
	    }
	}
    }
    
  if(getNDF() != 0)
    {
      os << "chisq/dof: " << getChisq()/getNDF() << " ";
    }
  os << std::endl;


  for (prop_map_t::const_iterator i = prop_map.begin(); i!= prop_map.end(); ++i)
    {
      PROPERTY prop_id = static_cast<PROPERTY>(i->first);
      pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << "\t" << i->first << ":\t" << property_info.first << " = \t";
      switch(property_info.second)
  {
  case type_int:
    cout << get_property_int(prop_id);
    break;
  case type_uint:
    cout << get_property_uint(prop_id);
    break;
  case type_float:
    cout << get_property_float(prop_id);
    break;
  default:
    cout << " unknown type ";
  }
      cout <<endl;
    }

  return;
}

void SvtxTrack::Reset()
{
  for(int i=0;i<100;i++)
  {
    clusterID[i]=-9999;
    for(int j=0;j<3;j++){
      position[i][j]=NAN;
    }
  }

  trackID = -1;
  momentum=NAN;
  for(int j=0;j<3;j++){
    mom3[j]=NAN;
  }

  charge=1;
  isprimary=false;
  ispositive=false;
  quality=NAN;
  
  DCA=NAN;
  DCA2D=NAN;

  chisq = NAN;
  chisqv = NAN;
  ndf = 0;
  
  for(int i=0;i<9;i++){
    scatter[i]=NAN;
  }

  for(int i=0;i<4;++i){
    cal_dphi[i] = NAN;
    cal_deta[i] = NAN;
    cal_energy_3x3[i] = NAN;
    cal_cluster_id[i] = -9999;
    cal_cluster_e[i] = NAN;
  }

  prop_map.clear();

  return;
}

int SvtxTrack::isValid() const
{
  return 1;
}


float SvtxTrack::getInnerMostHitPosition(int coor) const
{
  int layer=-1;
  for(int i=0;i<100;i++)
  {
    if(clusterID[i]>=0){layer=i;break;}
  }
  if(layer==-1){return NAN;}
  return position[layer][coor];
}


short SvtxTrack::getNhits() const
{
  int count = 100;
  for(unsigned int i = 0; i < 100; i++){	
    if(clusterID[i] == -9999)
      count--;
  }
  return count;
}





std::pair<const std::string,SvtxTrack::PROPERTY_TYPE>
SvtxTrack::get_property_info(const PROPERTY prop_id)
{
  switch (prop_id)
  {
  case  prop_FastSim_TruthID:
    return make_pair("Truth ID (available for fast sim track only)",SvtxTrack::type_int);
  default:
    cout << "unknown index " << prop_id << endl;
    exit(1);
  }
}


bool
SvtxTrack::check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type)
{
  pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
  if (property_info.second != prop_type)
    {
      return false;
    }
  return true;
}

string
SvtxTrack::get_property_type(const PROPERTY_TYPE prop_type)
{
  switch(prop_type)
    {
    case type_int:
      return "int";
    case type_uint:
      return "unsigned int";
    case type_float:
      return "float";
    default:
      return "unkown";
    }
}

bool
SvtxTrack::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i!=prop_map.end();
}

float
SvtxTrack::get_property_float(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_float))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_float) << endl;
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return (i==prop_map.end())? NAN : u_property(i->second).fdata ;
}

int
SvtxTrack::get_property_int(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_int) << endl;
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return (i==prop_map.end())? NAN : u_property(i->second).idata ;
}

unsigned int
SvtxTrack::get_property_uint(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_uint) << endl;
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return (i==prop_map.end())? NAN : u_property(i->second).uidata ;
}

void
SvtxTrack::set_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id,type_float))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_float) << endl;
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

void
SvtxTrack::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_int) << endl;
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

void
SvtxTrack::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second)
     << " not " << get_property_type(type_uint) << endl;
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

unsigned int
SvtxTrack::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
    {
      return iter->second;
    }
  return prop_MAX_NUMBER;
}


