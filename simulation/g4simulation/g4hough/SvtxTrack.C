#include "SvtxTrack.h"

#include <math.h>

ClassImp(SvtxTrack)

using namespace std;

SvtxTrack::SvtxTrack()
  : _track_id(-1),
    _phi(0.0),
    _d(0.0),
    _kappa(0.0),
    _z0(0.0),
    _dzdl(0.0),
    _clusterID(),
    _x(0.0),
    _y(0.0),
    _z(0.0),
    _covariance(6,6) {
  Reset();
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
  return;
}

void SvtxTrack::Reset()
{
  for(int i=0;i<100;i++)
  {
    _clusterID[i]=-9999;
    for(int j=0;j<3;j++){
      _position[i][j]=NAN;
    }
  }

  _track_id = -1;
  _momentum=NAN;
  for(int j=0;j<3;j++){
    _mom3[j]=NAN;
  }

  _charge=1;
  _isprimary=false;
  _ispositive=false;
  _quality=NAN;
  
  _DCA=NAN;
  _DCA2D=NAN;

  _chisq = NAN;
  _chisqv = NAN;
  _ndf = 0;
  
  for(int i=0;i<9;i++){
    _scatter[i]=NAN;
  }

  for(int i=0;i<4;++i){
    _cal_dphi[i] = NAN;
    _cal_deta[i] = NAN;
    _cal_energy_3x3[i] = NAN;
    _cal_cluster_id[i] = -9999;
    _cal_cluster_e[i] = NAN;
  }

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
    if(_clusterID[i]>=0){layer=i;break;}
  }
  if(layer==-1){return NAN;}
  return _position[layer][coor];
}


short SvtxTrack::getNhits() const
{
  int count = 100;
  for(unsigned int i = 0; i < 100; i++){	
    if(_clusterID[i] == -9999)
      count--;
  }
  return count;
}
