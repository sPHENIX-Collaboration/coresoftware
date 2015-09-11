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
    _x(0.0),
    _y(0.0),
    _z(0.0),
    _covariance(6,6),
    _cluster_ids() {
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
	  if (hasCluster(i)) {
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

void SvtxTrack::Reset() {

  _track_id = -1;

  _cluster_ids.clear();
  for (int i=0;i<100;i++) {
    for(int j=0;j<3;j++){
      _position[i][j]=NAN;
    }
  }

  _momentum=NAN;
  for(int j=0;j<3;j++){
    _mom3[j]=NAN;
  }

  _charge=1;
  _ispositive=false;
  _quality=NAN;
  
  _DCA=NAN;
  _DCA2D=NAN;

  _chisq = NAN;
  _chisqv = NAN;
  _ndf = 0;
  
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


float SvtxTrack::getInnerMostHitPosition(int coor) const {
  if (_cluster_ids.empty()) return NAN;  
  int layer = _cluster_ids.begin()->first; 
  return _position[layer][coor];
}

short SvtxTrack::getNhits() const {
  return _cluster_ids.size();
}
