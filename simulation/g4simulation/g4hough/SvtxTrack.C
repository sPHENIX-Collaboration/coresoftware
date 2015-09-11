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
    x(0.0),
    y(0.0),
    z(0.0),
    covariance(6,6) {
  Reset();
}

SvtxTrack::SvtxTrack(SvtxTrack *track) : covariance( *(track->getCovariance()) )
{
  _phi   = track->get_phi();
  _d     = track->get_d();
  _kappa = track->get_kappa();
  _z0    = track->get_z0();
  _dzdl  = track->get_dzdl();
  setDCA2Dsigma(track->getDCA2Dsigma());
  
  for(int i=0;i<100;i++)
  {
    clusterID[i] = track->getClusterID(i);
    for(int j=0;j<3;j++){
      position[i][j] = track->getHitPosition(i,j);
    }
  }
  
  _track_id = track->getTrackID();
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
}

SvtxTrack::SvtxTrack(const SvtxTrack& track) : covariance( *(track.getCovariance()) )
{
  _phi   = track.get_phi();
  _d     = track.get_d();
  _kappa = track.get_kappa();
  _z0    = track.get_z0();
  _dzdl  = track.get_dzdl();
  setDCA2Dsigma(track.getDCA2Dsigma());
  
  for(int i=0;i<100;i++)
  {
    clusterID[i] = track.getClusterID(i);
    for(int j=0;j<3;j++){
      position[i][j] = track.getHitPosition(i,j);
    }
  }
  
  _track_id = track.getTrackID();
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
    clusterID[i]=-9999;
    for(int j=0;j<3;j++){
      position[i][j]=NAN;
    }
  }

  _track_id = -1;
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
