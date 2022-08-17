#include "ClusterErrorPara.h"
#include <trackbase_historic/TrackSeed.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>

#include <TF1.h>
#include <string>
#include <iostream>
#include <cmath>

namespace
{

  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

ClusterErrorPara::ClusterErrorPara()
{
  f0 = new TF1("f0","[0] + x*[1]",0,10);
  f0->SetParameter(0,0.0158);
  f0->SetParameter(1,0.00772013);

  f1 = new TF1("f1","[0] + x*[1]",0,10);
  f1->SetParameter(0,0.0102342);
  f1->SetParameter(1,0.0320397);

  f2 = new TF1("f2","[0] + x*[1]",0,10);
  f2->SetParameter(0,0.010688);
  f2->SetParameter(1,0.0300981);

  fz0 = new TF1("fz0","pol1",-2,2);
  fz0->SetParameter(0,0.0520295);
  fz0->SetParameter(1,0.0129029);

  fz = new TF1("fz","pol2",-2,2);
  fz->SetParameter(0,0.0409945);
  fz->SetParameter(1,-0.00530857);
  fz->SetParameter(2,0.0573064);

  fmm_55_2 = new TF1("fmm_55_2","pol2",-2,2);
  fmm_55_2->SetParameter(0,0.0430592);
  fmm_55_2->SetParameter(1,-0.000177174);
  fmm_55_2->SetParameter(2,0.0914288);

  fmm_56_2 = new TF1("fmm_56_2","pol2",-2,2);
  fmm_56_2->SetParameter(0,0.00363897);
  fmm_56_2->SetParameter(1,0.0109713);
  fmm_56_2->SetParameter(2,0.032354);

  fmm_3 = new TF1("fmm_3","pol2",-2,2);
  fmm_3->SetParameter(0,0.00305396);
  fmm_3->SetParameter(1,0.00505814);
  fmm_3->SetParameter(2,0.0395137);

  static const double invsqrt12 = 1./std::sqrt(12);

  pitcherr_phi_mvtx = 0.002688 * invsqrt12;
  pitcherr_phi_intt = 0.0078 * invsqrt12;
  pitcherr_phi_mm1 = 0.1* invsqrt12 ;
  pitcherr_phi_mm2 = 31.6* invsqrt12 ;
  pitcherr_z_mvtx = 0.002924* invsqrt12 ;
  pitcherr_z_intt = 1.6* invsqrt12 ;
  pitcherr_z_mm1 = 54.2* invsqrt12 ;
  pitcherr_z_mm2 = 0.2* invsqrt12 ;

}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_cluster_error(TrackSeed *seed, TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key)
{
  float r = cluster_r;
  float R = TMath::Abs(1.0/seed->get_qOverR());
  double alpha = (r*r) /(2*r*R);
  double beta = atan(seed->get_slope());
  return get_cluster_error(cluster, key, alpha, beta);
}
//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_cluster_error(TrkrCluster* cluster,  TrkrDefs::cluskey key, double alpha, double beta)
{

  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror   = 0;

  switch( TrkrDefs::getTrkrId( key ) )
    {
   
    default: break;

    case TrkrDefs::micromegasId:
      if(layer==55){
	zerror = pitcherr_z_mm1;
	if(cluster->getPhiSize()==1){
	  phierror = pitcherr_phi_mm1;
	}else if(cluster->getPhiSize()==2){
	  phierror = fmm_55_2->Eval(alpha);
	}else if(cluster->getPhiSize()>=3){
	  phierror = fmm_3->Eval(alpha);
	}
	phierror *=scale_mm_0;
      }else if(layer==56){
	phierror =pitcherr_phi_mm2;
	if(cluster->getZSize()==1){
	  zerror = pitcherr_z_mm2;
	}else if(cluster->getZSize()==2){
	  zerror = fmm_56_2->Eval(beta);
	}else if(cluster->getZSize()>=3){
	  zerror = fmm_3->Eval(beta);
	}
	zerror *=scale_mm_1;
      }
      break;
    case TrkrDefs::mvtxId:
      static constexpr std::array<double, 7> scalefactors_mvtx_phi = {{ 0.36, 0.6,0.37,0.49,0.4,0.37,0.33 }};
      phierror = pitcherr_phi_mvtx;
      if(cluster->getPhiSize() == 1 && cluster->getZSize() == 1) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[0];
      else if(cluster->getPhiSize() == 2 && cluster->getZSize() == 1) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[1];
      else if(cluster->getPhiSize() == 1 && cluster->getZSize() == 2) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[2];
      else if( cluster->getPhiSize() == 2 && cluster->getZSize() == 2 ) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[0];
      else if( cluster->getPhiSize() == 2 && cluster->getZSize() == 3 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[1];
      else if( cluster->getPhiSize() == 3 && cluster->getZSize() == 2 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[2];
      else if( cluster->getPhiSize() == 3 && cluster->getZSize() == 3 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[3];
      else phierror = pitcherr_phi_mvtx*cluster->getPhiSize();

      phierror *= scale_mvtx;
      zerror = pitcherr_z_mvtx;
      static constexpr std::array<double, 4> scalefactors_z = {{ 0.47, 0.48, 0.71, 0.55 }};
      if( cluster->getZSize() == 2 && cluster->getPhiSize() == 2 )       zerror = pitcherr_z_mvtx*scalefactors_z[0];
      else if( cluster->getZSize() == 2 && cluster->getPhiSize() == 3 )  zerror = pitcherr_z_mvtx*scalefactors_z[1];
      else if( cluster->getZSize() == 3 && cluster->getPhiSize() == 2 )  zerror = pitcherr_z_mvtx*scalefactors_z[2];
      else if( cluster->getZSize() == 3 && cluster->getPhiSize() == 3 )  zerror = pitcherr_z_mvtx*scalefactors_z[3];
      zerror *= scale_mvtx_z;
      // else zerror = pitcherr_z_mvtx*cluster->getZSize();
      break;
      
    case TrkrDefs::inttId:
      static constexpr std::array<double, 3> scalefactors_intt_phi = {{ 0.85, 0.4, 0.33 }};
      
      if( cluster->getPhiSize() == 1 && layer < 5) phierror = pitcherr_phi_intt*scalefactors_intt_phi[0];
      else if( cluster->getPhiSize() == 2 && layer < 5) phierror = pitcherr_phi_intt*scalefactors_intt_phi[1];
      else if( cluster->getPhiSize() == 2 && layer > 4) phierror = pitcherr_phi_intt*scalefactors_intt_phi[2];
      else phierror = pitcherr_phi_intt*cluster->getPhiSize();
      if( layer == 3) phierror *= scale_intt_3;
      if( layer == 4) phierror *= scale_intt_4;
      if( layer == 5) phierror *= scale_intt_5;
      if( layer == 6) phierror *= scale_intt_6;

      zerror = pitcherr_z_intt*cluster->getZSize();
      break;
      
    case TrkrDefs::tpcId:
      
      int sector = -1;
      if(layer >=7 && layer < 23){
	sector = 0;
      }else if(layer>=23 && layer <39){
	sector = 1;
      }else if(layer>=39 && layer <55){
	sector = 2;
      }
      
      phierror = 0.0005;
      if(sector!=-1)
	phierror = 0.0005;
      if(sector==0){
	phierror = f0->Eval(alpha);
	phierror *=  scale_tpc_0;
	zerror = fz0->Eval(beta);
	zerror *=  scale_tpc_0_z;
      }
      if(sector==1){
	phierror = f1->Eval(alpha);
	phierror *=  scale_tpc_1;
	zerror = fz->Eval(beta);
	zerror *=  scale_tpc_1_z;
      }
      if(sector==2){
	phierror = f2->Eval(alpha);
	phierror *=  scale_tpc_2;
	zerror = fz->Eval(beta);
	zerror *=  scale_tpc_2_z;
      }


      TrkrClusterv4 *clusterv4 = dynamic_cast<TrkrClusterv4 *>(cluster);
      if(clusterv4->getEdge()>=1)
	phierror *= 2;
      break;
     
    }
  if(phierror==0) phierror = 100;
  if(zerror==0) zerror = 100;
  
  return std::make_pair(square(phierror),square(zerror));  
  //  return std::make_pair(phierror,zerror);  

}
//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_simple_cluster_error(TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key)
{

  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror   = 0;
  if(cluster_r>100){
    phierror *=1;
  }
  switch( TrkrDefs::getTrkrId( key ) )
    {
   
    default: break;

    case TrkrDefs::micromegasId:
      if(layer==55){
	phierror = pitcherr_phi_mm1;
	zerror = pitcherr_z_mm1;
      }else if(layer==56){
	phierror =pitcherr_phi_mm2;// = 31.6;
	zerror = pitcherr_z_mm2;// = 0.2;
      }
      break;
    case TrkrDefs::mvtxId:
      static constexpr std::array<double, 7> scalefactors_mvtx_phi = {{ 0.36, 0.6,0.37,0.49,0.4,0.37,0.33 }};
      
      if(cluster->getPhiSize() == 1 && cluster->getZSize() == 1) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[0];
      else if(cluster->getPhiSize() == 2 && cluster->getZSize() == 1) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[1];
      else if(cluster->getPhiSize() == 1 && cluster->getZSize() == 2) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[2];
      else if( cluster->getPhiSize() == 2 && cluster->getZSize() == 2 ) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[0];
      else if( cluster->getPhiSize() == 2 && cluster->getZSize() == 3 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[1];
      else if( cluster->getPhiSize() == 3 && cluster->getZSize() == 2 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[2];
      else if( cluster->getPhiSize() == 3 && cluster->getZSize() == 3 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[3];
      else phierror = pitcherr_phi_mvtx*cluster->getPhiSize();
      static constexpr std::array<double, 4> scalefactors_z = {{ 0.47, 0.48, 0.71, 0.55 }};
      if( cluster->getZSize() == 2 && cluster->getPhiSize() == 2 ) zerror = pitcherr_z_mvtx*scalefactors_z[0];
      else if( cluster->getZSize() == 2 && cluster->getPhiSize() == 3 )  zerror = pitcherr_z_mvtx*scalefactors_z[1];
      else if( cluster->getZSize() == 3 && cluster->getPhiSize() == 2 )  zerror = pitcherr_z_mvtx*scalefactors_z[2];
      else if( cluster->getZSize() == 3 && cluster->getPhiSize() == 3 )  zerror = pitcherr_z_mvtx*scalefactors_z[3];
      else zerror = pitcherr_z_mvtx*cluster->getZSize();
      break;
      
    case TrkrDefs::inttId:
      static constexpr std::array<double, 3> scalefactors_intt_phi = {{ 0.85, 0.4, 0.33 }};
      
      if( cluster->getPhiSize() == 1 && layer < 5) phierror = pitcherr_phi_intt*scalefactors_intt_phi[0];
      else if( cluster->getPhiSize() == 2 && layer < 5) phierror = pitcherr_phi_intt*scalefactors_intt_phi[1];
      else if( cluster->getPhiSize() == 2 && layer > 4) phierror = pitcherr_phi_intt*scalefactors_intt_phi[2];
      else phierror = pitcherr_phi_intt*cluster->getPhiSize();

      zerror = pitcherr_z_intt*cluster->getZSize();
      break;
      
    case TrkrDefs::tpcId:
      
      int sector = -1;
      if(layer >=7 && layer < 23){
	sector = 0;
      }else if(layer>=23 && layer <39){
	sector = 1;
      }else if(layer>=39 && layer <55){
	sector = 2;
      }
      
      float alpha = 0.1;
      double beta = 0.2;
      phierror = 0.0005;
      if(sector!=-1)
	phierror = 0.0005;
      if(sector==0){
	phierror = f0->Eval(alpha);
      }
      if(sector==1){
	phierror = f1->Eval(alpha);
      }
      if(sector==2){
	phierror = f2->Eval(alpha);
      }
      zerror = fz->Eval(beta);
      //      phierror*=cluster_r;
      TrkrClusterv4 *clusterv4 = dynamic_cast<TrkrClusterv4 *>(cluster);
      if(clusterv4->getEdge()>=1)
	phierror *= 2;
      
      break;
     
    }
  if(phierror==0) phierror = 100;
  if(zerror==0) zerror = 100;
  
  return std::make_pair(square(phierror),square(zerror));  
  //  return std::make_pair(phierror,zerror);  

}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_fix_tpc_cluster_error(TrkrCluster* cluster, TrkrDefs::cluskey key)
{

  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror   = 0;
  
  int sector = -1;
  if(layer >=7 && layer < 23){
    sector = 0;
  }else if(layer>=23 && layer <39){
    sector = 1;
  }else if(layer>=39 && layer <55){
    sector = 2;
  }
  
  float alpha = 0.1;
  
  phierror = 0.0005;
  if(sector!=-1)
    phierror = 0.0005;
  if(sector==0){
    phierror = f0->Eval(alpha);
  }
  if(sector==1){
    phierror = f1->Eval(alpha);
  }
  if(sector==2){
    phierror = f2->Eval(alpha);
  }
  double beta = 0.1;
  zerror = fz->Eval(beta);

  TrkrClusterv4 *clusterv4 = dynamic_cast<TrkrClusterv4 *>(cluster);
  if(clusterv4->getEdge()>=1)
    phierror *= 2;
  if(phierror==0) phierror = 100;
  if(zerror==0) zerror = 100;

  return std::make_pair(square(phierror),square(zerror));  
  //  return std::make_pair(phierror,zerror);  

}

//_________________________________________________________________________________
ClusterErrorPara::error_t ClusterErrorPara::get_si_cluster_error(const TrkrCluster* cluster, TrkrDefs::cluskey key)
{

  int layer = TrkrDefs::getLayer(key);
  double phierror = 0;
  double zerror   = 0;

  switch( TrkrDefs::getTrkrId( key ) )
    {
   
    default: break;

    case TrkrDefs::mvtxId:
      static constexpr std::array<double, 7> scalefactors_mvtx_phi = {{ 0.36, 0.6,0.37,0.49,0.4,0.37,0.33 }};
      
      if(cluster->getPhiSize() == 1 && cluster->getZSize() == 1) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[0];
      else if(cluster->getPhiSize() == 2 && cluster->getZSize() == 1) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[1];
      else if(cluster->getPhiSize() == 1 && cluster->getZSize() == 2) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[2];
      else if( cluster->getPhiSize() == 2 && cluster->getZSize() == 2 ) phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[0];
      else if( cluster->getPhiSize() == 2 && cluster->getZSize() == 3 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[1];
      else if( cluster->getPhiSize() == 3 && cluster->getZSize() == 2 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[2];
      else if( cluster->getPhiSize() == 3 && cluster->getZSize() == 3 )  phierror = pitcherr_phi_mvtx*scalefactors_mvtx_phi[3];
      else phierror = pitcherr_phi_mvtx*cluster->getPhiSize();
      phierror *= scale_mvtx;
      static constexpr std::array<double, 4> scalefactors_z = {{ 0.47, 0.48, 0.71, 0.55 }};
      if( cluster->getZSize() == 2 && cluster->getPhiSize() == 2 )       zerror = pitcherr_z_mvtx*scalefactors_z[0];
      else if( cluster->getZSize() == 2 && cluster->getPhiSize() == 3 )  zerror = pitcherr_z_mvtx*scalefactors_z[1];
      else if( cluster->getZSize() == 3 && cluster->getPhiSize() == 2 )  zerror = pitcherr_z_mvtx*scalefactors_z[2];
      else if( cluster->getZSize() == 3 && cluster->getPhiSize() == 3 )  zerror = pitcherr_z_mvtx*scalefactors_z[3];
      else zerror = pitcherr_z_mvtx*cluster->getZSize();
      zerror *= scale_mvtx_z;
      break;
      
    case TrkrDefs::inttId:
      static constexpr std::array<double, 3> scalefactors_intt_phi = {{ 0.85, 0.4, 0.33 }};
      
      if( cluster->getPhiSize() == 1 && layer < 5) phierror = pitcherr_phi_intt*scalefactors_intt_phi[0];
      else if( cluster->getPhiSize() == 2 && layer < 5) phierror = pitcherr_phi_intt*scalefactors_intt_phi[1];
      else if( cluster->getPhiSize() == 2 && layer > 4) phierror = pitcherr_phi_intt*scalefactors_intt_phi[2];
      else phierror = pitcherr_phi_intt*cluster->getPhiSize();
      if( layer == 3) phierror *= scale_intt_3;
      if( layer == 4) phierror *= scale_intt_4;
      if( layer == 5) phierror *= scale_intt_5;
      if( layer == 6) phierror *= scale_intt_6;

      zerror = pitcherr_z_intt*cluster->getZSize();
      break;
      
    }
  if(phierror==0) phierror = 100;
  if(zerror==0) zerror = 100;

  return std::make_pair(square(phierror),square(zerror));  
  //  return std::make_pair(phierror,zerror);  

}
