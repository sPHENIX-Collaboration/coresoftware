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
  f0 = new TF1("f0","[0] + x*[1] + x*x*[2]*[2]",0,10);
  f0->SetParameter(0,0.0169157);
  f0->SetParameter(1,0.0187465);
  f0->SetParameter(2,0.0199696);

  f1 = new TF1("f1","[0] + x*[1] + x*x*[2]*[2]",0,10);
  f1->SetParameter(0,0.0122044);
  f1->SetParameter(1,0.0159068);
  f1->SetParameter(2,0.367468);

  f2 = new TF1("f2","[0] + x*[1] + x*x*[2]*[2] ",0,10);
  f2->SetParameter(0,0.0118987);
  f2->SetParameter(1,0.0302875);
  f2->SetParameter(2,0.0271744);

  fz0 = new TF1("fz0","pol1",-2,2);
  fz0->SetParameter(0,0.0531544);
  fz0->SetParameter(1,0.0181959);

  fz = new TF1("fz","pol2",-2,2);
  fz->SetParameter(0,0.0404929);
  fz->SetParameter(1,-0.00159808);
  fz->SetParameter(2,0.0617323);

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

  fadcz0 = new TF1("fadcz0","[0]+([1]/pow(x-[2],2))",0,20000);
  fadcz0->SetParameter(0,0.479638);
  fadcz0->SetParameter(1,28289.9);
  fadcz0->SetParameter(2,-126.307);

  fadcz1 = new TF1("fadcz1","[0]+([1]/pow(x-[2],2))",0,20000);
  fadcz1->SetParameter(0,0.508903);
  fadcz1->SetParameter(1,94583.5);
  fadcz1->SetParameter(2,-238.725);

  fadcz2 = new TF1("fadcz2","[0]+([1]/pow(x-[2],2))",0,20000);
  fadcz2->SetParameter(0,0.52688);
  fadcz2->SetParameter(1,57381.8);
  fadcz2->SetParameter(2,-173.369);

  fadcphi0 = new TF1("fadcphi0","[0]+([1]/pow(x-[2],2))",0,20000);
  fadcphi0->SetParameter(0,0.409897);
  fadcphi0->SetParameter(1,80395.38);
  fadcphi0->SetParameter(2,-229.33);

  fadcphi1 = new TF1("fadcphi1","[0]+([1]/pow(x-[2],2))",0,20000);
  fadcphi1->SetParameter(0,0.392476);
  fadcphi1->SetParameter(1,243084);
  fadcphi1->SetParameter(2,-380.402);

  fadcphi2 = new TF1("fadcphi2","[0]+([1]/pow(x-[2],2))",0,20000);
  fadcphi2->SetParameter(0,0.366911);
  fadcphi2->SetParameter(1,128396);
  fadcphi2->SetParameter(2,-243.454);

  fadcphi2fine = new TF1("fadcphi2fine","pol1",0,20000);
  fadcphi2fine->SetParameter(0,1.04);
  fadcphi2fine->SetParameter(1,-0.000325796);

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
  double beta = TMath::Abs(atan(seed->get_slope()));
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
      

      TrkrClusterv4 *clusterv4 = dynamic_cast<TrkrClusterv4 *>(cluster);

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
	//	phierror *=  scale_tpc_0;
	zerror = fz0->Eval(beta);
	//zerror *=  scale_tpc_0_z;
	zerror *= fadcz0->Eval(clusterv4->getAdc());
	phierror *= fadcphi0->Eval(clusterv4->getAdc());
      }
      if(sector==1){
	phierror = f1->Eval(alpha);
	//	phierror *=  scale_tpc_1;
	zerror = fz->Eval(beta);
	//zerror *=  scale_tpc_1_z;
	zerror *= fadcz1->Eval(clusterv4->getAdc());
	phierror *= fadcphi1->Eval(clusterv4->getAdc());

      }
      if(sector==2){
	phierror = f2->Eval(alpha);
	//	phierror *=  scale_tpc_2;
	zerror = fz->Eval(beta);
	zerror *= fadcz2->Eval(clusterv4->getAdc());
	phierror *= fadcphi2->Eval(clusterv4->getAdc());
	phierror *= fadcphi2fine->Eval(clusterv4->getAdc());
	//	zerror *=  scale_tpc_2_z;

      }

      if(clusterv4->getEdge()>=1)
	phierror *= 2;
      if(clusterv4->getPhiSize()==1)
	phierror *= 2;
      if(clusterv4->getPhiSize()==4)
	phierror *= 2;
      if(clusterv4->getPhiSize()>=5)
	phierror *= 3;
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
