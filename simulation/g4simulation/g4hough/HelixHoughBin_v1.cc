#include "HelixHoughBin_v1.h"

#include <phool/phool.h>  // for PHWHERE

#include <cstdlib>       // for exit, NULL
#include <iostream>

using namespace std;

HelixHoughBin_v1::HelixHoughBin_v1(unsigned int bin)
  : _cluster_IDs(),
    _global_bin(),
    _bin(),
    _kappa_bins(),
    _phi_bins(),
    _phi_high_bins(),
    _phi_low_bins(),
    _d_bins(),
    _dzdl_bins(),
    _dzdl_high_bins(),
    _dzdl_low_bins(),
    _z0_bins(),
    _zoomlevel(0),
    _hough_space(NULL) {

  set_bin(_zoomlevel,bin);
//  cout<<"HoughBin:: bin "<<bin<<endl;
}

void HelixHoughBin_v1::identify(std::ostream& os) const {
  os << "---HelixHoughBin_v1-------------" << endl;
  os << "zoom level: "<< _zoomlevel<<" bin :" <<_bin<<endl;
  os << " helix hough bins (kappa, phi, d, dzdl, z0) = (" 
     << get_kappa_bin(_zoomlevel) << ","
     << get_phi_bin(_zoomlevel) << ","
     << get_d_bin(_zoomlevel) << ","
     << get_dzdl_bin(_zoomlevel) << ","
     << get_z0_bin(_zoomlevel) << ")" << endl;
  os << "--------------------------------" << endl;
}

void HelixHoughBin_v1::init() {

//  cout<<"HoughBin:: zoom level "<<_zoomlevel<<" bin "<< _bin[_zoomlevel]<<endl;
  set_bins(_zoomlevel,_bin[_zoomlevel]);

}

void HelixHoughBin_v1::set_hough_space(HelixHoughSpace* hough_space) {
  
  _hough_space = hough_space->Clone();
}


void HelixHoughBin_v1::set_bins(unsigned int zoomlevel, unsigned int bin) {

 if (!_hough_space) 
 {
   cerr << PHWHERE << "::Error - HelixHoughSpace is not set!! "<< endl;
   exit(1); 
 }

 unsigned int n_kappa_bins = _hough_space->get_n_kappa_bins(zoomlevel); 
 unsigned int n_phi_bins = _hough_space->get_n_phi_bins(zoomlevel); 
 unsigned int n_d_bins = _hough_space->get_n_d_bins(zoomlevel);
 unsigned int n_dzdl_bins = _hough_space->get_n_dzdl_bins(zoomlevel);

 unsigned int temp1 = n_kappa_bins * n_phi_bins * n_d_bins * n_dzdl_bins;
 unsigned int z0_bin = bin/temp1;

 unsigned int temp2 = bin - z0_bin * temp1;
  	      temp1 /= n_dzdl_bins;
 unsigned int dzdl_bin = temp2/temp1;

	      temp2 -= dzdl_bin * temp1;
	      temp1 /= n_d_bins;	      
 unsigned int d_bin = temp2/temp1;
 	      
	      temp2 -= d_bin * temp1;
	      temp1 /= n_phi_bins;
 unsigned int phi_bin = temp2/temp1;
	
	      temp2 -= phi_bin * temp1;
	      temp1 /= n_kappa_bins;
 unsigned int kappa_bin = temp2/temp1;


//   cout<<"HoughBin:: bin "<<bin<<endl;
//   cout<<"k "<<kappa_bin<<" phi " <<phi_bin<<" d "<<d_bin<<" dzdl "<< dzdl_bin<<" z0 "<< z0_bin<<endl;

 set_kappa_bin(zoomlevel, kappa_bin);
 set_phi_bin(zoomlevel, phi_bin); 
 set_d_bin(zoomlevel, d_bin);
 set_dzdl_bin(zoomlevel, dzdl_bin);
 set_z0_bin(zoomlevel, z0_bin);
 set_global_bin(zoomlevel);

// bin = kappabin + nkappabin* (phibin + nphibin* (dbin + ndbin *( dzdlbin + ndzdlbin*z0bin) ) )
} 

unsigned int HelixHoughBin_v1::get_global_bin(unsigned int zoomlevel)
{
	return _global_bin;
}

void HelixHoughBin_v1::set_global_bin(unsigned int zoomlevel)
{

// global_bin = zoom2bin + nzoom2bin*(zoom1bin + nzoom1bin*zoom0bin)
	_global_bin = _bin[0];
//	cout<<"global bin init"<<_global_bin<<endl;
	for (unsigned int izoom=1; izoom<zoomlevel+1; izoom++){
	unsigned int ntotalbins = _hough_space->get_n_kappa_bins(izoom);
	ntotalbins *= _hough_space->get_n_phi_bins(izoom);
	ntotalbins *= _hough_space->get_n_d_bins(izoom);
	ntotalbins *= _hough_space->get_n_dzdl_bins(izoom);
	ntotalbins *= _hough_space->get_n_z0_bins(izoom);
	_global_bin += _bin[izoom] + ntotalbins * _global_bin;
//	cout<<"global bin first "<<_global_bin<<endl;
	}
}


unsigned int HelixHoughBin_v1::get_neighbors_global_bin(unsigned int zoomlevel, unsigned int var, unsigned int bit_sign )
{

	unsigned int nbins[5]= {_hough_space->get_n_kappa_bins(zoomlevel),_hough_space->get_n_phi_bins(zoomlevel),_hough_space->get_n_d_bins(zoomlevel),_hough_space->get_n_dzdl_bins(zoomlevel),_hough_space->get_n_z0_bins(zoomlevel)};
	unsigned int bins[5] = {get_kappa_bin(zoomlevel),get_phi_bin(zoomlevel),get_d_bin(zoomlevel),get_dzdl_bin(zoomlevel),get_z0_bin(zoomlevel)};
//	cout<<"bins[] "<<bins[var]<<endl;
//	unsigned int bins_var[5] = {get_kappa_bin(zoomlevel),get_phi_bin(zoomlevel),get_d_bin(zoomlevel),get_dzdl_bin(zoomlevel),get_z0_bin(zoomlevel)};
	for (unsigned int i=0;i<5;i++)
	{
	int sign = 2*((bit_sign>>i)%2)-1;
	unsigned int bit_var = (var>>i)%2;
	if (!bit_var || (bins[i]==(nbins[i]-1)&& sign>0) || (bins[i]==0 && sign<0) ) continue;
	bins[i] += sign;
	}
	unsigned int _neighbor_bin = _hough_space->get_bin(zoomlevel,bins);
//	cout<<"sign "<<sign<< "bins[]"<<bins[var]<<endl;

	unsigned int _global_neighbor_bin;
	if (zoomlevel==0) _global_neighbor_bin = _neighbor_bin;
	else _global_neighbor_bin = _bin[0];
	for (unsigned int izoom=1; izoom<zoomlevel+1; izoom++){
	unsigned int ntotalbins = _hough_space->get_n_kappa_bins(izoom);
	ntotalbins *= _hough_space->get_n_phi_bins(izoom);
	ntotalbins *= _hough_space->get_n_d_bins(izoom);
	ntotalbins *= _hough_space->get_n_dzdl_bins(izoom);
	ntotalbins *= _hough_space->get_n_z0_bins(izoom);
	if (izoom == zoomlevel){
	_global_neighbor_bin += _neighbor_bin + ntotalbins * _global_neighbor_bin;
	}else{
	_global_neighbor_bin += _bin[izoom] + ntotalbins * _global_neighbor_bin;
	}
	}
	return _global_neighbor_bin;
}

float HelixHoughBin_v1::get_kappa_center(unsigned int zoomlevel)
{
        float center = _hough_space->get_kappa_min();
        for (unsigned int i=0; i<zoomlevel+1; i++){
        center += _hough_space->get_kappa_bin_size(i)*_kappa_bins[i];        
        }
        center += 0.5*_hough_space->get_kappa_bin_size(zoomlevel);
        return center;
}

float HelixHoughBin_v1::get_phi_center(unsigned int zoomlevel)
{
        float center = _hough_space->get_phi_min();
        for (unsigned int i=0; i<zoomlevel+1; i++){
        center += _hough_space->get_phi_bin_size(i)*_phi_bins[i];
        }
        center += 0.5*_hough_space->get_phi_bin_size(zoomlevel);
        return center;
}

float HelixHoughBin_v1::get_d_center(unsigned int zoomlevel)
{
        float center = _hough_space->get_d_min();
        for (unsigned int i=0; i<zoomlevel+1; i++){
        center += _hough_space->get_d_bin_size(i)*_d_bins[i];
        }
        center += 0.5*_hough_space->get_d_bin_size(zoomlevel);
        return center;
}

float HelixHoughBin_v1::get_dzdl_center(unsigned int zoomlevel)
{
        float center = _hough_space->get_dzdl_min();
        for (unsigned int i=0; i<zoomlevel+1; i++){
        center += _hough_space->get_dzdl_bin_size(i)*_dzdl_bins[i];
        }
        center += 0.5*_hough_space->get_dzdl_bin_size(zoomlevel);
        return center;
}

float HelixHoughBin_v1::get_z0_center(unsigned int zoomlevel)
{
        float center = _hough_space->get_z0_min();
        for (unsigned int i=0; i<zoomlevel+1; i++){
        center += _hough_space->get_z0_bin_size(i)*_z0_bins[i];
        }
        center += 0.5*_hough_space->get_z0_bin_size(zoomlevel);
        return center;
}


