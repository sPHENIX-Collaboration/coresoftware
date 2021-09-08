#include "HelixHoughSpace_v1.h"

using namespace std;

HelixHoughSpace_v1::HelixHoughSpace_v1()
  : _para_min(),
    _para_max(),
    _zoom_profile(),
    _max_zoom() {

  for (unsigned int i=0; i<5; i++) 
  {    
  _para_min[i] = -9999;
  _para_max[i] = -9999;

    for (unsigned int j=0; j<ZOOMLEVEL_MAX; j++)
    {
    _zoom_profile[j][i]= 1;  
    }
  }
  _max_zoom = 0;
}

HelixHoughSpace_v1::HelixHoughSpace_v1(const HelixHoughSpace_v1& hough_space) {
  *this = hough_space;
  return;
}

void HelixHoughSpace_v1::add_one_zoom(std::vector<unsigned int>& one_zoom) {

//  cout<<"one_zoom added "<< one_zoom[0]<<endl;
  _zoom_profile[_max_zoom][0] = one_zoom[0];
  _zoom_profile[_max_zoom][1] = one_zoom[1];
  _zoom_profile[_max_zoom][2] = one_zoom[2];
  _zoom_profile[_max_zoom][3] = one_zoom[3];
  _zoom_profile[_max_zoom][4] = one_zoom[4];
//  cout<<"zoom_profile "<<_zoom_profile[_max_zoom][0]<<endl;
  _max_zoom++;
//  cout<<"max_zoom "<< _max_zoom<<endl;
}

unsigned int HelixHoughSpace_v1::get_max_zoom() {

  return _max_zoom;

}

void HelixHoughSpace_v1::print_zoom_profile() {

  cout << "HelixHoughSpace::print_zoom_profile:" << endl;

  for (unsigned int i=0; i<5; i++) {
    for (unsigned int j=0; j<_max_zoom; j++) {
	cout << "[" << j << "][" << i << "] = " << _zoom_profile[j][i] << endl;
    }
  }
}

void HelixHoughSpace_v1::print_para_range() {

    cout << "HelixHoughSpace::print_para_range:" <<endl;

    cout << "kappa max = " << _para_max[0] << endl;
    cout << "kappa min = " << _para_min[0] << endl;
    cout << "phi max = " << _para_max[1] << endl;
    cout << "phi min = " << _para_min[1] << endl;
    cout << "d max = " << _para_max[2] << endl;
    cout << "d min = " << _para_min[2] << endl;
    cout << "dzdl max = " << _para_max[3] << endl;
    cout << "dzdl min = " << _para_min[3] << endl;
    cout << "z0 max = " << _para_max[4] << endl;
    cout << "z0 min = " << _para_min[4] << endl;

}

float HelixHoughSpace_v1::get_kappa_bin_size(unsigned int zoomlevel) const {
	
	float binsize = _para_max[0]-_para_min[0];
	for (unsigned int izoom = 0; izoom<zoomlevel+1; izoom++){
	binsize /= (float) _zoom_profile[izoom][0];
	}
	return binsize;
}

float HelixHoughSpace_v1::get_phi_bin_size(unsigned int zoomlevel) const {

        float binsize = _para_max[1]-_para_min[1];
        for (unsigned int izoom = 0; izoom<zoomlevel+1; izoom++){
        binsize /= (float) _zoom_profile[izoom][1];
        }
        return binsize;
}

float HelixHoughSpace_v1::get_d_bin_size(unsigned int zoomlevel) const {

        float binsize = _para_max[2]-_para_min[2];
        for (unsigned int izoom = 0; izoom<zoomlevel+1; izoom++){
        binsize /= (float) _zoom_profile[izoom][2];
        }
        return binsize;
}

float HelixHoughSpace_v1::get_dzdl_bin_size(unsigned int zoomlevel) const {

        float binsize = _para_max[3]-_para_min[3];
        for (unsigned int izoom = 0; izoom<zoomlevel+1; izoom++){
        binsize /= (float) _zoom_profile[izoom][3];
        }
        return binsize;
}

float HelixHoughSpace_v1::get_z0_bin_size(unsigned int zoomlevel) const {

        float binsize = _para_max[4]-_para_min[4];
//	cout<<"z0_bin_size "<<binsize<<endl;
        for (unsigned int izoom = 0; izoom<zoomlevel+1; izoom++){
        binsize /= (float) _zoom_profile[izoom][4];
        }
//	cout<<"z0_bin_size after "<<binsize<<endl;
        return binsize;
}
/*
float HelixHoughSpace_v1::get_kappa_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ik) const {
	float center = get_kappa_min();
	for (unsigned int i=0; i<zoomlevel+1; i++){
	center += get_kappa_bin_size(i)*v_ik[i]; 	
	}
	center += 0.5*get_kappa_bin_size(zoomlevel);
	return center;
}

float HelixHoughSpace_v1::get_phi_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ip) const {
	float center = get_phi_min();
	for (unsigned int i=0; i<zoomlevel+1; i++){
	center += get_phi_bin_size(i)*v_ip[i];
	}
	center += 0.5*get_phi_bin_size(zoomlevel);
	return center;
}
float HelixHoughSpace_v1::get_d_center(unsigned int zoomlevel, std::vector<unsigned int>& v_id) const {
	float center = get_d_min();
	for (unsigned int i=0; i<zoomlevel+1; i++){
	center += get_d_bin_size(i)*v_id[i];
	}
	center += 0.5*get_d_bin_size(zoomlevel);
	return center;
}

float HelixHoughSpace_v1::get_dzdl_center(unsigned int zoomlevel, std::vector<unsigned int>& v_il) const {
	float center = get_dzdl_min();
	for (unsigned int i=0; i<zoomlevel+1; i++) {
	center += get_dzdl_bin_size(i)*v_il[i];
	}
	center += 0.5*get_dzdl_bin_size(zoomlevel);
	return center;
}
float HelixHoughSpace_v1::get_z0_center(unsigned int zoomlevel, std::vector<unsigned int>& v_iz) const {
	float center = get_z0_min();
	for (unsigned int i=0; i<zoomlevel+1; i++){
	center += get_z0_bin_size(i)*v_iz[i];
	}
	center += 0.5*get_z0_bin_size(zoomlevel);
	return center;
}
*/
unsigned int HelixHoughSpace_v1::get_kappa_bin(unsigned int zoomlevel, float kappa) const {
	unsigned int kbin=0;
	kappa -= get_kappa_min();
	for (unsigned int i=0; i<zoomlevel+1; i++ ) {
	kbin = kappa/get_kappa_bin_size(i);
        kappa -= (kbin)*get_kappa_bin_size(i);
	}
	return kbin;
}

unsigned int HelixHoughSpace_v1::get_phi_bin(unsigned int zoomlevel, float phi) const {
	unsigned int pbin =0;
//	cout<<"phi "<< phi<<endl;
	phi -= get_phi_min();
//	cout<<"phi-phi_min "<<phi<<endl;
	for (unsigned int i=0; i<zoomlevel+1; i++) {
//	cout<<"phi-phi_min-phi_prev "<< phi<<endl;
//	cout<<"phi_bin_size "<<get_phi_bin_size(i)<<endl;
	pbin = phi/get_phi_bin_size(i);
        phi -= (pbin)*get_phi_bin_size(i);
	}

	return pbin;
}

unsigned int HelixHoughSpace_v1::get_d_bin(unsigned int zoomlevel, float d) const {
	unsigned int dbin =0;
	d -= get_d_min();
	for (unsigned int i=0; i<zoomlevel+1; i++) {
	dbin = d/get_d_bin_size(i);
	d -= (dbin)*get_d_bin_size(zoomlevel);
	}
	return dbin;
}

unsigned int HelixHoughSpace_v1::get_dzdl_bin(unsigned int zoomlevel, float dzdl) const {
	unsigned int lbin =0;
//	cout<<"dzdl "<<dzdl<<endl;
//	cout<<"dzdl_min "<<get_dzdl_min()<<endl;
	dzdl -= get_dzdl_min();
//	cout<<"dzdl-dzdl_min "<<dzdl<<endl;
	for (unsigned int i=0; i<zoomlevel+1; i++) {
//	dzdl -= (lbin)*get_dzdl_bin_size(i);
	lbin = dzdl/get_dzdl_bin_size(i);
	dzdl -= (lbin)*get_dzdl_bin_size(i);
//	cout<<"dzdl_bin_size "<<get_dzdl_bin_size(i)<<endl;
//	cout<<"dzdl-prevzoom "<<dzdl<<endl;
// 	if (zoomlevel==1) cout<<"i "<<i<<" lbin "<<lbin<<endl;
	if (lbin >= get_n_dzdl_bins(i)) return 999;
	if (lbin<0) return 9999;
	}	
//	cout<<"lbin "<<lbin<<endl;
	return lbin;
}

unsigned int HelixHoughSpace_v1::get_z0_bin(unsigned int zoomlevel, float z0) const {
	unsigned int zbin =0;
	for (unsigned int i =0; i<zoomlevel+1; i++) {
	zbin = z0/get_z0_bin_size(i);
	z0 -= (zbin)*get_z0_bin_size(i);
	}
	return zbin;
}

unsigned int HelixHoughSpace_v1::get_bin(unsigned int zoomlevel, unsigned int* bins) const {

//    cout<<"HoughSpace:: zoom "<<zoomlevel<<" kappa "<<bins[0]<<endl;
    unsigned int bin = 0;
    for (unsigned int izoom = 0; izoom<zoomlevel; izoom++) bin += get_n_kappa_bins(izoom)*get_n_phi_bins(izoom)*get_n_d_bins(izoom)*get_n_dzdl_bins(izoom)*get_n_z0_bins(izoom);  
    bin = bins[0] + _zoom_profile[zoomlevel][0] * (bins[1] + _zoom_profile[zoomlevel][1] * (bins[2] + _zoom_profile[zoomlevel][2] *( bins[3] + _zoom_profile[zoomlevel][3] * bins[4])));
    return bin;

// bin = kappabin + nkappabin* (phibin + nphibin* (dbin + ndbin *( dzdlbin + ndzdlbin*z0bin) ) )

}
/*
unsigned int HelixHoughSpace_v1::get_bin() {
  cout<<"HoughSpace:: zoom "<<zoomlevel<<" kappa "<<kappa_bin<<endl;
  unsigned int bin =0; 
 //= kappa_bin + _zoom_profile[zoomlevel][0] * (phi_bin + _zoom_profile[zoomlevel][1] * (d_bin + _zoom_profile[zoomlevel][2] *( dzdl_bin + _zoom_profile[zoomlevel][3] * z0_bin)));
  cout<<"HoughSpace:: bin "<<bin<<endl;
  return bin;
}
*/
