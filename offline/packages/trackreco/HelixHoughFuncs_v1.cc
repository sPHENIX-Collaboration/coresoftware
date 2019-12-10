#include "HelixHoughFuncs_v1.h"

#include <cmath>
#include <vector>

using namespace std;

HelixHoughFuncs_v1::HelixHoughFuncs_v1()
 : _hough_space(nullptr),_cur_zoom(0)
{

}

HelixHoughFuncs_v1::HelixHoughFuncs_v1(const HelixHoughFuncs_v1& hough_funcs) {
  *this = hough_funcs;
  return;
}

void HelixHoughFuncs_v1::set_hough_space(HelixHoughSpace* hough_space) {

  _hough_space = hough_space;
}


void HelixHoughFuncs_v1::calculate_dzdl_range(float* hitpos3d, std::vector<float>& z0_range, std::vector<float>& kappa_phi_d_ranges, float* dzdl_range) {

	float x = hitpos3d[0]; 
	float y = hitpos3d[1]; 
	float z = hitpos3d[2];
//	_hough_space->print_para_range();
 	// iz has to be an array of iz[zoomlevel], sot that it contains history
	float z0_min = z0_range[0]; 
	float z0_max = z0_range[1];
	float kappa_min = kappa_phi_d_ranges[0];
	float kappa_max = kappa_phi_d_ranges[1];
	float phi_min = kappa_phi_d_ranges[2];
	float phi_max = kappa_phi_d_ranges[3];

	float cosphi_min = cos(phi_min); 
	float cosphi_max = cos(phi_max);
	float sinphi_min = sin(phi_min);
	float sinphi_max = sin(phi_max);

	if (phi_max -phi_min > M_PI) { // for voting z only
	cosphi_min = 1.;
	cosphi_max = -1;
	sinphi_min = 0.;
	sinphi_max = 0.;
	}

	float d_min = kappa_phi_d_ranges[4]; 
	float d_max = kappa_phi_d_ranges[5];  

//        cout<<"x "<<x <<" z0_min "<<z0_min<<" cosphi_min "<<cosphi_min<<" cosphi_max "<<cosphi_max<<endl;
//	cout<<"d_min "<<d_min<<" d_max "<<d_max<<endl;
//	cout<<"kappa_min "<<kappa_min<<" kappa_max "<<kappa_max<<endl;
//	cout<<"phi_min "<<phi_min<<" phi_max "<<phi_max<<endl; 
	float d, k, dx, dy;
	float dzdl[4]={-999,-999,-999,-999};
	for (int icase=0; icase<2; icase++){
		switch (icase){
			case 0:
			d= d_min; dx=d*cosphi_min; dy=d*sinphi_min;k=kappa_min; 
			break;
			case 1:
			d= d_max; dx=d*cosphi_max; dy=d*sinphi_max;k=kappa_max;
			break;
		}

		float D = sqrt(pow(x-dx,2)+pow(y-dy,2));
		float a = k*D/2.;float a2 = a*a;
		float s = -999.;;
		if (a > .999) {a = 0.999; s = 2*asin(a)/k; }
		else if (a < 0.01){ s = D + a2*D/6. + 3/40.*a2*a2*D + 5/112.*a2*a2*a2;}
		else { s = 2*asin(a)/k;} 

		switch (icase){
			case 0:
                        dzdl[0] = pow(z-z0_max,2)/(pow(z-z0_max,2)+s*s);
			dzdl[0] = sqrt(dzdl[0]);
                        if (z < z0_max) dzdl[0] *= -1;
                        dzdl[1] = pow(z-z0_min,2)/(pow(z-z0_min,2)+s*s);
			dzdl[1] = sqrt(dzdl[1]);
                        if (z < z0_min) dzdl[1] *= -1;
//			cout<<"dzdl[0] "<< dzdl[0]<<endl;
//			cout<<"dzdl[1] "<<dzdl[1]<<endl;
			break;

			case 1:
                        dzdl[2] = pow(z-z0_max,2)/(pow(z-z0_max,2)+s*s);
			dzdl[2] = sqrt(dzdl[2]);
                        if (z < z0_max) dzdl[2] *= -1;
                        dzdl[3] = pow(z-z0_min,2)/(pow(z-z0_min,2)+s*s);
			dzdl[3] = sqrt(dzdl[3]);
                        if (z < z0_min) dzdl[3] *= -1;
//			cout<<"dzdl[2] "<<dzdl[2]<<endl;
//			cout<<"dzdl[3] "<<dzdl[3]<<endl;
			break;
		}
//		if (z0_min<5. && z0_max>5.){
//		cout<<"x "<<x<<" y "<<y<< " z "<<z<<endl; 
//		cout<<"d "<<d<<" k "<<k<<" dx "<<dx<<" dy "<<dy<<" z-z0_max "<<z-z0_max<<" z-z0_min "<<z-z0_min<<endl;
//		cout<<"D "<<D<<" s "<<s<<endl;
//		}
	}
	float dzdl_min = -999;
	float dzdl_max = -999;

//	if (z0_min<5. && z0_max>5.){
//	cout<<"dzdl (dmin,kmin, phimin) max "<<dzdl[0]<<" dzdl (dmin,kmin,phimin) min "<<dzdl[1]<<endl;
//	cout<<"dzdl (dmax,kmax, phimax) max "<<dzdl[2]<<" dzdl (dmax,kmax,phimax) min "<<dzdl[3]<<endl;
//	}

	dzdl_min = dzdl[0];
	for (int i = 1; i<4; i++) {
	if (dzdl_min > dzdl[i])
	dzdl_min = dzdl[i];
	}
	dzdl_max = dzdl[0];
	for (int j = 1; j<4; j++) {
	if (dzdl_max < dzdl[j])
	dzdl_max = dzdl[j];
	}
//	if (z0_min<5. && z0_max>5.)
//	if (_cur_zoom==1)
//	cout<<"dzdl_min "<<dzdl_min<<" dzdl_max "<<dzdl_max<<endl;

//	dzdl_range[0]= _hough_space->get_dzdl_bin(_cur_zoom,dzdl_min); // min bin
//	dzdl_range[1]= _hough_space->get_dzdl_bin(_cur_zoom,dzdl_max); // max bin
	dzdl_range[0] = dzdl_min;
	dzdl_range[1] = dzdl_max;

}

// if not separate by helicity
void HelixHoughFuncs_v1::calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges,  float* phi_r_range, float* phi_l_range){

	float x = hitpos2d[0];
	float y = hitpos2d[1];
	float hitphi= atan2(y,x); 
	if (hitphi < 0.) hitphi += 2.*M_PI;
	
	float kmin = kappa_d_ranges[0];
	float kmax = kappa_d_ranges[1];
	float dmin = kappa_d_ranges[2];
	float dmax = kappa_d_ranges[3];

//	cout<<"kmin "<<kmin<<" kmax "<<kmax<<" dmin "<< dmin<<" dmax "<<dmax<<endl;	

	float d, k; 
	float D, Dinv, ak, hk, hksq, xk1, xk2, yk1,yk2;
	float phi[8] = {-999,-999,-999,-999,-999,-999,-999,-999};
	
	for (int icase=0; icase<4; icase++ ){
		switch (icase){
			case 0:
			d= dmin; k=kmin;
			break;
			case 1:
			d= dmin; k=kmax;
			break;
			case 2:
			d= dmax; k=kmin;
			break;
			case 3:
			d= dmax; k=kmax;
			break;
		}

		// calc phi business
		D = sqrt(x*x+y*y);
		Dinv = 1./D;	
		ak = (2.*d + d*d*k + D*D*k)/2.*Dinv;
		hksq = pow(d*k+1.,2) - ak*ak;
		if (hksq>=0.) {
		hk=sqrt(hksq); 
		xk1 = (ak*x + hk*y)*Dinv;
		xk2 = (ak*x - hk*y)*Dinv;
		yk1 = (ak*y - hk*x)*Dinv;
		yk2 = (ak*y + hk*x)*Dinv;
		}

		switch (icase){
			case 0:
			if (hksq<0.){
			phi[0]= hitphi;//right
			phi[1]= hitphi;//left
			}else {
			phi[0]= atan2(yk1,xk1);// right
			if (phi[0]<0.) phi[0]+= 2.*M_PI;
			phi[1]= atan2(yk2,xk2);// left
			if (phi[1]<0.) phi[1]+= 2.*M_PI;
			}
			break;
		
			case 1:
			if (hksq<0.){
			phi[2]= hitphi;
			phi[3]= hitphi;
			}else {
			phi[2]=atan2(yk1,xk1);// right
			if (phi[2]<0.) phi[2]+= 2.*M_PI;
			phi[3]= atan2(yk2,xk2);// left
			if (phi[3]<0.) phi[3]+= 2.*M_PI;
			}
			break;
	
			case 2:
			if (hksq<0.){
			phi[4]= hitphi;
			phi[5]= hitphi;
			}else{
			phi[4]=atan2(yk1,xk1);//right
			if (phi[4]<0.) phi[4]+= 2.*M_PI;
			phi[5]=atan2(yk2,xk2);//left
			if (phi[5]<0.) phi[5]+= 2.*M_PI;
			}
			break;

			case 3:
			if (hksq<0.){
			phi[6]= hitphi;
			phi[7]= hitphi;
			} else{
			phi[6]=atan2(yk1,xk1);//right
			if (phi[6]<0.) phi[6]+= 2.*M_PI;
			phi[7]=atan2(yk2,xk2);//left
			if (phi[7]<0.) phi[7]+= 2.*M_PI; 
			}
			break;
		}	

	}	
//	for (int i = 0; i<8;i++)
//	cout<<"phi["<<i<<"]="<<phi[i]<<endl;
	float twopi = 2.*M_PI;
	float piover2=M_PI/2.;
	float threepiover2=3.*M_PI/2.;
	bool first = (phi[0]<piover2) || (phi[2]<piover2) || (phi[4]<piover2 || (phi[6]<piover2));
	bool fourth = (phi[0]>threepiover2) || (phi[2]>threepiover2) || (phi[4]>threepiover2) || (phi[6]>threepiover2);

	if (first && fourth){
		if (phi[0]>threepiover2) phi[0] -= twopi;
		if (phi[2]>threepiover2) phi[2] -= twopi;
		if (phi[4]>threepiover2) phi[4] -= twopi;
		if (phi[6]>threepiover2) phi[6] -= twopi;
	}

	float phi_r_min = phi[0]; 
	for (int i =1; i<4; i++){
		if(phi[2*i]<phi_r_min) phi_r_min= phi[2*i];
	}
	float phi_r_max	= phi[0];
	for (int i =1; i<4; i++){
		if(phi[2*i]>phi_r_max) phi_r_max=phi[2*i];
	}

	phi_r_range[0] = phi_r_min;
	phi_r_range[1] = phi_r_max;


        first = (phi[1]<piover2) || (phi[3]<piover2) || (phi[5]<piover2 || (phi[7]<piover2));
        fourth = (phi[1]>threepiover2) || (phi[3]>threepiover2) || (phi[5]>threepiover2) || (phi[7]>threepiover2);

        if (first && fourth){
                if (phi[1]>threepiover2) phi[1] -= twopi;
                if (phi[3]>threepiover2) phi[3] -= twopi;
                if (phi[5]>threepiover2) phi[5] -= twopi;
                if (phi[7]>threepiover2) phi[7] -= twopi;
        }

        float phi_l_min = phi[1];
        for (int i =1; i<4; i++){
                if(phi[2*i+1]<phi_l_min) phi_l_min= phi[2*i+1];
        }
        float phi_l_max = phi[1];
        for (int i =1; i<4; i++){
                if(phi[2*i+1]>phi_l_max) phi_l_max=phi[2*i+1];
        }

	phi_l_range[0] = phi_l_min;
	phi_l_range[1] = phi_l_max;	

//        if (dmin<.02 && dmax>0.02)
//	cout<<"phi_l_min "<<phi_l_min<<" phi_l_max "<<phi_l_max<<" phi_r_min "<<phi_r_min<<" phi_r_max "<<phi_r_max<<endl;

}

// if separage by helicity
// initial k
void HelixHoughFuncs_v1::calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_next_range){

        float x = hitpos2d[0];
        float y = hitpos2d[1];
        float hitphi= atan2(y,x);
        if (hitphi < 0.) hitphi += 2.*M_PI;

        float kmin = kappa_d_ranges[0];
        float kmax = kappa_d_ranges[1];
        float dmin = kappa_d_ranges[2];
        float dmax = kappa_d_ranges[3];

        float d, k;
        float D, Dinv, ak, hk=0., hksq, xk1, xk2, xk, yk1,yk2, yk;
	float cross, phi[4];

        for (int icase=0; icase<4; icase++ ){
                switch (icase){
                        case 0:
                        d= dmin; k=kmin;
                        break;
                        case 1:
                        d= dmin; k=kmax;
                        break;
                        case 2:
                        d= dmax; k=kmin;
                        break;
                        case 3:
                        d= dmax; k=kmax;
                        break;
                }

                // calc phi business
                D = sqrt(x*x+y*y);
                Dinv = 1./D;
                ak = (2.*d + d*d*k + D*D*k)/(2.*D);
                hksq = pow(d*k+1.,2) - ak*ak;
                if (hksq>=0.)
                hk=sqrt(hksq);
                xk1 = (ak*x + hk*y)*Dinv;
                xk2 = (ak*x - hk*y)*Dinv;
                yk1 = (ak*y - hk*x)*Dinv;
                yk2 = (ak*y + hk*x)*Dinv;
                

		cross = x*yk1 - y*xk1;
		if (cross*helicity>0){
		xk=xk1;
		yk=yk1;
		} else {
		xk=xk2;
		yk=yk2;
		}

                switch (icase){
                        case 0:
                        if (hksq<0.){
                        phi[0]= hitphi;//left
                        }else {
                        phi[0]= atan2(yk,xk);// left
                        if (phi[0]<0.) phi[0]+= 2.*M_PI;
                        }
                        break;

                        case 1: // min d, max k next_range[0]
                        if (hksq<0.){
                        phi[1]= hitphi;
                        }else {
                        phi[1]= atan2(yk,xk);// left
                        if (phi[1]<0.) phi[1]+= 2.*M_PI;
                        }
                        break;

                        case 2:
                        if (hksq<0.){
                        phi[2]= hitphi;
                        }else{
                        phi[2]=atan2(yk,xk);
                        if (phi[2]<0.) phi[2]+= 2.*M_PI;
                        }
                        break;

                        case 3: // max d, max k phi_next_range[1]
                        if (hksq<0.){
                        phi[3]= hitphi;
                        } else{
                        phi[3]=atan2(yk,xk);
                        if (phi[3]<0.) phi[3]+= 2.*M_PI;
                        }
                        break;
                }

        }// cases

	phi_next_range[0] = phi[1];
	phi_next_range[1] = phi[3];

        float twopi = 2.*M_PI;
        float piover2=M_PI/2.;
        float threepiover2=3.*M_PI/2.;
        bool first = (phi[0]<piover2) || (phi[1]<piover2) || (phi[2]<piover2 || (phi[3]<piover2));
        bool fourth = (phi[0]>threepiover2) || (phi[1]>threepiover2) || (phi[2]>threepiover2) || (phi[3]>threepiover2);

        if (first && fourth){
                if (phi[0]>threepiover2) phi[0] -= twopi;
                if (phi[1]>threepiover2) phi[1] -= twopi;
                if (phi[2]>threepiover2) phi[2] -= twopi;
                if (phi[3]>threepiover2) phi[3] -= twopi;
        }

	float phi_min = phi[0];
        for (int i =0; i<4; i++){
                if(phi[i]<phi_min) phi_min= phi[i];
        }
        float phi_max = phi[0];
        for (int i =0; i<4; i++){
                if(phi[i]>phi_max) phi_max=phi[i];
        }

        phi_range[0] = phi_min;
        phi_range[1] = phi_max;

}


void HelixHoughFuncs_v1::calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_prev_range, float* phi_next_range){

        float x = hitpos2d[0];
        float y = hitpos2d[1];
        float hitphi= atan2(y,x);
        if (hitphi < 0.) hitphi += 2.*M_PI;

//        float kmin = kappa_d_ranges[0];
        float kmax = kappa_d_ranges[1];
        float dmin = kappa_d_ranges[2];
        float dmax = kappa_d_ranges[3];

        float d, k;
        float D, Dinv, ak, hk=0., hksq, xk1, xk2, xk, yk1,yk2, yk;
        float cross, phi[4];


        for (int icase=0; icase<2; icase++ ){
                switch (icase){
                        case 0:
                        d= dmin; k=kmax;
                        break;
                        case 1:
                        d= dmax; k=kmax;
                        break;
                }

                // calc phi business
                                   D = sqrt(x*x+y*y);
                Dinv = 1./D;
                ak = (2.*d + d*d*k + D*D*k)/(2.*D);
                hksq = pow(d*k+1.,2) - ak*ak;
                if (hksq>=0.)
                hk=sqrt(hksq);
                xk1 = (ak*x + hk*y)*Dinv;
                xk2 = (ak*x - hk*y)*Dinv;
                yk1 = (ak*y - hk*x)*Dinv;
                yk2 = (ak*y + hk*x)*Dinv;


                cross = x*yk1 - y*xk1;
                if (cross*helicity>0){
                xk=xk1;
                yk=yk1;
                } else {
                xk=xk2;
                yk=yk2;
                }

                switch (icase){
                        case 0:
                        if (hksq<0.){
                        phi[0]= hitphi;//left
                        }else {
                        phi[0]= atan2(yk,xk);// left
                        if (phi[0]<0.) phi[0]+= 2.*M_PI;
                        }
                        break;

                        case 1: // min d, max k next_range[0]
                        if (hksq<0.){
                        phi[1]= hitphi;
                        }else {
                        phi[1]= atan2(yk,xk);// left
                        if (phi[1]<0.) phi[1]+= 2.*M_PI;
                        }
                        break;

                }

        }// cases

        phi_next_range[0] = phi[0];
        phi_next_range[1] = phi[1];
	phi[2] = phi_prev_range[0];
	phi[3] = phi_prev_range[1];

        float twopi = 2.*M_PI;
        float piover2=M_PI/2.;
        float threepiover2=3.*M_PI/2.;
        bool first = (phi[0]<piover2) || (phi[1]<piover2) || (phi[2]<piover2 || (phi[3]<piover2));
        bool fourth = (phi[0]>threepiover2) || (phi[1]>threepiover2) || (phi[2]>threepiover2) || (phi[3]>threepiover2);

        if (first && fourth){
                if (phi[0]>threepiover2) phi[0] -= twopi;
                if (phi[1]>threepiover2) phi[1] -= twopi;
                if (phi[2]>threepiover2) phi[2] -= twopi;
                if (phi[3]>threepiover2) phi[3] -= twopi;
        }

        float phi_min = phi[0];
        for (int i =0; i<4; i++){
                if(phi[i]<phi_min) phi_min= phi[i];
        }
        float phi_max = phi[0];
        for (int i =0; i<4; i++){
                if(phi[i]>phi_max) phi_max=phi[i];
        }

        phi_range[0] = phi_min;
        phi_range[1] = phi_max;

}
/*

float HelixHoughFuncs_v1::calc_dzdl_error() {
float dzdl_err = ;

return dzdl_err;
}

float HelixHoughFuncs_v1::calc_phi_error() {
float phi_err =;
 
return phi_err;
}

*/
