//-------------------------------------------
// Read in a HepMC file with HIJING event
// output, and impose azimuthal modulation
// on individual tracks based on their
// transverse momentum and event centrality,
// writing out a new HepMC file
//
// JDOK
// 06-06-18
//-------------------------------------------

#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TStyle.h"

using namespace std;

//---------------------------------
// Variables
//---------------------------------

struct particle {
	int barcode;
	int pdg_id;
	float px;
	float py;
	float pz;
	float energy;
	float mass;
	int status;
};

//Parameterization of v2(pT) for each centrality
TF1 *f_v2_0_20;
TF1 *f_v2_20_40;
TF1 *f_v2_40_60;

//Values of impact parameter that define centrality categories
const float B_0_20 = 6.9;
const float B_20_40 = 9.7;
const float B_40_60 = 11.9;

//File to write out
ofstream myfile;

//Reaction plane angle and impact parameter for the event at hand
float current_psi;
float current_b;

//---------------------------------
// Functions
//---------------------------------

/*
 * Initialize v2(pT) parameterization for each centrality class
 */
void initializeV2Parameterization()
{
	f_v2_0_20 = new TF1("f_v2_0_20", "(x < 1.346590909)*(0.0540357*x + -0.00531966) + (pol2 + [3]/TMath::Power(x,[4]))*(x >= 1.346590909 && x < 5.0) + (-0.0121271*x + 0.141288)*(x>=5.0)", 0, 10);
	f_v2_0_20->SetParameter(0, -0.556786);
	f_v2_0_20->SetParameter(1, -0.113381);
	f_v2_0_20->SetParameter(2, 0.00408226);
	f_v2_0_20->SetParameter(3, 0.709262);
	f_v2_0_20->SetParameter(4, -0.273957);

	f_v2_20_40 = new TF1("f_v2_20_40", "(x < 1.266311167)*(0.0801973*x + 0.0398895) + (pol2 + [3]/TMath::Power(x,[4]))*(x >= 1.266311167 && x < 5.0) + (-0.0231989*x + 0.252772)*(x>=5.0)", 0, 10);
	f_v2_20_40->SetParameter(0, -0.855345);
	f_v2_20_40->SetParameter(1, -0.182768);
	f_v2_20_40->SetParameter(2, 0.00715066);
	f_v2_20_40->SetParameter(3, 1.14565);
	f_v2_20_40->SetParameter(4, -0.255077);

	f_v2_40_60 = new TF1("f_v2_40_60", "(x < 1.304806202)*(0.0554162*x + 0.0979441) + (pol2 + [3]/TMath::Power(x,[4]))*(x >= 1.304806202 && x < 5.0) + (-0.0173517*x + 0.250273)*(x>=5.0)", 0, 10);
	f_v2_40_60->SetParameter(0, -0.137921);
	f_v2_40_60->SetParameter(1, -0.413409);
	f_v2_40_60->SetParameter(2, 0.0102499);
	f_v2_40_60->SetParameter(3, 0.689962);
	f_v2_40_60->SetParameter(4, -0.695189);
}


/*
 * Transcendental equation from solving a differential equation, as described in the documentation
 * This equation implicitly defines the desired transformation on phi for a given v2
 */
float f(float x, float phi, float v2)
{
	float f = x + v2 * TMath::Sin(2 * x) - phi - 0.5;
	return f;
}


/*
 * Implementation of the secant root-finding method to solve the transcendental equation defined in function f(x, phi, v2)
 */
float secant(float phi, float v2)
{
	//Code from https://www.geeksforgeeks.org/program-to-find-root-of-an-equations-using-secant-method/

	//x1 and x2 define the range where the solution is said to exist
	//E is a small error tolerance parameter
	float x1 = 0;
	float x2 = 20;
	float E = 1E-3;

	float n = 0, xm, x0, c;
	if (f(x1, phi, v2) * f(x2, phi, v2) < 0) {
		do {
			// calculate the intermediate value
			x0 = (x1 * f(x2, phi, v2) - x2 * f(x1, phi, v2)) / (f(x2, phi, v2) - f(x1, phi, v2));

			// check if x0 is root of equation or not
			c = f(x1, phi, v2) * f(x0, phi, v2);

			// update the value of interval
			x1 = x2;
			x2 = x0;

			// update number of iteration
			n++;

			// if x0 is the root of equation then break the loop
			if (c == 0)
				break;
			xm = (x1 * f(x2, phi, v2) - x2 * f(x1, phi, v2)) / (f(x2, phi, v2) - f(x1, phi, v2));
		} while (fabs(xm - x0) >= E); // repeat the loop until convergence

		return x0;
	}
	else
	{
		return -9999;
	}
}


/*
 * For a given particle, take the event impact parameter,
 * and the event plane angle and add flow modulations based on
 * particle pT and event centrality
 */
particle processParticle(particle p)
{
	float pT = TMath::Sqrt(p.px * p.px + p.py * p.py);
	float phi = TMath::ATan2(p.py, p.px);

	//If it's a very low-pT particle, don't add any flow
	if (pT < 0.05) return p;

	//The range of ATan2 is [-pi, pi], so add pi to distribute the particles within [0, 2pi]
	phi = phi + TMath::Pi();

	//Get the desired v2 for the particle depending on its pT and event centrality
	//Since we don't have a parameterization for the 60-100% class based on published data, use that for 40-60%
	float v2 = -999;
	if (current_b >= 0.0 && current_b < B_0_20)
	{
		v2 = f_v2_0_20->Eval(pT);
	}
	else if (current_b >= B_0_20 && current_b < B_20_40)
	{
		v2 = f_v2_20_40->Eval(pT);
	}
	else if(current_b >= B_20_40)
	{
		v2 = f_v2_40_60->Eval(pT);
	}

	//Apply mapping by numerically solving transcendental eq. with parameter v2
	//
	float phi_prime = secant(phi, v2);

	//Wrap phi_prime angle around
	if (phi_prime > 2 * TMath::Pi())
	{
		phi_prime = phi_prime - 2 * TMath::Pi();
	}
	else if (phi_prime < 0)
	{
		phi_prime = phi_prime + 2 * TMath::Pi();
	}

	p.px = pT * TMath::Cos(phi_prime);
	p.py = pT * TMath::Sin(phi_prime);

	return p;
}


/*
 * Take a HepMC file, from HIJING output, and parse it line by line.
 * For every event, the final-state particles are stored in a vector and processed to add flow modulations
 */
void hepMCFlowAfterburner(string input_file="sHijing.dat", string output_file="sHijing_wflow.dat")
{
	//Initialize v2 parameterizations for each centrality class
	initializeV2Parameterization();

	//Load in HepMC file and parse line by line
	ifstream infile(input_file.c_str());
	string line;

	//Open file to write out modified HepMC file
	myfile.open (output_file.c_str());

	int evtnumber = 0;

	while (getline(infile, line))
	{
		//If this is not a particle line, write it out to the modified file
		if (line.c_str()[0] != 'P')	myfile << line << endl;

		//Have we found a new event?
		if (line.c_str()[0] == 'E')
		{
			evtnumber++;
			
			//Look two lines down for impact parameter and event plane information
			getline(infile, line);
			myfile << line << endl;
			getline(infile, line);
			myfile << line << endl;

			string junk1;
			float n_hardscat;
			float n_proj_part;
			float n_targ_part;
			float n_coll;
			float n_spect_n;
			float n_spect_p;
			float n_wound1;
			float n_wound2;
			float n_wound3;
			float impact_par;
			float ep_angle;
			float eccen;
			float cross_sec;

			stringstream s(line);
			if (!(s >> junk1 >> n_hardscat >> n_proj_part >> n_targ_part >> n_coll >> n_spect_n >> n_spect_p >> n_wound1 >> n_wound2 >> n_wound3 >> impact_par >> ep_angle >> eccen >> cross_sec)) break;

			current_b = impact_par;
			current_psi = ep_angle;
		}

		//Have we found a particle?
		if (line.c_str()[0] == 'P')
		{
			string junk1;
			int barcode;
			int pid;
			float px, py, pz;
			float energy;
			float mass;
			int status;
			float pol_theta;
			float pol_phi;
			int barcode_vtx;
			int num_entries;

			stringstream s(line);
			if (!(s >> junk1 >> barcode >> pid >> px >> py >> pz >> energy >> mass >> status >> pol_theta >> pol_phi >> barcode_vtx >> num_entries)) break;

			//Store it in the vector containing the particles of the event at hand
			particle p;
			p.px = px;
			p.py = py;
			p.pz = pz;
			p.barcode = barcode;
			p.pdg_id = pid;
			p.energy = energy;
			p.mass = mass;
			p.status = status;

			particle p_modif = processParticle(p);
			myfile << "P " << p_modif.barcode << " " << p_modif.pdg_id << " " << p_modif.px << " " << p_modif.py << " " << p_modif.pz << " " << p_modif.energy << " " << p_modif.mass << " " << " " << p_modif.status << " 0 0 0 0" << endl;
		}
	}
}