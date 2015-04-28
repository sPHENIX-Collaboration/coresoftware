//
// Inspired by code from ATLAS.  Thanks!
//
#include <iostream>
#include <cstdlib>
#include <string>

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenRanges.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "gsl/gsl_histogram.h"

int
main ()
{
  using boost::property_tree::ptree;
  ptree pt;

  std::ifstream config_file("test.xml");

  if (config_file)
    {
      // Read XML configuration file.
      read_xml (config_file, pt);
    }

  std::string input = pt.get("TEST.INPUT", "test.dat");

  // Try to open input file.
  std::ifstream istr(input.c_str());
  if (!istr)
    {
      std::cout << __PRETTY_FUNCTION__ << ": input file \"" 
		<< input << "\" not found!" << std::endl;
      return 1;
    }

  // Book a GSL histogram
  size_t n = 32;
  double a, b;
  double w = M_PI / n;
  a = -M_PI / 2.0 - w / 2.0;
  b = M_PI / 2.0 - w / 2.0;
  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h, a, b);  

  size_t N = 20;
  double A = 0.2, B = 5.0;
  gsl_histogram *v2h = gsl_histogram_alloc (N);
  gsl_histogram_set_ranges_uniform (v2h, A, B);  
  gsl_histogram *v2hn = gsl_histogram_alloc (N);
  gsl_histogram_set_ranges_uniform (v2hn, A, B);  

  HepMC::IO_GenEvent ascii_in (istr);
  HepMC::GenEvent *evt;

  while (ascii_in >> evt)
    {
      HepMC::HeavyIon *hi = evt->heavy_ion();
      HepMC::GenVertex *primary_vtx = evt->barcode_to_vertex(-1);
      double phi0 = hi->event_plane_angle();
	
      HepMC::GenVertexParticleRange r(*primary_vtx, HepMC::children);
      for (HepMC::GenVertex::particle_iterator it = r.begin (); it != r.end (); it++)
	{
	  HepMC::GenParticle *p1 = (*it);
	  if (p1->status() == 103)
	    {
	      continue;
	    }
	  double eta1 = p1->momentum().eta();
	  if (fabs(eta1) > 1.5)
	    {
	      continue;
	    }
	  double pt = p1->momentum().perp();
	  if (pt < 0.2 || pt > 10.0)
	    {
	      continue;
	    }
	  double phi1 = p1->momentum().phi();
	  double dphi = phi1 - phi0;
	  dphi = atan2(sin(dphi),cos(dphi));
	  gsl_histogram_accumulate(v2h, pt, cos(2*dphi));
	  gsl_histogram_increment(v2hn, pt);
	  gsl_histogram_increment(h, dphi);
	}
      
      delete evt;
      ascii_in >> evt;
    }
  
  std::cout << "pt,v2,v2e" << std::endl;
  gsl_histogram_div(v2h,v2hn);
  double lower, upper;
  for (size_t i = 0; i < N; i++)
    {
      gsl_histogram_get_range(v2h, i, &lower, &upper);
      double pt = (lower + upper) / 2.0;
      double counts = gsl_histogram_get(v2hn, i);
      double val = gsl_histogram_get(v2h, i);
      double err = val/sqrt(counts);
      std::cout << pt << ", " << val << ", " << err << std::endl;
    }
 
  return 0;
}

  
