//
// Inspired by code from ATLAS.  Thanks!
//
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenRanges.h>
#include <HepMC/GenVertex.h>
#include <HepMC/HeavyIon.h>
#include <HepMC/IO_BaseClass.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/IteratorRange.h>
#include <HepMC/SimpleVector.h>

// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp>  // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 8 && (BOOST_VERSION == 106000 || BOOST_VERSION == 106700 || BOOST_VERSION == 107000))
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma message "ignoring bogus gcc warning in boost header ptree.hpp"
#include <boost/property_tree/ptree.hpp>
#pragma GCC diagnostic warning "-Wshadow"
#pragma GCC diagnostic warning "-Wunused-parameter"
#else
#include <boost/property_tree/ptree.hpp>
#endif

#include <boost/operators.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <gsl/gsl_histogram.h>

#include <algorithm>  // for max
#include <cmath>      // for cos
#include <cstdlib>
#include <iostream>
#include <string>

int main()
{
  using boost::property_tree::ptree;
  ptree proptree;

  std::ifstream config_file("test.xml");

  if (config_file)
  {
    // Read XML configuration file.
    read_xml(config_file, proptree);
  }

  std::string input = proptree.get("TEST.INPUT", "test.dat");

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
  gsl_histogram *h = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(h, a, b);

  size_t N = 20;
  double A = 0.2, B = 5.0;
  gsl_histogram *v2h = gsl_histogram_alloc(N);
  gsl_histogram_set_ranges_uniform(v2h, A, B);
  gsl_histogram *v2hn = gsl_histogram_alloc(N);
  gsl_histogram_set_ranges_uniform(v2hn, A, B);

  HepMC::IO_GenEvent ascii_in(istr);
  HepMC::GenEvent *evt;

  while (ascii_in >> evt)
  {
    HepMC::HeavyIon *hi = evt->heavy_ion();
    HepMC::GenVertex *primary_vtx = evt->barcode_to_vertex(-1);
    double phi0 = hi->event_plane_angle();

    HepMC::GenVertexParticleRange r(*primary_vtx, HepMC::children);
    for (HepMC::GenVertex::particle_iterator it = r.begin(); it != r.end(); it++)
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
      dphi = atan2(sin(dphi), cos(dphi));
      gsl_histogram_accumulate(v2h, pt, cos(2 * dphi));
      gsl_histogram_increment(v2hn, pt);
      gsl_histogram_increment(h, dphi);
    }

    delete evt;
    ascii_in >> evt;
  }

  std::cout << "pt,v2,v2e" << std::endl;
  gsl_histogram_div(v2h, v2hn);
  double lower, upper;
  for (size_t i = 0; i < N; i++)
  {
    gsl_histogram_get_range(v2h, i, &lower, &upper);
    double pt = (lower + upper) / 2.0;
    double counts = gsl_histogram_get(v2hn, i);
    double val = gsl_histogram_get(v2h, i);
    double err = val / sqrt(counts);
    std::cout << pt << ", " << val << ", " << err << std::endl;
  }

  return 0;
}
