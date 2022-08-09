#include <HepMC3/LHEF.h>
#include <functional>
#include <ios>
#include <iterator>
#include <locale>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <functional>
#include <string>
#include <HepMC3/Version.h>
#include <HepMC3/Reader.h>
#include <HepMC3/Writer.h>
#include <HepMC3/Print.h>
#include <src/stl_binders.hpp>
#include <src/binders.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

void bind_pyHepMC3_14(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // LHEF::OAttr file:HepMC3/LHEF.h line:45
		pybind11::class_<LHEF::OAttr<std::string>, std::shared_ptr<LHEF::OAttr<std::string>>> cl(M("LHEF"), "OAttr_std_string_t", "");
		cl.def( pybind11::init<std::string, const std::string &>(), pybind11::arg("n"), pybind11::arg("v") );

		cl.def( pybind11::init( [](LHEF::OAttr<std::string> const &o){ return new LHEF::OAttr<std::string>(o); } ) );
		cl.def_readwrite("name", &LHEF::OAttr<std::string>::name);
		cl.def_readwrite("val", &LHEF::OAttr<std::string>::val);
		cl.def("assign", (struct LHEF::OAttr<std::string > & (LHEF::OAttr<std::string>::*)(const struct LHEF::OAttr<std::string > &)) &LHEF::OAttr<std::string>::operator=, "C++: LHEF::OAttr<std::string>::operator=(const struct LHEF::OAttr<std::string > &) --> struct LHEF::OAttr<std::string > &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](LHEF::OAttr<std::string> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // LHEF::OAttr file:HepMC3/LHEF.h line:45
		pybind11::class_<LHEF::OAttr<long>, std::shared_ptr<LHEF::OAttr<long>>> cl(M("LHEF"), "OAttr_long_t", "");
		cl.def( pybind11::init<std::string, const long &>(), pybind11::arg("n"), pybind11::arg("v") );

		cl.def( pybind11::init( [](LHEF::OAttr<long> const &o){ return new LHEF::OAttr<long>(o); } ) );
		cl.def_readwrite("name", &LHEF::OAttr<long>::name);
		cl.def_readwrite("val", &LHEF::OAttr<long>::val);
		cl.def("assign", (struct LHEF::OAttr<long> & (LHEF::OAttr<long>::*)(const struct LHEF::OAttr<long> &)) &LHEF::OAttr<long>::operator=, "C++: LHEF::OAttr<long>::operator=(const struct LHEF::OAttr<long> &) --> struct LHEF::OAttr<long> &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](LHEF::OAttr<long> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // LHEF::OAttr file:HepMC3/LHEF.h line:45
		pybind11::class_<LHEF::OAttr<double>, std::shared_ptr<LHEF::OAttr<double>>> cl(M("LHEF"), "OAttr_double_t", "");
		cl.def( pybind11::init<std::string, const double &>(), pybind11::arg("n"), pybind11::arg("v") );

		cl.def( pybind11::init( [](LHEF::OAttr<double> const &o){ return new LHEF::OAttr<double>(o); } ) );
		cl.def_readwrite("name", &LHEF::OAttr<double>::name);
		cl.def_readwrite("val", &LHEF::OAttr<double>::val);
		cl.def("assign", (struct LHEF::OAttr<double> & (LHEF::OAttr<double>::*)(const struct LHEF::OAttr<double> &)) &LHEF::OAttr<double>::operator=, "C++: LHEF::OAttr<double>::operator=(const struct LHEF::OAttr<double> &) --> struct LHEF::OAttr<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](LHEF::OAttr<double> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // LHEF::OAttr file:HepMC3/LHEF.h line:45
		pybind11::class_<LHEF::OAttr<int>, std::shared_ptr<LHEF::OAttr<int>>> cl(M("LHEF"), "OAttr_int_t", "");
		cl.def( pybind11::init<std::string, const int &>(), pybind11::arg("n"), pybind11::arg("v") );

		cl.def( pybind11::init( [](LHEF::OAttr<int> const &o){ return new LHEF::OAttr<int>(o); } ) );
		cl.def_readwrite("name", &LHEF::OAttr<int>::name);
		cl.def_readwrite("val", &LHEF::OAttr<int>::val);
		cl.def("assign", (struct LHEF::OAttr<int> & (LHEF::OAttr<int>::*)(const struct LHEF::OAttr<int> &)) &LHEF::OAttr<int>::operator=, "C++: LHEF::OAttr<int>::operator=(const struct LHEF::OAttr<int> &) --> struct LHEF::OAttr<int> &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](LHEF::OAttr<int> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	// LHEF::oattr(std::string, const std::string &) file:HepMC3/LHEF.h line:68
	M("LHEF").def("oattr", (struct LHEF::OAttr<std::string > (*)(std::string, const std::string &)) &LHEF::oattr<std::string>, "C++: LHEF::oattr(std::string, const std::string &) --> struct LHEF::OAttr<std::string >", pybind11::arg("name"), pybind11::arg("value"));

	// LHEF::oattr(std::string, const long &) file:HepMC3/LHEF.h line:68
	M("LHEF").def("oattr", (struct LHEF::OAttr<long> (*)(std::string, const long &)) &LHEF::oattr<long>, "C++: LHEF::oattr(std::string, const long &) --> struct LHEF::OAttr<long>", pybind11::arg("name"), pybind11::arg("value"));

	// LHEF::oattr(std::string, const double &) file:HepMC3/LHEF.h line:68
	M("LHEF").def("oattr", (struct LHEF::OAttr<double> (*)(std::string, const double &)) &LHEF::oattr<double>, "C++: LHEF::oattr(std::string, const double &) --> struct LHEF::OAttr<double>", pybind11::arg("name"), pybind11::arg("value"));

	// LHEF::oattr(std::string, const int &) file:HepMC3/LHEF.h line:68
	M("LHEF").def("oattr", (struct LHEF::OAttr<int> (*)(std::string, const int &)) &LHEF::oattr<int>, "C++: LHEF::oattr(std::string, const int &) --> struct LHEF::OAttr<int>", pybind11::arg("name"), pybind11::arg("value"));

	{ // LHEF::XMLTag file:HepMC3/LHEF.h line:87
		pybind11::class_<LHEF::XMLTag, std::shared_ptr<LHEF::XMLTag>> cl(M("LHEF"), "XMLTag", "The XMLTag struct is used to represent all information within an\n XML tag. It contains the attributes as a map, any sub-tags as a\n vector of pointers to other XMLTag objects, and any other\n information as a single string.");
		cl.def( pybind11::init( [](){ return new LHEF::XMLTag(); } ) );
		cl.def( pybind11::init( [](LHEF::XMLTag const &o){ return new LHEF::XMLTag(o); } ) );
		cl.def_readwrite("name", &LHEF::XMLTag::name);
		cl.def_readwrite("attr", &LHEF::XMLTag::attr);
		cl.def_readwrite("tags", &LHEF::XMLTag::tags);
		cl.def_readwrite("contents", &LHEF::XMLTag::contents);
		cl.def("getattr", (bool (LHEF::XMLTag::*)(std::string, double &) const) &LHEF::XMLTag::getattr, "Find an attribute named  and set the double variable  to\n the corresponding value. \n\n false if no attribute was found.\n\nC++: LHEF::XMLTag::getattr(std::string, double &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::XMLTag::*)(std::string, bool &) const) &LHEF::XMLTag::getattr, "Find an attribute named  and set the bool variable  to\n true if the corresponding value is \"yes\". \n\n false if no\n attribute was found.\n\nC++: LHEF::XMLTag::getattr(std::string, bool &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::XMLTag::*)(std::string, long &) const) &LHEF::XMLTag::getattr, "Find an attribute named  and set the long variable  to\n the corresponding value. \n\n false if no attribute was found.\n\nC++: LHEF::XMLTag::getattr(std::string, long &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::XMLTag::*)(std::string, int &) const) &LHEF::XMLTag::getattr, "Find an attribute named  and set the long variable  to\n the corresponding value. \n\n false if no attribute was found.\n\nC++: LHEF::XMLTag::getattr(std::string, int &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::XMLTag::*)(std::string, std::string &) const) &LHEF::XMLTag::getattr, "Find an attribute named  and set the string variable  to\n the corresponding value. \n\n false if no attribute was found.\n\nC++: LHEF::XMLTag::getattr(std::string, std::string &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def_static("findXMLTags", [](std::string const & a0) -> std::vector<struct LHEF::XMLTag *> { return LHEF::XMLTag::findXMLTags(a0); }, "", pybind11::arg("str"));
		cl.def_static("findXMLTags", (class std::vector<struct LHEF::XMLTag *> (*)(std::string, std::string *)) &LHEF::XMLTag::findXMLTags, "Scan the given string and return all XML tags found as a vector\n of pointers to XMLTag objects. Text which does not belong to any\n tag is stored in tags without name and in the string pointed to\n by leftover (if not null).\n\nC++: LHEF::XMLTag::findXMLTags(std::string, std::string *) --> class std::vector<struct LHEF::XMLTag *>", pybind11::arg("str"), pybind11::arg("leftover"));
		cl.def_static("deleteAll", (void (*)(class std::vector<struct LHEF::XMLTag *> &)) &LHEF::XMLTag::deleteAll, "Delete all tags in a vector.\n\nC++: LHEF::XMLTag::deleteAll(class std::vector<struct LHEF::XMLTag *> &) --> void", pybind11::arg("tags"));
		cl.def("assign", (struct LHEF::XMLTag & (LHEF::XMLTag::*)(const struct LHEF::XMLTag &)) &LHEF::XMLTag::operator=, "C++: LHEF::XMLTag::operator=(const struct LHEF::XMLTag &) --> struct LHEF::XMLTag &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::XMLTag>(cl);
	}
	// LHEF::hashline(std::string) file:HepMC3/LHEF.h line:328
	M("LHEF").def("hashline", (std::string (*)(std::string)) &LHEF::hashline, "Helper function to make sure that each line in the string  starts with a\n #-character and that the string ends with a new-line.\n\nC++: LHEF::hashline(std::string) --> std::string", pybind11::arg("s"));

	{ // LHEF::TagBase file:HepMC3/LHEF.h line:345
		pybind11::class_<LHEF::TagBase, std::shared_ptr<LHEF::TagBase>> cl(M("LHEF"), "TagBase", "This is the base class of all classes representing xml tags.");
		cl.def( pybind11::init( [](){ return new LHEF::TagBase(); } ) );
		cl.def( pybind11::init( [](const class std::map<std::string, std::string > & a0){ return new LHEF::TagBase(a0); } ), "doc" , pybind11::arg("attr"));
		cl.def( pybind11::init<const class std::map<std::string, std::string > &, std::string>(), pybind11::arg("attr"), pybind11::arg("conts") );

		cl.def( pybind11::init( [](LHEF::TagBase const &o){ return new LHEF::TagBase(o); } ) );
		cl.def_readwrite("attributes", &LHEF::TagBase::attributes);
		cl.def_readwrite("contents", &LHEF::TagBase::contents);
		cl.def("getattr", [](LHEF::TagBase &o, std::string const & a0, double & a1) -> bool { return o.getattr(a0, a1); }, "", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::TagBase::*)(std::string, double &, bool)) &LHEF::TagBase::getattr, "Find an attribute named  and set the double variable  to\n the corresponding value. Remove the correspondig attribute from\n the list if found and  is true. \n\n false if no\n attribute was found.\n\nC++: LHEF::TagBase::getattr(std::string, double &, bool) --> bool", pybind11::arg("n"), pybind11::arg("v"), pybind11::arg("erase"));
		cl.def("getattr", [](LHEF::TagBase &o, std::string const & a0, bool & a1) -> bool { return o.getattr(a0, a1); }, "", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::TagBase::*)(std::string, bool &, bool)) &LHEF::TagBase::getattr, "Find an attribute named  and set the bool variable  to\n true if the corresponding value is \"yes\". Remove the correspondig\n attribute from the list if found and  is true. \n\n\n false if no attribute was found.\n\nC++: LHEF::TagBase::getattr(std::string, bool &, bool) --> bool", pybind11::arg("n"), pybind11::arg("v"), pybind11::arg("erase"));
		cl.def("getattr", [](LHEF::TagBase &o, std::string const & a0, long & a1) -> bool { return o.getattr(a0, a1); }, "", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::TagBase::*)(std::string, long &, bool)) &LHEF::TagBase::getattr, "Find an attribute named  and set the long variable  to\n the corresponding value. Remove the correspondig attribute from\n the list if found and  is true. \n\n false if no\n attribute was found.\n\nC++: LHEF::TagBase::getattr(std::string, long &, bool) --> bool", pybind11::arg("n"), pybind11::arg("v"), pybind11::arg("erase"));
		cl.def("getattr", [](LHEF::TagBase &o, std::string const & a0, int & a1) -> bool { return o.getattr(a0, a1); }, "", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::TagBase::*)(std::string, int &, bool)) &LHEF::TagBase::getattr, "Find an attribute named  and set the long variable  to\n the corresponding value. Remove the correspondig attribute from\n the list if found and  is true. \n\n false if no\n attribute was found.\n\nC++: LHEF::TagBase::getattr(std::string, int &, bool) --> bool", pybind11::arg("n"), pybind11::arg("v"), pybind11::arg("erase"));
		cl.def("getattr", [](LHEF::TagBase &o, std::string const & a0, std::string & a1) -> bool { return o.getattr(a0, a1); }, "", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (LHEF::TagBase::*)(std::string, std::string &, bool)) &LHEF::TagBase::getattr, "Find an attribute named  and set the string variable  to\n the corresponding value. Remove the correspondig attribute from\n the list if found and  is true. \n\n false if no\n attribute was found.\n\nC++: LHEF::TagBase::getattr(std::string, std::string &, bool) --> bool", pybind11::arg("n"), pybind11::arg("v"), pybind11::arg("erase"));
		cl.def_static("yes", (std::string (*)()) &LHEF::TagBase::yes, "Static string token for truth values.\n\nC++: LHEF::TagBase::yes() --> std::string");
		cl.def("assign", (struct LHEF::TagBase & (LHEF::TagBase::*)(const struct LHEF::TagBase &)) &LHEF::TagBase::operator=, "C++: LHEF::TagBase::operator=(const struct LHEF::TagBase &) --> struct LHEF::TagBase &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_LHEFTagBase_binder(cl);
	}
	{ // LHEF::Generator file:HepMC3/LHEF.h line:474
		pybind11::class_<LHEF::Generator, std::shared_ptr<LHEF::Generator>, LHEF::TagBase> cl(M("LHEF"), "Generator", "The Generator class contains information about a generator used in a run.");
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::Generator const &o){ return new LHEF::Generator(o); } ) );
		cl.def_readwrite("name", &LHEF::Generator::name);
		cl.def_readwrite("version", &LHEF::Generator::version);
		cl.def("assign", (struct LHEF::Generator & (LHEF::Generator::*)(const struct LHEF::Generator &)) &LHEF::Generator::operator=, "C++: LHEF::Generator::operator=(const struct LHEF::Generator &) --> struct LHEF::Generator &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		  binder::custom_T_binder<LHEF::Generator>(cl);
	}
	{ // LHEF::XSecInfo file:HepMC3/LHEF.h line:511
		pybind11::class_<LHEF::XSecInfo, std::shared_ptr<LHEF::XSecInfo>, LHEF::TagBase> cl(M("LHEF"), "XSecInfo", "The XSecInfo class contains information given in the xsecinfo tag.");
		cl.def( pybind11::init( [](){ return new LHEF::XSecInfo(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::XSecInfo const &o){ return new LHEF::XSecInfo(o); } ) );
		cl.def_readwrite("neve", &LHEF::XSecInfo::neve);
		cl.def_readwrite("ntries", &LHEF::XSecInfo::ntries);
		cl.def_readwrite("totxsec", &LHEF::XSecInfo::totxsec);
		cl.def_readwrite("xsecerr", &LHEF::XSecInfo::xsecerr);
		cl.def_readwrite("maxweight", &LHEF::XSecInfo::maxweight);
		cl.def_readwrite("meanweight", &LHEF::XSecInfo::meanweight);
		cl.def_readwrite("negweights", &LHEF::XSecInfo::negweights);
		cl.def_readwrite("varweights", &LHEF::XSecInfo::varweights);
		cl.def_readwrite("weightname", &LHEF::XSecInfo::weightname);
		cl.def("assign", (struct LHEF::XSecInfo & (LHEF::XSecInfo::*)(const struct LHEF::XSecInfo &)) &LHEF::XSecInfo::operator=, "C++: LHEF::XSecInfo::operator=(const struct LHEF::XSecInfo &) --> struct LHEF::XSecInfo &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::XSecInfo>(cl);
	}
	{ // LHEF::EventFile file:HepMC3/LHEF.h line:617
		pybind11::class_<LHEF::EventFile, std::shared_ptr<LHEF::EventFile>, LHEF::TagBase> cl(M("LHEF"), "EventFile", "Simple struct to store information about separate eventfiles to be\n loaded.");
		cl.def( pybind11::init( [](){ return new LHEF::EventFile(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::EventFile const &o){ return new LHEF::EventFile(o); } ) );
		cl.def_readwrite("filename", &LHEF::EventFile::filename);
		cl.def_readwrite("neve", &LHEF::EventFile::neve);
		cl.def_readwrite("ntries", &LHEF::EventFile::ntries);
		cl.def("assign", (struct LHEF::EventFile & (LHEF::EventFile::*)(const struct LHEF::EventFile &)) &LHEF::EventFile::operator=, "C++: LHEF::EventFile::operator=(const struct LHEF::EventFile &) --> struct LHEF::EventFile &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::EventFile>(cl);
	}
	{ // LHEF::Cut file:HepMC3/LHEF.h line:669
		pybind11::class_<LHEF::Cut, std::shared_ptr<LHEF::Cut>, LHEF::TagBase> cl(M("LHEF"), "Cut", "The Cut class represents a cut used by the Matrix Element generator.");
		cl.def( pybind11::init( [](){ return new LHEF::Cut(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &, const class std::map<std::string, class std::set<long> > &>(), pybind11::arg("tag"), pybind11::arg("ptypes") );

		cl.def( pybind11::init( [](LHEF::Cut const &o){ return new LHEF::Cut(o); } ) );
		cl.def_readwrite("type", &LHEF::Cut::type);
		cl.def_readwrite("p1", &LHEF::Cut::p1);
		cl.def_readwrite("np1", &LHEF::Cut::np1);
		cl.def_readwrite("p2", &LHEF::Cut::p2);
		cl.def_readwrite("np2", &LHEF::Cut::np2);
		cl.def_readwrite("min", &LHEF::Cut::min);
		cl.def_readwrite("max", &LHEF::Cut::max);
		cl.def("match", [](LHEF::Cut const &o, long const & a0) -> bool { return o.match(a0); }, "", pybind11::arg("id1"));
		cl.def("match", (bool (LHEF::Cut::*)(long, long) const) &LHEF::Cut::match, "Check if a  matches p1 and  matches p2. Only non-zero\n values are considered.\n\nC++: LHEF::Cut::match(long, long) const --> bool", pybind11::arg("id1"), pybind11::arg("id2"));
		cl.def("passCuts", (bool (LHEF::Cut::*)(const class std::vector<long> &, const class std::vector<class std::vector<double> > &) const) &LHEF::Cut::passCuts, "Check if the particles given as a vector of PDG  numbers,\n and a vector of vectors of momentum components,  will pass\n the cut defined in this event.\n\nC++: LHEF::Cut::passCuts(const class std::vector<long> &, const class std::vector<class std::vector<double> > &) const --> bool", pybind11::arg("id"), pybind11::arg("p"));
		cl.def_static("eta", (double (*)(const class std::vector<double> &)) &LHEF::Cut::eta, "Return the pseudorapidity of a particle with momentum \n   \n\nC++: LHEF::Cut::eta(const class std::vector<double> &) --> double", pybind11::arg("p"));
		cl.def_static("rap", (double (*)(const class std::vector<double> &)) &LHEF::Cut::rap, "Return the true rapidity of a particle with momentum \n   \n\nC++: LHEF::Cut::rap(const class std::vector<double> &) --> double", pybind11::arg("p"));
		cl.def_static("deltaR", (double (*)(const class std::vector<double> &, const class std::vector<double> &)) &LHEF::Cut::deltaR, "Return the delta-R of a particle pair with momenta  and \n   \n\nC++: LHEF::Cut::deltaR(const class std::vector<double> &, const class std::vector<double> &) --> double", pybind11::arg("p1"), pybind11::arg("p2"));
		cl.def("outside", (bool (LHEF::Cut::*)(double) const) &LHEF::Cut::outside, "Return true if the given  is outside limits.\n\nC++: LHEF::Cut::outside(double) const --> bool", pybind11::arg("value"));
		cl.def("assign", (struct LHEF::Cut & (LHEF::Cut::*)(const struct LHEF::Cut &)) &LHEF::Cut::operator=, "C++: LHEF::Cut::operator=(const struct LHEF::Cut &) --> struct LHEF::Cut &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::Cut>(cl);
	}
	{ // LHEF::ProcInfo file:HepMC3/LHEF.h line:915
		pybind11::class_<LHEF::ProcInfo, std::shared_ptr<LHEF::ProcInfo>, LHEF::TagBase> cl(M("LHEF"), "ProcInfo", "The ProcInfo class represents the information in a procinfo tag.");
		cl.def( pybind11::init( [](){ return new LHEF::ProcInfo(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::ProcInfo const &o){ return new LHEF::ProcInfo(o); } ) );
		cl.def_readwrite("iproc", &LHEF::ProcInfo::iproc);
		cl.def_readwrite("loops", &LHEF::ProcInfo::loops);
		cl.def_readwrite("qcdorder", &LHEF::ProcInfo::qcdorder);
		cl.def_readwrite("eworder", &LHEF::ProcInfo::eworder);
		cl.def_readwrite("fscheme", &LHEF::ProcInfo::fscheme);
		cl.def_readwrite("rscheme", &LHEF::ProcInfo::rscheme);
		cl.def_readwrite("scheme", &LHEF::ProcInfo::scheme);
		cl.def("assign", (struct LHEF::ProcInfo & (LHEF::ProcInfo::*)(const struct LHEF::ProcInfo &)) &LHEF::ProcInfo::operator=, "C++: LHEF::ProcInfo::operator=(const struct LHEF::ProcInfo &) --> struct LHEF::ProcInfo &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::ProcInfo>(cl);
	}
}
