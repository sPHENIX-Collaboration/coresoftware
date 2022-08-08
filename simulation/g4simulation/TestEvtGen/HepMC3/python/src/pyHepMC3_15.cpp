#include <HepMC3/LHEF.h>
#include <ios>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
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

void bind_pyHepMC3_15(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // LHEF::MergeInfo file:HepMC3/LHEF.h line:992
		pybind11::class_<LHEF::MergeInfo, std::shared_ptr<LHEF::MergeInfo>, LHEF::TagBase> cl(M("LHEF"), "MergeInfo", "The MergeInfo class represents the information in a mergeinfo tag.");
		cl.def( pybind11::init( [](){ return new LHEF::MergeInfo(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::MergeInfo const &o){ return new LHEF::MergeInfo(o); } ) );
		cl.def_readwrite("iproc", &LHEF::MergeInfo::iproc);
		cl.def_readwrite("mergingscale", &LHEF::MergeInfo::mergingscale);
		cl.def_readwrite("maxmult", &LHEF::MergeInfo::maxmult);
		cl.def("assign", (struct LHEF::MergeInfo & (LHEF::MergeInfo::*)(const struct LHEF::MergeInfo &)) &LHEF::MergeInfo::operator=, "C++: LHEF::MergeInfo::operator=(const struct LHEF::MergeInfo &) --> struct LHEF::MergeInfo &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::MergeInfo>(cl);
	}
	{ // LHEF::WeightInfo file:HepMC3/LHEF.h line:1042
		pybind11::class_<LHEF::WeightInfo, std::shared_ptr<LHEF::WeightInfo>, LHEF::TagBase> cl(M("LHEF"), "WeightInfo", "The WeightInfo class encodes the description of a given weight\n present for all events.");
		cl.def( pybind11::init( [](){ return new LHEF::WeightInfo(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::WeightInfo const &o){ return new LHEF::WeightInfo(o); } ) );
		cl.def_readwrite("inGroup", &LHEF::WeightInfo::inGroup);
		cl.def_readwrite("isrwgt", &LHEF::WeightInfo::isrwgt);
		cl.def_readwrite("name", &LHEF::WeightInfo::name);
		cl.def_readwrite("muf", &LHEF::WeightInfo::muf);
		cl.def_readwrite("mur", &LHEF::WeightInfo::mur);
		cl.def_readwrite("pdf", &LHEF::WeightInfo::pdf);
		cl.def_readwrite("pdf2", &LHEF::WeightInfo::pdf2);
		cl.def("assign", (struct LHEF::WeightInfo & (LHEF::WeightInfo::*)(const struct LHEF::WeightInfo &)) &LHEF::WeightInfo::operator=, "C++: LHEF::WeightInfo::operator=(const struct LHEF::WeightInfo &) --> struct LHEF::WeightInfo &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::WeightInfo>(cl);
	}
	{ // LHEF::WeightGroup file:HepMC3/LHEF.h line:1128
		pybind11::class_<LHEF::WeightGroup, std::shared_ptr<LHEF::WeightGroup>, LHEF::TagBase> cl(M("LHEF"), "WeightGroup", "The WeightGroup assigns a group-name to a set of WeightInfo objects.");
		cl.def( pybind11::init( [](){ return new LHEF::WeightGroup(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &, int, class std::vector<struct LHEF::WeightInfo> &>(), pybind11::arg("tag"), pybind11::arg("groupIndex"), pybind11::arg("wiv") );

		cl.def( pybind11::init( [](LHEF::WeightGroup const &o){ return new LHEF::WeightGroup(o); } ) );
		cl.def_readwrite("type", &LHEF::WeightGroup::type);
		cl.def_readwrite("combine", &LHEF::WeightGroup::combine);
		cl.def("assign", (struct LHEF::WeightGroup & (LHEF::WeightGroup::*)(const struct LHEF::WeightGroup &)) &LHEF::WeightGroup::operator=, "C++: LHEF::WeightGroup::operator=(const struct LHEF::WeightGroup &) --> struct LHEF::WeightGroup &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // LHEF::Weight file:HepMC3/LHEF.h line:1169
		pybind11::class_<LHEF::Weight, std::shared_ptr<LHEF::Weight>, LHEF::TagBase> cl(M("LHEF"), "Weight", "The Weight class represents the information in a weight tag.");
		cl.def( pybind11::init( [](){ return new LHEF::Weight(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::Weight const &o){ return new LHEF::Weight(o); } ) );
		cl.def_readwrite("name", &LHEF::Weight::name);
		cl.def_readwrite("iswgt", &LHEF::Weight::iswgt);
		cl.def_readwrite("born", &LHEF::Weight::born);
		cl.def_readwrite("sudakov", &LHEF::Weight::sudakov);
		cl.def_readwrite("weights", &LHEF::Weight::weights);
		cl.def_readwrite("indices", &LHEF::Weight::indices);
		cl.def("assign", (struct LHEF::Weight & (LHEF::Weight::*)(const struct LHEF::Weight &)) &LHEF::Weight::operator=, "C++: LHEF::Weight::operator=(const struct LHEF::Weight &) --> struct LHEF::Weight &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::Weight>(cl);
	}
	{ // LHEF::Clus file:HepMC3/LHEF.h line:1250
		pybind11::class_<LHEF::Clus, std::shared_ptr<LHEF::Clus>, LHEF::TagBase> cl(M("LHEF"), "Clus", "The Clus class represents a clustering of two particle entries into\n one as defined in a clustering tag.");
		cl.def( pybind11::init( [](){ return new LHEF::Clus(); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::Clus const &o){ return new LHEF::Clus(o); } ) );
		cl.def_readwrite("p1", &LHEF::Clus::p1);
		cl.def_readwrite("p2", &LHEF::Clus::p2);
		cl.def_readwrite("p0", &LHEF::Clus::p0);
		cl.def_readwrite("scale", &LHEF::Clus::scale);
		cl.def_readwrite("alphas", &LHEF::Clus::alphas);
		cl.def("assign", (struct LHEF::Clus & (LHEF::Clus::*)(const struct LHEF::Clus &)) &LHEF::Clus::operator=, "C++: LHEF::Clus::operator=(const struct LHEF::Clus &) --> struct LHEF::Clus &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::Clus>(cl);
	}
	{ // LHEF::Scale file:HepMC3/LHEF.h line:1313
		pybind11::class_<LHEF::Scale, std::shared_ptr<LHEF::Scale>, LHEF::TagBase> cl(M("LHEF"), "Scale", "Store special scales from within a scales tag.");
		cl.def( pybind11::init( [](){ return new LHEF::Scale(); } ), "doc" );
		cl.def( pybind11::init( [](std::string const & a0){ return new LHEF::Scale(a0); } ), "doc" , pybind11::arg("st"));
		cl.def( pybind11::init( [](std::string const & a0, int const & a1){ return new LHEF::Scale(a0, a1); } ), "doc" , pybind11::arg("st"), pybind11::arg("emtr"));
		cl.def( pybind11::init<std::string, int, double>(), pybind11::arg("st"), pybind11::arg("emtr"), pybind11::arg("sc") );

		cl.def( pybind11::init<const struct LHEF::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](LHEF::Scale const &o){ return new LHEF::Scale(o); } ) );
		cl.def_readwrite("stype", &LHEF::Scale::stype);
		cl.def_readwrite("emitter", &LHEF::Scale::emitter);
		cl.def_readwrite("recoilers", &LHEF::Scale::recoilers);
		cl.def_readwrite("emitted", &LHEF::Scale::emitted);
		cl.def_readwrite("scale", &LHEF::Scale::scale);
		cl.def("assign", (struct LHEF::Scale & (LHEF::Scale::*)(const struct LHEF::Scale &)) &LHEF::Scale::operator=, "C++: LHEF::Scale::operator=(const struct LHEF::Scale &) --> struct LHEF::Scale &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::Scale>(cl);
	}
	{ // LHEF::Scales file:HepMC3/LHEF.h line:1416
		pybind11::class_<LHEF::Scales, std::shared_ptr<LHEF::Scales>, LHEF::TagBase> cl(M("LHEF"), "Scales", "Collect different scales relevant for an event.");
		cl.def( pybind11::init( [](){ return new LHEF::Scales(); } ), "doc" );
		cl.def( pybind11::init( [](double const & a0){ return new LHEF::Scales(a0); } ), "doc" , pybind11::arg("defscale"));
		cl.def( pybind11::init<double, int>(), pybind11::arg("defscale"), pybind11::arg("npart") );

		cl.def( pybind11::init( [](const struct LHEF::XMLTag & a0){ return new LHEF::Scales(a0); } ), "doc" , pybind11::arg("tag"));
		cl.def( pybind11::init( [](const struct LHEF::XMLTag & a0, double const & a1){ return new LHEF::Scales(a0, a1); } ), "doc" , pybind11::arg("tag"), pybind11::arg("defscale"));
		cl.def( pybind11::init<const struct LHEF::XMLTag &, double, int>(), pybind11::arg("tag"), pybind11::arg("defscale"), pybind11::arg("npart") );

		cl.def( pybind11::init( [](LHEF::Scales const &o){ return new LHEF::Scales(o); } ) );
		cl.def_readwrite("muf", &LHEF::Scales::muf);
		cl.def_readwrite("mur", &LHEF::Scales::mur);
		cl.def_readwrite("mups", &LHEF::Scales::mups);
		cl.def_readwrite("SCALUP", &LHEF::Scales::SCALUP);
		cl.def_readwrite("scales", &LHEF::Scales::scales);
		cl.def("hasInfo", (bool (LHEF::Scales::*)() const) &LHEF::Scales::hasInfo, "Check if this object contains useful information besides SCALUP.\n\nC++: LHEF::Scales::hasInfo() const --> bool");
		cl.def("getScale", (double (LHEF::Scales::*)(std::string, int, int, int) const) &LHEF::Scales::getScale, "Return the scale of type st for a given emission of particle type\n pdgem from the emitter with number emr and a recoiler rec. (Note\n that the indices for emr and rec starts at 1 and 0 is interpreted\n as any particle.) First it will check for Scale object with an\n exact match. If not found, it will search for an exact match for\n the emitter and recoiler with an undefined emitted particle. If\n not found, it will look for a match for only emitter and emitted,\n of if not found, a match for only the emitter. Finally a general\n Scale object will be used, or if nothing matches, the mups will\n be returned.\n\nC++: LHEF::Scales::getScale(std::string, int, int, int) const --> double", pybind11::arg("st"), pybind11::arg("pdgem"), pybind11::arg("emr"), pybind11::arg("rec"));
		cl.def("assign", (struct LHEF::Scales & (LHEF::Scales::*)(const struct LHEF::Scales &)) &LHEF::Scales::operator=, "C++: LHEF::Scales::operator=(const struct LHEF::Scales &) --> struct LHEF::Scales &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::Scales>(cl);
	}
	{ // LHEF::PDFInfo file:HepMC3/LHEF.h line:1540
		pybind11::class_<LHEF::PDFInfo, std::shared_ptr<LHEF::PDFInfo>, LHEF::TagBase> cl(M("LHEF"), "PDFInfo", "The PDFInfo class represents the information in a pdfinto tag.");
		cl.def( pybind11::init( [](){ return new LHEF::PDFInfo(); } ), "doc" );
		cl.def( pybind11::init<double>(), pybind11::arg("defscale") );

		cl.def( pybind11::init( [](const struct LHEF::XMLTag & a0){ return new LHEF::PDFInfo(a0); } ), "doc" , pybind11::arg("tag"));
		cl.def( pybind11::init<const struct LHEF::XMLTag &, double>(), pybind11::arg("tag"), pybind11::arg("defscale") );

		cl.def( pybind11::init( [](LHEF::PDFInfo const &o){ return new LHEF::PDFInfo(o); } ) );
		cl.def_readwrite("p1", &LHEF::PDFInfo::p1);
		cl.def_readwrite("p2", &LHEF::PDFInfo::p2);
		cl.def_readwrite("x1", &LHEF::PDFInfo::x1);
		cl.def_readwrite("x2", &LHEF::PDFInfo::x2);
		cl.def_readwrite("xf1", &LHEF::PDFInfo::xf1);
		cl.def_readwrite("xf2", &LHEF::PDFInfo::xf2);
		cl.def_readwrite("scale", &LHEF::PDFInfo::scale);
		cl.def_readwrite("SCALUP", &LHEF::PDFInfo::SCALUP);
		cl.def("assign", (struct LHEF::PDFInfo & (LHEF::PDFInfo::*)(const struct LHEF::PDFInfo &)) &LHEF::PDFInfo::operator=, "C++: LHEF::PDFInfo::operator=(const struct LHEF::PDFInfo &) --> struct LHEF::PDFInfo &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		 binder::custom_T_binder<LHEF::PDFInfo>(cl);
	}
	{ // LHEF::HEPRUP file:HepMC3/LHEF.h line:1627
		pybind11::class_<LHEF::HEPRUP, std::shared_ptr<LHEF::HEPRUP>, LHEF::TagBase> cl(M("LHEF"), "HEPRUP", "The HEPRUP class is a simple container corresponding to the Les Houches\n accord (<A HREF=\"http://arxiv.org/abs/hep-ph/0109068\">hep-ph/0109068</A>)\n common block with the same name. The members are named in the same\n way as in the common block. However, fortran arrays are represented\n by vectors, except for the arrays of length two which are\n represented by pair objects.");
		cl.def( pybind11::init( [](){ return new LHEF::HEPRUP(); } ) );
		cl.def( pybind11::init( [](LHEF::HEPRUP const &o){ return new LHEF::HEPRUP(o); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &, int>(), pybind11::arg("tagin"), pybind11::arg("versin") );

		cl.def_readwrite("IDBMUP", &LHEF::HEPRUP::IDBMUP);
		cl.def_readwrite("EBMUP", &LHEF::HEPRUP::EBMUP);
		cl.def_readwrite("PDFGUP", &LHEF::HEPRUP::PDFGUP);
		cl.def_readwrite("PDFSUP", &LHEF::HEPRUP::PDFSUP);
		cl.def_readwrite("IDWTUP", &LHEF::HEPRUP::IDWTUP);
		cl.def_readwrite("NPRUP", &LHEF::HEPRUP::NPRUP);
		cl.def_readwrite("XSECUP", &LHEF::HEPRUP::XSECUP);
		cl.def_readwrite("XERRUP", &LHEF::HEPRUP::XERRUP);
		cl.def_readwrite("XMAXUP", &LHEF::HEPRUP::XMAXUP);
		cl.def_readwrite("LPRUP", &LHEF::HEPRUP::LPRUP);
		cl.def_readwrite("xsecinfos", &LHEF::HEPRUP::xsecinfos);
		cl.def_readwrite("eventfiles", &LHEF::HEPRUP::eventfiles);
		cl.def_readwrite("cuts", &LHEF::HEPRUP::cuts);
		cl.def_readwrite("ptypes", &LHEF::HEPRUP::ptypes);
		cl.def_readwrite("procinfo", &LHEF::HEPRUP::procinfo);
		cl.def_readwrite("mergeinfo", &LHEF::HEPRUP::mergeinfo);
		cl.def_readwrite("generators", &LHEF::HEPRUP::generators);
		cl.def_readwrite("weightinfo", &LHEF::HEPRUP::weightinfo);
		cl.def_readwrite("weightmap", &LHEF::HEPRUP::weightmap);
		cl.def_readwrite("weightgroup", &LHEF::HEPRUP::weightgroup);
		cl.def_readwrite("junk", &LHEF::HEPRUP::junk);
		cl.def_readwrite("version", &LHEF::HEPRUP::version);
		cl.def_readwrite("dprec", &LHEF::HEPRUP::dprec);
		cl.def("assign", (class LHEF::HEPRUP & (LHEF::HEPRUP::*)(const class LHEF::HEPRUP &)) &LHEF::HEPRUP::operator=, "Assignment operator.\n\nC++: LHEF::HEPRUP::operator=(const class LHEF::HEPRUP &) --> class LHEF::HEPRUP &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("weightNameHepMC", (std::string (LHEF::HEPRUP::*)(int) const) &LHEF::HEPRUP::weightNameHepMC, "Return the name of the weight with given index suitable to ne\n used for HepMC3 output.\n\nC++: LHEF::HEPRUP::weightNameHepMC(int) const --> std::string", pybind11::arg("i"));
		cl.def("clear", (void (LHEF::HEPRUP::*)()) &LHEF::HEPRUP::clear, "Clear all information.\n\nC++: LHEF::HEPRUP::clear() --> void");
		cl.def("resize", (void (LHEF::HEPRUP::*)(int)) &LHEF::HEPRUP::resize, "Set the NPRUP variable, corresponding to the number of\n sub-processes, to  and resize all relevant vectors\n accordingly.\n\nC++: LHEF::HEPRUP::resize(int) --> void", pybind11::arg("nrup"));
		cl.def("resize", (void (LHEF::HEPRUP::*)()) &LHEF::HEPRUP::resize, "Assuming the NPRUP variable, corresponding to the number of\n sub-processes, is correctly set, resize the relevant vectors\n accordingly.\n\nC++: LHEF::HEPRUP::resize() --> void");
		cl.def("weightIndex", (int (LHEF::HEPRUP::*)(std::string) const) &LHEF::HEPRUP::weightIndex, "the index of the weight with the given \n   \n\nC++: LHEF::HEPRUP::weightIndex(std::string) const --> int", pybind11::arg("name"));
		cl.def("nWeights", (int (LHEF::HEPRUP::*)() const) &LHEF::HEPRUP::nWeights, "the number of weights (including the nominial one).\n\nC++: LHEF::HEPRUP::nWeights() const --> int");
		cl.def("getXSecInfo", [](LHEF::HEPRUP &o) -> LHEF::XSecInfo & { return o.getXSecInfo(); }, "", pybind11::return_value_policy::automatic);
		cl.def("getXSecInfo", (struct LHEF::XSecInfo & (LHEF::HEPRUP::*)(std::string)) &LHEF::HEPRUP::getXSecInfo, "the XSecInfo object corresponding to the named weight \n If no such object exists, it will be created.\n\nC++: LHEF::HEPRUP::getXSecInfo(std::string) --> struct LHEF::XSecInfo &", pybind11::return_value_policy::automatic, pybind11::arg("weightname"));

		 binder::custom_T_binder<LHEF::HEPRUP>(cl);
	}
	{ // LHEF::EventGroup file:HepMC3/LHEF.h line:2069
		pybind11::class_<LHEF::EventGroup, std::shared_ptr<LHEF::EventGroup>, std::vector<LHEF::HEPEUP *>> cl(M("LHEF"), "EventGroup", "The EventGroup represents a set of events which are to be\n considered together.");
		cl.def( pybind11::init( [](){ return new LHEF::EventGroup(); } ) );
		cl.def( pybind11::init( [](LHEF::EventGroup const &o){ return new LHEF::EventGroup(o); } ) );
		cl.def_readwrite("nreal", &LHEF::EventGroup::nreal);
		cl.def_readwrite("ncounter", &LHEF::EventGroup::ncounter);
		cl.def("assign", (struct LHEF::EventGroup & (LHEF::EventGroup::*)(const struct LHEF::EventGroup &)) &LHEF::EventGroup::operator=, "The assignment also copies the included HEPEUP object.\n\nC++: LHEF::EventGroup::operator=(const struct LHEF::EventGroup &) --> struct LHEF::EventGroup &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("clear", (void (LHEF::EventGroup::*)()) &LHEF::EventGroup::clear, "Remove all subevents.\n\nC++: LHEF::EventGroup::clear() --> void");
	}
	{ // LHEF::HEPEUP file:HepMC3/LHEF.h line:2117
		pybind11::class_<LHEF::HEPEUP, std::shared_ptr<LHEF::HEPEUP>, LHEF::TagBase> cl(M("LHEF"), "HEPEUP", "The HEPEUP class is a simple container corresponding to the Les Houches accord\n (<A HREF=\"http://arxiv.org/abs/hep-ph/0109068\">hep-ph/0109068</A>)\n common block with the same name. The members are named in the same\n way as in the common block. However, fortran arrays are represented\n by vectors, except for the arrays of length two which are\n represented by pair objects.");
		cl.def( pybind11::init( [](){ return new LHEF::HEPEUP(); } ) );
		cl.def( pybind11::init( [](LHEF::HEPEUP const &o){ return new LHEF::HEPEUP(o); } ) );
		cl.def( pybind11::init<const struct LHEF::XMLTag &, class LHEF::HEPRUP &>(), pybind11::arg("tagin"), pybind11::arg("heprupin") );

		cl.def_readwrite("NUP", &LHEF::HEPEUP::NUP);
		cl.def_readwrite("IDPRUP", &LHEF::HEPEUP::IDPRUP);
		cl.def_readwrite("XWGTUP", &LHEF::HEPEUP::XWGTUP);
		cl.def_readwrite("XPDWUP", &LHEF::HEPEUP::XPDWUP);
		cl.def_readwrite("SCALUP", &LHEF::HEPEUP::SCALUP);
		cl.def_readwrite("AQEDUP", &LHEF::HEPEUP::AQEDUP);
		cl.def_readwrite("AQCDUP", &LHEF::HEPEUP::AQCDUP);
		cl.def_readwrite("IDUP", &LHEF::HEPEUP::IDUP);
		cl.def_readwrite("ISTUP", &LHEF::HEPEUP::ISTUP);
		cl.def_readwrite("MOTHUP", &LHEF::HEPEUP::MOTHUP);
		cl.def_readwrite("ICOLUP", &LHEF::HEPEUP::ICOLUP);
		cl.def_readwrite("PUP", &LHEF::HEPEUP::PUP);
		cl.def_readwrite("VTIMUP", &LHEF::HEPEUP::VTIMUP);
		cl.def_readwrite("SPINUP", &LHEF::HEPEUP::SPINUP);
		cl.def_readwrite("namedweights", &LHEF::HEPEUP::namedweights);
		cl.def_readwrite("weights", &LHEF::HEPEUP::weights);
		cl.def_readwrite("clustering", &LHEF::HEPEUP::clustering);
		cl.def_readwrite("pdfinfo", &LHEF::HEPEUP::pdfinfo);
		cl.def_readwrite("PDFGUPsave", &LHEF::HEPEUP::PDFGUPsave);
		cl.def_readwrite("PDFSUPsave", &LHEF::HEPEUP::PDFSUPsave);
		cl.def_readwrite("scales", &LHEF::HEPEUP::scales);
		cl.def_readwrite("ntries", &LHEF::HEPEUP::ntries);
		cl.def_readwrite("isGroup", &LHEF::HEPEUP::isGroup);
		cl.def_readwrite("subevents", &LHEF::HEPEUP::subevents);
		cl.def_readwrite("junk", &LHEF::HEPEUP::junk);
		cl.def("setEvent", (class LHEF::HEPEUP & (LHEF::HEPEUP::*)(const class LHEF::HEPEUP &)) &LHEF::HEPEUP::setEvent, "Copy information from the given HEPEUP. Sub event information is\n left untouched.\n\nC++: LHEF::HEPEUP::setEvent(const class LHEF::HEPEUP &) --> class LHEF::HEPEUP &", pybind11::return_value_policy::automatic, pybind11::arg("x"));
		cl.def("assign", (class LHEF::HEPEUP & (LHEF::HEPEUP::*)(const class LHEF::HEPEUP &)) &LHEF::HEPEUP::operator=, "Assignment operator.\n\nC++: LHEF::HEPEUP::operator=(const class LHEF::HEPEUP &) --> class LHEF::HEPEUP &", pybind11::return_value_policy::automatic, pybind11::arg("x"));
		cl.def("reset", (void (LHEF::HEPEUP::*)()) &LHEF::HEPEUP::reset, "Reset the HEPEUP object (does not touch the sub events).\n\nC++: LHEF::HEPEUP::reset() --> void");
		cl.def("clear", (void (LHEF::HEPEUP::*)()) &LHEF::HEPEUP::clear, "Clear the HEPEUP object.\n\nC++: LHEF::HEPEUP::clear() --> void");
		cl.def("resize", (void (LHEF::HEPEUP::*)(int)) &LHEF::HEPEUP::resize, "Set the NUP variable, corresponding to the number of particles in\n the current event, to  and resize all relevant vectors\n accordingly.\n\nC++: LHEF::HEPEUP::resize(int) --> void", pybind11::arg("nup"));
		cl.def("totalWeight", [](LHEF::HEPEUP const &o) -> double { return o.totalWeight(); }, "");
		cl.def("totalWeight", (double (LHEF::HEPEUP::*)(int) const) &LHEF::HEPEUP::totalWeight, "Return the total weight for this event (including all sub\n evenets) for the given index.\n\nC++: LHEF::HEPEUP::totalWeight(int) const --> double", pybind11::arg("i"));
		cl.def("totalWeight", (double (LHEF::HEPEUP::*)(std::string) const) &LHEF::HEPEUP::totalWeight, "Return the total weight for this event (including all sub\n evenets) for the given weight name.\n\nC++: LHEF::HEPEUP::totalWeight(std::string) const --> double", pybind11::arg("name"));
		cl.def("weight", [](LHEF::HEPEUP const &o) -> double { return o.weight(); }, "");
		cl.def("weight", (double (LHEF::HEPEUP::*)(int) const) &LHEF::HEPEUP::weight, "Return the weight for the given index.\n\nC++: LHEF::HEPEUP::weight(int) const --> double", pybind11::arg("i"));
		cl.def("weight", (double (LHEF::HEPEUP::*)(std::string) const) &LHEF::HEPEUP::weight, "Return the weight for the given weight name.\n\nC++: LHEF::HEPEUP::weight(std::string) const --> double", pybind11::arg("name"));
		cl.def("setWeight", (void (LHEF::HEPEUP::*)(int, double)) &LHEF::HEPEUP::setWeight, "Set the weight with the given index.\n\nC++: LHEF::HEPEUP::setWeight(int, double) --> void", pybind11::arg("i"), pybind11::arg("w"));
		cl.def("setWeight", (bool (LHEF::HEPEUP::*)(std::string, double)) &LHEF::HEPEUP::setWeight, "Set the weight with the given name.\n\nC++: LHEF::HEPEUP::setWeight(std::string, double) --> bool", pybind11::arg("name"), pybind11::arg("w"));
		cl.def("resize", (void (LHEF::HEPEUP::*)()) &LHEF::HEPEUP::resize, "Assuming the NUP variable, corresponding to the number of\n particles in the current event, is correctly set, resize the\n relevant vectors accordingly.\n\nC++: LHEF::HEPEUP::resize() --> void");
		cl.def("setWeightInfo", (bool (LHEF::HEPEUP::*)(unsigned int)) &LHEF::HEPEUP::setWeightInfo, "Setup the current event to use weight i. If zero, the default\n weight will be used.\n\nC++: LHEF::HEPEUP::setWeightInfo(unsigned int) --> bool", pybind11::arg("i"));
		cl.def("setSubEvent", (bool (LHEF::HEPEUP::*)(unsigned int)) &LHEF::HEPEUP::setSubEvent, "Setup the current event to use sub event i. If zero, no sub event\n will be chsen.\n\nC++: LHEF::HEPEUP::setSubEvent(unsigned int) --> bool", pybind11::arg("i"));

		 binder::custom_T_binder<LHEF::HEPEUP>(cl);
	}
}
