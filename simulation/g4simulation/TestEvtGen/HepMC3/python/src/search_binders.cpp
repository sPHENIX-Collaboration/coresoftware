#include "search_binders.h"
namespace binder {
		void custom_Selector_binder( pybind11::class_<HepMC3::StandardSelector, std::shared_ptr<HepMC3::StandardSelector>> cl)
		{
		cl.def_static("STATUS", (class HepMC3::SelectorWrapper<int> (*)( )) &HepMC3::StandardSelector::STATUS, "C++: HepMC3::StandardSelector::STATUS()");
		cl.def_static("PDG_ID", (class HepMC3::SelectorWrapper<int> (*)( )) &HepMC3::StandardSelector::PDG_ID, "C++: HepMC3::StandardSelector::PDG_ID()");
		cl.def_static("PT", (class HepMC3::SelectorWrapper<double> (*)( )) &HepMC3::StandardSelector::PT, "C++: HepMC3::StandardSelector::PT()");
		cl.def_static("ENERGY", (class HepMC3::SelectorWrapper<double> (*)( )) &HepMC3::StandardSelector::ENERGY, "C++: HepMC3::StandardSelector::ENERGY()");
		cl.def_static("RAPIDITY", (class HepMC3::SelectorWrapper<double> (*)( )) &HepMC3::StandardSelector::RAPIDITY, "C++: HepMC3::StandardSelector::RAPIDITY()");
        cl.def_static("ETA", (class HepMC3::SelectorWrapper<double> (*)( )) &HepMC3::StandardSelector::ETA, "C++: HepMC3::StandardSelector::ETA()");
        cl.def_static("PHI", (class HepMC3::SelectorWrapper<double> (*)( )) &HepMC3::StandardSelector::PHI, "C++: HepMC3::StandardSelector::PHI()");
        cl.def_static("ET", (class HepMC3::SelectorWrapper<double> (*)( )) &HepMC3::StandardSelector::ET, "C++: HepMC3::StandardSelector::ET()");
        cl.def_static("MASS", (class HepMC3::SelectorWrapper<double> (*)( )) &HepMC3::StandardSelector::MASS, "C++: HepMC3::StandardSelector::MASS()");
	}

	void	search_binder(pybind11::module &M)
	{
	M.def("applyFilter", (class std::vector<class std::shared_ptr<class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<class HepMC3::GenParticle> > > (*)(const class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> &, const class std::vector<class std::shared_ptr<class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<class HepMC3::GenParticle> > > &)) &HepMC3::applyFilter, "Apply a Filter to a list of GenParticles\n Returns a vector of GenParticles that satisfy the Filter\n\nC++: HepMC3::applyFilter(const class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> &, const class std::vector<class std::shared_ptr<class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<class HepMC3::GenParticle> > > &) --> class std::vector<class std::shared_ptr<class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<class HepMC3::GenParticle> > >", pybind11::arg("filter"), pybind11::arg("particles"));

	// HepMC3::applyFilter(const class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> &, const class std::vector<class std::shared_ptr<const class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<const class HepMC3::GenParticle> > > &) file: line:31
	M.def("applyFilter", (class std::vector<class std::shared_ptr<const class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<const class HepMC3::GenParticle> > > (*)(const class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> &, const class std::vector<class std::shared_ptr<const class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<const class HepMC3::GenParticle> > > &)) &HepMC3::applyFilter, "Apply a Filter to a list of ConstGenParticles\n Returns a vector of ConstGenParticles that satisfy the Filter\n\nC++: HepMC3::applyFilter(const class std::function<bool (class std::shared_ptr<const class HepMC3::GenParticle>)> &, const class std::vector<class std::shared_ptr<const class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<const class HepMC3::GenParticle> > > &) --> class std::vector<class std::shared_ptr<const class HepMC3::GenParticle>, class std::allocator<class std::shared_ptr<const class HepMC3::GenParticle> > >", pybind11::arg("filter"), pybind11::arg("particles"));

	// HepMC3::ACCEPT_ALL(class std::shared_ptr<const class HepMC3::GenParticle>) file: line:41
	M.def("ACCEPT_ALL", (bool (*)(class std::shared_ptr<const class HepMC3::GenParticle>)) &HepMC3::ACCEPT_ALL, "A Filter that will accept all particles\n This might be needed if a signature requires a default Filter\n\nC++: HepMC3::ACCEPT_ALL(class std::shared_ptr<const class HepMC3::GenParticle>) --> bool", pybind11::arg("dummy"));


  M.def("children_particles",   (std::vector<HepMC3::GenParticlePtr>(*)(HepMC3::GenVertexPtr)      ) &HepMC3::children_particles,  "See documentation");
  M.def("children_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*)(HepMC3::ConstGenVertexPtr) ) &HepMC3::children_particles,  "See documentation");
  M.def("children_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)(HepMC3::GenParticlePtr)        ) &HepMC3::children_vertices,  "See documentation");
  M.def("children_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*)(HepMC3::ConstGenParticlePtr)   ) &HepMC3::children_vertices,  "See documentation");

  M.def("grandchildren_particles",   (std::vector<HepMC3::GenParticlePtr>(*)(HepMC3::GenParticlePtr)       ) &HepMC3::grandchildren_particles,  "See documentation");
  M.def("grandchildren_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*) (HepMC3::ConstGenParticlePtr)) &HepMC3::grandchildren_particles,  "See documentation");
  M.def("grandchildren_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)(HepMC3::GenVertexPtr)          ) &HepMC3::grandchildren_vertices,  "See documentation");
  M.def("grandchildren_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*)(HepMC3::ConstGenVertexPtr)    ) &HepMC3::grandchildren_vertices,  "See documentation");

  M.def("parent_particles",   (std::vector<HepMC3::GenParticlePtr>(*)     (HepMC3::GenVertexPtr)       ) &HepMC3::parent_particles,  "See documentation");
  M.def("parent_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*)(HepMC3::ConstGenVertexPtr)  ) &HepMC3::parent_particles,  "See documentation");
  M.def("parent_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)        (HepMC3::GenParticlePtr)         ) &HepMC3::parent_vertices,  "See documentation");
  M.def("parent_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*)   (HepMC3::ConstGenParticlePtr)    ) &HepMC3::parent_vertices,  "See documentation");

  M.def("grandparent_particles",   (std::vector<HepMC3::GenParticlePtr>(*)     (HepMC3::GenParticlePtr)       ) &HepMC3::grandparent_particles,  "See documentation");
  M.def("grandparent_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*)(HepMC3::ConstGenParticlePtr)  ) &HepMC3::grandparent_particles,  "See documentation");
  M.def("grandparent_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)        (HepMC3::GenVertexPtr)         ) &HepMC3::grandparent_vertices,  "See documentation");
  M.def("grandparent_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*)   (HepMC3::ConstGenVertexPtr)    ) &HepMC3::grandparent_vertices,  "See documentation");

  M.def("descendant_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*)(HepMC3::ConstGenVertexPtr) ) &HepMC3::descendant_particles,  "See documentation");
  M.def("descendant_particles",   (std::vector<HepMC3::GenParticlePtr>(*)     (HepMC3::GenVertexPtr)     ) &HepMC3::descendant_particles,  "See documentation");
  M.def("descendant_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*)(HepMC3::ConstGenParticlePtr) ) &HepMC3::descendant_particles,  "See documentation");
  M.def("descendant_particles",   (std::vector<HepMC3::GenParticlePtr>(*)     (HepMC3::GenParticlePtr)      ) &HepMC3::descendant_particles,  "See documentation");

  M.def("descendant_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*) (HepMC3::ConstGenParticlePtr)   ) &HepMC3::descendant_vertices,  "See documentation");
  M.def("descendant_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)      (HepMC3::GenParticlePtr)        ) &HepMC3::descendant_vertices,  "See documentation");
  M.def("descendant_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*) (HepMC3::ConstGenVertexPtr)  ) &HepMC3::descendant_vertices,  "See documentation");
  M.def("descendant_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)      (HepMC3::GenVertexPtr)         ) &HepMC3::descendant_vertices,  "See documentation");

  M.def("ancestor_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*)(HepMC3::ConstGenVertexPtr) ) &HepMC3::ancestor_particles,  "See documentation");
  M.def("ancestor_particles",   (std::vector<HepMC3::GenParticlePtr>(*)     (HepMC3::GenVertexPtr)      ) &HepMC3::ancestor_particles,  "See documentation");
  M.def("ancestor_particles",   (std::vector<HepMC3::ConstGenParticlePtr>(*)(HepMC3::ConstGenParticlePtr) ) &HepMC3::ancestor_particles,  "See documentation");
  M.def("ancestor_particles",   (std::vector<HepMC3::GenParticlePtr>(*)     (HepMC3::GenParticlePtr)      ) &HepMC3::ancestor_particles,  "See documentation");


  M.def("ancestor_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*) (HepMC3::ConstGenParticlePtr)   ) &HepMC3::ancestor_vertices,  "See documentation");
  M.def("ancestor_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)      (HepMC3::GenParticlePtr)           ) &HepMC3::ancestor_vertices,  "See documentation");
  M.def("ancestor_vertices",   (std::vector<HepMC3::ConstGenVertexPtr>(*) (HepMC3::ConstGenVertexPtr)  ) &HepMC3::ancestor_vertices,  "See documentation");
  M.def("ancestor_vertices",   (std::vector<HepMC3::GenVertexPtr>(*)      (HepMC3::GenVertexPtr)            ) &HepMC3::ancestor_vertices,  "See documentation");

       }

} // namespace binder
