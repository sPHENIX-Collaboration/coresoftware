#include <HepMC3/Data/GenParticleData.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/LHEF.h>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream> // __str__
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

void bind_pyHepMC3_3(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<std::shared_ptr<HepMC3::GenParticle>,std::allocator<std::shared_ptr<HepMC3::GenParticle> >>(M("std"), "std_shared_ptr_HepMC3_GenParticle_t", "std_allocator_std_shared_ptr_HepMC3_GenParticle_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<char,std::allocator<char>>(M("std"), "char", "std_allocator_char_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<float,std::allocator<float>>(M("std"), "float", "std_allocator_float_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<long double,std::allocator<long double>>(M("std"), "long_double", "std_allocator_long_double_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<long long,std::allocator<long long>>(M("std"), "long_long", "std_allocator_long_long_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<unsigned int,std::allocator<unsigned int>>(M("std"), "unsigned_int", "std_allocator_unsigned_int_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<unsigned long,std::allocator<unsigned long>>(M("std"), "unsigned_long", "std_allocator_unsigned_long_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<unsigned long long,std::allocator<unsigned long long>>(M("std"), "unsigned_long_long", "std_allocator_unsigned_long_long_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<int,std::allocator<int>>(M("std"), "int", "std_allocator_int_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<long,std::allocator<long>>(M("std"), "long", "std_allocator_long_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<double,std::allocator<double>>(M("std"), "double", "std_allocator_double_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<std::string,std::allocator<std::string >>(M("std"), "std_string", "std_allocator_std_string_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<std::shared_ptr<HepMC3::GenVertex>,std::allocator<std::shared_ptr<HepMC3::GenVertex> >>(M("std"), "std_shared_ptr_HepMC3_GenVertex_t", "std_allocator_std_shared_ptr_HepMC3_GenVertex_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<LHEF::XMLTag *,std::allocator<LHEF::XMLTag *>>(M("std"), "LHEF_XMLTag_*", "std_allocator_LHEF_XMLTag_*_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<std::vector<double>,std::allocator<std::vector<double> >>(M("std"), "std_vector_double_t", "std_allocator_std_vector_double_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<LHEF::WeightInfo,std::allocator<LHEF::WeightInfo>>(M("std"), "LHEF_WeightInfo", "std_allocator_LHEF_WeightInfo_t");

	// std::vector file:bits/stl_vector.h line:339
	binder::vector_binder<LHEF::HEPEUP *,std::allocator<LHEF::HEPEUP *>>(M("std"), "LHEF_HEPEUP_*", "std_allocator_LHEF_HEPEUP_*_t");

}
