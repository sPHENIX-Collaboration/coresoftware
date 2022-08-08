// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_HEPEVT_WRAPPER_TEMPLATE_H
#define HEPMC3_HEPEVT_WRAPPER_TEMPLATE_H
#include <iostream>
#include <cstdio>
#include <set>
#include <map>
#include <cstring> // memset
#include <cassert>
#include <algorithm> //min max for VS2017
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/HEPEVT_Helpers.h"

/**
 *  @file HEPEVT_Wrapper_Template.h
 *  @brief Definition of \b class HEPEVT_Wrapper_Template
 *
 *  @class HepMC3::HEPEVT_Wrapper_Template
 *  @brief An interface to HEPEVT common block implemented as template class
 */
namespace HepMC3
{

template <int max_particles, typename momentum_type = double>
class HEPEVT_Wrapper_Template
{
//
// Functions
//
public:
    /** @brief Default constructor */
    HEPEVT_Wrapper_Template() { m_hepevtptr=nullptr;};
    /** @brief Default destructor */
    ~HEPEVT_Wrapper_Template() {};
    /** @brief Print information from HEPEVT common block */
    void print_hepevt( std::ostream& ostr = std::cout )  const;
    /** @brief Print particle information */
    void print_hepevt_particle( int index, std::ostream& ostr = std::cout )  const;
    /** @brief Set all entries in HEPEVT to zero */
    void zero_everything();
    /** @brief Convert GenEvent to HEPEVT*/
    bool GenEvent_to_HEPEVT( const GenEvent* evt )  { return GenEvent_to_HEPEVT_nonstatic(evt, this);};
    /** @brief Convert HEPEVT to GenEvent*/
    bool HEPEVT_to_GenEvent( GenEvent* evt )  const { return HEPEVT_to_GenEvent_nonstatic(evt, this);};
    /** @brief Tries to fix list of daughters */
    bool fix_daughters();
    /** @brief  Fortran common block HEPEVT */
    struct HEPEVT_Templated<max_particles, momentum_type>* m_hepevtptr;
private:
    /** @brief  Internalstorage storage. Optional.*/
    std::shared_ptr<struct HEPEVT_Templated<max_particles, momentum_type> > m_internal_storage;
//
// Accessors
//
public:
    void   allocate_internal_storage(); //!< Allocates m_internal_storage storage in smart pointer to hold HEPEVT of fixed size
    void   copy_to_internal_storage( char *c, int N ); //!< Copies the content of foreign common block into the internal storage
    void   set_max_number_entries( unsigned int size ) { if (size != max_particles) printf("This implementation does not support change of the block size.\n"); assert(size == max_particles); }//!< Set block size
    void   set_hepevt_address(char *c) { m_hepevtptr = (struct HEPEVT_Templated<max_particles, momentum_type>*)c;          } //!< Set Fortran block address
    int    max_number_entries()   const    { return max_particles;                              } //!< Block size
    int    event_number()       const      { return m_hepevtptr->nevhep;             } //!< Get event number
    int    number_entries()     const      { return m_hepevtptr->nhep;               } //!< Get number of entries
    int    status(const int index )  const      { return m_hepevtptr->isthep[index-1];    } //!< Get status code
    int    id(const int index )     const       { return m_hepevtptr->idhep[index-1];     } //!< Get PDG particle id
    int    first_parent(const int index )  const { return m_hepevtptr->jmohep[index-1][0]; } //!< Get index of 1st mother
    int    last_parent(const int index )  const { return m_hepevtptr->jmohep[index-1][1]; } //!< Get index of last mother
    int    first_child(const int index )   const { return m_hepevtptr->jdahep[index-1][0]; } //!< Get index of 1st daughter
    int    last_child(const int index )   const { return m_hepevtptr->jdahep[index-1][1]; } //!< Get index of last daughter
    double px(const int index )    const        { return m_hepevtptr->phep[index-1][0];   } //!< Get X momentum
    double py(const int index )    const        { return m_hepevtptr->phep[index-1][1];   } //!< Get Y momentum
    double pz(const int index )    const        { return m_hepevtptr->phep[index-1][2];   } //!< Get Z momentum
    double e(const int index )     const        { return m_hepevtptr->phep[index-1][3];   } //!< Get Energy
    double m(const int index )     const        { return m_hepevtptr->phep[index-1][4];   } //!< Get generated mass
    double x(const int index )     const        { return m_hepevtptr->vhep[index-1][0];   } //!< Get X Production vertex
    double y(const int index )     const        { return m_hepevtptr->vhep[index-1][1];   } //!< Get Y Production vertex
    double z(const int index )      const       { return m_hepevtptr->vhep[index-1][2];   } //!< Get Z Production vertex
    double t(const int index )      const       { return m_hepevtptr->vhep[index-1][3];   } //!< Get production time
    int    number_parents(const int index ) const;                                   //!< Get number of parents
    int    number_children(const int index ) const;                                  //!< Get number of children from the range of daughters
    int    number_children_exact(const int index ) const;                            //!< Get number of children by counting
    void set_event_number( const int evtno )       { m_hepevtptr->nevhep = evtno;         } //!< Set event number
    void set_number_entries( const int noentries ) { m_hepevtptr->nhep = noentries;       } //!< Set number of entries
    void set_status( const int index, const int status ) { m_hepevtptr->isthep[index-1] = status; } //!< Set status code
    void set_id(const int index, const int id )         { m_hepevtptr->idhep[index-1] = id;      } //!< Set PDG particle id
    void set_parents( const int index, const int firstparent, const int lastparent );              //!< Set parents
    void set_children( const int index, const int firstchild, const int lastchild );               //!< Set children
    void set_momentum( const int index, const double px, const double py, const double pz, const double e );   //!< Set 4-momentum
    void set_mass( const int index, double mass );                                     //!< Set mass
    void set_position( const int index, const double x, const double y, const double z, const double t );      //!< Set position in time-space
};

//
// inline definitions
//
template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::print_hepevt( std::ostream& ostr )  const
{
    ostr << " Event No.: " << m_hepevtptr->nevhep << std::endl;
    ostr << "  Nr   Type   Parent(s)  Daughter(s)      Px       Py       Pz       E    Inv. M." << std::endl;
    for ( int i = 1; i <= m_hepevtptr->nhep; ++i )
    {
        print_hepevt_particle( i, ostr );
    }
}
template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::print_hepevt_particle( int index, std::ostream& ostr )  const
{
    char buf[255];//Note: the format is fixed, so no reason for complicated treatment

    sprintf(buf, "%5i %6i", index, m_hepevtptr->idhep[index-1]);
    ostr << buf;
    sprintf(buf, "%4i - %4i  ", m_hepevtptr->jmohep[index-1][0], m_hepevtptr->jmohep[index-1][1]);
    ostr << buf;
    sprintf(buf, "%4i - %4i ", m_hepevtptr->jdahep[index-1][0], m_hepevtptr->jdahep[index-1][1]);
    ostr << buf;
    sprintf(buf, "%8.2f %8.2f %8.2f %8.2f %8.2f", m_hepevtptr->phep[index-1][0], m_hepevtptr->phep[index-1][1], m_hepevtptr->phep[index-1][2], m_hepevtptr->phep[index-1][3], m_hepevtptr->phep[index-1][4]);
    ostr << buf << std::endl;
}

template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::allocate_internal_storage()
{
    m_internal_storage = std::make_shared<struct HEPEVT_Templated<max_particles, momentum_type>>();
    m_hepevtptr = m_internal_storage.get();
}

template <int max_particles, typename momentum_type>
void HEPEVT_Wrapper_Template<max_particles, momentum_type>::copy_to_internal_storage(char *c, int N)
{
    if ( N < 1 || N > max_particles) return;
    m_internal_storage = std::make_shared<struct HEPEVT_Templated<max_particles, momentum_type>>();
    m_hepevtptr = m_internal_storage.get();
    char* x = c;
    m_hepevtptr->nevhep = *((int*)x);
    x += sizeof(int);
    m_hepevtptr->nhep = *((int*)x);
    x += sizeof(int);
    memcpy(m_hepevtptr->isthep, x, N*sizeof(int));
    x += sizeof(int)*N;
    memcpy(m_hepevtptr->idhep, x, N*sizeof(int));
    x += sizeof(int)*N;
    memcpy(m_hepevtptr->jmohep, x, 2*N*sizeof(int));
    x += sizeof(int)*N*2;
    memcpy(m_hepevtptr->jdahep, x, 2*N*sizeof(int));
    x += sizeof(int)*N*2;
    memcpy(m_hepevtptr->phep, x, 5*N*sizeof(momentum_type));
    x += sizeof(momentum_type)*N*5;
    memcpy(m_hepevtptr->vhep, x, 4*N*sizeof(momentum_type));
}

template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::zero_everything()
{
    memset(m_hepevtptr, 0, sizeof(struct HEPEVT_Templated<max_particles, momentum_type>));
}

template <int max_particles, typename momentum_type>
inline int HEPEVT_Wrapper_Template<max_particles, momentum_type>::number_parents( const int index )  const
{
    return (m_hepevtptr->jmohep[index-1][0]) ? (m_hepevtptr->jmohep[index-1][1]) ? m_hepevtptr->jmohep[index-1][1]-m_hepevtptr->jmohep[index-1][0] : 1 : 0;
}

template <int max_particles, typename momentum_type>
inline int HEPEVT_Wrapper_Template<max_particles, momentum_type>::number_children( const int index )  const
{
    return (m_hepevtptr->jdahep[index-1][0]) ? (m_hepevtptr->jdahep[index-1][1]) ? m_hepevtptr->jdahep[index-1][1]-m_hepevtptr->jdahep[index-1][0] : 1 : 0;
}

template <int max_particles, typename momentum_type>
inline int HEPEVT_Wrapper_Template<max_particles, momentum_type>::number_children_exact( const int index )  const
{
    int nc = 0;
    for ( int i = 1; i <= m_hepevtptr->nhep; ++i )
        if (((m_hepevtptr->jmohep[i-1][0] <= index && m_hepevtptr->jmohep[i-1][1] >= index)) || (m_hepevtptr->jmohep[i-1][0] == index) || (m_hepevtptr->jmohep[i-1][1]==index)) nc++;
    return nc;
}

template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::set_parents( const int index,  const int firstparent, const int lastparent )
{
    m_hepevtptr->jmohep[index-1][0] = firstparent;
    m_hepevtptr->jmohep[index-1][1] = lastparent;
}

template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::set_children(  const int index,  const int  firstchild,  const int lastchild )
{
    m_hepevtptr->jdahep[index-1][0] = firstchild;
    m_hepevtptr->jdahep[index-1][1] = lastchild;
}

template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::set_momentum( const int index, const double px, const double py, const double pz, const double e )
{
    m_hepevtptr->phep[index-1][0] = px;
    m_hepevtptr->phep[index-1][1] = py;
    m_hepevtptr->phep[index-1][2] = pz;
    m_hepevtptr->phep[index-1][3] = e;
}

template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::set_mass( const int index, double mass )
{
    m_hepevtptr->phep[index-1][4] = mass;
}

template <int max_particles, typename momentum_type>
inline void HEPEVT_Wrapper_Template<max_particles, momentum_type>::set_position( const int index, const double x, const double y, const double z, const double t )
{
    m_hepevtptr->vhep[index-1][0] = x;
    m_hepevtptr->vhep[index-1][1] = y;
    m_hepevtptr->vhep[index-1][2] = z;
    m_hepevtptr->vhep[index-1][3] = t;
}

template <int max_particles, typename momentum_type>
inline bool HEPEVT_Wrapper_Template<max_particles, momentum_type>::fix_daughters()
{
    /*AV The function should be called  for a record that has correct particle ordering and mother ids.
    As a result it produces a record with ranges where the daughters can be found.
    Not every particle in the range will be a daughter. It is true only for proper events.
    The return tells if the record was fixed succesfully.
    */
    for ( int i = 1; i <= number_entries(); i++ )
        for ( int k=1; k <= number_entries(); k++ ) if (i != k)
                if ((first_parent(k) <= i) && (i <= last_parent(k)))
                    set_children(i, (first_child(i) == 0 ? k : std::min(first_child(i), k)), (last_child(i) == 0 ? k : std::max(last_child(i), k)));
    bool is_fixed = true;
    for ( int i = 1; i <= number_entries(); i++ )
        is_fixed = (is_fixed && (number_children_exact(i) == number_children(i)));
    return is_fixed;
}

} // namespace HepMC3
#endif
