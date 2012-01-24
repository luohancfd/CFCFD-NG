// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#ifndef CHEMICAL_SPECIES_LIBRARY_HH
#define CHEMICAL_SPECIES_LIBRARY_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "chemical-species.hh"

bool chemical_species_library_initialised();

int initialise_chemical_species_library( double min_massf, lua_State * L );

int clear_chemical_species_library();

int get_library_nsp();

Chemical_species * get_library_species_pointer( int isp );

Chemical_species * get_library_species_pointer_from_name( std::string name );

int get_library_natoms();

Atomic_species * get_library_atom_pointer( int ia );

int get_library_nfully_coupled_diatoms();

Fully_coupled_diatomic_species * get_library_fully_coupled_diatom_pointer( int id );

Fully_coupled_diatomic_species * get_library_fully_coupled_diatom_pointer_from_name( std::string name );

int get_library_ndiatoms();

Diatomic_species * get_library_diatom_pointer( int id );

Diatomic_species * get_library_diatom_pointer_from_name( std::string name );

int get_library_npolyatoms();

Polyatomic_species * get_library_polyatom_pointer( int ip );

Polyatomic_species * get_library_polyatom_pointer_from_name( std::string name );

int get_library_nelecs();

Free_electron_species * get_library_electron_pointer();

int get_library_index_from_name( std::string name );

std::string get_library_species_name( int isp );

#endif
