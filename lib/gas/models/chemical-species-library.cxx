// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#include <iostream>
#include <sstream>
#include <cstdlib>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "chemical-species-library.hh"

using namespace std;

static bool initialised_flag = false;

static vector<Atomic_species*> atoms;
static vector<Fully_coupled_diatomic_species*> fully_coupled_diatoms;
static vector<Diatomic_species*> diatoms;
static vector<Polyatomic_species*> polyatoms;
static Free_electron_species* electron = 0;

static vector<Chemical_species*> species;

bool chemical_species_library_initialised()
{
    return initialised_flag;
}

int initialise_chemical_species_library( double min_massf, lua_State * L )
{
    // 1. Create species
    lua_getglobal(L, "species");

    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "initialise_chemical_species_library():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }
    
    int nsp = lua_objlen(L, -1);
    
    cout << "initialise_chemical_species_library() - creating "
         << nsp << " chemical species" << endl;

    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1);
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "initialise_chemical_species_library():\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	string species_type = get_string( L, -1, "species_type" );
	cout << "- Creating " << sp << " as a new " << species_type << " species" << endl;
	if ( species_type.find("monatomic")!=string::npos ) {
	    atoms.push_back( new Atomic_species( string(sp), species_type, isp, min_massf, L ) );
	    species.push_back( atoms.back() );
	}
	else if ( species_type.find("fully coupled diatomic")!=string::npos ) {
	    fully_coupled_diatoms.push_back( new Fully_coupled_diatomic_species( string(sp), species_type, isp, min_massf, L ) );
	    species.push_back( fully_coupled_diatoms.back() );
	}
	else if ( species_type.find("diatomic")!=string::npos ) {
	    diatoms.push_back( new Diatomic_species( string(sp), species_type, isp, min_massf, L ) );
	    species.push_back( diatoms.back() );
	}
	else if ( species_type.find("polyatomic")!=string::npos ) {
	    polyatoms.push_back( new Polyatomic_species( string(sp), species_type, isp, min_massf, L ) );
	    species.push_back( polyatoms.back() );
	}
	else if ( species_type.find("free electron")!=string::npos ) {
	    electron = new Free_electron_species( string(sp), species_type, isp, min_massf, L );
	    species.push_back( electron );
	}
	else {
	    ostringstream ost;
	    ost << "initialise_chemical_species_library():\n";
	    ost << "Could not decode type label for species " << sp << ": " << species_type << endl;
	    input_error(ost);
	}
	// else if ( species_type.find("diatomic")!=string::npos )
	// else if ( species_type.find("polyatomic")!=string::npos )

	lua_pop(L, 1);	// pop 'sp'
    }
    
    lua_pop(L, 1);	// pop 'species'
    
    // set the initialised flag
    initialised_flag = true;
    
    return SUCCESS;
}

int clear_chemical_species_library()
{
    for ( size_t ia=0; ia<atoms.size(); ++ia )
    	delete atoms[ia];
    
    for ( size_t ifcd=0; ifcd<fully_coupled_diatoms.size(); ++ifcd )
    	delete fully_coupled_diatoms[ifcd];
    
    for ( size_t id=0; id<diatoms.size(); ++id )
    	delete diatoms[id];
    
    for ( size_t ip=0; ip<polyatoms.size(); ++ip )
    	delete polyatoms[ip];
    
    if ( electron ) delete electron;
    
    atoms.clear();
    fully_coupled_diatoms.clear();
    diatoms.clear();
    polyatoms.clear();
    electron = 0;
    species.clear();
    
    initialised_flag = false;
    
    return SUCCESS;
}

int get_library_nsp()
{
    return (int) species.size();
}

Chemical_species * get_library_species_pointer( int isp )
{
    return species[isp];
}

Chemical_species * get_library_species_pointer_from_name( string name )
{
    for ( size_t isp=0; isp<species.size(); ++isp ) {
    	if ( species[isp]->get_name()==name )
    	    return species[isp];
    }
    
    cout << "get_library_species_pointer_from_name()\n";
    cout << "species: " << name << " is not in the species library.\n";
    cout << "Bailing out!\n";
    exit(BAD_INPUT_ERROR);
}

int get_library_natoms()
{
    return (int) atoms.size();
}

Atomic_species * get_library_atom_pointer( int ia )
{
    return atoms[ia];
}

int get_library_nfully_coupled_diatoms()
{
    return (int) fully_coupled_diatoms.size();
}

Fully_coupled_diatomic_species * get_library_fully_coupled_diatom_pointer( int ifcd )
{
    return fully_coupled_diatoms[ifcd];
}

Fully_coupled_diatomic_species * get_library_fully_coupled_diatom_pointer_from_name( string name )
{
    
    for ( size_t ifcd=0; ifcd<fully_coupled_diatoms.size(); ++ifcd ) {
    	if ( fully_coupled_diatoms[ifcd]->get_name()==name )
    	    return fully_coupled_diatoms[ifcd];
    }
    
    cout << "get_library_fully_coupled_diatom_pointer_from_name()\n";
    cout << "fully coupled diatom: " << name << " is not in the species library.\n";
    cout << "Bailing out!\n";
    exit(BAD_INPUT_ERROR);
}

int get_library_ndiatoms()
{
    return (int) diatoms.size();
}

Diatomic_species * get_library_diatom_pointer( int id )
{
    return diatoms[id];
}

Diatomic_species * get_library_diatom_pointer_from_name( string name )
{
    
    for ( size_t id=0; id<diatoms.size(); ++id ) {
    	if ( diatoms[id]->get_name()==name )
    	    return diatoms[id];
    }
    
    cout << "get_library_diatom_pointer_from_name()\n";
    cout << "diatom: " << name << " is not in the species library.\n";
    cout << "Bailing out!\n";
    exit(BAD_INPUT_ERROR);
}

int get_library_npolyatoms()
{
    return (int) polyatoms.size();
}

Polyatomic_species * get_library_polyatom_pointer( int ip )
{
    return polyatoms[ip];
}

Polyatomic_species * get_library_polyatom_pointer_from_name( string name )
{
    
    for ( size_t ip=0; ip<polyatoms.size(); ++ip ) {
    	if ( polyatoms[ip]->get_name()==name )
    	    return polyatoms[ip];
    }
    
    cout << "get_library_polyatom_pointer_from_name()\n";
    cout << "polyatom: " << name << " is not in the species library.\n";
    cout << "Bailing out!\n";
    exit(BAD_INPUT_ERROR);
}

int get_library_nelecs()
{
    if ( electron ) return 1;
    else return 0;
}

Free_electron_species * get_library_electron_pointer()
{
    return electron;
}

int get_library_index_from_name( string name )
{
    int index = -1;
    
    for ( size_t isp=0; isp<species.size(); ++isp ) {
    	if ( species[isp]->get_name()==name ) {
    	    index = isp;
    	    break;
    	}
    }
    
    if ( index < 0 ) {
	ostringstream ost;
	ost << "get_library_index_from_name() (see chemical-species-library.cxx)\n";
	ost << "Species with name: " << name << " could not be found.\n";
	input_error(ost);
    }
    
    return index;
}

string get_library_species_name( int isp )
{
    return species[isp]->get_name();
}

