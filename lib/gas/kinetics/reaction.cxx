// Author: Rowan J. Gollan
// Date: 12-Sep-2008
// Place: Hampton, Virginia, USA

#include <map>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "reaction.hh"
#include "normal-reaction.hh"
#include "third-body-reaction.hh"

using namespace std;

Reaction::
Reaction(lua_State *L, Gas_model &g)
{
    // Setup nu first.
    map<int, int> f_coeffs;
    read_table_as_map(L, -1, "f_coeffs", f_coeffs);

    map<int, int> b_coeffs;
    read_table_as_map(L, -1, "b_coeffs", b_coeffs);

    nu_ = b_coeffs;

    map<int, int>::const_iterator it;
    for ( it = f_coeffs.begin(); it != f_coeffs.end(); ++it ) {
	if ( b_coeffs.find(it->first) != b_coeffs.end() ) {
	    nu_[it->first] = nu_[it->first] - it->second;
	}
	else {
	    nu_.insert(pair<int, int>(it->first, -it->second));
	}
    }

    lua_getfield(L, -1, "frc");
    if ( lua_isnil(L, -1) ) {
	frc_ = 0;
    }
    else {
	frc_ = create_Reaction_rate_coefficient(L, g);
    }
    lua_pop(L, 1);

    lua_getfield(L, -1, "brc");
    if ( lua_isnil(L, -1) ) {
	brc_ = 0;
    }
    else {
	brc_ = create_Reaction_rate_coefficient(L, g);
    }
    lua_pop(L, 1);

    lua_getfield(L, -1, "ec");
    if ( lua_isnil(L, -1) ) {
	ec_ = 0;
    }
    else {
	ec_ = create_Equilibrium_constant(L, nu_, g);
    }
    lua_pop(L, 1);

    // Test that at least one of forward or backward rate coefficients are specified.
    if ( (frc_ == 0) && (brc_ == 0) ) {
	ostringstream ost;
	ost << "Reaction::Reaction():\n";
	ost << "Error in specification of Reaction.\n";
	ost << "At least one of: forward reaction rate coefficient\n";
	ost << "                 backward reaction rate coefficient\n";
	ost << "must be specified.  None were found.\n";
	input_error(ost);
    }

    if ( (frc_ != 0) && ( (ec_ != 0) || (brc_ != 0 ) ) ) {
	compute_kf_first_ = true;
    }
    else {
	compute_kf_first_ = false;
    }
}

Reaction::
~Reaction()
{
    delete frc_;
    delete brc_;
    delete ec_;
}

double
Reaction::
production(int isp)
{
    if ( nu_[isp] > 0 ) {
	return nu_[isp]*w_f_;
    }
    else {
	return -nu_[isp]*w_b_;
    }
}

double
Reaction::
loss(int isp)
{
    if ( nu_[isp] > 0 ) {
	return nu_[isp]*w_b_;
    }
    else {
	return  -nu_[isp]*w_f_;
    }
}

int
Reaction::
get_nu(int isp)
{
    if ( nu_.find(isp) == nu_.end() )
	return 0;
    else
	return nu_[isp];
}

void
Reaction::
get_rate_controlling_temperatures( const Gas_data &Q, double &T_f, double &T_b )
{
    if ( frc_ ) {
    	T_f = -1.0;
    	if ( ec_ ) {
    	    T_b = Q.T[ec_->get_iT()];
    	}
    	else {
    	    T_b = -1.0;
    	}
    }
    else if ( brc_ ) {
    	T_b = -1.0;
    	if ( ec_ ) {
    	    T_f = Q.T[ec_->get_iT()];
    	}
    	else {
    	   T_f = -1.0;
    	}
    }
    else {
    	T_b = -1.0;
    	T_f = -1.0;
    }
    
    return;
}

double
Reaction::
s_compute_k_f(const Gas_data &Q)
{
    if ( frc_ != 0 ) {
	if ( frc_->eval(Q) == SUCCESS ) {
	    return frc_->k();
	}
    }
    else {
	if ( ec_ != 0 ) {
	    double K_eq = ec_->eval(Q);
	    return k_b_ * K_eq;
	}
	else {
	    return 0.0;
	}
    }
    // If we got here, then the computation
    // was NOT succesful.
    cout << "Error computing forward rate coefficient.\n";
    cout << "Bailing out!\n";
    exit(BAD_REACTION_RATE_ERROR);
}

const double K_SMALL = 1.0e-50;

double
Reaction::
s_compute_k_b(const Gas_data &Q)
{
    if ( brc_ != 0 ) {
	if ( brc_->eval(Q) == SUCCESS ) {
	    return brc_->k();
	}
    }
    else {
	if ( ec_ != 0 ) {
	    double K_eq = ec_->eval(Q);
	    double k_b = K_eq < K_SMALL ? 0.0 : k_f_ / K_eq;
	    return k_b;
	}
	else {
	    return 0.0;
	}
    }
    // If we got here, then the computation
    // was NOT succesful.
    cout << "Error computing backward rate coefficient.\n";
    cout << "Bailing out!\n";
    exit(BAD_REACTION_RATE_ERROR);
}

Reaction* create_Reaction(lua_State *L, Gas_model &g)
{
    // Later implement as an object factory.
    map<string, Reaction* (*)(lua_State *, Gas_model &g)> r_types;
    r_types.insert(pair<string, Reaction* (*)(lua_State *, Gas_model &g)>("normal reaction",
									  create_Normal_reaction));
    r_types.insert(pair<string, Reaction* (*)(lua_State *, Gas_model &g)>("third body reaction",
									  create_Third_body_reaction));

    string type = get_string(L, -1, "type");

    if ( r_types.find(type) == r_types.end() ) {
	ostringstream ost;
	ost << "create_Reaction():\n";
	ost << "Error in specification of reaction type.\n";
	ost << "The selected type: " << type << " is unknown.\n";
	ost << "The available types are: " << endl;
	map<string, Reaction* (*)(lua_State*, Gas_model&)>::const_iterator it;
	for ( it = r_types.begin(); it != r_types.end(); ++it ) {
	    ost << "   " << it->first << endl;
	}
	input_error(ost);
    }
    
    Reaction* r = r_types[type](L, g);

    if ( r == 0 ) {
	ostringstream ost;
	ost << "create_Reaction():\n";
	ost << "Error trying to create reaction type: " << type << endl;
	input_error(ost);
    }

    return r;
}


Reaction* get_reaction_from_file(int ir, string cfile, Gas_model &g)
{
    // Setup lua_State for parsing.
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    // Set up species table
    lua_newtable(L);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	// This species table maps to C++ indices, because
	// it is used to setup the integer maps for
	// the reaction coefficients.
	lua_pushinteger(L, isp);
	lua_setfield(L, -2, g.species_name(isp).c_str());
    }
    // Plus add a field 'size': no of species
    lua_pushinteger(L, g.get_number_of_species());
    lua_setfield(L, -2, "size");
    lua_setglobal(L, "species");

    // Path to reaction parsing script
    char *e3bin = getenv("E3BIN");
    string home;
    if ( e3bin == NULL ) {
	// Assume default location of $HOME/e3bin
	home.append(getenv("HOME")); home.append("/e3bin");
    }
    else {
	home.append(e3bin);
    }
    string script_file(home);
    script_file.append("/reaction_parser.lua");

    if ( luaL_dofile(L, script_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error in loading script file: " << script_file << endl;
	input_error(ost);
    }
    
    // Parse the input file...
    lua_getglobal(L, "main");
    lua_pushstring(L, cfile.c_str());
    if ( lua_pcall(L, 1, 0, 0) != 0 ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error trying to load reaction scheme file: " << cfile << endl;
	ost << "Lua error message: " << lua_tostring(L, -1) << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "reaction-rate-coeff.cxx::\n";
	ost << "Error interpreting 'reactions'; a table of reactions is expected.\n";
	input_error(ost);
    }

    lua_rawgeti(L, -1, ir);

    Reaction *reac = create_Reaction(L, g);

    lua_close(L);
    return reac;
}
    
    
    
