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
Reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower)
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
	frc_ = create_Reaction_rate_coefficient(L, g, T_upper, T_lower);
    }
    lua_pop(L, 1);

    lua_getfield(L, -1, "brc");
    if ( lua_isnil(L, -1) ) {
	brc_ = 0;
    }
    else {
	brc_ = create_Reaction_rate_coefficient(L, g, T_upper, T_lower);
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

    // This flag indicates if the forward and backward rate controlling temperatures
    // are different (true) or the same (false)
    different_rcts_=false;
    if ( frc_ != 0 )
        if ( frc_->get_type().find("dissociation")!=string::npos)
            different_rcts_=true;
    if ( brc_ != 0 )
        if ( brc_->get_type().find("dissociation")!=string::npos)
            different_rcts_=true;

    // If different we will need a gas data instance to recompute the other rate
    if ( different_rcts_ )
        Q_ = new Gas_data(g.get_number_of_species(),g.get_number_of_modes());
}

Reaction::
~Reaction()
{
    delete frc_;
    delete brc_;
    delete ec_;
    if ( different_rcts_ )
        delete Q_;
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
	    if (different_rcts_)
	    {
		// The reaction rate coefficient models that have different
		// rate controlling temperatures should be functions only of
		// temperatures. Here all temperatures are set to the rate
		// controlling temperature.
		double T = ec_->get_rate_controlling_temperature(Q);
		for (size_t itm=0; itm<Q_->T.size(); itm++)
		    Q_->T[itm] = T;
		if ( brc_->eval(*Q_) == SUCCESS )
		    return brc_->k() * K_eq;
	    }
	    else return k_b_ * K_eq;
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
	    if (K_eq < K_SMALL) return 0.0;
	    if (different_rcts_)
	    {
		// The reaction rate coefficient models that have different
		// rate controlling temperatures should be functions only of
		// temperatures. Here all temperatures are set to the rate
		// controlling temperature.
		double T = ec_->get_rate_controlling_temperature(Q);
		for (size_t itm=0; itm<Q_->T.size(); itm++)
		    Q_->T[itm] = T;
		if ( frc_->eval(*Q_) == SUCCESS )
		    return frc_->k() / K_eq;
	    }
	    else return k_f_ / K_eq;
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

Reaction* create_Reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower)
{
    // Later implement as an object factory.
    map<string, Reaction* (*)(lua_State *, Gas_model &g, double, double)> r_types;
    r_types.insert(pair<string, Reaction* (*)(lua_State *, Gas_model &g, double, double)>("normal reaction",
									  create_Normal_reaction));
    r_types.insert(pair<string, Reaction* (*)(lua_State *, Gas_model &g, double, double)>("third body reaction",
									  create_Third_body_reaction));

    string type = get_string(L, -1, "type");

    if ( r_types.find(type) == r_types.end() ) {
	ostringstream ost;
	ost << "create_Reaction():\n";
	ost << "Error in specification of reaction type.\n";
	ost << "The selected type: " << type << " is unknown.\n";
	ost << "The available types are: " << endl;
	map<string, Reaction* (*)(lua_State*, Gas_model&, double, double)>::const_iterator it;
	for ( it = r_types.begin(); it != r_types.end(); ++it ) {
	    ost << "   " << it->first << endl;
	}
	input_error(ost);
    }
    
    Reaction* r = r_types[type](L, g, T_upper, T_lower);

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

    lua_getglobal(L, "scheme_t");
    lua_getfield(L, -1, "temperature_limits");
    double T_lower = get_positive_number(L, -1, "lower");
    double T_upper = get_positive_number(L, -1, "upper");
    lua_pop(L, 1);
    lua_pop(L, 1); // pop scheme_t

    lua_getglobal(L, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "reaction-rate-coeff.cxx::\n";
	ost << "Error interpreting 'reactions'; a table of reactions is expected.\n";
	input_error(ost);
    }

    lua_rawgeti(L, -1, ir);

    Reaction *reac = create_Reaction(L, g, T_upper, T_lower);

    lua_close(L);
    return reac;
}
    
    
    
