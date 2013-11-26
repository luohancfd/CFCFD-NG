// bc_catalytic.cxx
//
// Catalytic wall BC
// Applies the wall mass-fractions due to surface catalyticity 

#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "kernel.hh"
#include "diffusion.hh"

#include "bc.hh"
#include "bc_catalytic.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"

CatalyticWallBC::CatalyticWallBC( int type_code )
: type_code( type_code ) {}

CatalyticWallBC::CatalyticWallBC( CatalyticWallBC &cw )
: type_code( cw.type_code ) {}

CatalyticWallBC::~CatalyticWallBC() {}

SuperCatalyticWallBC::SuperCatalyticWallBC( vector<double> massf_wall )
: CatalyticWallBC( SUPER_CATALYTIC ), massf_wall( massf_wall )
{
    Gas_model * gm = get_gas_model_ptr();
    size_t nsp = gm->get_number_of_species();
    if ( massf_wall.size() != nsp ) {
    	ostringstream oss;
    	oss << "SuperCatalyticWallBC::SuperCatalyticWallBC()" << endl
    	    << "The given vector of mass-fractions is not the correct length." << endl;
    	input_error( oss );
    }
}

SuperCatalyticWallBC::SuperCatalyticWallBC( SuperCatalyticWallBC &cw )
: CatalyticWallBC( SUPER_CATALYTIC ), massf_wall( cw.massf_wall ) {}

SuperCatalyticWallBC::~SuperCatalyticWallBC() {}

int
SuperCatalyticWallBC::
apply( Gas_data &Q, vector<double> &massf )
{
    for( size_t isp = 0; isp < massf.size(); ++isp ) {
    	massf[isp] = massf_wall[isp];
    }
    
    return SUCCESS;
}

PartiallyCatalyticWallBC::PartiallyCatalyticWallBC( string input_file )
: CatalyticWallBC( PARTIALLY_CATALYTIC )
{
	    gmodel = get_gas_model_ptr();
	    rupdate_wc = create_Reaction_update(input_file, *gmodel);
}

PartiallyCatalyticWallBC::PartiallyCatalyticWallBC( PartiallyCatalyticWallBC &cw )
: CatalyticWallBC( PARTIALLY_CATALYTIC ) {}

PartiallyCatalyticWallBC::~PartiallyCatalyticWallBC()
{
    delete rupdate_wc;
}

int
PartiallyCatalyticWallBC::
apply( Gas_data &Q, vector<double> &massf )
{
       Gas_data Q_copy(Q);
       global_data &G = *get_global_data_ptr();
       double dt = G.dt_global;

       double dt_suggest;

        int flag = rupdate_wc->update_state(Q_copy, dt, dt_suggest, gmodel);
        gmodel->eval_thermo_state_rhoe(Q_copy);

        // Species densities: mass of species isp per unit volume.
        for ( size_t isp = 0; isp < Q_copy.massf.size(); ++isp )
            massf[isp] = Q_copy.massf[isp];
        return flag;
}

EquilibriumCatalyticWallBC::EquilibriumCatalyticWallBC( string fname )
: CatalyticWallBC( EQUIL_CATALYTIC )
{
    string line, buffer;
    istringstream ss;
    vector<string> tokens;
    size_t n_entries;
    
    Gas_model * gm = get_gas_model_ptr();
    size_t nsp = gm->get_number_of_species();

    ifstream infile( fname.c_str() );

    /* Read information about table rows, etc. */
    getline( infile, line, '\n' );
    ss.str(line);
    while( ss >> buffer )
    	tokens.push_back( buffer );
    if( tokens.size() != 3 ) {
    	cerr << "EquilibriumCatalyticWallBC::EquilibriumCatalyticWallBC()\n"
    	     << "The look-up table file for the catalytic wall b.c.\n"
 	     << "does NOT begin with 3 entries: lpmin, dlp, n_entries\n"
 	     << "This occured while trying to read file: " << fname << endl;
 	exit(BAD_INPUT_ERROR);
    }
    ss.str(tokens[0]); ss >> lpmin;
    ss.str(tokens[1]); ss >> dlp;
    ss.str(tokens[2]); ss >> n_entries;

    if( n_entries > MAX_EQ_WC_TABLE_ENTRIES ) {
	cout << "EquilibriumCatalyticWallBC::EquilibriumCatalyticWallBC()\n"
    	     << "WARNING: only the first " << MAX_EQ_WC_TABLE_ENTRIES << " will be read in\n"
 	     << "for the catalytic wall b.c. from file: " << fname << endl
 	     << "If you require more entries, please amaned the value MAX_TABLE_ENTRIES at the\n"
 	     << "top of cns_bc.cxx and re-compile.\n";
 	ipmax = MAX_EQ_WC_TABLE_ENTRIES - 1;
    }
    else {
    	ipmax = n_entries - 1;
    }
    
    for( size_t i = 0; i < MAX_EQ_WC_TABLE_ENTRIES; ++i ) {
   	if( i > ipmax ) break;
	// Read in logp and f[] information
	getline( infile, line, '\n' );
	tokens.clear();
	ss.str(line);
	while( ss >> buffer )
	    tokens.push_back( buffer );
	if( tokens.size() != nsp ) {
	    cerr << "EquilibriumCatalyticWallBC::EquilibriumCatalyticWallBC()\n"
    	         << "There was a problem reading entry " << i << " in file:\n"
		 << fname << endl
		 << "The correct number of entries was not present.\n";
	    exit(BAD_INPUT_ERROR);
	}
	fC[i].resize(nsp);
	for( size_t isp = 0; isp < nsp; ++isp ) {
	    ss.str(tokens[isp]);
	    ss >> fC[i][isp];
	}
    }
}

EquilibriumCatalyticWallBC::EquilibriumCatalyticWallBC( EquilibriumCatalyticWallBC &cw )
: CatalyticWallBC( EQUIL_CATALYTIC ), lpmin( cw.lpmin ), 
  dlp( cw.dlp ), ipmax( cw.ipmax )
{
    for ( size_t i=0; i<MAX_EQ_WC_TABLE_ENTRIES; ++i )
    	fC[i] = cw.fC[i];
}

EquilibriumCatalyticWallBC::~EquilibriumCatalyticWallBC() {}

int
EquilibriumCatalyticWallBC::
apply( Gas_data &Q, vector<double> &massf )
{
    double logp, lpfrac;
    size_t ip;

    logp = log10( Q.p );
    ip = static_cast<int>((logp - lpmin) / dlp);
    
    /*
     * Ensure that index is in bounds of array.
     */
    // if( ip < 0 ) { // Unnecessary test for unsigned values.
    //	ip = 0;
    //}
    if( ip > (ipmax - 1) ) {
	ip = ipmax - 1;
    }

    /*
     * Find the interpolation fraction
     */
    lpfrac = (logp - (lpmin + ip * dlp)) / dlp;

    /*
     * Interpolate for each of the species mass fractions
     */
    double sum = 0.0;
    for( size_t isp = 0; isp < massf.size(); ++isp ) {
	massf[isp] = (1.0 - lpfrac) * fC[ip][isp] + lpfrac * fC[ip+1][isp];
	sum += massf[isp];
    }

    /*
     * Normalise the mass fractions as they have probably have not
     * been conserved by the interpolation.
     */
    for( size_t isp = 0; isp < massf.size(); ++isp ) massf[isp] /= sum;

    return 0;
}
