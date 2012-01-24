// bc_catalytic.cxx
//
// Catalytic wall BC
// Applies the wall mass-fractions due to surface catalyticity 

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

CatalyticWallBC::CatalyticWallBC( int type_code )
: type_code( type_code ) {}

CatalyticWallBC::CatalyticWallBC( CatalyticWallBC &cw )
: type_code( cw.type_code ) {}

CatalyticWallBC::~CatalyticWallBC() {}

SuperCatalyticWallBC::SuperCatalyticWallBC( vector<double> massf_wall )
: CatalyticWallBC( SUPER_CATALYTIC ), massf_wall( massf_wall )
{
    Gas_model * gm = get_gas_model_ptr();
    int nsp = gm->get_number_of_species();
    if ( (int) massf_wall.size() != nsp ) {
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

EquilibriumCatalyticWallBC::EquilibriumCatalyticWallBC( string fname )
: CatalyticWallBC( EQUIL_CATALYTIC )
{
    string line, buffer;
    istringstream ss;
    vector<string> tokens;
    int n_entries;
    
    Gas_model * gm = get_gas_model_ptr();
    int nsp = gm->get_number_of_species();

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
    
    for( int i = 0; i < MAX_EQ_WC_TABLE_ENTRIES; ++i ) {
   	if( i > ipmax ) break;
	// Read in logp and f[] information
	getline( infile, line, '\n' );
	tokens.clear();
	ss.str(line);
	while( ss >> buffer )
	    tokens.push_back( buffer );
	if( (int) tokens.size() != nsp ) {
	    cerr << "EquilibriumCatalyticWallBC::EquilibriumCatalyticWallBC()\n"
    	         << "There was a problem reading entry " << i << " in file:\n"
		 << fname << endl
		 << "The correct number of entries was not present.\n";
	    exit(BAD_INPUT_ERROR);
	}
	fC[i].resize(nsp);
	for( int isp = 0; isp < nsp; ++isp ) {
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
    int ip;

    logp = log10( Q.p );
    ip = (int) ((logp - lpmin) / dlp );
    
    /*
     * Ensure that index is in bounds of array
     */
    if( ip < 0 ) {
	ip = 0;
    }
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
