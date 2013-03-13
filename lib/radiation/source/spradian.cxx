/** \file spradian.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 8-Mar-12 : initial framework
 *
 **/
 
#include <cstring>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

#include "../../util/source/lua_service.hh"

#include "spradian.hh"
#include "radiation_constants.hh"


using namespace std;

void mystrcpy(char* fstring, size_t fstring_len, const char* cstring);

SpradianParams::SpradianParams( vector<SpradianRadiator*> &rad_vec, double lambda_min, double lambda_max, int spectral_points )
{
    // 0. Initialise the parameters that are not used or always the same
    method = 1;
    depth = 0.0;
    stand_off = 0.0;
    nose_radius = 0.0;
    nnode = 1;
    z[0] = 0.0;
    tblack[0] = 0.0;
    tblack[1] = 0.0;

    // 1. Initialise the char arrays that are constant for all calculations
    wavmin = lambda_min * 10;
    wavmax = lambda_max * 10;
    nwav = spectral_points;
    avg_num = 1;

    set_char_arrays_to_white_space();
    // NOTE: We are using ions to indicate continuum systems are requested.
    //       This is consistent with photaura and makes sense (i.e. if ions
    //       are not included in the CFD calculation then bound-free radiation
    //       should not be calculated).
    // NOTE: Also setting concentration pointers while we are at it
    for ( size_t irad=0; irad<rad_vec.size(); ++irad ) {
        if ( rad_vec[irad]->name=="C" ) {
            mystrcpy( atom_rads[0][0], 2, "C " );
            mystrcpy( atom_rads[0][1], 2, "bb" );
            rad_vec[irad]->set_conc_pointer( &(concC[0]) );
        }
        else if ( rad_vec[irad]->name=="C_plus" ) {
            mystrcpy( atom_rads[1][0], 2, "C " );
            mystrcpy( atom_rads[1][1], 2, "bf" );
            mystrcpy( atom_rads[2][0], 2, "C " );
            mystrcpy( atom_rads[2][1], 2, "bf" );
            rad_vec[irad]->set_conc_pointer( &(concCp[0]) );
        }
        else if ( rad_vec[irad]->name=="H" ) {
            mystrcpy( atom_rads[3][0], 2, "H " );
            mystrcpy( atom_rads[3][1], 2, "bb" );
            rad_vec[irad]->set_conc_pointer( &(concH[0]) );
        }
        else if ( rad_vec[irad]->name=="H_plus" ) {
            mystrcpy( atom_rads[4][0], 2, "H " );
            mystrcpy( atom_rads[4][1], 2, "bf" );
            mystrcpy( atom_rads[5][0], 2, "H " );
            mystrcpy( atom_rads[5][1], 2, "bf" );
            rad_vec[irad]->set_conc_pointer( &(concHp[0]) );
        }
        else if ( rad_vec[irad]->name=="N" ) {
            mystrcpy( atom_rads[6][0], 2, "N " );
            mystrcpy( atom_rads[6][1], 2, "bb" );
            rad_vec[irad]->set_conc_pointer( &(concN[0]) );
        }
        else if ( rad_vec[irad]->name=="N_plus" ) {
            mystrcpy( atom_rads[7][0], 2, "N " );
            mystrcpy( atom_rads[7][1], 2, "bf" );
            mystrcpy( atom_rads[8][0], 2, "N " );
            mystrcpy( atom_rads[8][1], 2, "bf" );
            rad_vec[irad]->set_conc_pointer( &(concNp[0]) );
        }
        else if ( rad_vec[irad]->name=="O" ) {
            mystrcpy( atom_rads[9][0], 2, "O " );
            mystrcpy( atom_rads[9][1], 2, "bb" );
            rad_vec[irad]->set_conc_pointer( &(concO[0]) );
        }
        else if ( rad_vec[irad]->name=="O_plus" ) {
            mystrcpy( atom_rads[10][0], 2, "O " );
            mystrcpy( atom_rads[10][1], 2, "bf" );
            mystrcpy( atom_rads[11][0], 2, "O " );
            mystrcpy( atom_rads[11][1], 2, "bf" );
            rad_vec[irad]->set_conc_pointer( &(concOp[0]) );
        }
        else if ( rad_vec[irad]->name=="C2" ) {
            mystrcpy( diatom_bands[0][0], 4, "C2  " );
            mystrcpy( diatom_bands[0][1], 4, "Swan" );
            mystrcpy( diatom_bands[1][0], 4, "C2  " );
            mystrcpy( diatom_bands[1][1], 4, "Phil" );
            mystrcpy( diatom_bands[2][0], 4, "C2  " );
            mystrcpy( diatom_bands[2][1], 4, "BL  " );
            mystrcpy( diatom_bands[3][0], 4, "C2  " );
            mystrcpy( diatom_bands[3][1], 4, "Frey" );
            mystrcpy( diatom_bands[4][0], 4, "C2  " );
            mystrcpy( diatom_bands[4][1], 4, "FH  " );
            mystrcpy( diatom_bands[5][0], 4, "C2  " );
            mystrcpy( diatom_bands[5][1], 4, "Mull" );
            mystrcpy( diatom_bands[6][0], 4, "C2  " );
            mystrcpy( diatom_bands[6][1], 4, "Db  " );
            rad_vec[irad]->set_conc_pointer( &(concC2[0]) );
        }
        else if ( rad_vec[irad]->name=="C2_plus" ) {
            mystrcpy( diatom_bands[7][0], 4, "C2nr" );
            mystrcpy( diatom_bands[7][1], 4, "cont" );
            mystrcpy( diatom_bands[8][0], 4, "C2  " );
            mystrcpy( diatom_bands[8][1], 4, "cont" );
            rad_vec[irad]->set_conc_pointer( &(concC2p[0]) );
        }
        else if ( rad_vec[irad]->name=="CN" ) {
            mystrcpy( diatom_bands[9][0], 4, "CN  " );
            mystrcpy( diatom_bands[9][1], 4, "Viol" );
            mystrcpy( diatom_bands[10][0], 4, "CN  " );
            mystrcpy( diatom_bands[10][1], 4, "Red " );
            rad_vec[irad]->set_conc_pointer( &(concCN[0]) );
        }
        else if ( rad_vec[irad]->name=="CO" ) {
            mystrcpy( diatom_bands[11][0], 4, "CO  " );
            mystrcpy( diatom_bands[11][1], 4, "4+  " );
            mystrcpy( diatom_bands[12][0], 4, "CO  " );
            mystrcpy( diatom_bands[12][1], 4, "BX  " );
            mystrcpy( diatom_bands[13][0], 4, "CO  " );
            mystrcpy( diatom_bands[13][1], 4, "CX  " );
            mystrcpy( diatom_bands[14][0], 4, "CO  " );
            mystrcpy( diatom_bands[14][1], 4, "EX  " );
            mystrcpy( diatom_bands[15][0], 4, "CO  " );
            mystrcpy( diatom_bands[15][1], 4, "FX  " );
            mystrcpy( diatom_bands[16][0], 4, "CO  " );
            mystrcpy( diatom_bands[16][1], 4, "GX  " );
            rad_vec[irad]->set_conc_pointer( &(concCO[0]) );
        }
        else if ( rad_vec[irad]->name=="CO_plus" ) {
            mystrcpy( diatom_bands[17][0], 4, "CO  " );
            mystrcpy( diatom_bands[17][1], 4, "cont" );
            rad_vec[irad]->set_conc_pointer( &(concCOp[0]) );
        }
        else if ( rad_vec[irad]->name=="H2" ) {
            mystrcpy( diatom_bands[18][0], 4, "H2  " );
            mystrcpy( diatom_bands[18][1], 4, "BX  " );
            mystrcpy( diatom_bands[19][0], 4, "H2  " );
            mystrcpy( diatom_bands[19][1], 4, "CX  " );
            rad_vec[irad]->set_conc_pointer( &(concH2[0]) );
        }
        else if ( rad_vec[irad]->name=="H_plus" ) {
            mystrcpy( diatom_bands[20][0], 4, "H2  " );
            mystrcpy( diatom_bands[20][1], 4, "cont" );
            rad_vec[irad]->set_conc_pointer( &(concHp[0]) );
        }
        else if ( rad_vec[irad]->name=="N2" ) {
            mystrcpy( diatom_bands[21][0], 4, "N2  " );
            mystrcpy( diatom_bands[21][1], 4, "1+  " );
            mystrcpy( diatom_bands[22][0], 4, "N2  " );
            mystrcpy( diatom_bands[22][1], 4, "2+  " );
            mystrcpy( diatom_bands[23][0], 4, "N2  " );
            mystrcpy( diatom_bands[23][1], 4, "W   " );
            mystrcpy( diatom_bands[24][0], 4, "N2  " );
            mystrcpy( diatom_bands[24][1], 4, "BH1 " );
            mystrcpy( diatom_bands[25][0], 4, "N2  " );
            mystrcpy( diatom_bands[25][1], 4, "BH2 " );
            mystrcpy( diatom_bands[26][0], 4, "N2  " );
            mystrcpy( diatom_bands[26][1], 4, "WJ  " );
            mystrcpy( diatom_bands[27][0], 4, "N2  " );
            mystrcpy( diatom_bands[27][1], 4, "CY  " );
            mystrcpy( diatom_bands[28][0], 4, "N2  " );
            mystrcpy( diatom_bands[28][1], 4, "LBH " );
            mystrcpy( diatom_bands[29][0], 4, "N2  " );
            mystrcpy( diatom_bands[29][1], 4, "VK  " );
            rad_vec[irad]->set_conc_pointer( &(concN2[0]) );
        }
        else if ( rad_vec[irad]->name=="N2_plus" ) {
            // NOTE: the vuv_bf routine fails for N2
            // mystrcpy( diatom_bands[30][0], 4, "N2  " );
            // mystrcpy( diatom_bands[30][1], 4, "cont" );
            mystrcpy( diatom_bands[31][0], 4, "N2+ " );
            mystrcpy( diatom_bands[31][1], 4, "1-  " );
            mystrcpy( diatom_bands[32][0], 4, "N2+ " );
            mystrcpy( diatom_bands[32][1], 4, "Mein" );
            mystrcpy( diatom_bands[33][0], 4, "N2+ " );
            mystrcpy( diatom_bands[33][1], 4, "2-  " );
            rad_vec[irad]->set_conc_pointer( &(concN2p[0]) );
        }
        else if ( rad_vec[irad]->name=="NO" ) {
            mystrcpy( diatom_bands[34][0], 4, "NO  " );
            mystrcpy( diatom_bands[34][1], 4, "Gamm" );
            mystrcpy( diatom_bands[35][0], 4, "NO  " );
            mystrcpy( diatom_bands[35][1], 4, "Beta" );
            mystrcpy( diatom_bands[36][0], 4, "NO  " );
            mystrcpy( diatom_bands[36][1], 4, "Delt" );
            mystrcpy( diatom_bands[37][0], 4, "NO  " );
            mystrcpy( diatom_bands[37][1], 4, "Epsi" );
            rad_vec[irad]->set_conc_pointer( &(concNO[0]) );
        }
        else if ( rad_vec[irad]->name=="O2" ) {
            mystrcpy( diatom_bands[38][0], 4, "O2  " );
            mystrcpy( diatom_bands[38][1], 4, "SR  " );
            rad_vec[irad]->set_conc_pointer( &(concO2[0]) );
        }
        else if ( rad_vec[irad]->name=="O2_plus" ) {
            mystrcpy( diatom_bands[39][0], 4, "O2  " );
            mystrcpy( diatom_bands[39][1], 4, "cont" );
            rad_vec[irad]->set_conc_pointer( &(concO2p[0]) );
        }
        else if ( rad_vec[irad]->name=="C3" ) {
            mystrcpy( triatom_bands[0][0], 4, "C3  " );
            mystrcpy( triatom_bands[0][1], 4, "Swin" );
            mystrcpy( triatom_bands[1][0], 4, "C3  " );
            mystrcpy( triatom_bands[1][1], 4, "vuv " );
            mystrcpy( triatom_bands[2][0], 4, "C3  " );
            mystrcpy( triatom_bands[2][1], 4, "phot" );
            rad_vec[irad]->set_conc_pointer( &(concC3[0]) );
        }
        else if ( rad_vec[irad]->name=="C2H" ) {
            mystrcpy( triatom_bands[3][0], 4, "C2H " );
            mystrcpy( triatom_bands[3][1], 4, "uv  " );
            rad_vec[irad]->set_conc_pointer( &(concC2H[0]) );
        }
        else if ( rad_vec[irad]->name=="CH" ) {
            rad_vec[irad]->set_conc_pointer( &(concCH[0]) );
        }
        else if ( rad_vec[irad]->name=="e_minus" ) {
            rad_vec[irad]->set_conc_pointer( &(concNE[0]) );
        }
        else if ( rad_vec[irad]->name=="OH" ) {
            rad_vec[irad]->set_conc_pointer( &(concOH[0]) );
        }
    }

    // 2. Initialise the parameters that will change for each calculation
    this->reset_variable_parameters();
}

SpradianParams::SpradianParams( string fname )
{
    this->read_from_file( fname );
}

SpradianParams::~SpradianParams() {}

void SpradianParams::reset_variable_parameters()
{
    conchvy[0] = 0.0;
    concC[0] = 0.0;
    concC2[0] = 0.0;
    concC2H[0] = 0.0;
    concC3[0] = 0.0;
    concCH[0] = 0.0;
    concCN[0] = 0.0;
    concCO[0] = 0.0;
    concCp[0] = 0.0;
    concH[0] = 0.0;
    concH2[0] = 0.0;
    concHp[0] = 0.0;
    concN[0] = 0.0;
    concN2[0] = 0.0;
    concN2p[0] = 0.0;
    concNE[0] = 0.0;
    concNO[0] = 0.0;
    concNp[0] = 0.0;
    concO[0] = 0.0;
    concO2[0] = 0.0;
    concOH[0] = 0.0;
    concOp[0] = 0.0;
    concatm[0] = 0.0;
    concmol[0] = 0.0;

    concC2p[0] = 0.0;
    concCOp[0] = 0.0;
    concO2p[0] = 0.0;

    tran[0] = 0.0;
    trot[0] = 0.0;
    tvib[0] = 0.0;
    telec[0] = 0.0;

    avg_molwt[0] = 0.0;
}

void SpradianParams::set_char_arrays_to_white_space()
{
    // Nonequilibrium atomic species
    for ( size_t i=0; i<SPR_N_ATM_NEQS; i++ ) {
        mystrcpy( atom_noneqs[i], 2, "  ");
    }

    // Atomic radiators
    for ( size_t i=0; i<SPR_N_ATM_RADS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            mystrcpy( atom_rads[i][j], 2, "  ");
        }
    }

    // Nonequilibrium diatomic species
    for ( size_t i=0; i<SPR_N_DTM_NEQS; i++ )
        mystrcpy( diatom_noneqs[i], 4, "    ");

    // Diatomic bands
    for ( size_t i=0; i<SPR_N_DTM_BNDS; i++ ) {
        for ( size_t j=0; j<2; j++ )
            mystrcpy( diatom_bands[i][j], 4, "    ");
    }

    // Triatomic bands
    for ( size_t i=0; i<SPR_N_TRM_BNDS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            mystrcpy( triatom_bands[i][j], 4, "    ");
        }
    }
}

void SpradianParams::write_to_file( string fname )
{
    ofstream pfile;
    pfile.open(fname.c_str());
    pfile << method << endl
          << depth << endl
          << stand_off << endl
          << nose_radius << endl
          << nnode << endl
          << z[0] << endl
          << wavmin << endl
          << wavmax << endl
          << nwav << endl
          << avg_num << endl
          << conchvy[0] << endl
          << concC[0] << endl
          << concC2[0] << endl
          << concC2H[0] << endl
          << concC3[0] << endl
          << concCH[0] << endl
          << concCN[0] << endl
          << concCO[0] << endl
          << concCp[0] << endl
          << concH[0] << endl
          << concH2[0] << endl
          << concHp[0] << endl
          << concN[0] << endl
          << concN2[0] << endl
          << concN2p[0] << endl
          << concNE[0] << endl
          << concNO[0] << endl
          << concNp[0] << endl
          << concO[0] << endl
          << concO2[0] << endl
          << concOH[0] << endl
          << concOp[0] << endl
          << concatm[0] << endl
          << concmol[0] << endl
          << concC2p[0] << endl
          << concCOp[0] << endl
          << concO2p[0] << endl
          << tran[0] << endl
          << trot[0] << endl
          << tvib[0] << endl
          << telec[0] << endl
          << avg_molwt[0] << endl
          << tblack[0] << endl
          << tblack[1] << endl;

    char three_chars[3] = { "  " };
    char five_chars[5] = { "    " };

    // Nonequilibrium atomic species
    for ( size_t i=0; i<SPR_N_ATM_NEQS; i++ ) {
        mystrcpy( three_chars, 2, atom_noneqs[i]);
        pfile << three_chars << endl;
    }

    // Atomic radiators
    for ( size_t i=0; i<SPR_N_ATM_RADS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            mystrcpy( three_chars, 2, atom_rads[i][j]);
            pfile << three_chars << endl;
        }
    }

    // Nonequilibrium diatomic species
    for ( size_t i=0; i<SPR_N_DTM_NEQS; i++ ) {
        mystrcpy( five_chars, 4, diatom_noneqs[i]);
        pfile << five_chars << endl;
    }

    // Diatomic bands
    for ( size_t i=0; i<SPR_N_DTM_BNDS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            mystrcpy( five_chars, 4, diatom_bands[i][j] );
            pfile << five_chars << endl;
        }
    }

    // Triatomic bands
    for ( size_t i=0; i<SPR_N_TRM_BNDS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            mystrcpy( five_chars, 4, triatom_bands[i][j]);
            pfile << five_chars << endl;
        }
    }

    pfile.close();
}

void SpradianParams::read_from_file( string fname )
{
    ifstream pfile(fname.c_str());
    if ( !pfile.is_open() ) {
        cout << "SpradianParams::read_from_file()" << endl
             << "Could not open spradian parameter file: " << fname << endl
             << "Exiting program." << endl;
        exit( FAILURE );
    }

    pfile >> method >> depth >> stand_off >> nose_radius >> nnode;
    pfile >> z[0] >> wavmin >> wavmax >> nwav >> avg_num;
    pfile >> conchvy[0] >> concC[0] >> concC2[0] >> concC2H[0];
    pfile >> concC3[0] >> concCH[0] >> concCN[0] >> concCO[0];
    pfile >> concCp[0] >> concH[0] >> concH2[0] >> concHp[0];
    pfile >> concN[0] >> concN2[0] >> concN2p[0] >> concNE[0];
    pfile >> concNO[0] >> concNp[0] >> concO[0] >> concO2[0];
    pfile >> concOH[0] >> concOp[0] >> concatm[0] >> concmol[0];
    pfile >> concC2p[0] >> concCOp[0] >> concO2p[0];
    pfile >> tran[0] >> trot[0] >> tvib[0] >> telec[0];
    pfile >> avg_molwt[0] >> tblack[0] >> tblack[1];

    string buff;

    // Get to the end of the current line
    getline(pfile,buff);

    // Nonequilibrium atomic species
    for ( size_t i=0; i<SPR_N_ATM_NEQS; i++ ) {
        getline(pfile,buff);
        mystrcpy( atom_noneqs[i], 2, buff.c_str());
    }

    // Atomic radiators
    for ( size_t i=0; i<SPR_N_ATM_RADS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            getline(pfile,buff);
            mystrcpy( atom_rads[i][j], 2, buff.c_str());
        }
    }

    // Nonequilibrium diatomic species
    for ( size_t i=0; i<SPR_N_DTM_NEQS; i++ ) {
        getline(pfile,buff);
        mystrcpy( diatom_noneqs[i], 4, buff.c_str());
    }

    // Diatomic bands
    for ( size_t i=0; i<SPR_N_DTM_BNDS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            getline(pfile,buff);
            mystrcpy( diatom_bands[i][j], 4, buff.c_str());
        }
    }

    // Triatomic bands
    for ( size_t i=0; i<SPR_N_TRM_BNDS; i++ ) {
        for ( size_t j=0; j<2; j++ ) {
            getline(pfile,buff);
            mystrcpy( triatom_bands[i][j], 4, buff.c_str());
        }
    }

    pfile.close();
}

int SpradianParams::call_radipac()
{
    /* Firstly test for atom.dat, diatom.dat and triatom.dat */
    const char *args[] = { "atom.dat", "diatom.dat", "triatom.dat" };
    vector<string> files(args, args+3);
    for ( size_t i=0; i<files.size(); ++i ) {
        string file = files[i];
	ifstream infile(file.c_str());
	if ( !infile.good() ) { 
	    cout << "SpradianParams::call_radipac()" << endl
	         << "The file: " << file << " is missing from the working directory." << endl
		 << "Exiting program." << endl;
	    exit( FAILURE );
        }
    }

    return radipac_( z, &nnode, &depth, &method, tblack, &stand_off, &nose_radius, \
                 &wavmin, &wavmax, &nwav, &avg_num, atom_noneqs, atom_rads, diatom_noneqs, diatom_bands, triatom_bands, \
                  tran, trot, tvib, telec, concC, concC2, concC2H, concC3, concCH, concCN, concCO, concCp, concH, concH2, \
                  concHp, concN, concN2, concN2p, concNE, concNO, concNp, concO, concO2, concOH, \
                  concOp, conchvy, concatm, concmol, avg_molwt, norm_int, dibydx, flux, dqbydx, \
                  divq );

}

Spradian::
Spradian( lua_State * L )
 : RadiationSpectralModel( L )
{
    this->initialise(L);
}

Spradian::
Spradian( const string input_file )
 : RadiationSpectralModel( input_file )
{
    // 1. Get spectral_model string from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, input_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Spradian::Spradian()\n";
	ost << "Error in input file: " << input_file << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "spectral_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Spradian::Spradian():\n";
	ost << "Error in the declaration of spectral_data: a table is expected.\n";
	input_error(ost);
    }
    
    this->initialise(L);

    lua_pop(L,1);	// pop spectral_data
    lua_close(L);
}

void Spradian::initialise( lua_State * L )
{
    iT  = get_int(L,-1,"iT");
    iTr = get_int(L,-1,"iTr");
    iTv = get_int(L,-1,"iTv");
    iTe = get_int(L,-1,"iTe");

    lua_getfield(L, -1, "radiators" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Spradian::Spradian():\n";
	ost << "Error in the declaration of radiators: a table is expected.\n";
	input_error(ost);
    }

    nrad = lua_objlen(L, -1);

    if ( nrad==0 ) {
	cout << "Spradian::Spradian()" << endl
             << "No radiators have been requested, exiting program" << endl;
        exit( BAD_INPUT_ERROR );
    }

    for ( int irad = 0; irad < nrad; ++irad ) {
	lua_rawgeti(L, -1, irad+1); // A Lua list is offset one from the C++ vector index
	const char* rad = luaL_checkstring(L, -1);
	rad_names.push_back(string(rad));
	lua_pop(L, 1);
    }

    lua_pop(L,1);	// pop radiators

    // Now construct the radiators
    for ( int irad = 0; irad < nrad; ++irad ) {
	radiators.push_back( create_new_spradian_radiator( L, rad_names[irad] ) );
	if ( radiators.back()->type == "electron_radiator" ) {
	    e_index = radiators.back()->isp;
	}
    }

    // Finally create the parameter class
    params = new SpradianParams( radiators, lambda_min, lambda_max, spectral_points );

    // create working directories for each thread
    if ( omp_get_thread_num()==0 ) {
        for ( int i=0; i<omp_get_max_threads(); ++i ) {
            ostringstream oss;
            oss << "rm -R -f spradian_working_dir_" << i;
            system(oss.str().c_str());
            oss.clear(); oss.str("");
            oss << "mkdir spradian_working_dir_" << i;
            system(oss.str().c_str());
            oss.clear(); oss.str("");
            oss << "cp $HOME/e3bin/atom.dat $HOME/e3bin/diatom.dat $HOME/e3bin/triatom.dat spradian_working_dir_" << i << "/";
            system(oss.str().c_str());
        }
    }
}

Spradian::~Spradian()
{
    // 1. Delete radiators
    for ( size_t irad=0; irad<radiators.size(); ++irad )
	delete radiators[irad];

    // 2. Delete SpradianParams class
    delete params;
}

string Spradian::str() const
{
    return "Spradian";
}

double
Spradian::
integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved )
{
    double j_total = 0.0;
    
    return j_total;
}

double
Spradian::
variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved )
{
    cout << "Spradian::variably_integrated_emission_for_gas_state()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X )
{
    // Set the adaptive flag
    X.adaptive = adaptive_spectral_grid;

    // 0. Make sure the vectors in CoeffSpectra are sized to zero.
    //    This is required for adaptive spectral distributions.
    X.nu.clear();
    X.j_nu.clear();
    X.kappa_nu.clear();

    /* 1. Initiliase radipac parameters */
    this->setup_parameters(Q);

    // 2. Move to the working directory
    ostringstream oss;
    oss << "spradian_working_dir_" << omp_get_thread_num();
    chdir(oss.str().c_str());

    /* 3. Call radipac fortran subroutine */
#   if 0
    params->write_to_file("spradian.params");
    params->call_radipac();
#   else
    params->write_to_file("spradian.params");
    oss.clear(); oss.str("");
    oss << "run_spradian.x spradian.params";
    system(oss.str().c_str());
#   endif

    // 3. Pick up the solution and insert it into CoeffSpectra
#   if 1
    ifstream specfile( "fort.10" );
    if ( !specfile.is_open() ) {
        cout << "Spradian::spectra_for_gas_state()" << endl
             << "Could not open spradian spectra file 'fort.10'." << endl
             << "Exiting program." << endl;
        exit( FAILURE );
    }

    // First line is a header
    string line;
    getline(specfile,line);
    // cout << "header: " << line << endl;

    // Remaining lines should be spectral data
    // Note that the spradian data starts from the lower wavelength whereas the
    // CoeffSpectra class starts from the highest wavelength, and spradian outputs
    // j_lambda whereas we need j_nu (hence the conversion)
    int k, m;
    double wavel, emis, absb, intens;
    double nu;
    while ( getline (specfile,line) ) {
        // cout << line << endl;
        istringstream iss(line);
        iss >> k >> m >> wavel >> emis >> absb >> intens;
        // cout << "k = " << k << ", m = " << m << ", wavel = " << wavel << ", emis = " << emis << ", absb = " << absb << ", intens = " << intens << endl;
        nu = lambda2nu( wavel/10.0 );   // Ang -> Hz
        X.nu.push_back( nu );
        X.j_nu.push_back( emis * 1.0e12 * RC_c_SI / nu / nu );  // W/cm3-um-sr -> W/m3-Hz-sr
        X.kappa_nu.push_back( absb*1.0e2 );  // cm-1 -> m-1
    }
    specfile.close();
#   else
    char line[100];
    int k, m;
    float wavel, emis, absb, intens;
    FILE *specFILE = fopen ("fort.10","r");
    if ( specFILE==NULL ) {
        cout << "Spradian::spectra_for_gas_state()" << endl
             << "Could not open spradian spectra file 'fort.10'." << endl
             << "Exiting program." << endl;
        exit( FAILURE );
    }
    // Discard the first line
    if ( fgets( line, sizeof line, specFILE ) == NULL ) {
        cout << "Spradian::spectra_for_gas_state()" << endl
             << "Problem reading data from file 'fort.10'." << endl
             << "Exiting program." << endl;
        exit( FAILURE );
    }
    // Get the spectral data
    while ( fgets( line, sizeof line, specFILE ) != NULL ) {
        if ( sscanf( line, "%d %d %f %e %e %e", &k, &m, &wavel, &emis, &absb, &intens ) != 6 ) {
            cout << "Spradian::spectra_for_gas_state()" << endl
                 << "Problem reading data from file 'fort.10'." << endl
                 << "Exiting program." << endl;
            exit( FAILURE );
        }
        cout << "k = " << k << ", m = " << m << ", wavel = " << wavel << ", emis = " << emis << ", absb = " << absb << ", intens = " << intens << endl;
        double nu = lambda2nu( wavel/10.0 );   // Ang -> Hz
        X.nu.push_back( nu );
        X.j_nu.push_back( emis * 1.0e12 * RC_c_SI / nu / nu );  // W/cm3-um-sr -> W/m3-Hz-sr
        X.kappa_nu.push_back( absb*1.0e2 );  // cm-1 -> m-1
    }
#   endif

    // Check that the spectra is the correct size
    if ( (int)X.nu.size()!=spectral_points ) {
        cout << "Spradian::spectra_for_gas_state()" << endl
             << "Only " << X.nu.size() << " spectral points were read, but " << spectral_points << " were expected." << endl
             << "Exiting program." << endl;
        exit( FAILURE );
    }

    // We want ascending frequencies for consistency with photaura model
    if ( X.nu.front() > X.nu.back() ) {
        reverse( X.nu.begin(), X.nu.end() );
        reverse( X.j_nu.begin(), X.j_nu.end() );
        reverse( X.kappa_nu.begin(), X.kappa_nu.end() );
    }

    // Move out of the working directory
    chdir("..");

    return;
}

void
Spradian::
spectral_distribution_for_gas_state(Gas_data &Q, vector<double> &nus)
{
    cout << "Spradian::spectral_distribution_for_gas_state()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
write_line_widths( Gas_data &Q )
{
    cout << "Spradian::write_line_widths()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
prep_rad_pop_files()
{
    cout << "Spradian::prep_rad_pop_files()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
append_current_rad_pops( double x )
{
    cout << "Spradian::append_current_rad_pops()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
write_QSS_analysis_files( Gas_data &Q, int index )
{
    cout << "Spradian::write_QSS_analysis_files()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void Spradian::setup_parameters( Gas_data &Q )
{
    // Temperatures
    params->tran[0]  = Q.T[iT];
    params->trot[0]  = Q.T[iTr];
    params->tvib[0]  = Q.T[iTv];
    params->telec[0] = Q.T[iTe];

    // Species concentrations
    params->concatm[0] = 0.0;
    params->concmol[0] = 0.0;
    for ( size_t irad=0; irad<radiators.size(); ++irad ) {
        // cout << "irad = " << irad << ", name = " << radiators[irad]->name << endl;
        int isp = radiators[irad]->isp;
        radiators[irad]->set_concentration( Q.rho*Q.massf[isp] );
        if ( radiators[irad]->type=="atomic radiator" ) {
            params->concatm[0] += radiators[irad]->get_concentration();
        }
        else if ( radiators[irad]->type=="diatomic radiator" ) {
            params->concmol[0] += radiators[irad]->get_concentration();
        }
        else if ( radiators[irad]->type=="triatomic radiator" ) {
            params->concmol[0] += radiators[irad]->get_concentration();
        }
    }
    params->conchvy[0] = ( Q.p - Q.p_e ) / Q.T[iT] / RC_k_SI;                                           // particles / m3
    params->avg_molwt[0] = Q.rho / ( params->conchvy[0] + params->concNE[0] ) * RC_Na * 1000.0;         // g / mol

    // cout << "concatm[0] = " << params->concatm[0] << endl;
    // cout << "concmol[0] = " << params->concmol[0] << endl;
    // cout << "conchvy[0] = " << params->conchvy[0] << endl;
    // cout << " concNE[0] = " <<  params->concNE[0] << endl;
    // cout << "avg_molwt[0] = " << params->avg_molwt[0] << endl;
}

void mystrcpy(char* fstring, size_t fstring_len, const char* cstring)
{
    size_t inlen = strlen(cstring);
    size_t cpylen = min(inlen, fstring_len);

    copy(cstring, cstring + cpylen, fstring);
    fill(fstring + cpylen, fstring + fstring_len, ' ');
}

