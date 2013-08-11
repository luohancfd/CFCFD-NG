/** \file poshax.cxx
 *
 *  \brief A program for computing the nonequilibrium flow behind a normal shock.
 *
 *  \author Rowan J. Gollan, DFP
 *  \version 04-Jan-07 : a port of the Python version
 *  \version 29-Mar-10 : a port of the original version to use new gas and radiation models
 *  \version 30-Jun-10 : now integrating the full ODE system in a fully coupled manner
 *
 **/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/radiation/source/LOS_pieces.hh"
#include "../../../lib/radiation/source/photaura.hh"

#include "post_shock_flow.hh"
#include "poshax_radiation_transport.hh"

using namespace std;

void print_usage()
{

    cout << "Usage poshax.x:\n";
    cout << " > poshax.x input.cfg\n";
    cout << " where input.cfg is a configuration file specifying the problem.\n";
}

int main(int argc, char *argv[])
{
    if( argc != 2 ) {
	print_usage();
	exit(1);
    }

    string input(argv[1]);
    ConfigParser cfg(input);

    string gas_model_file;
    Gas_model * gmodel = 0;
    if( ! cfg.parse_string("models", "gas_model_file", gas_model_file, 
    	                   "no-gas-model" ) ||
    	gas_model_file == "no-gas-model" ) {
	cout << "Error reading gas_model_file in [models] section of " << input
	     << endl;
	cout << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    else {
    	gmodel = create_gas_model( gas_model_file );
    }
    
    string reaction_file;
    Reaction_update * rupdate = 0;
    if( ! cfg.parse_string("models", "reaction_file", reaction_file,
    	                   "no-reactions" ) ) {
	cout << "Error reading reaction_file in [models] section of " << input
	     << endl;
	cout << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    else if ( reaction_file == "no-reactions" ) {
    	cout << "No reaction_file specified, chemical reactions inactive\n";
    }
    else {
    	rupdate = create_Reaction_update(reaction_file,*gmodel);
    }
    
    string energy_exchange_file;
    Energy_exchange_update * eeupdate = 0;
    if( ! cfg.parse_string("models", "energy_exchange_file", 
    	                    energy_exchange_file, "no-energy-exchange" ) ) {
	cout << "Error reading energy_exchange_file in [models] section of " 
	     << input << endl;
	cout << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    else if ( energy_exchange_file == "no-energy-exchange" ) {
    	cout << "No energy_exchange_file specified, thermal energy exchange"
    	     << " inactive" << endl;
    }
    else {
    	eeupdate = create_Energy_exchange_update(energy_exchange_file,*gmodel);
    }
    
    string radiation_file;
    PoshaxRadiationTransportModel * rtmodel = 0;
    if( ! cfg.parse_string("models", "radiation_file", radiation_file, 
    	                   "no-radiation" ) ) {
	cout << "Error reading radiation_file in [models] section of " << input
	     << endl;
	cout << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    else if ( radiation_file == "no-radiation" ) {
    	cout << "No radiation_file specified, radiative cooling inactive\n";
    }
    else {
    	rtmodel = create_poshax_radiation_transport_model(radiation_file);
    }

    // Default to loose coupling
    string coupling_default = "loose";
    
    string coupling_str;
    if( ! cfg.parse_string("controls", "source_term_coupling", coupling_str, 
    	                    coupling_default) ) {
	cout << "Error read source_term_coupling in [controls] section of " 
	     << input << endl
	     << "This parameter controls how the source terms are applied.\n"
	     << "The available options are 'full' or 'loose'\n"
	     << "Exiting program!\n";
	exit(BAD_INPUT_ERROR);
    }

    double dx;
    if( ! cfg.parse_double("controls", "dx", dx, 0.0) ||
	dx <= 0.0 ) {
	cout << "Error reading dx in [controls] section of " << input << endl
	     << "The value is less than or equal to zero OR it has not been"
	     << " specified." << endl
	     << "It should be a positive value that will be used as the first\n"
	     << "estimate for the stepping distance behind the shock."
	     << "Exiting program!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    bool adaptive_dx;
    if( ! cfg.parse_boolean("controls", "adaptive_dx", adaptive_dx, true) ) {
	cout << "Error reading adaptive_dx in [controls] section of " 
	     << input << endl
	     << "Exiting program!\n";
	exit(BAD_INPUT_ERROR);
    }

    double final_x;
    if( ! cfg.parse_double("controls", "final_x", final_x, 1.0e-3) ||
	final_x <= 0.0 ) {
	cout << "Error reading final_x in [controls] section of " << input 
	     << endl
	     << "The value is less than or equal to zero OR it has not been"
	     << " specified." << endl
	     << "It should be a positive value that will be used as the final\n"
	     << " distance behind the shock for computing the solution."
	     << "Exiting program!\n";
	exit(BAD_INPUT_ERROR);
    }

    double plot_dx;
    if( ! cfg.parse_double("controls", "plot_dx", plot_dx, final_x*1.0e-3) ||
	plot_dx <= 0 ) {
	cout << "Error reading plot_dx in [controls] section of " << input 
	     << endl
	     << "This value should be a positive double and indicates hown"
	     << "often a solution is written to file."
	     << "Exiting program!\n"
	     << endl;
	exit(BAD_INPUT_ERROR);
    }


    double dx_scale;
    if( ! cfg.parse_double("controls", "dx_scale", dx_scale, 1.0) ||
	plot_dx <= 0 ) {
	cout << "Error reading dx_scale in [controls] section of " << input 
	     << endl
	     << "This value should be a positive double and scales dx.\n"
	     << "Exiting program!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    string output_file_name;
    if( ! cfg.parse_string("controls", "output_file", output_file_name, 
    	                   "poshax_output.data" ) ) {
	cout << "Error reading output_file in [controls] section of " << input
	     << endl
	     << "Exiting program!\n";
	exit(BAD_INPUT_ERROR);
    }

    string species_output_type;
    if ( !cfg.parse_string("controls", "species_output", species_output_type,
			   "massf") ) {
	cout << "Error reading species_output in [controls] section of " << input
	     << endl
	     << "Exiting program!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    bool write_reaction_rates;
    if ( !cfg.parse_boolean("controls", "write_reaction_rates", write_reaction_rates,
                           false) ) {
        cout << "Error reading write_reaction_rates in [controls] section of " << input
             << endl
             << "Exiting program!\n";
        exit(BAD_INPUT_ERROR);
    }

    // Parse the radiation output options
    double rad_dx;
    cfg.parse_double("radiation", "rad_dx", rad_dx, 0.0);
    double TS_dx;
    cfg.parse_double("radiation", "TS_dx", TS_dx, 0.0);
    double path_length;
    cfg.parse_double("radiation", "path_length", path_length, 0.0);
    bool write_rad_level_pops;
    cfg.parse_boolean("radiation", "write_rad_level_pops", 
    	               write_rad_level_pops, false );
    bool write_rad_emissions;
    cfg.parse_boolean("radiation", "write_rad_emissions", 
    	               write_rad_emissions, false );
    bool write_spectra;
    cfg.parse_boolean("radiation", "write_spectra", write_spectra, false );
    vector<double> vdnf;
    vector<double> lambda_min;
    cfg.parse_vector_of_doubles("radiation", "lambda_min", lambda_min, vdnf);
    vector<double> lambda_max;
    cfg.parse_vector_of_doubles("radiation", "lambda_max", lambda_max, vdnf);
    if ( lambda_max.size()!=lambda_min.size() ) {
    	cout << "lambda_max vector is the wrong size!" << endl
    	     << "Bailing out!" << endl;
    	exit(BAD_INPUT_ERROR);
    }
    vector<double> dx_smear;
    cfg.parse_vector_of_doubles("radiation", "dx_smear", dx_smear, vdnf);
    if ( dx_smear.size()!=lambda_min.size() ) {
    	cout << "dx_smear vector is the wrong size!" << endl
    	     << "Bailing out!" << endl;
    	exit(BAD_INPUT_ERROR);
    }
    if ( !rtmodel && ( rad_dx > 0.0 || path_length > 0.0 || write_rad_level_pops 
    	|| write_rad_emissions || lambda_min.size()
    	|| lambda_max.size() || dx_smear.size() ) ) {
    	cout << "Cannot compute radiation without a radiation model!" << endl
    	     << "Bailing out!" << endl;
    	exit(BAD_INPUT_ERROR);
    }

    // Parse the pre-shock conditions.
    
    double rho_inf;
    cfg.parse_double("initial-conditions", "rho_inf", rho_inf, -1.0);
    
    double p_inf;
    cfg.parse_double("initial-conditions", "p_inf", p_inf, -1.0);
    
    if ( rho_inf < 0.0 && p_inf < 0.0 ) {
    	cout << "Error reading p_inf or rho_inf in [initial-conditions] section"
    	     << " of " << input << endl;
    	cout << "Bailing out!\n";
    	exit(BAD_INPUT_ERROR);
    }

    double u_inf;
    cfg.parse_double("initial-conditions", "u_inf", u_inf, -1.0);

    double M_inf;
    cfg.parse_double("initial-conditions", "M_inf", M_inf, -1.0);
    
    if ( u_inf < 0.0 && M_inf < 0.0 ) {
    	cout << "Error reading u_inf or M_inf in [initial-conditions] section"
    	     << " of " << input << endl;
    	cout << "Bailing out!\n";
    	exit(BAD_INPUT_ERROR);
    }
    
    vector<double> massf_inf;
    if ( !cfg.parse_vector_of_doubles("initial-conditions", "massf_inf", massf_inf, vdnf) ||
	 ( massf_inf.size() > 0 && massf_inf.size() != (size_t) gmodel->get_number_of_species() ) ) {
    	cout << "Error reading massf_inf in [initial-conditions] section of " 
    	     << input << endl;
    	cout << "Bailing out!\n";
    	exit(BAD_INPUT_ERROR);
    }
    
    vector<double> molef_inf;
    if ( !cfg.parse_vector_of_doubles("initial-conditions", "molef_inf", molef_inf, vdnf) ||
	 ( molef_inf.size() > 0 && molef_inf.size() != (size_t) gmodel->get_number_of_species() ) ) {
    	cout << "Error reading molef_inf in [initial-conditions] section of "
    	     << input << endl;
    	cout << "Bailing out!\n";
    	exit(BAD_INPUT_ERROR);
    }
    
    if ( massf_inf.size()==0 && molef_inf.size()==0 ) {
    	cout << "Neither 'molef_inf' or 'massf_inf' is present in the "
    	     << "[initial-conditions] section of " << input << endl;
    	cout << "Bailing out!\n";
    	exit(BAD_INPUT_ERROR);
    }
    else if ( molef_inf.size() > 0 ) {
    	// convert mole fractions to mass fractions
    	double MW = 0.0;
    	for ( int isp=0; isp<gmodel->get_number_of_species(); ++isp  )
    	    MW += molef_inf[isp] * gmodel->molecular_weight(isp);
    	for ( int isp=0; isp<gmodel->get_number_of_species(); ++isp  )
    	    massf_inf.push_back( molef_inf[isp]*gmodel->molecular_weight(isp) 
    	    	                 / MW );
    }

    // Make sure the mass-fractions sum to one
    scale_mass_fractions(massf_inf);

    vector<double> T_inf;
    if ( !cfg.parse_vector_of_doubles("initial-conditions", "T_inf", T_inf, vdnf)
	 || T_inf.size() != (size_t) gmodel->get_number_of_modes() ) {
    	cout << "Error reading T_inf in [initial-conditions] section of " 
    	     << input << endl;
	cout << "Incorrect number of entries for number of thermal modes in gas model.\n";
    	cout << "Bailing out!\n";
    	exit(BAD_INPUT_ERROR);
    }
    
    // Calculate the freestream density, and make sure rho_inf and u_inf is known
    
    Gas_data Q(gmodel);
    int nsp = gmodel->get_number_of_species();
    int ntm = gmodel->get_number_of_modes();
    Q.T = T_inf;
    Q.massf = massf_inf;
    if ( rho_inf < 0.0 ) {
    	Q.p = p_inf;
    	gmodel->eval_thermo_state_pT(Q);
    }
    else {
    	Q.rho = rho_inf;
    	gmodel->eval_thermo_state_rhoT(Q);
    }
    if ( u_inf < 0.0 ) {
    	u_inf = M_inf * Q.a;
    }
    
    // Create the initial flow state
    
    Flow_state initial_condition( Q, u_inf );
    
    // Initialise the post shock flow solver
    Post_shock_flow * psr = 0;
    if ( coupling_str=="loose" )
    	psr = new Loosely_coupled_post_shock_flow( initial_condition, gmodel,
    	   	   	   	   	   	   rupdate, eeupdate, 
    	   	   	   	   	   	   rtmodel );
    else if ( coupling_str=="full" )
    	psr = new Fully_coupled_post_shock_flow( initial_condition, gmodel,
    	   		   	   	   	 rupdate, eeupdate, 
    	   	   	   	   	   	 rtmodel );
    else {
    	cout << "Source term coupling model: " << coupling_str
    	     << " not recognised." << endl
    	     << "The available options are 'full' or 'loose' coupling\n"
    	     << "Bailing out!" << endl;
    	exit(BAD_INPUT_ERROR);
    }
    
    // Print out pre-shock flow conditions
    
    cout << setprecision(6) << showpoint;
    cout << "Pre-shock conditions:\n";
    cout << "---------------------\n";
    cout << "p_inf     = " << psr->icflow.Q->p << endl;
    for ( size_t iT=0; iT<psr->icflow.Q->T.size(); ++iT )
    	cout << "T_s[" << iT << "]    = " << psr->icflow.Q->T[iT] << endl;
    cout << "rho_inf   = " << psr->icflow.Q->rho << endl;
    cout << "u_inf     = " << psr->icflow.u << endl;
    
    cout << "Post-shock conditions:\n";
    cout << "----------------------\n";
    cout << "p_s       = " << psr->psflow.Q->p << endl;
    for ( size_t iT=0; iT<psr->psflow.Q->T.size(); ++iT )
    	cout << "T_s[" << iT << "]    = " << psr->psflow.Q->T[iT] << endl;
    cout << "rho_s     = " << psr->psflow.Q->rho << endl;
    cout << "u_s       = " << psr->psflow.u << endl;

    ofstream outfile;
    outfile.open(output_file_name.c_str());
    if( outfile.fail() ) {
	cout << "Error opening file: " << output_file_name << endl;
	cout << "Bailing Out!\n";
	exit(FILE_ERROR);
    }

    outfile << setprecision(12) << showpoint;

    outfile << "# " << output_file_name << endl;
    outfile << "# Columns:\n";
    int col = 1;
    outfile << "# " << col << ": x (m)\n";
    ++col;
    for ( int itm=0; itm<ntm; ++itm ) {
    	outfile << "# " << col << ": T[" << itm << "] (K)\n";
    	++col;
    }
    outfile << "# " << col << ": p (Pa)\n";
    ++col;
    outfile << "# " << col << ": rho (kg/m^3)\n";
    ++col;
    outfile << "# " << col << ": u (m/s)\n";
    ++col;
    if ( species_output_type == "molef" ) {
	for ( int isp = 0; isp < nsp; ++isp ) {
	    outfile << "# " << col << ": molef[" << isp << "]-" << gmodel->species_name(isp) << "\n";
	    ++col;
	}
    }
    else if ( species_output_type == "moles" ) {
	for ( int isp = 0; isp < nsp; ++isp ) {
	    outfile << "# " << col << ": moles[" << isp << "]-" << gmodel->species_name(isp) << "\n";
	    ++col;
	}
    }
    else {
	for ( int isp=0; isp<nsp; ++isp ) {
	    outfile << "# " << col << ": massf[" << isp << "]-" << gmodel->species_name(isp) << "\n";
	    ++col;
	}
    }
    if ( rtmodel ) {
    	outfile << "# " << col << ": Q_rad (W/m**3)\n";
    	++col;
    }

    /* Pieces for radiation calculation */
    int rad_count = 0;
    double next_spectra_x = rad_dx;
    double next_TS_x = TS_dx;
    vector<double> rad_x;
    RadiationSpectralModel * rsm = 0;
    TS_data * TS = 0;
    vector<double> divq_rad;
    if ( rad_dx > 0.0 || TS_dx > 0.0 ) {
    	rsm = rtmodel->get_rsm_pointer();
    }
    if ( TS_dx > 0.0 ) {
    	int nTS = int(final_x/TS_dx)+1;
    	cout << "init TS_data with nn = " << nTS << endl;
    	TS = new TS_data( rsm, nTS );
    	TS->T_i_ = 0.0;
    	TS->T_f_ = 0.0;
    	divq_rad.resize( nTS );
    }
    if ( write_rad_level_pops ) {
    	rsm->prep_radiator_population_files();
    }
    ofstream emissions_outfile;
    Photaura * psm = 0;
    if ( write_rad_emissions ) {
    	psm = dynamic_cast<Photaura*>(rtmodel->get_rsm_pointer());
    	emissions_outfile.open( "rad_emissions.txt" );
    	emissions_outfile << "# Column 1: x location (m)" << endl;
    	for ( int irad=0; irad<psm->get_nrad(); ++irad ) {
    	    emissions_outfile << "# Column " << irad + 2 << ": Radiator " 
    	                      << psm->get_rad_name(irad) 
    	                      << " emission (W/m**3-sr)" << endl;
    	}
	emissions_outfile << setprecision(6) << showpoint;
    }
    vector<IntensityProfile> I_vec(lambda_min.size()+1);
    
    // Prep output file for reaction rates if requested
    ofstream reaction_rates_outfile;
    if ( write_reaction_rates ) {
        reaction_rates_outfile.open( "reaction_rates.txt" );
        reaction_rates_outfile << "# Column 1: x location (m)" << endl;
        reaction_rates_outfile << "# Columns 2 onwards: net reaction rate, one column for each reaction" << endl;
        reaction_rates_outfile << setprecision(6) << showpoint;
    }

    double x = 0.0, new_dx, next_plot_x = 0.0;
    int count = 0;
    int TS_count = 0;
    outfile << setw(20) << x << psr->psflow.str(bool(rtmodel), species_output_type, gmodel->M()) << endl;
    if ( TS_dx > 0 ) {
	// tangent slab problem - create first shock point 
	cout << "creating TS point " << TS_count << " at x = " << 0 << endl;
	TS->set_rad_point( TS_count, psr->psflow.Q, &(divq_rad[TS_count]), 0.0,
	                   TS_dx );
	TS_count += 1;
    }
    while( x < final_x ) {
	new_dx = psr->increment_in_space(x, dx);
	x = x + dx;
	if ( adaptive_dx ) dx = dx_scale*new_dx;
	count++;
	if( x > next_plot_x ) {
	    outfile << setw(20) << x << psr->psflow.str(bool(rtmodel), species_output_type, gmodel->M()) << endl;
	    cout << setprecision(6) << showpoint;
	    cout << "count = " << setw(6) << count << " :: ";
	    cout << "x = " << setw(12) << x << ' ';
	    cout << "dx = " << setw(12) << dx << ' ';
	    cout << "rho = " << setw(12) << psr->psflow.Q->rho << ' ';
	    for ( size_t iT=0; iT<T_inf.size(); ++iT )
	    	cout << "T[" << iT << "] = " << setw(12) << psr->psflow.Q->T[iT] << ' ';
	    cout << "u = " << setw(12) << psr->psflow.u << ' ';
	    if ( bool(rtmodel) ) {
		double e = accumulate(psr->psflow.Q->e.begin(), psr->psflow.Q->e.end(), 0.0);
		double u = psr->psflow.u;
		double rho = psr->psflow.Q->rho;
		double E_total = rho*u*(e + 0.5*u*u) + u*psr->psflow.Q->p;
		cout << "E_total = " << setw(12) << E_total << ' ';
	    	cout << "Q_rad = " << setw(12) << psr->psflow.Q_rad << endl;
            }
	    else cout << endl;
	    if ( write_rad_emissions ) {
	    	emissions_outfile << setw(20) << x;
	    	for ( int irad=0; irad<psm->get_nrad(); ++irad )
	    	    emissions_outfile << setw(20) << 
	    	    psm->integrate_emission_spectrum_for_radiator( *psr->psflow.Q,
	    	    	   	   	   	   	   	    irad );
	    	emissions_outfile << endl;
	    }
	    if ( write_reaction_rates )
	        psr->write_reaction_rates_to_file( x, reaction_rates_outfile );
	    next_plot_x += plot_dx;
	}
	if ( next_spectra_x > 0.0 && x > next_spectra_x ) {
	    cout << "- Computing radiation spectra " << rad_count << endl;
	    LOS_data LOS( rsm, 1, 0.0, 0.0 );
	    double * div_q = 0;
	    rad_x.push_back( x );
	    LOS.set_rad_point( 0, psr->psflow.Q, div_q, path_length / 2.0, 
	    	               path_length );
	    SpectralIntensity S(rsm);
	    I_vec[0].add_new_point( x, LOS.integrate_LOS( S ) );
	    for ( size_t i=0; i<lambda_min.size(); ++i ){
	    	I_vec[i+1].add_new_point( x, 
	    	    S.integrate_intensity_spectra( lambda_min[i], 
	    	    	                           lambda_max[i] ) );
	    }
	    if ( write_rad_level_pops ) {
	    	rsm->append_current_radiator_populations( x );
	    	rsm->write_QSS_population_analysis_files( *psr->psflow.Q, 
	    	                                           rad_count );
	    }
	    if ( write_spectra ) {
	        ostringstream oss;
	        oss << "coeff-spectra-" << int(x*1.0e3) << "mm.txt";
	        LOS.get_rpoint_pointer(0)->X_->write_to_file(oss.str());
	        oss.str("");
                oss << "intensity-spectra-" << int(x*1.0e3) << "mm.txt";
                S.write_to_file(oss.str());
	    }
	    ++rad_count;
	    next_spectra_x += rad_dx;
	}
	if ( next_TS_x > 0.0 && x > next_TS_x ) {
	    // tangent slab problem
	    cout << "creating TS point " << TS_count << " at x = " << x << endl;
	    TS->set_rad_point( TS_count, psr->psflow.Q, &(divq_rad[TS_count]), 
	    	               x, TS_dx );
	    next_TS_x += TS_dx;
	    TS_count += 1;
	}
    }

    cout << "Final conditions:\n";
    cout << "-----------------\n";
    cout << "p_f     = " << psr->psflow.Q->p << endl;
    cout << "T_f     = " << psr->psflow.Q->T[0] << endl;
    cout << "rho_f   = " << psr->psflow.Q->rho << endl;
    cout << "u_f     = " << psr->psflow.u << endl;
    
    // Radiation post-processing
    if ( rad_dx > 0.0 ) {
    	// unsmeared intensity profiles
    	for ( size_t i=0; i<I_vec.size(); ++i ) {
    	    ostringstream oss;
    	    if ( i==0 ) oss << "IvX-unsmeared_total.txt";
    	    else oss << "IvX-unsmeared_" << int(lambda_min[i-1]) << "-" 
    	    	     << int(lambda_max[i-1]) << "nm.txt";
    	    I_vec[i].write_to_file( oss.str() );
    	}
    	// smeared intensity profiles
    	for ( size_t i=0; i<I_vec.size(); ++i ) {
    	    ostringstream oss;
    	    double smear_dx = 0.0;
    	    if ( dx_smear.size()==0 )
    	        smear_dx = -1.0;
    	    else if ( i==0 ) {
    	    	oss << "IvX-smeared_total.txt";
    	    	smear_dx = dx_smear[0];
    	    }
    	    else {
    	    	oss << "IvX-smeared_" << int(lambda_min[i-1]) << "-" 
    	    	    << int(lambda_max[i-1]) << "nm.txt";
    	    	smear_dx = dx_smear[i-1];
    	    }
    	    if ( smear_dx > 0.0 ) {
    	        I_vec[i].spatially_smear( smear_dx );
    	        I_vec[i].write_to_file( oss.str() );
    	    }
    	}

    }
    if ( TS_dx > 0.0 ) {
    	double q_rad = TS->quick_solve_for_divq();
    	cout << "Tangent slab flux q_rad = " << q_rad*1.0e-4 << " W/cm**2\n";
    	ofstream TS_outfile;
    	TS_outfile.open( "TS_divq.txt" );
    	TS_outfile << "# Column 1: x location (m)" << endl;
    	TS_outfile << "# Column 2: divq_rad (W/m**3)" << endl;
    	TS_outfile << "# Column 3: - 4pi * j_total (W/m**3)" << endl;
	TS_outfile << setprecision(6) << showpoint;
	for ( size_t iTS=0; iTS<divq_rad.size(); ++iTS ) {
	    double j_total = TS->rpoints_[iTS]->X_->integrate_emission_spectra();
	    // ostringstream oss;
	    // oss << "coeff-spectra-" << iTS << ".txt";
	    // TS->rpoints_[iTS]->X_->write_to_file( oss.str() );
	    TS_outfile << setw(20) << TS->rpoints_[iTS]->s_ << setw(20) 
	               << divq_rad[iTS] << setw(20) << - 4.0 * M_PI * j_total
	               << endl;
	}
	TS_outfile.close();
    }

    outfile.close();
    delete gmodel;
    if ( rupdate ) delete rupdate;
    if ( eeupdate ) delete eeupdate;
    if ( rtmodel ) delete rtmodel;
    if ( write_rad_emissions ) {
    	emissions_outfile.close();
    }
    if ( TS ) delete TS;

    if ( write_reaction_rates )
        reaction_rates_outfile.close();

    cout << "Data created in: " << output_file_name << endl;
    cout << "Done.\n";

    return 0;
}
