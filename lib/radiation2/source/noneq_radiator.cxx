/** \file noneq_radiators.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 31-Mar-09
 *  \brief Definitions for nonequilibrium radiators
 *
 **/

#include "noneq_radiators.hh"
#include "amalgam.hh"

using namespace std;

static std::vector<NoneqRadiator*> noneq_radiators;

NoneqRadiator::NoneqRadiator( const string name, const string section, ConfigParser * cfg, bool radiators_present )
 : name_( name )
{
    // Check if boltzmann equilibriations have been requested
    cfg->parse_boolean( section, "boltzmann_eqs", boltz_eqs_, false );
    
    // Can only initialise rad_pointer_ if all radiators have been created
    if ( radiators_present ) {
    	irad_ = get_irad_from_name( name );
	rad_pointer_ = get_radiator_pointer( irad_ );
	// Initialise boltz_fractions_ if required
	// NOTE: -1's will remain in base positions for whole simulations
	if ( boltz_eqs_ ) boltz_fractions_.resize(rad_pointer_->ie_n, -1.0);
    }
    // else -> will be done manually in QSSNoneqSystem::complete_initialisation()
    
    cfg->parse_boolean( section, "all_levels_noneq", all_levels_noneq_, false );
    
    vector<int> noneq_elec_states;
    vector<string> noneq_elec_labels;
    
    if ( !all_levels_noneq_ ) {
	// explicitly defined noneq levels
	vector<int> nfvi;
	if ( !cfg->parse_vector_of_ints( section, "noneq_elec_states", noneq_elec_states, nfvi ) || noneq_elec_states.size()==0 ) {
	    cout << "QSSNoneqSystem::QSSNoneqSystem()" << endl
		 << "noneq_elec_states not present in section: " << section << endl
		 << "Bailing out!" << endl;
	    exit( BAD_INPUT_ERROR );
	}
	
	vector<string> nfvs;
	if ( !cfg->parse_vector_of_strings( section, "noneq_elec_labels", noneq_elec_labels, nfvs ) || noneq_elec_labels.size()==0 ) {
	    cout << "QSSNoneqSystem::QSSNoneqSystem()" << endl
		 << "noneq_elec_labels not present in section: " << section << endl
		 << "Bailing out!" << endl;
	    exit( BAD_INPUT_ERROR );
	}
	
	if ( noneq_elec_states.size() != noneq_elec_labels.size() ) {
	    cout << "QSSNoneqSystem::QSSNoneqSystem()" << endl
		 << "noneq_elec_labels and noneq_elec_states vectors in section: " << section << endl
		 << "are of different lengths." << endl
		 << "Bailing out!" << endl;
	    exit( BAD_INPUT_ERROR );
	}
    }
    else if ( radiators_present ) {
	// auto-create all electronic levels as nonequilibrium states if radiators are present
	for ( int ie=0; ie<rad_pointer_->ie_n; ++ie ) {
	    ostringstream label;
	    label << "e" << ie;
	    noneq_elec_states.push_back( ie );
	    noneq_elec_labels.push_back( label.str() );
	}
    }

    // create NoneqElecState's
    for ( size_t ins=0; ins<noneq_elec_states.size(); ++ins )
	elec_states_.push_back( new NoneqElecState( noneq_elec_states[ins], noneq_elec_labels[ins] ) );
	    
    if ( ECHO_INPUT > 0 ) {
	for ( size_t ilev=0; ilev<elec_states_.size(); ++ilev )
	    cout << " - Noneq State " << ilev << ": " << elec_states_[ilev]->label_ << endl;
    }
}

NoneqRadiator::
~NoneqRadiator()
{
    for ( size_t ins=0; ins<elec_states_.size(); ++ins )
	delete elec_states_[ins];
}

size_t NoneqRadiator::get_nelev_from_elev( int ie_index )
{
    for ( size_t ins=0; ins<elec_states_.size(); ++ins ) {
	if ( elec_states_[ins]->index_==ie_index ) return ins;
    }
    
    // if we get here the search failed
    cout << "NoneqRadiator::get_nelev_from_elev()" << endl
         << "ie_index = " << ie_index << " not found for radiator: " << rad_pointer_->rad_name << endl
         << "Bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}

NoneqRadiator * get_noneq_radiator_pointer( int ne_irad )
{
    if ( ne_irad >= int(noneq_radiators.size()) ) {
	cout << "get_noneq_radiator_pointer()" << endl
	     << "requested irad: " << ne_irad << " not in range." << endl          
	     << "Bailing out!" << endl;
	exit( BAD_INPUT_ERROR );
    }
    return noneq_radiators[ne_irad];
}

NoneqRadiator * get_noneq_radiator_pointer_from_name( string name )
{
    for ( size_t ne_irad = 0; ne_irad<noneq_radiators.size(); ++ne_irad ) {
	if ( noneq_radiators[ne_irad]->rad_pointer_->rad_name == name )
	    return noneq_radiators[ne_irad];
    }
    
    // If we get here, name search failed
    cout << "get_noneq_radiator_pointer_from_name()" << endl
         << "requested radiator: " << name << " not found." << endl          
         << "Bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}

NoneqRadiator * create_new_noneq_radiator( const string name, const string section, ConfigParser * cfg, bool radiators_present )
{
    noneq_radiators.push_back( new NoneqRadiator( name, section, cfg, radiators_present ) );
    
    // Return pointer to this newly created NoneqRadiator
    return noneq_radiators.back();
}

void clear_noneq_radiators( void ) 
{
    for ( size_t irad=0; irad<noneq_radiators.size(); ++irad )
	delete noneq_radiators[irad];
}

void decode_noneq_label( string label, int &irad, int &ielev, int &ivlev, int &ne_index )
{
    // NOTE: expecting a '_' character after the species name
    for ( size_t ne_irad=0; ne_irad<noneq_radiators.size(); ++ne_irad ) {
	string species = noneq_radiators[ne_irad]->rad_pointer_->rad_name;
	int index = -1;
	index = label.find( species + "_" );
	if ( index >= 0 ) {
	    irad = get_irad_from_name( species );
	    string level;
	    level.assign( label, species.size()+1, label.size()-species.size()-1 );
	    // string all '_' characters
	    while( 1 ) {
		int c_index = level.find('_');
		if ( c_index >= 0 ) level.erase( c_index, 1 );
		else break;
	    }
	    // determine vibrational state if present
	    int v_index = level.find('v');
	    if ( v_index>=0 ) {
		// pull out v index
		string v_str;
		v_str.assign( level, v_index+1, level.size()-v_index-1 );
		ivlev = atoi(v_str.c_str());
		// erase v index characters
		level.erase(v_index,v_str.size()+1);
	    }
	    else ivlev = -1;
	    // level should now be the e string
	    string e_str = level;
	    // search noneq elec state labels
	    for ( size_t ne_ilev=0; ne_ilev<noneq_radiators[ne_irad]->elec_states_.size(); ++ne_ilev ) {
		if ( e_str==noneq_radiators[ne_irad]->elec_states_[ne_ilev]->label_ ) {
		    ielev = noneq_radiators[ne_irad]->elec_states_[ne_ilev]->index_;
		    ne_index = noneq_radiators[ne_irad]->elec_states_[ne_ilev]->ne_index_;
		    return;
		}
	    }
	}
    }
    
    cout << "decode_noneq_label()" << endl
         << "label: " << label << " was not able to be decoded successfully." << endl
         << "Bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}

int get_number_of_noneq_radiators( void )
{
    return int(noneq_radiators.size());
}

