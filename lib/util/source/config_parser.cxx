/** \file config_parser.cxx
 *  \brief Definitions for the ConfigParser class.
 *  \author RJG
 *  \version 12-Feb-06
 *  \version 30-Jun-08 PJ zero-length and multiword strings handled; boolean values also.
 **/

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include "config_parser.hh"
using namespace std;

// Constructor
ConfigParser::ConfigParser( string fname )
    : file_name( fname )
{
    // Now we have to do some 'real' work parsing the input file
    // and setting up the config_map
    FILE* fp = NULL;
    while ( fp == NULL ) {
        fp = fopen( fname.c_str(), "r" );
    }      
    fclose(fp);      
    ifstream infile( fname.c_str() );    
    if( ! infile ) {
	cout << "ConfigParser - unable to open file: " << fname << endl
	     << "No configuration information has been read!" << endl;
	throw runtime_error("File not found.");
    }
    
    string line, buffer, key_name;
    string hash_sign( "#" );
    string open_brace( "[" );
    string equal_sign( "=" );
    string colon( ":" );
    string semi_colon( ";" );
    vector<string> tokens;
    vector<string> val_vec;
    string sec_name = "_no_section_found_yet_";
    string filt_elems( "[]" );
    map<string, vector<string> > sec_map;

    while( getline( infile, line, '\n') ) {
	istringstream ss(line);
	// Tokenize by white space
	while( ss >> buffer )
	    tokens.push_back( buffer );
	// Now we have the string tokenized (delimited by whitespace)
	
	if( tokens.empty() ) {
	    // We've picked up a blank line in all likely hood.
	    tokens.clear();
	    continue;
	}

	// ----------- A coment line marked by # or ; ----------------- //
       
	else if( (tokens[0].substr(0,1) == hash_sign) ||
	    (tokens[0].substr(0,1) == semi_colon) ) {
	    // We've found a comment marked by # or ;
	    tokens.clear();
	    continue;
	}
	
	// ---------- A new section beginning ------------------------- //

	else if( tokens[0].substr(0,1) == open_brace ) {
	    // We've found a section beginning. Let's filter off the
	    // braces to get the section name.

	    string::size_type pos = 0;
	    while( (pos = tokens[0].find_first_of(filt_elems, pos)) != string::npos)
		tokens[0].erase(pos, 1);
	    
	    if( sec_name == "_no_section_found_yet_" ) {
		sec_name = tokens[0];
		tokens.clear();
	    }
	    else {
		// We insert the old map previously created.
		config_map.insert(map< string, map<string, vector<string> > >::value_type( sec_name, sec_map ));
		// We set the new name.
		sec_name = tokens[0];
		// We clear the old map and tokens.
		sec_map.clear();
		tokens.clear();
	    }
	}
	
	// --------- An entry within a section ----------------------- //

	else if( tokens[1].substr(0,1) == equal_sign ||
		 tokens[1].substr(0,1) == colon ) {
	    // We've found a key/val pair within a section (hopefully)
	    if( sec_name == "_no_section_found_yet_" ) {
		tokens.clear();
	    }
	    else { // We have something to do.
		key_name = tokens[0];
		for (int i = 2; i < int(tokens.size()); ++i)
		    val_vec.push_back(tokens[i]);
		sec_map.insert( map<string, vector<string> >::value_type( key_name, val_vec ));
		val_vec.clear();
		tokens.clear();
	    }
	}
	
	// -------- A line with which we know not what to do --------- //
	
	else {
	    // All we can to is clear the tokens and keep pushing on.
	    // Perhaps the user delimited the key/value pairs with a different
	    // token.  Maybe a != as a comment type declaration.
	    tokens.clear();
	}
    }
    
    if( sec_name != "_no_section_found_yet_" ) {
	// We need to inster the final section map found.
	config_map.insert(map< string, map<string, vector<string> > >::value_type( sec_name, sec_map ));
    }
    // Now we should have a filled config_map.
    return;
}

// Destructor
ConfigParser::~ConfigParser() {}

// String representation
// They say C++ can allow you to write some cryptic code.

string ConfigParser::str() const
{
    ostringstream ost;
    typedef map<string, map<string, vector<string> > > cfgmap;
    typedef map<string, vector<string> > secmap;
    
    cfgmap::const_iterator iter, iter_end;
    secmap::const_iterator mini_iter, mini_iter_end;

    vector<string>::const_iterator vec_iter, vec_iter_end;
    iter = config_map.begin();
    iter_end = config_map.end();


    while( iter != iter_end ) {
	mini_iter = (*iter).second.begin();
	mini_iter_end = (*iter).second.end();
	ost << "[" << (*iter).first << "]\n";
	while( mini_iter != mini_iter_end ) {
	    vec_iter = (*mini_iter).second.begin();
	    vec_iter_end = (*mini_iter).second.end();
	    ost << (*mini_iter).first << " = ";
	    while( vec_iter != vec_iter_end ) {
		ost << (*vec_iter) << " ";
		++vec_iter;
	    }
	    ost << "\n";
	    ++mini_iter;
	}
	++iter;
    }
    return ost.str();
}
	    

bool ConfigParser::has_section_and_key( const string section, const string key )
{
    if( config_map.count(section) == 0 ) {
	return false;
    }
    else {
	if( config_map[section].count(key) == 0 ) {
	    return false;
	}
    }
    // If we've made it this far..
    return true;
}
	

// Public functions
bool ConfigParser::parse_string( const string section, const string key,
				 string &val, const string notfound )
{
    if ( ! has_section_and_key(section, key) ) {
	val = notfound;
	// Though we didn't find anything we were able to return
	// a string in val.
	return true;
    }
    // It may be that the value is a zero-length string.
    // In that case, no value tokens were previously stored.
    // It may otherwise be the case that the value is a multiword string
    // and the tokenizer would have stripped the white-space characters.
    val = "";
    for ( size_t i = 0; i < config_map[section][key].size(); ++i ) {
	val += config_map[section][key][i];
	if ( i < config_map[section][key].size()-1 ) val += " ";
    }
    return true;
}


bool
ConfigParser::
parse_vector_of_strings( const string section, const string key,
			 vector<string> &val, const vector<string> notfound )
{

    if( ! has_section_and_key(section, key) ) {
	val = notfound;
	// Though the return value is still true as we are able
	// to return a vector of strings in val
	return true;
    }
    val = config_map[section][key];
    return true;
}

bool ConfigParser::parse_size_t( const string section, const string key,
				 size_t &val, const size_t notfound )
{
    if( ! has_section_and_key(section, key) ) {
	val = notfound;
	return true;
    }
    
    istringstream ss( config_map[section][key][0] );
    
    return ( ss >> val ? true : false );

}

bool ConfigParser::parse_int( const string section, const string key,
			      int &val, const int notfound )
{
    if( ! has_section_and_key(section, key) ) {
	val = notfound;
	return true;
    }
    
    istringstream ss( config_map[section][key][0] );
    
    return ( ss >> val ? true : false );

}

bool
ConfigParser::
parse_vector_of_ints( const string section, const string key,
		      vector<int> &val, const vector<int> notfound )
{
    if( ! has_section_and_key(section, key) ) {
	val = notfound;
	// Though we couldn't find the key/value you pair
	// we were able to put an int in val.
	return true;
    }
	
    val.clear();
    int ival;
    
    for( vector<string>::iterator it = config_map[section][key].begin();
	 it != config_map[section][key].end(); ++it ) {
	istringstream ss( (*it) );
	if( ss >> ival ) {
	    val.push_back(ival);
	}
	else {
	    // We may have some integers but there was a problem converting
	    // one of the values.
	    return false;
	}
    }
    // If we've made it this far..
    return true;

}

bool ConfigParser::parse_double( const string section, const string key,
				 double &val, const double notfound )
{
    if( ! has_section_and_key(section, key) ) {
	val = notfound;
	return true;
    }
    
    istringstream ss( config_map[section][key][0] );
    
    return ( ss >> val ? true : false );

}

bool
ConfigParser::
parse_vector_of_doubles( const string section, const string key,
			 vector<double> &val, const vector<double> notfound )
{
    if( ! has_section_and_key(section, key) ) {
	val = notfound;
	// Though we couldn't find the key/value pair in the collection,
	// we were able to put an int in val.
	return true;
    }
	
    val.clear();
    double dval;
    
    for( vector<string>::iterator it = config_map[section][key].begin();
	 it != config_map[section][key].end(); ++it ) {
	istringstream ss( (*it) );
	if( ss >> dval ) {
	    val.push_back(dval);
	}
	else {
	    // We may have some doubles but there was a problem converting
	    // at least one of the values.
	    return false;
	}
    }
    // If we've made it this far..
    return true;

}


bool ConfigParser::parse_boolean( const string section, const string key,
				  bool &val, const bool notfound )
{
    if( ! has_section_and_key(section, key) ) {
	val = notfound;
	return true;
    }
    // Python writes True/False
    // C++ needs true/false
    string value_string = config_map[section][key][0];
    for ( size_t i = 0; i < value_string.length(); ++i )
	value_string[i] = tolower(value_string[i]);
    // cout << "Value_string=" << value_string << endl;

    // istringstream ss( value_string );
    // return ( ss >> val ? true : false );
    if ( value_string == "0" || value_string == "false" ||
	 value_string == "f" ) {
	val = false;
    } else {
	// Defaults to 'true' for all other cases
	val = true;
    }
    return true;
}


ostream& operator<<( ostream &os, const ConfigParser &cfg )
{
    os << cfg.str();
    return os;
}
