/** \file config_parser.hh
 *  \brief A C++ class for config parsing from .ini-type files.
 *  \author RJG
 *  \version 11-Feb-06
 *
 **/

#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

#ifndef CFG_PARSER_HH
#define CFG_PARSER_HH

class ConfigParser{
public:
    // Public data members
    string file_name;
    map< string, map<string, vector<string> > > config_map;

    // Constructors
    ConfigParser( string fname );

    // Destructors
    virtual ~ConfigParser();

    // Public member functions
    // These function provide the useful
    // behaviour of the class

    string str() const;

    bool parse_string( const string section, const string key,
		       string &val, const string notfound );

    bool parse_vector_of_strings( const string section, const string key,
				  vector<string> &val, const vector<string> notfound );

    bool parse_size_t( const string section, const string key,
		       size_t &val, const size_t notfound );

    bool parse_int( const string section, const string key,
		    int &val, const int notfound );

    bool parse_vector_of_ints( const string section, const string key,
			       vector<int> &val, const vector<int> notfound );

    bool parse_double( const string section, const string key,
		       double &val, const double notfound );
    bool parse_vector_of_doubles( const string section, const string key,
				  vector<double> &val, const vector<double> notfound );

    bool parse_boolean( const string section, const string key,
			bool &val, const bool notfound );

    bool has_section_and_key( const string section, const string key );

};

ostream& operator<<( ostream &os, const ConfigParser &cfg );
    
#endif
