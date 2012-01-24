/** \file config_parser_test.cxx
 *  \brief A C++ program to test the config_parser.
 *  \author RJG
 *  \version 11-Feb-06
 *
 **/

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "config_parser.hh"
using namespace std;

int main() {
    bool truth_val;
    ConfigParser cfg( string("test.inp") );

    cout << setprecision(3) << showpoint;
    
    cout << "---------------------------------------------------\n";
    cout << " Test 0: Check the string representation           \n";
    cout << "---------------------------------------------------\n";

    cout << cfg;

    cout << "---------------------------------------------------\n";
    cout << " Test 1: Retrieve a string using parse_string()    \n";
    cout << "---------------------------------------------------\n";
    
    string test_str;
    string nf_str("notfound");
    cfg.parse_string( string("global"), string("test_str"),
		      test_str, nf_str );
    cout << "Value in test_str:\n";
    cout << test_str << endl;


    cout << "---------------------------------------------------\n";
    cout << " Test 2: Retrieve a vector of string using         \n";
    cout << "         parse_vector_of_strings() with valid key  \n";
    cout << "---------------------------------------------------\n";
    
    vector<string> test_vector;
    vector<string> notfound;
    notfound.push_back( string("not") );
    notfound.push_back( string("found") );
    cfg.parse_vector_of_strings( string("global"), string("test_str_vector"),
				 test_vector, notfound );
    
    cout << "Values in test_vector:\n";
    for( vector<string>::iterator it = test_vector.begin();
	 it != test_vector.end(); ++it ) {
	cout << (*it)  << " ";
    }
    cout << endl;

    
    cout << "---------------------------------------------------\n";
    cout << " Test 3: Retrieve a vector of string using         \n";
    cout << "         parse_vector_of_strings() but invalid key \n";
    cout << "---------------------------------------------------\n";
    
    test_vector.clear();
    cfg.parse_vector_of_strings( string("global"), string("invalid_key"),
				 test_vector, notfound );
    
    cout << "Values in test_vector:\n";
    for( vector<string>::iterator it = test_vector.begin();
	 it != test_vector.end(); ++it ) {
	cout << (*it)  << " ";
    }
    cout << endl;

    cout << "---------------------------------------------------\n";
    cout << " Test 4: Retrieve a single int using               \n";
    cout << "         parse_int() with a valid key              \n";
    cout << "---------------------------------------------------\n";
    
    int ival;
    truth_val = cfg.parse_int( string("global"), string("test_int"),
			       ival, 100 );
    cout << "Value in ival:\n";
    cout << ival << endl;
    cout << "truth_val= " << truth_val << endl;

    cout << "---------------------------------------------------\n";
    cout << " Test 5: Retrieve a single int using               \n";
    cout << "         parse_int() but double is found           \n";
    cout << "---------------------------------------------------\n";

    truth_val = cfg.parse_int( string("global"), string("test_dbl"),
			       ival, 100 );
    cout << "Value in ival:\n";
    cout << ival << endl;
    cout << "truth_val= " << truth_val << endl;

    cout << "---------------------------------------------------\n";
    cout << " Test 6: Retrieve a single int using               \n";
    cout << "         parse_int() but string is found           \n";
    cout << "         ival = 0 before call                      \n";
    cout << "---------------------------------------------------\n";

    ival = 0;
    truth_val = cfg.parse_int( string("global"), string("test_str"),
			       ival, 100 );
    cout << "Value in ival:\n";
    cout << ival << endl;
    cout << "truth_val= " << truth_val << endl;


    cout << "---------------------------------------------------\n";
    cout << " Test 7: Retrieve a vector of ints using           \n";
    cout << "         parse_vector_of_ints() with valid key     \n";
    cout << "---------------------------------------------------\n";
    
    vector<int> int_vec;
    vector<int> nf_int_vec;
    nf_int_vec.push_back( 0 );
    nf_int_vec.push_back( 1 );
    truth_val = cfg.parse_vector_of_ints( string("global"), string("test_int_vector"),
					       int_vec, nf_int_vec );
    
    cout << "Values in int_vec:\n";
    for( vector<int>::iterator it = int_vec.begin();
	 it != int_vec.end(); ++it ) {
	cout << (*it)  << " ";
    }
    cout << endl;
    cout << "truth_val= " << truth_val << endl;
    
    cout << "---------------------------------------------------\n";
    cout << " Test 8: Retrieve a vector of ints using           \n";
    cout << "         parse_vector_of_ints() but no ints exist  \n";
    cout << "---------------------------------------------------\n";
    
    nf_int_vec.push_back( 0 );
    nf_int_vec.push_back( 1 );
    truth_val = cfg.parse_vector_of_ints( string("global"), string("test_str_vector"),
					  int_vec, nf_int_vec );
    
    cout << "Values in int_vec:\n";
    for( vector<int>::iterator it = int_vec.begin();
	 it != int_vec.end(); ++it ) {
	cout << (*it)  << " ";
    }
    cout << endl;
    cout << "truth_val= " << truth_val << endl;

    cout << "---------------------------------------------------\n";
    cout << " Test 9: Retrieve a single double using            \n";
    cout << "         parse_double() with a valid key           \n";
    cout << "---------------------------------------------------\n";
    
    double dval;
    truth_val = cfg.parse_double( string("global"), string("test_dbl"),
				  dval, 100.0 );
    cout << "Value in dval:\n";
    cout << dval << endl;
    cout << "truth_val= " << truth_val << endl;

    cout << "---------------------------------------------------\n";
    cout << " Test 10: Retrieve a single double using           \n";
    cout << "          parse_double() but int is found          \n";
    cout << "          dval = 0 before call.                    \n";
    cout << "---------------------------------------------------\n";

    dval = 0;
    truth_val = cfg.parse_double( string("global"), string("test_int"),
				  dval, 100.0 );
    cout << "Value in dval:\n";
    cout << dval << endl;
    cout << "truth_val= " << truth_val << endl;

    cout << "---------------------------------------------------\n";
    cout << " Test 11 Retrieve a single double  using           \n";
    cout << "         parse_double() but string is found        \n";
    cout << "         dval = 0 before call                      \n";
    cout << "---------------------------------------------------\n";

    dval = 0;
    truth_val = cfg.parse_double( string("global"), string("test_str"),
				  dval, 100 );
    cout << "Value in ival:\n";
    cout << dval << endl;
    cout << "truth_val= " << truth_val << endl;


    cout << "---------------------------------------------------\n";
    cout << " Test 12: Retrieve a vector of doubles using       \n";
    cout << "          parse_vector_of_doublea() with valid key \n";
    cout << "---------------------------------------------------\n";
    
    vector<double> dbl_vec;
    vector<double> nf_dbl_vec;
    nf_dbl_vec.push_back( 0.0 );
    nf_dbl_vec.push_back( 1.0 );
    truth_val = cfg.parse_vector_of_doubles( string("global"), string("test_dbl_vector"),
					     dbl_vec, nf_dbl_vec );
    
    cout << "Values in dbl_vec:\n";
    for( vector<double>::iterator it = dbl_vec.begin();
	 it != dbl_vec.end(); ++it ) {
	cout << (*it)  << " ";
    }
    cout << endl;
    cout << "truth_val= " << truth_val << endl;
    
    cout << "---------------------------------------------------------\n";
    cout << " Test 13: Retrieve a vector of doubles using             \n";
    cout << "          parse_vector_of_doubles() but no doubles exist \n";
    cout << "---------------------------------------------------------\n";
    
    nf_dbl_vec.push_back( 0.0 );
    nf_dbl_vec.push_back( 1.0 );
    truth_val = cfg.parse_vector_of_doubles( string("global"), string("test_str_vector"),
					     dbl_vec, nf_dbl_vec );
    
    cout << "Values in dbl_vec:\n";
    for( vector<double>::iterator it = dbl_vec.begin();
	 it != dbl_vec.end(); ++it ) {
	cout << (*it)  << " ";
    }
    cout << endl;
    cout << "truth_val= " << truth_val << endl;
    
    cout << "---------------------------------------------------------\n";
    cout << " Test 14: Retrieve a vector of doubles using             \n";
    cout << "          parse_vector_of_doubles() but invalid key      \n";
    cout << "---------------------------------------------------------\n";
    
    nf_dbl_vec.push_back( 0.0 );
    nf_dbl_vec.push_back( 1.0 );
    truth_val = cfg.parse_vector_of_doubles( string("global"), string("dummy"),
					     dbl_vec, nf_dbl_vec );
    
    cout << "Values in dbl_vec:\n";
    for( vector<double>::iterator it = dbl_vec.begin();
	 it != dbl_vec.end(); ++it ) {
	cout << (*it)  << " ";
    }
    cout << endl;
    cout << "truth_val= " << truth_val << endl;
    
    

    cout << "\nDone.\n";

    return 0;
}
