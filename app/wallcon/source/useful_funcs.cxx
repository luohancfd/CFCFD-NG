#include "useful_funcs.hh"

#include <iostream>
#include <sstream>
#include <iterator>
#include <iostream>
#include <stdlib.h>
#include <fstream> //class to read from files
#include <string>
#include <cstdlib>



std::vector<std::string> strip_string(std::string &line){
    /* strip_string(string) -> vector<string>
       Takes a line read by getline and strips white space returning a vector.
     */

    if (line == ""){
    	std::vector<std::string> tokens;
    	tokens.push_back("");
    	return tokens;
    }
    std::istringstream iss(line);
    std::istream_iterator<std::string> being(iss), end;
    std::vector<std::string> tokens(being, end);
    return tokens;
}

std::vector<double> read_time_varying_bc(std::string boundary, std::string dir, int iteration){
    std::vector<double> terms;
    std::string filename;
    std::vector<double> user_q;
    std::string line;

    //Set the path
    filename = dir + "/" + boundary + "/" +std::to_string(static_cast<long long>(iteration)) ;

    //Read file with terms
	std::ifstream user_file( filename.c_str() );
    if (!(user_file.is_open())){
	std::cout << "wallcon read_time_varying_bc in bc_userdef_t failed to open file at pat "<< filename << std::endl;
	exit(1);
    }

    //Create vector to pass to boundary condition.
    while ( !(user_file.eof()) ){
	getline(user_file,line);
	terms.push_back( atof(line.c_str()));
    }
    return terms;
}

std::vector<double> read_time_varying_bc2(std::string bc_type, std::string boundary, std::string dir, int iteration){
    std::vector<double> terms;
    std::string filename;
    std::vector<double> user_q;
    std::string line;

    //Set the path
    filename = dir + "/" + boundary + "/" +std::to_string(static_cast<long long>(iteration)) + "/" + bc_type ;

    //Read file with terms
	std::ifstream user_file( filename.c_str() );
    if (!(user_file.is_open())){
	std::cout << "wallcon read_time_varying_bc in bc_userdef_t failed to open file at pat "<< filename << std::endl;
	exit(1);
    }

    //Create vector to pass to boundary condition.
    while ( !(user_file.eof()) ){
	getline(user_file,line);
	terms.push_back( atof(line.c_str()));
    }
    return terms;
}

