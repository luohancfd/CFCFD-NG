#ifndef USEFUL_FUNCS_HH
#define USEFUL_FUNCS_HH

#include <string>
#include <vector>

std::vector<std::string> strip_string(std::string &line);
std::vector<double> read_time_varying_bc(std::string boundary, std::string dir, int iteration);
std::vector<double> read_time_varying_bc2(std::string bc_type, std::string boundary, std::string dir, int iteration);


#endif
