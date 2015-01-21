#ifndef PLOTTING_HH
#define PLOTTING_HH

#include <string>
#include "solid_block.hh"

int write_grid_to_file(SolidBlock &blk);
int write_temperatures_to_file(SolidBlock &blk, double time, int print_number);
int write_temperatures_to_file2(SolidBlock &blk);
int create_directory(std::string path);
int write_solid_times(double time, int print_number, int first_flag);

#endif
